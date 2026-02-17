-- | Pure generation forest for fullerene enumeration.
--
-- The search forest has three roots (C20, C28, C30). Each node's children
-- are computed by: bounding → expansion enumeration → Rule 2 filtering →
-- pure graph surgery → canonical test → automorphism group computation.
--
-- The forest is lazy: children are computed on demand. This enables
-- different traversal strategies (DFS, BFS, parallel) without changing
-- the generation logic.
--
-- For parallel evaluation: each subtree is independent (no shared state).
-- Apply strategy annotations (par, parList, parBuffer) to 'gnChildren'
-- at the desired fork depth.

module GenForest
    ( -- * Configuration
      GenConfig(..)
    , mkConfig
      -- * Seeds
    , seedsWithAuts
      -- * Pure child generation
    , generateChildren
      -- * Lazy forest
    , GenTree(..)
    , genForest
      -- * Traversal schemes
    , dfsCountBySize
    , parCountBySize
    , parMutCountBySize
      -- ** Lazy DFS stream
    , allGraphs
      -- ** BFS (breadth-first)
    , bfsByLevel
    , bfsCountBySize
      -- ** Generic fold
    , foldForest
      -- ** Parallel map
    , parMapForest
      -- ** Target-size generation
    , graphsOfSize
      -- ** Work-queue (dynamic load balance)
    , workQueueCountBySize
      -- * Analysis
    , depthProfile
    , subtreeSizes
    ) where

import Seeds (DualGraph(..), c20, c28, c30)
import Expansion
import Canonical
import MutGraph (newMutGraph, loadGraph, freezeGraph, unsafeFreezeGraph,
                 applyExpansionM, applyRingM, undoMutation)
import Data.Array.Unboxed (UArray)
import Data.Map.Strict (Map)
import qualified Data.Map.Strict as Map
import Data.List (foldl', sortBy)
import Data.Ord (Down(..))
import Control.Monad (foldM, forM_, replicateM)
import Control.Monad.ST (runST)
import Control.Parallel.Strategies (using, parList, rdeepseq)
import Control.DeepSeq (NFData(..))
import Control.Concurrent (forkIO, MVar, newEmptyMVar, putMVar, takeMVar)
import Control.Concurrent.STM
import Control.Exception (SomeException, try)
import Data.IORef

---------------------------------------------------------------------
-- Configuration
---------------------------------------------------------------------

data GenConfig = GenConfig
    { gcMaxDV    :: !Int              -- max dual vertices (C_n/2 + 2)
    , gcMslTable :: !(UArray Int Int) -- precomputed geometric bounds
    }

-- | Build a GenConfig for a given maximum dual vertex count.
mkConfig :: Int -> GenConfig
mkConfig maxDV = GenConfig maxDV (buildMaxStraightLengths maxDV)

---------------------------------------------------------------------
-- Seeds
---------------------------------------------------------------------

-- | The three seed graphs with their automorphism groups.
-- These are the roots of the generation forest.
seedsWithAuts :: [(DualGraph, [Automorphism])]
seedsWithAuts =
    [ (g, canonAuts (canonicalBFSAndGroup g))
    | g <- [c20, c28, c30]
    ]

---------------------------------------------------------------------
-- Pure child generation
---------------------------------------------------------------------

-- | Compute all canonical children of a graph with their automorphism groups.
--
-- This is the core generation step, completely pure:
--   1. Enumerate reductions (for bounding lemmas)
--   2. Compute maximum expansion length from bounds
--   3. Enumerate expansion sites within bound
--   4. Filter by Rule 2 (automorphism orbit equivalence)
--   5. Apply each expansion (pure graph surgery)
--   6. Canonical test: accept only if inverse is the canonical reduction
--   7. Compute automorphism group for accepted children
--
-- The nanotube ring (F) expansion is handled separately — it is always
-- canonical for the (5,0) nanotube family from C30.
generateChildren :: GenConfig -> DualGraph -> [Automorphism]
                 -> [(DualGraph, [Automorphism])]
generateChildren cfg g auts = regular ++ ring
  where
    maxDV = gcMaxDV cfg
    nv    = numVertices g

    -- Reductions for bounding (only need length ≤ 2)
    allReds = allReductionsUpTo 2 g
    maxLen  = maxExpansionLength maxDV (gcMslTable cfg) g allReds

    -- Expansion sites: enumerate, then filter by Rule 2
    expsR2 = filterByRule2 g auts (expansions maxLen g)

    -- Apply each expansion, keep those that pass the canonical test
    regular =
        [ (child, canonAuts (canonicalBFSAndGroup child))
        | e <- expsR2
        , nv + numNewVertices (expKind e) <= maxDV
        , let (child, _) = applyExpansion e g
        , isCanonical e nv child
        ]

    -- Nanotube ring expansion (always canonical, only for (5,0) nanotubes)
    ring = case findNanotubeRing g of
        Just (ringVerts, outer) | nv + 5 <= maxDV ->
            let child = applyRing g ringVerts outer
            in [(child, canonAuts (canonicalBFSAndGroup child))]
        _ -> []

---------------------------------------------------------------------
-- Lazy forest
---------------------------------------------------------------------

-- | A node in the generation forest: a graph with its automorphism group
-- and a lazily-computed list of children.
data GenTree = GenTree
    { gnGraph    :: !DualGraph
    , gnAuts     :: ![Automorphism]
    , gnChildren :: [GenTree]       -- lazy
    }

-- | Build the complete generation forest from the three seeds.
-- The tree is lazy: each node's children are computed only when demanded.
genForest :: GenConfig -> [GenTree]
genForest cfg = map mkTree seedsWithAuts
  where
    mkTree (g, auts) = GenTree g auts
        [mkTree (c, ca) | (c, ca) <- generateChildren cfg g auts]

---------------------------------------------------------------------
-- Traversal
---------------------------------------------------------------------

-- | Count isomers by dual vertex count via sequential DFS.
dfsCountBySize :: GenConfig -> Map Int Int
dfsCountBySize cfg = foldl' countTree Map.empty (genForest cfg)
  where
    maxDV = gcMaxDV cfg
    countTree acc (GenTree g _ children)
        | numVertices g > maxDV = acc
        | otherwise =
            let acc' = Map.insertWith (+) (numVertices g) 1 acc
            in foldl' countTree acc' children

-- | Count isomers using parallel evaluation with a flat fork at the given depth.
--
-- Forces the tree sequentially to the fork depth, collecting all subtree roots
-- into a single flat list. Then sparks each subtree's DFS count in parallel via
-- 'parList rdeepseq'. This avoids recursive sparking where the parent thread
-- races workers — instead, all work units are sparked at once.
--
-- Typical usage:
--   parCountBySize (mkConfig 32) 6  +RTS -N4
parCountBySize :: GenConfig -> Int -> Map Int Int
parCountBySize cfg depth =
    let (topCount, subtrees) = collectAtDepth depth (genForest cfg)
        subtreeCounts = map countSeqTree subtrees `using` parList rdeepseq
    in mergeAll (topCount : subtreeCounts)
  where
    maxDV = gcMaxDV cfg

    mergeAll = foldl' (Map.unionWith (+)) Map.empty

    -- Walk the forest to the fork depth, counting nodes above the fork
    -- and collecting subtree roots at the fork depth.
    collectAtDepth :: Int -> [GenTree] -> (Map Int Int, [GenTree])
    collectAtDepth d trees = foldl' go (Map.empty, []) trees
      where
        go (acc, subs) t@(GenTree g _ children)
            | numVertices g > maxDV = (acc, subs)
            | d <= 0 = (acc, t : subs)
            | otherwise =
                let acc' = Map.insertWith (+) (numVertices g) 1 acc
                    (childAcc, childSubs) = collectAtDepth (d - 1) children
                in (Map.unionWith (+) acc' childAcc, subs ++ childSubs)

    -- Count an entire subtree sequentially, returning a Map.
    countSeqTree :: GenTree -> Map Int Int
    countSeqTree (GenTree g _ children)
        | numVertices g > maxDV = Map.empty
        | otherwise = foldl' countSeq (Map.singleton (numVertices g) 1) children

    -- Sequential DFS accumulator (no sparks, no intermediate maps)
    countSeq :: Map Int Int -> GenTree -> Map Int Int
    countSeq acc (GenTree g _ children)
        | numVertices g > maxDV = acc
        | otherwise =
            let acc' = Map.insertWith (+) (numVertices g) 1 acc
            in foldl' countSeq acc' children

-- | Hybrid parallel: pure forest at the top for parallelism, MutGraph DFS
-- within each subtree for sequential speed.
--
-- Forces the tree to the fork depth using pure 'generateChildren', collecting
-- subtree roots. Each subtree gets its own 'MutGraph' allocation in a separate
-- 'runST' block, enabling independent parallel evaluation with zero shared state.
-- Within each subtree, uses the apply/unsafeFreeze/undo pattern for minimal
-- allocation.
parMutCountBySize :: GenConfig -> Int -> Map Int Int
parMutCountBySize cfg depth =
    let (topCount, subtrees) = collectSubtreeRoots depth (genForest cfg)
        subtreeCounts = map mutCountSubtree subtrees `using` parList rdeepseq
    in mergeAll (topCount : subtreeCounts)
  where
    maxDV = gcMaxDV cfg

    mergeAll = foldl' (Map.unionWith (+)) Map.empty

    -- Walk the forest to the fork depth using pure generation,
    -- collecting (graph, automorphisms) pairs at the fork depth.
    collectSubtreeRoots :: Int -> [GenTree] -> (Map Int Int, [(DualGraph, [Automorphism])])
    collectSubtreeRoots d trees = foldl' go (Map.empty, []) trees
      where
        go (acc, subs) (GenTree g auts children)
            | numVertices g > maxDV = (acc, subs)
            | d <= 0 = (acc, (g, auts) : subs)
            | otherwise =
                let acc' = Map.insertWith (+) (numVertices g) 1 acc
                    (childAcc, childSubs) = collectSubtreeRoots (d - 1) children
                in (Map.unionWith (+) acc' childAcc, subs ++ childSubs)

    -- Count an entire subtree using MutGraph in its own runST block.
    mutCountSubtree :: (DualGraph, [Automorphism]) -> Map Int Int
    mutCountSubtree (g, auts) = runST $ do
        mg <- newMutGraph maxDV
        loadGraph mg g
        mutDFS mg g auts

    -- MutGraph-based DFS: apply expansion → unsafeFreeze for canonical test →
    -- safe freeze only for accepted children → recurse → undo.
    mutDFS mg g auts
        | numVertices g > maxDV = return Map.empty
        | otherwise = do
            let nv = numVertices g
                allReds = allReductionsUpTo 2 g
                maxLen = maxExpansionLength maxDV (gcMslTable cfg) g allReds
                expsR2 = filterByRule2 g auts (expansions maxLen g)

            -- Process regular expansions
            count1 <- foldM (\acc e -> do
                let newVerts = numNewVertices (expKind e)
                if nv + newVerts > maxDV
                    then return acc
                    else do
                        (undo, _) <- applyExpansionM mg e g
                        childUnsafe <- unsafeFreezeGraph mg
                        let !verdict = isCanonicalV e nv childUnsafe
                        acc' <- case verdict of
                            CanonAccept -> do
                                child <- freezeGraph mg
                                let childAuts = canonAuts (canonicalBFSAndGroup child)
                                result <- mutDFS mg child childAuts
                                return (Map.unionWith (+) acc result)
                            _ -> return acc
                        undoMutation mg undo
                        return acc'
                ) (Map.singleton nv 1) expsR2

            -- F expansion for (5,0) nanotubes
            case findNanotubeRing g of
                Just (ring, outer) | nv + 5 <= maxDV -> do
                    undo <- applyRingM mg ring outer
                    child <- freezeGraph mg
                    let childAuts = canonAuts (canonicalBFSAndGroup child)
                    result <- mutDFS mg child childAuts
                    undoMutation mg undo
                    return (Map.unionWith (+) count1 result)
                _ -> return count1

---------------------------------------------------------------------
-- Lazy DFS stream
---------------------------------------------------------------------

-- | Lazy stream of all graphs in DFS pre-order with their automorphism groups.
-- Graphs are produced on demand — useful for pipelines, testing, or taking
-- a prefix without forcing the entire tree.
allGraphs :: GenConfig -> [(DualGraph, [Automorphism])]
allGraphs cfg = concatMap go (genForest cfg)
  where
    maxDV = gcMaxDV cfg
    go (GenTree g auts children)
        | numVertices g > maxDV = []
        | otherwise = (g, auts) : concatMap go children

---------------------------------------------------------------------
-- BFS (breadth-first)
---------------------------------------------------------------------

-- | BFS traversal: returns graphs grouped by tree depth (generation level).
-- Level 0 = seeds, level 1 = children of seeds, etc.
-- Each inner list contains all graphs at that tree depth.
-- The tree is materialized one level at a time — all nodes at depth d
-- are produced before any at depth d+1.
bfsByLevel :: GenConfig -> [[(DualGraph, [Automorphism])]]
bfsByLevel cfg = go (genForest cfg)
  where
    maxDV = gcMaxDV cfg
    go [] = []
    go level =
        let valid = [ (g, auts, children)
                    | GenTree g auts children <- level
                    , numVertices g <= maxDV ]
            graphs = [(g, auts) | (g, auts, _) <- valid]
            nextLevel = concatMap (\(_, _, cs) -> cs) valid
        in graphs : go nextLevel

-- | Count isomers by dual vertex count using BFS traversal.
bfsCountBySize :: GenConfig -> Map Int Int
bfsCountBySize cfg = foldl' countLevel Map.empty (bfsByLevel cfg)
  where
    countLevel acc graphs = foldl' (\m (g, _) ->
        Map.insertWith (+) (numVertices g) 1 m) acc graphs

---------------------------------------------------------------------
-- Generic fold
---------------------------------------------------------------------

-- | Fold over all graphs in the forest with an arbitrary accumulator.
-- Traversal is DFS pre-order. The function receives the accumulator,
-- graph, and automorphism group at each node.
--
-- This generalises all counting schemes:
--   dfsCountBySize = foldForest (\m g _ -> Map.insertWith (+) (numVertices g) 1 m) Map.empty
foldForest :: (a -> DualGraph -> [Automorphism] -> a) -> a -> GenConfig -> a
foldForest f z cfg = foldl' go z (genForest cfg)
  where
    maxDV = gcMaxDV cfg
    go acc (GenTree g auts children)
        | numVertices g > maxDV = acc
        | otherwise = foldl' go (f acc g auts) children

---------------------------------------------------------------------
-- Parallel map
---------------------------------------------------------------------

-- | Apply a function to every graph in the forest, in parallel.
-- Uses flat-fork at the given depth: pure DFS collects results above
-- the fork, parList sparks subtrees below.
parMapForest :: NFData b => GenConfig -> Int -> (DualGraph -> [Automorphism] -> b) -> [b]
parMapForest cfg depth f =
    let (topResults, subtrees) = collectResults depth (genForest cfg)
        subtreeResults = map (concatMap go) subtrees `using` parList rdeepseq
    in topResults ++ concat subtreeResults
  where
    maxDV = gcMaxDV cfg

    go (GenTree g auts children)
        | numVertices g > maxDV = []
        | otherwise = f g auts : concatMap go children

    collectResults d trees = foldl' step ([], []) trees
      where
        step (rs, subs) t@(GenTree g auts children)
            | numVertices g > maxDV = (rs, subs)
            | d <= 0 = (rs, [t] : subs)
            | otherwise =
                let r = f g auts
                    (childRs, childSubs) = collectResults (d - 1) children
                in (rs ++ [r] ++ childRs, subs ++ childSubs)

---------------------------------------------------------------------
-- Target-size generation
---------------------------------------------------------------------

-- | Generate only graphs of exactly the target dual vertex count.
-- Prunes subtrees where all descendants would exceed the target.
-- More efficient than generating everything and filtering when the
-- target is small relative to maxDV.
graphsOfSize :: GenConfig -> Int -> [DualGraph]
graphsOfSize cfg target = concatMap go (genForest cfg)
  where
    go (GenTree g auts children)
        | nv > target = []  -- prune: all descendants are larger
        | nv == target = [g]  -- emit, don't recurse (children are larger)
        | otherwise = concatMap go children  -- nv < target, keep going
      where nv = numVertices g

---------------------------------------------------------------------
-- Work-queue (dynamic load balance)
---------------------------------------------------------------------

-- | Dynamic work-queue parallel generation.
--
-- Separates three concerns:
--   * Generation: pure 'generateChildren' computes subtree roots
--   * Distribution: 'TQueue' — workers pull the next chunk
--   * Consumption: per-worker 'IORef' accumulators, merged at the end
--
-- The tree is split at a fork depth (like flat-fork), but instead of
-- statically assigning subtrees to sparks, the chunks go into a TQueue.
-- Workers pull chunks and DFS each one. When a fast worker finishes,
-- it immediately grabs the next chunk — automatic load balancing.
-- This handles uneven subtree sizes better than parList, where the
-- overall time is dominated by the largest spark.
--
-- Termination: an 'inFlight' counter tracks chunks not yet fully
-- processed. Workers block (STM retry) when queue empty but inFlight > 0.
workQueueCountBySize :: GenConfig -> Int -> IO (Map Int Int)
workQueueCountBySize cfg nWorkers = do
    -- Split the tree at a sensible fork depth
    let (topCount, subtrees) = collectSubtreeRoots' (genForest cfg)
    queue <- newTQueueIO
    inFlight <- newTVarIO (0 :: Int)

    -- Seed the queue with subtree roots
    atomically $ do
        mapM_ (writeTQueue queue) subtrees
        writeTVar inFlight (length subtrees)

    -- Per-worker accumulators
    refs <- replicateM nWorkers (newIORef Map.empty)
    dones <- replicateM nWorkers newEmptyMVar

    -- Launch workers
    forM_ (zip refs dones) $ \(ref, done) ->
        forkIO (wqWorker cfg queue inFlight ref >> putMVar done ())

    -- Wait for all workers to finish
    mapM_ takeMVar dones

    -- Merge per-worker results with the top-level counts
    maps <- mapM readIORef refs
    return $! foldl' (Map.unionWith (+)) topCount maps
  where
    maxDV = gcMaxDV cfg
    forkDepth = max 4 (ceiling (logBase 2 (fromIntegral (nWorkers * 8)) :: Double) :: Int)

    -- Collect (graph, auts) pairs at fork depth, counting nodes above.
    collectSubtreeRoots' :: [GenTree] -> (Map Int Int, [(DualGraph, [Automorphism])])
    collectSubtreeRoots' trees = go forkDepth trees
      where
        go d ts = foldl' step (Map.empty, []) ts
          where
            step (acc, subs) (GenTree g auts children)
                | numVertices g > maxDV = (acc, subs)
                | d <= 0 = (acc, (g, auts) : subs)
                | otherwise =
                    let acc' = Map.insertWith (+) (numVertices g) 1 acc
                        (childAcc, childSubs) = go (d - 1) children
                    in (Map.unionWith (+) acc' childAcc, subs ++ childSubs)

-- | Worker loop: pull a subtree root from the queue, DFS it, repeat.
wqWorker :: GenConfig -> TQueue (DualGraph, [Automorphism])
         -> TVar Int -> IORef (Map Int Int) -> IO ()
wqWorker cfg queue inFlight countRef = loop
  where
    maxDV = gcMaxDV cfg

    loop = do
        mItem <- atomically $ do
            mi <- tryReadTQueue queue
            case mi of
                Just item -> return (Just item)
                Nothing -> do
                    nf <- readTVar inFlight
                    if nf == 0
                        then return Nothing
                        else retry

        case mItem of
            Nothing -> return ()
            Just (g, auts) -> do
                let !count = localDFS g auts
                modifyIORef' countRef (Map.unionWith (+) count)
                atomically $ modifyTVar' inFlight (subtract 1)
                loop

    localDFS :: DualGraph -> [Automorphism] -> Map Int Int
    localDFS g auts
        | numVertices g > maxDV = Map.empty
        | otherwise =
            let children = generateChildren cfg g auts
            in foldl' (\acc (c, ca) -> Map.unionWith (+) acc (localDFS c ca))
                      (Map.singleton (numVertices g) 1) children

---------------------------------------------------------------------
-- Analysis
---------------------------------------------------------------------

-- | Count nodes at each depth in the generation forest.
-- Returns [(depth, nodeCount)] for understanding tree shape and
-- choosing the right parallel fork depth.
depthProfile :: GenConfig -> [(Int, Int)]
depthProfile cfg = Map.toAscList $ foldl' (countAt 0) Map.empty (genForest cfg)
  where
    maxDV = gcMaxDV cfg
    countAt d acc (GenTree g _ children)
        | numVertices g > maxDV = acc
        | otherwise =
            let acc' = Map.insertWith (+) d 1 acc
            in foldl' (countAt (d + 1)) acc' children

-- | Compute the size (node count) of each subtree at the fork depth.
-- Returns sizes sorted descending — reveals work imbalance.
subtreeSizes :: GenConfig -> Int -> [Int]
subtreeSizes cfg depth = sortBy (compare `on` Down) sizes
  where
    (_, subtrees) = collectAtDepth' depth (genForest cfg)
    sizes = map treeSize subtrees
    maxDV = gcMaxDV cfg

    on f g x y = f (g x) (g y)

    collectAtDepth' :: Int -> [GenTree] -> (Int, [GenTree])
    collectAtDepth' d trees = foldl' go (0, []) trees
      where
        go (n, subs) t@(GenTree g _ children)
            | numVertices g > maxDV = (n, subs)
            | d <= 0 = (n, t : subs)
            | otherwise =
                let (cn, csubs) = collectAtDepth' (d - 1) children
                in (n + 1 + cn, subs ++ csubs)

    treeSize :: GenTree -> Int
    treeSize (GenTree g _ children)
        | numVertices g > maxDV = 0
        | otherwise = 1 + sum (map treeSize children)
