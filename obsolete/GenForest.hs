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
      -- ** Subtree collection (for Search module)
    , collectSubtreeRoots
      -- ** Work-queue (dynamic load balance)
    , workQueueCountBySize
    , wqMutCountBySize
      -- ** Hierarchical parallel
    , hierarchicalMutCount
      -- ** Parallel BFS + DFS
    , parBfsDfsCount
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
import Control.DeepSeq (NFData(..), rwhnf)
import Control.Concurrent (forkIO, MVar, newEmptyMVar, putMVar, takeMVar)
import Control.Concurrent.STM
import Data.IORef
import qualified Data.IntMap.Strict as IM

-- NFData instances for parallel evaluation of split results.
-- All fields are strict/unboxed, so WHNF suffices for the arrays.
instance NFData DualGraph where
    rnf (DG nv nbrs d5 el af df) =
        rnf nv `seq` rnf nbrs `seq` rnf d5 `seq` rnf el `seq`
        rwhnf af `seq` rwhnf df

instance NFData Automorphism where
    rnf (Aut p d) = rwhnf p `seq` rnf d

instance NFData Orientation where
    rnf Preserving = ()
    rnf Reversing  = ()

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
    let (topCount, subtrees) = collectSubtreeRoots cfg depth (genForest cfg)
        subtreeCounts = map (mutCountSubtree cfg) subtrees `using` parList rdeepseq
    in foldl' (Map.unionWith (+)) Map.empty (topCount : subtreeCounts)

-- | Walk the forest to a fork depth, counting nodes above and collecting
-- (graph, automorphisms) pairs at the fork depth.
collectSubtreeRoots :: GenConfig -> Int -> [GenTree]
                    -> (Map Int Int, [(DualGraph, [Automorphism])])
collectSubtreeRoots cfg d trees = foldl' go (Map.empty, []) trees
  where
    maxDV = gcMaxDV cfg
    go (acc, subs) (GenTree g auts children)
        | numVertices g > maxDV = (acc, subs)
        | d <= 0 = (acc, (g, auts) : subs)
        | otherwise =
            let acc' = Map.insertWith (+) (numVertices g) 1 acc
                (childAcc, childSubs) = collectSubtreeRoots cfg (d - 1) children
            in (Map.unionWith (+) acc' childAcc, subs ++ childSubs)

-- | Count an entire subtree using MutGraph DFS in its own 'runST' block.
-- Pre-allocates a single 'MutGraph', loads the root graph, then does DFS
-- with apply/unsafeFreeze/undo. Used by both flat-fork and work-queue
-- parallel schemes.
mutCountSubtree :: GenConfig -> (DualGraph, [Automorphism]) -> Map Int Int
mutCountSubtree cfg (g, auts) = runST $ do
    mg <- newMutGraph maxDV
    loadGraph mg g
    mutDFS mg g auts
  where
    maxDV = gcMaxDV cfg

    -- MutGraph-based DFS: apply expansion → unsafeFreeze for canonical test →
    -- safe freeze only for accepted children → recurse → undo.
    mutDFS mg g' auts'
        | numVertices g' > maxDV = return Map.empty
        | otherwise = do
            let nv = numVertices g'
                allReds = allReductionsUpTo 2 g'
                maxLen = maxExpansionLength maxDV (gcMslTable cfg) g' allReds
                expsR2 = filterByRule2 g' auts' (expansions maxLen g')

            -- Process regular expansions
            count1 <- foldM (\acc e -> do
                let newVerts = numNewVertices (expKind e)
                if nv + newVerts > maxDV
                    then return acc
                    else do
                        (undo, _) <- applyExpansionM mg e g'
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
            case findNanotubeRing g' of
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

-- | Run a work-queue over pre-built chunks.
-- Shared plumbing for all work-queue variants.
runWorkQueueOn :: Int
               -> ((DualGraph, [Automorphism]) -> Map Int Int)
               -> Map Int Int                              -- counts above the chunks
               -> [(DualGraph, [Automorphism])]            -- chunks to distribute
               -> IO (Map Int Int)
runWorkQueueOn nWorkers countChunk topCount chunks = do
    queue <- newTQueueIO
    inFlight <- newTVarIO (0 :: Int)

    -- Seed the queue
    atomically $ do
        mapM_ (writeTQueue queue) chunks
        writeTVar inFlight (length chunks)

    -- Per-worker accumulators
    refs <- replicateM nWorkers (newIORef Map.empty)
    dones <- replicateM nWorkers newEmptyMVar

    -- Launch workers
    forM_ (zip refs dones) $ \(ref, done) ->
        forkIO (wqWorkerLoop countChunk queue inFlight ref >> putMVar done ())

    -- Wait for all workers to finish
    mapM_ takeMVar dones

    -- Merge per-worker results with the above-fork counts
    maps <- mapM readIORef refs
    return $! foldl' (Map.unionWith (+)) topCount maps

-- | Worker loop: pull a subtree root from the queue, count it, repeat.
wqWorkerLoop :: ((DualGraph, [Automorphism]) -> Map Int Int)
             -> TQueue (DualGraph, [Automorphism])
             -> TVar Int -> IORef (Map Int Int) -> IO ()
wqWorkerLoop countChunk queue inFlight countRef = loop
  where
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
            Just chunk -> do
                let !count = countChunk chunk
                modifyIORef' countRef (Map.unionWith (+) count)
                atomically $ modifyTVar' inFlight (subtract 1)
                loop

-- | Dynamic work-queue with pure DFS backend.
--
-- Pre-splits the tree at a fork depth, puts chunks into a 'TQueue'.
-- Workers pull chunks and DFS each one using pure 'generateChildren'.
-- A fast worker that finishes a small subtree immediately grabs the
-- next chunk — automatic load balancing.
workQueueCountBySize :: GenConfig -> Int -> Int -> IO (Map Int Int)
workQueueCountBySize cfg nWorkers depth =
    let (topCount, chunks) = collectSubtreeRoots cfg depth (genForest cfg)
    in runWorkQueueOn nWorkers pureCount topCount chunks
  where
    maxDV = gcMaxDV cfg

    pureCount :: (DualGraph, [Automorphism]) -> Map Int Int
    pureCount (g, auts) = localDFS g auts

    localDFS :: DualGraph -> [Automorphism] -> Map Int Int
    localDFS g auts
        | numVertices g > maxDV = Map.empty
        | otherwise =
            let children = generateChildren cfg g auts
            in foldl' (\acc (c, ca) -> Map.unionWith (+) acc (localDFS c ca))
                      (Map.singleton (numVertices g) 1) children

-- | Dynamic work-queue with MutGraph DFS backend.
--
-- Same work distribution as 'workQueueCountBySize', but each chunk is
-- processed using a fresh 'MutGraph' with in-place surgery and
-- unsafeFreeze. Combines the work-queue's dynamic load balancing with
-- MutGraph's low-allocation inner loop.
wqMutCountBySize :: GenConfig -> Int -> Int -> IO (Map Int Int)
wqMutCountBySize cfg nWorkers depth =
    let (topCount, chunks) = collectSubtreeRoots cfg depth (genForest cfg)
    in runWorkQueueOn nWorkers (mutCountSubtree cfg) topCount chunks

---------------------------------------------------------------------
-- Hierarchical parallel (parallel work generation)
---------------------------------------------------------------------

-- | Hierarchical parallel generation with MutGraph.
--
-- Three phases:
--   1. Sequential split to depth d1 (tiny: ~44 nodes at d1=4).
--   2. Parallel split: each level-1 subtree is expanded to depth d2
--      using MutGraph, producing fine-grained chunks. All level-1
--      subtrees are split simultaneously via parList.
--   3. Work-queue distributes all chunks to workers for MutGraph DFS.
--
-- This avoids the sequential bottleneck of single-level splitting:
-- the expensive generation work between d1 and d1+d2 is parallelized.
-- At C400+ the nodes above the fork depth number in the millions;
-- a single-level split would serialize all that work.
--
-- The two depth parameters control the granularity:
--   d1: sequential bootstrap depth (keep small, e.g. 3-5)
--   d2: parallel split depth per level-1 subtree
--   effective total depth ≈ d1 + d2
hierarchicalMutCount :: GenConfig -> Int -> Int -> Int -> IO (Map Int Int)
hierarchicalMutCount cfg nWorkers d1 d2 = do
    -- Phase 1: Sequential split to d1 (fast, few nodes)
    let (topCount, level1) = collectSubtreeRoots cfg d1 (genForest cfg)

    -- Phase 2: Parallel split each level-1 subtree to d2
        splitResults = map (mutSplitSubtree cfg d2) level1
                       `using` parList rdeepseq
        splitCounts = map fst splitResults
        allChunks   = concatMap snd splitResults
        aboveCount  = foldl' (Map.unionWith (+)) topCount splitCounts

    -- Phase 3: Work-queue distributes chunks to workers
    runWorkQueueOn nWorkers (mutCountSubtree cfg) aboveCount allChunks

-- | Split a subtree to depth d using MutGraph, returning:
--   (counts of nodes above the split, chunk roots at the split level)
--
-- Uses the same apply/freeze/undo pattern as 'mutCountSubtree', but
-- instead of counting leaf nodes, collects them as work chunks.
-- Each call gets its own MutGraph in a fresh 'runST' block.
mutSplitSubtree :: GenConfig -> Int -> (DualGraph, [Automorphism])
                -> (Map Int Int, [(DualGraph, [Automorphism])])
mutSplitSubtree cfg depth (g, auts) = runST $ do
    mg <- newMutGraph maxDV
    loadGraph mg g
    go mg depth g auts
  where
    maxDV = gcMaxDV cfg

    go mg d g' auts'
        | numVertices g' > maxDV = return (Map.empty, [])
        | d <= 0 = return (Map.empty, [(g', auts')])
        | otherwise = do
            let nv = numVertices g'
                allReds = allReductionsUpTo 2 g'
                maxLen = maxExpansionLength maxDV (gcMslTable cfg) g' allReds
                expsR2 = filterByRule2 g' auts' (expansions maxLen g')

            -- Process expansions: apply → freeze → recurse → undo
            (count1, roots1) <- foldM (\(accC, accR) e -> do
                let newVerts = numNewVertices (expKind e)
                if nv + newVerts > maxDV
                    then return (accC, accR)
                    else do
                        (undo, _) <- applyExpansionM mg e g'
                        childUnsafe <- unsafeFreezeGraph mg
                        let !verdict = isCanonicalV e nv childUnsafe
                        result <- case verdict of
                            CanonAccept -> do
                                child <- freezeGraph mg
                                let childAuts = canonAuts (canonicalBFSAndGroup child)
                                go mg (d - 1) child childAuts
                            _ -> return (Map.empty, [])
                        undoMutation mg undo
                        let (cc, cr) = result
                        return (Map.unionWith (+) accC cc, accR ++ cr)
                ) (Map.singleton nv 1, []) expsR2

            -- Nanotube ring
            case findNanotubeRing g' of
                Just (ring, outer) | nv + 5 <= maxDV -> do
                    undo <- applyRingM mg ring outer
                    child <- freezeGraph mg
                    let childAuts = canonAuts (canonicalBFSAndGroup child)
                    (cc, cr) <- go mg (d - 1) child childAuts
                    undoMutation mg undo
                    return (Map.unionWith (+) count1 cc, roots1 ++ cr)
                _ -> return (count1, roots1)

---------------------------------------------------------------------
-- Parallel BFS + DFS (self-generating work queue)
---------------------------------------------------------------------

-- | Parallel BFS generates work for a parallel DFS.
--
-- Two queues, two worker pools:
--   * BFS queue: items are (graph, auts, depth). BFS workers call
--     'generateChildren', count the parent, and either re-enqueue
--     children (depth < forkDepth) or push them onto the DFS queue
--     (depth >= forkDepth).
--   * DFS queue: items are (graph, auts). DFS workers pull a subtree
--     root and count it via 'mutCountSubtree'.
--
-- The BFS pool starts with 3 items (seeds) and fills itself
-- exponentially — no sequential bottleneck. Once all BFS items reach
-- the fork depth, the BFS queue drains and BFS workers exit. Their
-- threads are freed for DFS workers.
--
-- nBfs: number of BFS worker threads (small, e.g. nWorkers)
-- nDfs: number of DFS worker threads (can overlap with BFS)
-- forkDepth: depth at which BFS items become DFS chunks
parBfsDfsCount :: GenConfig -> Int -> Int -> Int -> IO (Map Int Int)
parBfsDfsCount cfg nBfs nDfs forkDepth = do
    -- BFS queue: (graph, auts, current depth)
    bfsQ <- newTQueueIO
    bfsInFlight <- newTVarIO (0 :: Int)

    -- DFS queue: (graph, auts) — leaf chunks
    dfsQ <- newTQueueIO
    dfsInFlight <- newTVarIO (0 :: Int)

    -- Signal: BFS phase is complete (all BFS workers exited)
    bfsDone <- newTVarIO False

    -- Seed the BFS queue
    atomically $ do
        forM_ seedsWithAuts $ \(g, auts) ->
            writeTQueue bfsQ (g, auts, 0 :: Int)
        writeTVar bfsInFlight (length seedsWithAuts)

    -- BFS workers: per-worker count accumulators
    bfsRefs <- replicateM nBfs (newIORef Map.empty)
    bfsDones <- replicateM nBfs newEmptyMVar
    forM_ (zip bfsRefs bfsDones) $ \(ref, done) ->
        forkIO (bfsWorker cfg forkDepth bfsQ bfsInFlight dfsQ dfsInFlight ref
                >> putMVar done ())

    -- DFS workers: per-worker count accumulators
    dfsRefs <- replicateM nDfs (newIORef Map.empty)
    dfsDones <- replicateM nDfs newEmptyMVar
    forM_ (zip dfsRefs dfsDones) $ \(ref, done) ->
        forkIO (dfsWorker cfg dfsQ dfsInFlight bfsDone ref
                >> putMVar done ())

    -- Wait for BFS to complete, then signal DFS workers
    mapM_ takeMVar bfsDones
    atomically $ writeTVar bfsDone True

    -- Wait for DFS to complete
    mapM_ takeMVar dfsDones

    -- Merge all results
    bfsMaps <- mapM readIORef bfsRefs
    dfsMaps <- mapM readIORef dfsRefs
    return $! foldl' (Map.unionWith (+)) Map.empty (bfsMaps ++ dfsMaps)

-- | BFS worker: generate children, route to BFS or DFS queue.
bfsWorker :: GenConfig -> Int
          -> TQueue (DualGraph, [Automorphism], Int)   -- BFS queue
          -> TVar Int                                   -- BFS inFlight
          -> TQueue (DualGraph, [Automorphism])         -- DFS queue
          -> TVar Int                                   -- DFS inFlight
          -> IORef (Map Int Int)                        -- per-worker counts
          -> IO ()
bfsWorker cfg forkDepth bfsQ bfsInFlight dfsQ dfsInFlight countRef = loop
  where
    maxDV = gcMaxDV cfg

    loop = do
        mItem <- atomically $ do
            mi <- tryReadTQueue bfsQ
            case mi of
                Just item -> return (Just item)
                Nothing -> do
                    nf <- readTVar bfsInFlight
                    if nf == 0
                        then return Nothing
                        else retry

        case mItem of
            Nothing -> return ()
            Just (g, auts, d)
                | numVertices g > maxDV -> do
                    atomically $ modifyTVar' bfsInFlight (subtract 1)
                    loop
                | otherwise -> do
                    -- Count this node
                    modifyIORef' countRef (Map.insertWith (+) (numVertices g) 1)
                    -- Generate children
                    let !children = generateChildren cfg g auts
                    -- Route children
                    atomically $ do
                        forM_ children $ \(c, ca) ->
                            if d + 1 < forkDepth
                                then do
                                    writeTQueue bfsQ (c, ca, d + 1)
                                    modifyTVar' bfsInFlight (+ 1)
                                else do
                                    writeTQueue dfsQ (c, ca)
                                    modifyTVar' dfsInFlight (+ 1)
                        modifyTVar' bfsInFlight (subtract 1)
                    loop

-- | DFS worker: pull a chunk and count it with MutGraph.
-- Waits for either: a DFS item, or BFS completion + empty DFS queue.
dfsWorker :: GenConfig
          -> TQueue (DualGraph, [Automorphism])
          -> TVar Int -> TVar Bool
          -> IORef (Map Int Int) -> IO ()
dfsWorker cfg dfsQ dfsInFlight bfsDone countRef = loop
  where
    loop = do
        mItem <- atomically $ do
            mi <- tryReadTQueue dfsQ
            case mi of
                Just item -> return (Just item)
                Nothing -> do
                    nf <- readTVar dfsInFlight
                    done <- readTVar bfsDone
                    if nf == 0 && done
                        then return Nothing    -- all work complete
                        else retry             -- wait for BFS to produce more

        case mItem of
            Nothing -> return ()
            Just chunk -> do
                let !count = mutCountSubtree cfg chunk
                modifyIORef' countRef (Map.unionWith (+) count)
                atomically $ modifyTVar' dfsInFlight (subtract 1)
                loop

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
