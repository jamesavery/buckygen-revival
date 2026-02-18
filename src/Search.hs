{-# OPTIONS_GHC -Wno-orphans #-}
-- | Unified search abstraction for fullerene generation.
--
-- Separates /what/ to collect ('Fold') from /how/ to traverse ('Schedule').
-- The single entry point 'search' dispatches to per-schedule implementations.
--
-- Contains the complete generation forest: configuration, seeds, child
-- generation, lazy tree, and all traversal backends.  This module is
-- self-contained — it does not depend on GenForest.
--
-- This is the @SearchMonad@ abstraction envisioned in @bucky_pseudocode.tex@,
-- made concrete as a fold algebra + schedule selector rather than a type class.

module Search
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
      -- * Collection algebra
    , Fold(..)
    , countBySize, collectAll, collectOfSize, mapNodes
      -- * Scheduling strategy
    , Schedule(..)
      -- * Unified search
    , search, searchPure
      -- * Monadic interface (pure DFS)
    , SearchM, yield, searchForest, runSearchM
      -- * Convenience traversals
    , allGraphs
    , graphsOfSize
    , bfsByLevel
      -- * Analysis
    , depthProfile
    , subtreeSizes
    ) where

import Seeds (DualGraph(..), c20, c28, c30)
import Expansion
    ( Expansion(..), findNanotubeRing, numNewVertices, expansions
    , applyExpansion, applyRing
    )
import Canonical
    ( Automorphism(..), Orientation(..), CanonVerdict(..), isCanonicalV
    , isCanonical
    , allReductionsUpTo, maxExpansionLength, filterByRule2
    , canonicalBFSAndGroup, canonAuts
    , buildMaxStraightLengths
    )
import MutGraph
    ( newMutGraph, loadGraph, freezeGraph, unsafeFreezeGraph
    , applyExpansionM, applyRingM, undoMutation
    )

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

---------------------------------------------------------------------
-- NFData instances for parallel evaluation
---------------------------------------------------------------------

instance NFData DualGraph where
    rnf (DG nv nbrs d5 af df) =
        rnf nv `seq` rnf nbrs `seq` rnf d5 `seq`
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
generateChildren :: GenConfig -> DualGraph -> [Automorphism]
                 -> [(DualGraph, [Automorphism])]
generateChildren cfg g auts = regular ++ ring
  where
    maxDV = gcMaxDV cfg
    nv    = numVertices g

    -- Reductions for bounding (only need length <= 2)
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
-- Fold: the collection algebra
---------------------------------------------------------------------

-- | A fold over the generation forest: what to collect at each node.
--
-- The three fields form a monoid homomorphism:
--   * 'foldVisit' maps each node to a result
--   * 'foldMerge' combines results (must be associative)
--   * 'foldEmpty' is the identity for 'foldMerge'
data Fold r = Fold
    { foldVisit :: !(DualGraph -> [Automorphism] -> r)
    , foldMerge :: !(r -> r -> r)
    , foldEmpty :: !r
    }

-- | Count isomers by dual vertex count: @Map Int Int@.
countBySize :: Fold (Map Int Int)
countBySize = Fold
    { foldVisit = \g _auts -> Map.singleton (numVertices g) 1
    , foldMerge = Map.unionWith (+)
    , foldEmpty = Map.empty
    }

-- | Collect all graphs with their automorphism groups.
collectAll :: Fold [(DualGraph, [Automorphism])]
collectAll = Fold
    { foldVisit = \g auts -> [(g, auts)]
    , foldMerge = (++)
    , foldEmpty = []
    }

-- | Collect only graphs of a specific dual vertex count.
collectOfSize :: Int -> Fold [DualGraph]
collectOfSize target = Fold
    { foldVisit = \g _auts ->
        if numVertices g == target then [g] else []
    , foldMerge = (++)
    , foldEmpty = []
    }

-- | Map a function over every node, collecting results in a list.
mapNodes :: (DualGraph -> [Automorphism] -> b) -> Fold [b]
mapNodes f = Fold
    { foldVisit = \g auts -> [f g auts]
    , foldMerge = (++)
    , foldEmpty = []
    }

---------------------------------------------------------------------
-- Schedule: traversal strategies
---------------------------------------------------------------------

-- | How to traverse the generation forest.
data Schedule
    = SeqDFS                         -- ^ Sequential depth-first
    | SeqBFS                         -- ^ Sequential breadth-first
    | ParFlat  !Int                  -- ^ Parallel flat-fork at depth, pure DFS
    | ParMut   !Int                  -- ^ Parallel flat-fork at depth, MutGraph DFS
    | WorkQ    !Int !Int             -- ^ Work-queue: nWorkers, forkDepth
    | WorkQMut !Int !Int             -- ^ Work-queue + MutGraph: nWorkers, forkDepth
    | HierMut  !Int !Int !Int        -- ^ Hierarchical: nWorkers, d1, d2
    | BfsDfs   !Int !Int !Int        -- ^ BFS+DFS: nBfs, nDfs, forkDepth
    deriving (Eq, Show)

---------------------------------------------------------------------
-- Unified search entry points
---------------------------------------------------------------------

-- | Unified search: run a fold over the generation forest using the
-- given schedule. Returns the accumulated result.
--
-- All schedules produce identical results for a given fold and config;
-- they differ only in traversal order and parallelism.
search :: NFData r => Fold r -> Schedule -> GenConfig -> IO r
search fold sched cfg = case sched of
    SeqDFS           -> return $! searchDFS fold cfg
    SeqBFS           -> return $! searchBFS fold cfg
    ParFlat  d       -> return $! searchParFlat fold cfg d
    ParMut   d       -> return $! searchParMut fold cfg d
    WorkQ    nw d    -> searchWorkQ fold cfg nw d
    WorkQMut nw d    -> searchWorkQMut fold cfg nw d
    HierMut  nw d1 d2 -> searchHierMut fold cfg nw d1 d2
    BfsDfs   nb nd d -> searchBfsDfs fold cfg nb nd d

-- | Pure sequential DFS search. No IO needed.
searchPure :: Fold r -> GenConfig -> r
searchPure = searchDFS

---------------------------------------------------------------------
-- Backend: Sequential DFS
---------------------------------------------------------------------

searchDFS :: Fold r -> GenConfig -> r
searchDFS fold cfg = foldl' (foldMerge fold) (foldEmpty fold)
                            (map go (genForest cfg))
  where
    maxDV = gcMaxDV cfg
    go (GenTree g auts children)
        | numVertices g > maxDV = foldEmpty fold
        | otherwise =
            let here = foldVisit fold g auts
                below = foldl' (foldMerge fold) (foldEmpty fold)
                               (map go children)
            in foldMerge fold here below

---------------------------------------------------------------------
-- Backend: Sequential BFS
---------------------------------------------------------------------

searchBFS :: Fold r -> GenConfig -> r
searchBFS fold cfg = go (foldEmpty fold) (genForest cfg)
  where
    maxDV = gcMaxDV cfg
    go acc [] = acc
    go acc level =
        let valid = [ (g, auts, children)
                    | GenTree g auts children <- level
                    , numVertices g <= maxDV ]
            levelResult = foldl' (\a (g, auts, _) ->
                foldMerge fold a (foldVisit fold g auts)) (foldEmpty fold) valid
            nextLevel = concatMap (\(_, _, cs) -> cs) valid
        in go (foldMerge fold acc levelResult) nextLevel

---------------------------------------------------------------------
-- Backend: Parallel flat-fork (pure DFS per subtree)
---------------------------------------------------------------------

searchParFlat :: NFData r => Fold r -> GenConfig -> Int -> r
searchParFlat fold cfg depth =
    let (topResult, subtrees) = collectWithFold fold cfg depth (genForest cfg)
        subtreeResults = map (foldTree fold cfg) subtrees `using` parList rdeepseq
    in foldl' (foldMerge fold) topResult subtreeResults

-- | Walk the forest to fork depth, folding nodes above the fork and
-- collecting subtree roots at the fork depth.
collectWithFold :: Fold r -> GenConfig -> Int -> [GenTree]
                -> (r, [GenTree])
collectWithFold fold cfg d trees = foldl' go (foldEmpty fold, []) trees
  where
    maxDV = gcMaxDV cfg
    go (acc, subs) t@(GenTree g auts children)
        | numVertices g > maxDV = (acc, subs)
        | d <= 0 = (acc, t : subs)
        | otherwise =
            let acc' = foldMerge fold acc (foldVisit fold g auts)
                (childAcc, childSubs) = collectWithFold fold cfg (d - 1) children
            in (foldMerge fold acc' childAcc, subs ++ childSubs)

-- | Fold an entire subtree (pure DFS).
foldTree :: Fold r -> GenConfig -> GenTree -> r
foldTree fold cfg (GenTree g auts children)
    | numVertices g > maxDV = foldEmpty fold
    | otherwise =
        let here = foldVisit fold g auts
            below = foldl' (foldMerge fold) (foldEmpty fold)
                           (map (foldTree fold cfg) children)
        in foldMerge fold here below
  where maxDV = gcMaxDV cfg

---------------------------------------------------------------------
-- Backend: Parallel flat-fork (MutGraph DFS per subtree)
---------------------------------------------------------------------

searchParMut :: NFData r => Fold r -> GenConfig -> Int -> r
searchParMut fold cfg depth =
    let (topResult, subtrees) =
            collectAboveWithFold fold cfg depth (genForest cfg)
        subtreeResults = map (mutFoldSubtree fold cfg) subtrees
                         `using` parList rdeepseq
    in foldl' (foldMerge fold) topResult subtreeResults

-- | Walk the forest to fork depth, applying the fold to nodes above
-- and collecting (graph, auts) pairs at the fork depth for MutGraph DFS.
collectAboveWithFold :: Fold r -> GenConfig -> Int -> [GenTree]
                     -> (r, [(DualGraph, [Automorphism])])
collectAboveWithFold fold cfg d trees = foldl' go (foldEmpty fold, []) trees
  where
    maxDV = gcMaxDV cfg
    go (acc, subs) (GenTree g auts children)
        | numVertices g > maxDV = (acc, subs)
        | d <= 0 = (acc, (g, auts) : subs)
        | otherwise =
            let acc' = foldMerge fold acc (foldVisit fold g auts)
                (childAcc, childSubs) =
                    collectAboveWithFold fold cfg (d - 1) children
            in (foldMerge fold acc' childAcc, subs ++ childSubs)

-- | MutGraph-based DFS of a subtree, parameterized by a Fold.
--
-- Pre-allocates a single MutGraph, loads the root, then does DFS
-- with apply/unsafeFreeze/undo. The fold's pure foldVisit and foldMerge
-- are called from within ST (they are pure, so this is safe).
mutFoldSubtree :: Fold r -> GenConfig -> (DualGraph, [Automorphism]) -> r
mutFoldSubtree fold cfg (g, auts) = runST $ do
    mg <- newMutGraph maxDV
    loadGraph mg g
    mutDFS mg g auts
  where
    maxDV = gcMaxDV cfg

    mutDFS mg g' auts'
        | numVertices g' > maxDV = return (foldEmpty fold)
        | otherwise = do
            let nv = numVertices g'
                allReds = allReductionsUpTo 2 g'
                maxLen = maxExpansionLength maxDV (gcMslTable cfg) g' allReds
                expsR2 = filterByRule2 g' auts' (expansions maxLen g')
                here = foldVisit fold g' auts'

            -- Process regular expansions
            below <- foldM (\acc e -> do
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
                                return (foldMerge fold acc result)
                            _ -> return acc
                        undoMutation mg undo
                        return acc'
                ) (foldEmpty fold) expsR2

            -- F expansion for (5,0) nanotubes
            result <- case findNanotubeRing g' of
                Just (ring, outer) | nv + 5 <= maxDV -> do
                    undo <- applyRingM mg ring outer
                    child <- freezeGraph mg
                    let childAuts = canonAuts (canonicalBFSAndGroup child)
                    r <- mutDFS mg child childAuts
                    undoMutation mg undo
                    return (foldMerge fold below r)
                _ -> return below

            return (foldMerge fold here result)

---------------------------------------------------------------------
-- Backend: Work-queue (pure DFS)
---------------------------------------------------------------------

searchWorkQ :: NFData r => Fold r -> GenConfig -> Int -> Int -> IO r
searchWorkQ fold cfg nWorkers depth = do
    let (topResult, subtrees) =
            collectAboveWithFold fold cfg depth (genForest cfg)
    runWorkQueueOnFold fold nWorkers (pureFoldSubtree fold cfg) topResult subtrees

-- | Pure DFS of a subtree parameterized by Fold.
pureFoldSubtree :: Fold r -> GenConfig -> (DualGraph, [Automorphism]) -> r
pureFoldSubtree fold cfg (g, auts) = localDFS g auts
  where
    maxDV = gcMaxDV cfg
    localDFS g' auts'
        | numVertices g' > maxDV = foldEmpty fold
        | otherwise =
            let here = foldVisit fold g' auts'
                children = generateChildren cfg g' auts'
                below = foldl' (\acc (c, ca) ->
                    foldMerge fold acc (localDFS c ca)) (foldEmpty fold) children
            in foldMerge fold here below

---------------------------------------------------------------------
-- Backend: Work-queue (MutGraph DFS)
---------------------------------------------------------------------

searchWorkQMut :: NFData r => Fold r -> GenConfig -> Int -> Int -> IO r
searchWorkQMut fold cfg nWorkers depth = do
    let (topResult, subtrees) =
            collectAboveWithFold fold cfg depth (genForest cfg)
    runWorkQueueOnFold fold nWorkers (mutFoldSubtree fold cfg) topResult subtrees

---------------------------------------------------------------------
-- Backend: Hierarchical (parallel split + work-queue MutGraph DFS)
---------------------------------------------------------------------

searchHierMut :: NFData r => Fold r -> GenConfig -> Int -> Int -> Int -> IO r
searchHierMut fold cfg nWorkers d1 d2 = do
    -- Phase 1: Sequential split to d1
    let (topResult, level1) =
            collectAboveWithFold fold cfg d1 (genForest cfg)

    -- Phase 2: Parallel split each level-1 subtree to d2
        splitResults = map (mutSplitFoldSubtree fold cfg d2) level1
                       `using` parList rdeepseq
        splitAbove = map fst splitResults
        allChunks  = concatMap snd splitResults
        aboveResult = foldl' (foldMerge fold) topResult splitAbove

    -- Phase 3: Work-queue distributes chunks to workers
    runWorkQueueOnFold fold nWorkers (mutFoldSubtree fold cfg) aboveResult allChunks

-- | Split a subtree to depth d using MutGraph, folding nodes above the
-- split and collecting chunk roots at the split level.
mutSplitFoldSubtree :: Fold r -> GenConfig -> Int
                    -> (DualGraph, [Automorphism])
                    -> (r, [(DualGraph, [Automorphism])])
mutSplitFoldSubtree fold cfg depth (g, auts) = runST $ do
    mg <- newMutGraph maxDV
    loadGraph mg g
    go mg depth g auts
  where
    maxDV = gcMaxDV cfg

    go mg d g' auts'
        | numVertices g' > maxDV = return (foldEmpty fold, [])
        | d <= 0 = return (foldEmpty fold, [(g', auts')])
        | otherwise = do
            let nv = numVertices g'
                allReds = allReductionsUpTo 2 g'
                maxLen = maxExpansionLength maxDV (gcMslTable cfg) g' allReds
                expsR2 = filterByRule2 g' auts' (expansions maxLen g')
                here = foldVisit fold g' auts'

            -- Process expansions
            (belowResult, roots1) <- foldM (\(accR, accRoots) e -> do
                let newVerts = numNewVertices (expKind e)
                if nv + newVerts > maxDV
                    then return (accR, accRoots)
                    else do
                        (undo, _) <- applyExpansionM mg e g'
                        childUnsafe <- unsafeFreezeGraph mg
                        let !verdict = isCanonicalV e nv childUnsafe
                        result <- case verdict of
                            CanonAccept -> do
                                child <- freezeGraph mg
                                let childAuts = canonAuts (canonicalBFSAndGroup child)
                                go mg (d - 1) child childAuts
                            _ -> return (foldEmpty fold, [])
                        undoMutation mg undo
                        let (cr, croots) = result
                        return (foldMerge fold accR cr, accRoots ++ croots)
                ) (foldEmpty fold, []) expsR2

            -- Nanotube ring
            (finalResult, finalRoots) <- case findNanotubeRing g' of
                Just (ring, outer) | nv + 5 <= maxDV -> do
                    undo <- applyRingM mg ring outer
                    child <- freezeGraph mg
                    let childAuts = canonAuts (canonicalBFSAndGroup child)
                    (cr, croots) <- go mg (d - 1) child childAuts
                    undoMutation mg undo
                    return (foldMerge fold belowResult cr, roots1 ++ croots)
                _ -> return (belowResult, roots1)

            return (foldMerge fold here finalResult, finalRoots)

---------------------------------------------------------------------
-- Backend: Parallel BFS + DFS
---------------------------------------------------------------------

searchBfsDfs :: NFData r => Fold r -> GenConfig -> Int -> Int -> Int -> IO r
searchBfsDfs fold cfg nBfs nDfs forkDepth = do
    -- BFS queue: (graph, auts, current depth)
    bfsQ <- newTQueueIO
    bfsInFlight <- newTVarIO (0 :: Int)

    -- DFS queue: (graph, auts) — leaf chunks
    dfsQ <- newTQueueIO
    dfsInFlight <- newTVarIO (0 :: Int)

    -- Signal: BFS phase is complete
    bfsDone <- newTVarIO False

    -- Seed the BFS queue
    atomically $ do
        forM_ seedsWithAuts $ \(g, auts) ->
            writeTQueue bfsQ (g, auts, 0 :: Int)
        writeTVar bfsInFlight (length seedsWithAuts)

    -- BFS workers: per-worker fold accumulators
    bfsRefs <- replicateM nBfs (newIORef (foldEmpty fold))
    bfsDones <- replicateM nBfs newEmptyMVar
    forM_ (zip bfsRefs bfsDones) $ \(ref, done) ->
        forkIO (bfsFoldWorker fold cfg forkDepth bfsQ bfsInFlight
                              dfsQ dfsInFlight ref >> putMVar done ())

    -- DFS workers: per-worker fold accumulators
    dfsRefs <- replicateM nDfs (newIORef (foldEmpty fold))
    dfsDones <- replicateM nDfs newEmptyMVar
    forM_ (zip dfsRefs dfsDones) $ \(ref, done) ->
        forkIO (dfsFoldWorker fold cfg dfsQ dfsInFlight bfsDone ref
                >> putMVar done ())

    -- Wait for BFS to complete, then signal DFS workers
    mapM_ takeMVar bfsDones
    atomically $ writeTVar bfsDone True

    -- Wait for DFS to complete
    mapM_ takeMVar dfsDones

    -- Merge all results
    bfsResults <- mapM readIORef bfsRefs
    dfsResults <- mapM readIORef dfsRefs
    return $! foldl' (foldMerge fold) (foldEmpty fold) (bfsResults ++ dfsResults)

-- | BFS worker: generate children, fold the current node, route children.
bfsFoldWorker :: Fold r -> GenConfig -> Int
              -> TQueue (DualGraph, [Automorphism], Int)
              -> TVar Int
              -> TQueue (DualGraph, [Automorphism])
              -> TVar Int
              -> IORef r -> IO ()
bfsFoldWorker fold cfg forkDepth bfsQ bfsInFlight dfsQ dfsInFlight countRef = loop
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
                    -- Fold this node
                    let !nodeResult = foldVisit fold g auts
                    modifyIORef' countRef (foldMerge fold nodeResult)
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

-- | DFS worker: pull a chunk and fold it with MutGraph.
dfsFoldWorker :: Fold r -> GenConfig
              -> TQueue (DualGraph, [Automorphism])
              -> TVar Int -> TVar Bool
              -> IORef r -> IO ()
dfsFoldWorker fold cfg dfsQ dfsInFlight bfsDone countRef = loop
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
                        then return Nothing
                        else retry

        case mItem of
            Nothing -> return ()
            Just chunk -> do
                let !result = mutFoldSubtree fold cfg chunk
                modifyIORef' countRef (foldMerge fold result)
                atomically $ modifyTVar' dfsInFlight (subtract 1)
                loop

---------------------------------------------------------------------
-- Work-queue plumbing (generic over Fold)
---------------------------------------------------------------------

-- | Run a work-queue distributing chunks to workers, each processing
-- with the given function. Generic over the fold result type.
runWorkQueueOnFold :: Fold r -> Int
                   -> ((DualGraph, [Automorphism]) -> r)
                   -> r                                    -- above-fork result
                   -> [(DualGraph, [Automorphism])]        -- chunks
                   -> IO r
runWorkQueueOnFold fold nWorkers processChunk topResult chunks = do
    queue <- newTQueueIO
    inFlight <- newTVarIO (0 :: Int)

    -- Seed the queue
    atomically $ do
        mapM_ (writeTQueue queue) chunks
        writeTVar inFlight (length chunks)

    -- Per-worker accumulators
    refs <- replicateM nWorkers (newIORef (foldEmpty fold))
    dones <- replicateM nWorkers newEmptyMVar

    -- Launch workers
    forM_ (zip refs dones) $ \(ref, done) ->
        forkIO (workerLoop processChunk fold queue inFlight ref >> putMVar done ())

    -- Wait for all workers to finish
    mapM_ takeMVar dones

    -- Merge per-worker results with the above-fork result
    maps <- mapM readIORef refs
    return $! foldl' (foldMerge fold) topResult maps

-- | Worker loop: pull a chunk, process it, merge into accumulator, repeat.
workerLoop :: ((DualGraph, [Automorphism]) -> r)
           -> Fold r
           -> TQueue (DualGraph, [Automorphism])
           -> TVar Int -> IORef r -> IO ()
workerLoop processChunk fold queue inFlight countRef = loop
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
                let !result = processChunk chunk
                modifyIORef' countRef (foldMerge fold result)
                atomically $ modifyTVar' inFlight (subtract 1)
                loop

---------------------------------------------------------------------
-- SearchM: monadic interface for custom pure DFS traversals
---------------------------------------------------------------------

-- | A monad for custom search strategies over the generation forest.
-- Wraps a fold computation, allowing 'yield' to emit results at each
-- node. Runs as pure DFS.
newtype SearchM r a = SearchM { unSearchM :: Fold r -> GenConfig -> (a, r) }

instance Functor (SearchM r) where
    fmap f (SearchM g) = SearchM $ \fold cfg ->
        let (a, r) = g fold cfg in (f a, r)

instance Applicative (SearchM r) where
    pure a = SearchM $ \_fold _cfg -> (a, foldEmpty _fold)
    SearchM mf <*> SearchM ma = SearchM $ \fold cfg ->
        let (f, r1) = mf fold cfg
            (a, r2) = ma fold cfg
        in (f a, foldMerge fold r1 r2)

instance Monad (SearchM r) where
    SearchM ma >>= f = SearchM $ \fold cfg ->
        let (a, r1) = ma fold cfg
            SearchM mb = f a
            (b, r2) = mb fold cfg
        in (b, foldMerge fold r1 r2)

-- | Emit a graph/automorphism pair as a result.
yield :: DualGraph -> [Automorphism] -> SearchM r ()
yield g auts = SearchM $ \fold _cfg -> ((), foldVisit fold g auts)

-- | Run a DFS over the entire generation forest, calling the given
-- action at each node. The action receives the graph and automorphisms
-- and can call 'yield' to emit results.
searchForest :: (DualGraph -> [Automorphism] -> SearchM r ())
             -> SearchM r ()
searchForest action = SearchM $ \fold cfg ->
    let maxDV = gcMaxDV cfg
        go (GenTree g auts children)
            | numVertices g > maxDV = foldEmpty fold
            | otherwise =
                let ((), here) = unSearchM (action g auts) fold cfg
                    below = foldl' (foldMerge fold) (foldEmpty fold)
                                   (map go children)
                in foldMerge fold here below
        result = foldl' (foldMerge fold) (foldEmpty fold)
                        (map go (genForest cfg))
    in ((), result)

-- | Run a SearchM computation, returning the accumulated result.
runSearchM :: Fold r -> GenConfig -> SearchM r a -> r
runSearchM fold cfg (SearchM f) = snd (f fold cfg)

---------------------------------------------------------------------
-- Convenience traversals
---------------------------------------------------------------------

-- | Lazy stream of all graphs in DFS pre-order with their automorphism groups.
allGraphs :: GenConfig -> [(DualGraph, [Automorphism])]
allGraphs cfg = concatMap go (genForest cfg)
  where
    maxDV = gcMaxDV cfg
    go (GenTree g auts children)
        | numVertices g > maxDV = []
        | otherwise = (g, auts) : concatMap go children

-- | Generate only graphs of exactly the target dual vertex count.
graphsOfSize :: GenConfig -> Int -> [DualGraph]
graphsOfSize cfg target = concatMap go (genForest cfg)
  where
    go (GenTree g _ children)
        | nv > target = []
        | nv == target = [g]
        | otherwise = concatMap go children
      where nv = numVertices g

-- | BFS traversal: returns graphs grouped by tree depth (generation level).
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

---------------------------------------------------------------------
-- Analysis
---------------------------------------------------------------------

-- | Count nodes at each depth in the generation forest.
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
    collectAtDepth' d trees = foldl' step (0, []) trees
      where
        step (n, subs) t@(GenTree g _ children)
            | numVertices g > maxDV = (n, subs)
            | d <= 0 = (n, t : subs)
            | otherwise =
                let (cn, csubs) = collectAtDepth' (d - 1) children
                in (n + 1 + cn, subs ++ csubs)

    treeSize :: GenTree -> Int
    treeSize (GenTree g _ children)
        | numVertices g > maxDV = 0
        | otherwise = 1 + sum (map treeSize children)
