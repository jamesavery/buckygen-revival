{-# LANGUAGE StrictData #-}
-- | Mutable graph representation for the generation hot path.
--
-- Pre-allocated Nf×6 flat array layout in the ST monad.
-- Expansions mutate in-place, freeze for the pure canonical test,
-- then undo to restore parent state.

module MutGraph
    ( MutGraph
    , UndoInfo
      -- * Lifecycle
    , newMutGraph
    , loadGraph
    , freezeGraph
    , unsafeFreezeGraph
      -- * Mutable expansion surgery
    , applyExpansionM
    , applyRingM
      -- * Undo
    , undoMutation
    ) where

import Seeds (DualGraph(..))
import Expansion

import Control.Monad (forM_, forM, foldM)
import Control.Monad.ST (ST)
import Data.Array.MArray (newArray, readArray, writeArray)
import Data.Array.ST (STUArray, freeze)
import Data.Array.Unsafe (unsafeFreeze)
import Data.Array.Unboxed (UArray)
import qualified Data.Array.Unboxed as UA
import qualified Data.IntMap.Strict as IM
import Data.List (sort)
import Data.STRef (STRef, newSTRef, readSTRef, writeSTRef)

---------------------------------------------------------------------
-- Types
---------------------------------------------------------------------

data MutGraph s = MutGraph
    { mgAdj   :: !(STUArray s Int Int)  -- maxNf × 6 flat
    , mgDeg   :: !(STUArray s Int Int)  -- maxNf
    , mgNv    :: !(STRef s Int)
    , mgDeg5  :: !(STRef s [Int])
    , mgMaxNf :: !Int
    }

data UndoInfo = UndoInfo
    { uiSaved    :: ![(Int, [Int], Int)]  -- [(vertex, slots[0..5], degree)]
    , uiNewVerts :: !Int
    , uiOldNv    :: !Int
    , uiOldDeg5  :: ![Int]
    }

---------------------------------------------------------------------
-- Lifecycle
---------------------------------------------------------------------

-- | Allocate a mutable graph with capacity for maxNf vertices.
newMutGraph :: Int -> ST s (MutGraph s)
newMutGraph maxNf = do
    adj <- newArray (0, maxNf * 6 - 1) (-1)
    deg <- newArray (0, maxNf - 1) 0
    nv  <- newSTRef 0
    d5  <- newSTRef []
    return (MutGraph adj deg nv d5 maxNf)

-- | Load an immutable DualGraph into the mutable graph.
loadGraph :: MutGraph s -> DualGraph -> ST s ()
loadGraph mg g = do
    let nv = numVertices g
        af = adjFlat g
        df = degFlat g
    -- Copy adjacency and degree arrays
    forM_ [0 .. nv * 6 - 1] $ \i ->
        writeArray (mgAdj mg) i (af UA.! i)
    forM_ [0 .. nv - 1] $ \i ->
        writeArray (mgDeg mg) i (df UA.! i)
    writeSTRef (mgNv mg) nv
    writeSTRef (mgDeg5 mg) (degree5 g)

-- | Freeze the mutable graph into an immutable DualGraph.
-- Uses safe freeze (copies arrays) to prevent aliasing issues when
-- the frozen graph is read after the MutGraph has been mutated.
-- This is necessary because canonicalBFSAndGroup may be evaluated lazily
-- after undoMutation has restored the MutGraph to parent state.
freezeGraph :: MutGraph s -> ST s DualGraph
freezeGraph mg = do
    nv <- readSTRef (mgNv mg)
    d5 <- readSTRef (mgDeg5 mg)
    af <- freeze (mgAdj mg)
    df <- freeze (mgDeg mg)
    return (DG nv IM.empty d5 af df)
{-# INLINE freezeGraph #-}

-- | Zero-copy freeze: aliases the mutable arrays as immutable UArrays.
-- The returned DualGraph shares memory with the MutGraph — any subsequent
-- mutation (undoMutation, applyExpansionM) will corrupt the DualGraph.
--
-- SAFETY: The caller MUST fully evaluate all reads from the returned
-- DualGraph BEFORE any mutation of the MutGraph. Pattern-matching on
-- a strict CanonVerdict from isCanonicalV satisfies this requirement.
--
-- Used in the generation hot loop: 97% of children fail isCanonical
-- and only need this zero-copy view. The 3% that pass use freezeGraph.
unsafeFreezeGraph :: MutGraph s -> ST s DualGraph
unsafeFreezeGraph mg = do
    nv <- readSTRef (mgNv mg)
    d5 <- readSTRef (mgDeg5 mg)
    af <- unsafeFreeze (mgAdj mg)
    df <- unsafeFreeze (mgDeg mg)
    return (DG nv IM.empty d5 af df)
{-# INLINE unsafeFreezeGraph #-}

---------------------------------------------------------------------
-- Mutable navigation
---------------------------------------------------------------------

readNbrM :: MutGraph s -> Int -> Int -> ST s Int
readNbrM mg u i = readArray (mgAdj mg) (u * 6 + i)
{-# INLINE readNbrM #-}

readDegM :: MutGraph s -> Int -> ST s Int
readDegM mg u = readArray (mgDeg mg) u
{-# INLINE readDegM #-}

-- | Find index of v in u's CW neighbor list (linear scan over 5-6 entries).
indexOfM :: MutGraph s -> Int -> Int -> ST s Int
indexOfM mg u v = do
    d <- readDegM mg u
    go 0 d
  where
    base = u * 6
    go i d
        | i >= d    = error $ "indexOfM: " ++ show v ++ " not in " ++ show u ++ "'s neighbors"
        | otherwise = do
            n <- readArray (mgAdj mg) (base + i)
            if n == v then return i else go (i + 1) d

-- | Previous neighbor (counter-clockwise) from v in u's neighbor list.
prevCWM :: MutGraph s -> Int -> Int -> ST s Int
prevCWM mg u v = do
    d <- readDegM mg u
    i <- indexOfM mg u v
    readNbrM mg u ((i - 1 + d) `mod` d)

---------------------------------------------------------------------
-- Mutable surgery primitives
---------------------------------------------------------------------

-- | Save vertex state for undo (6 neighbor slots + degree).
saveVertex :: MutGraph s -> Int -> ST s (Int, [Int], Int)
saveVertex mg v = do
    d <- readDegM mg v
    slots <- forM [0..5] $ \i -> readArray (mgAdj mg) (v * 6 + i)
    return (v, slots, d)

-- | Set all neighbors for a vertex.
setVertexNbrs :: MutGraph s -> Int -> [Int] -> ST s ()
setVertexNbrs mg u ns = do
    let d = length ns
    forM_ (zip [0..] ns) $ \(i, n) ->
        writeArray (mgAdj mg) (u * 6 + i) n
    -- Pad remaining slots with -1
    forM_ [d..5] $ \i ->
        writeArray (mgAdj mg) (u * 6 + i) (-1)
    writeArray (mgDeg mg) u d

-- | Insert 'new' after 'after' in u's CW list (degree 5→6).
-- Shifts entries right by 1, writes 'new' at the vacated slot.
insertAfterM :: MutGraph s -> Int -> Int -> Int -> ST s ()
insertAfterM mg u after new = do
    d <- readDegM mg u
    pos <- indexOfM mg u after
    -- Shift entries at pos+1..d-1 right by 1 (from end to avoid overwrite)
    forM_ [d-1, d-2 .. pos+1] $ \i -> do
        v <- readArray (mgAdj mg) (u * 6 + i)
        writeArray (mgAdj mg) (u * 6 + i + 1) v
    -- Write new at pos+1
    writeArray (mgAdj mg) (u * 6 + pos + 1) new
    -- Increment degree
    writeArray (mgDeg mg) u (d + 1)

-- | Replace 'old' with 'new' in u's CW list (value-based scan).
replaceNbrM :: MutGraph s -> Int -> Int -> Int -> ST s ()
replaceNbrM mg u old new = do
    d <- readDegM mg u
    go 0 d
  where
    base = u * 6
    go i d
        | i >= d = error $ "replaceNbrM: " ++ show old ++ " not in " ++ show u ++ "'s neighbors"
        | otherwise = do
            v <- readArray (mgAdj mg) (base + i)
            if v == old
                then writeArray (mgAdj mg) (base + i) new
                else go (i + 1) d

---------------------------------------------------------------------
-- Undo
---------------------------------------------------------------------

-- | Restore saved state.
undoMutation :: MutGraph s -> UndoInfo -> ST s ()
undoMutation mg info = do
    -- Restore modified vertices
    forM_ (uiSaved info) $ \(v, slots, d) -> do
        forM_ (zip [0..5] slots) $ \(i, n) ->
            writeArray (mgAdj mg) (v * 6 + i) n
        writeArray (mgDeg mg) v d
    -- Clear new vertices
    let oldNv = uiOldNv info
    forM_ [oldNv .. oldNv + uiNewVerts info - 1] $ \v -> do
        forM_ [0..5] $ \i -> writeArray (mgAdj mg) (v * 6 + i) (-1)
        writeArray (mgDeg mg) v 0
    -- Restore nv and deg5
    writeSTRef (mgNv mg) (uiOldNv info)
    writeSTRef (mgDeg5 mg) (uiOldDeg5 info)

---------------------------------------------------------------------
-- Helpers for building undo info
---------------------------------------------------------------------

beginUndo :: MutGraph s -> ST s (Int, [Int])
beginUndo mg = do
    nv <- readSTRef (mgNv mg)
    d5 <- readSTRef (mgDeg5 mg)
    return (nv, d5)

mkUndo :: [(Int, [Int], Int)] -> Int -> Int -> [Int] -> UndoInfo
mkUndo saved newVerts oldNv oldDeg5 =
    UndoInfo saved newVerts oldNv oldDeg5

---------------------------------------------------------------------
-- Mutable expansion functions
---------------------------------------------------------------------

-- | Apply L0 expansion (adds 2 degree-5 vertices).
applyL0M :: MutGraph s -> PathInfo -> Dir -> ST s UndoInfo
applyL0M mg (PathInfo path par) dir = do
    (nv, oldDeg5) <- beginUndo mg
    let a = nv; b = nv + 1
        p0 = path !! 0; p1 = path !! 1; p2 = path !! 2
        q0 = par  !! 0; q1 = par  !! 1; q2 = par  !! 2

    -- Save modified vertices
    saved <- mapM (saveVertex mg) [p0, p1, p2, q0, q1, q2]

    -- Set new vertex neighbor lists
    let aNbrs = case dir of
            DRight -> [p0, q0, q1, b, p1]
            DLeft  -> [p0, p1, b, q1, q0]
        bNbrs = case dir of
            DRight -> [a, q1, q2, p2, p1]
            DLeft  -> [a, p1, p2, q2, q1]
    setVertexNbrs mg a aNbrs
    setVertexNbrs mg b bNbrs

    -- insertAfter at p0
    case dir of
        DRight -> insertAfterM mg p0 q0 a
        DLeft  -> insertAfterM mg p0 p1 a

    -- insertAfter at q2
    case dir of
        DRight -> insertAfterM mg q2 p2 b
        DLeft  -> do
            prev <- prevCWM mg q2 p2
            insertAfterM mg q2 prev b

    -- replaceNbr at path/par vertices
    replaceNbrM mg p1 q0 a
    replaceNbrM mg p1 q1 b
    replaceNbrM mg p2 q1 b
    replaceNbrM mg q0 p1 a
    replaceNbrM mg q1 p1 a
    replaceNbrM mg q1 p2 b

    -- Update nv and deg5
    let deg5' = sort $ a : b : filter (\v -> v /= p0 && v /= q2) oldDeg5
    writeSTRef (mgNv mg) (nv + 2)
    writeSTRef (mgDeg5 mg) deg5'

    return (mkUndo saved 2 nv oldDeg5)

-- | Apply straight expansion L_i (i >= 1, adds pathlength vertices).
applyStraightM :: MutGraph s -> PathInfo -> Dir -> Int -> ST s UndoInfo
applyStraightM mg (PathInfo path par) dir pathlength = do
    (nv, oldDeg5) <- beginUndo mg

    -- Vertices to save: path[0], par[pathlength] (endpoints, degree 5→6)
    --                   path[1..pathlength], par[0..pathlength-1]
    let endpointVerts = [path !! 0, par !! pathlength]
        pathVerts = [path !! i | i <- [1..pathlength]]
        parVerts  = [par !! i  | i <- [0..pathlength-1]]
    saved <- mapM (saveVertex mg) (endpointVerts ++ pathVerts ++ parVerts)

    -- Set new vertex neighbor lists
    forM_ [0..pathlength-1] $ \i -> do
        let nbrsI
              | i == 0 && i == pathlength - 1 = error "applyStraightM: pathlength must be >= 2"
              | i == 0 = case dir of
                  DRight -> [path!!0, par!!0, par!!1, nv+1, path!!1]
                  DLeft  -> [path!!0, path!!1, nv+1, par!!1, par!!0]
              | i == pathlength - 1 = case dir of
                  DRight -> [nv+i-1, par!!i, par!!(i+1), path!!(i+1), path!!i]
                  DLeft  -> [nv+i-1, path!!i, path!!(i+1), par!!(i+1), par!!i]
              | otherwise = case dir of
                  DRight -> [nv+i-1, par!!i, par!!(i+1), nv+i+1, path!!(i+1), path!!i]
                  DLeft  -> [nv+i-1, path!!i, path!!(i+1), nv+i+1, par!!(i+1), par!!i]
        setVertexNbrs mg (nv + i) nbrsI

    -- insertAfter at endpoints (degree 5→6)
    case dir of
        DRight -> insertAfterM mg (path!!0) (par!!0) nv
        DLeft  -> insertAfterM mg (path!!0) (path!!1) nv
    case dir of
        DRight -> insertAfterM mg (par!!pathlength) (path!!pathlength) (nv+pathlength-1)
        DLeft  -> do
            prev <- prevCWM mg (par!!pathlength) (path!!pathlength)
            insertAfterM mg (par!!pathlength) prev (nv+pathlength-1)

    -- Replace neighbors on path vertices
    forM_ [1..pathlength] $ \i -> do
        replaceNbrM mg (path!!i) (par!!(i-1)) (nv+i-1)
        if i < pathlength
            then replaceNbrM mg (path!!i) (par!!i) (nv+i)
            else return ()

    -- Replace neighbors on parallel path vertices
    forM_ [0..pathlength-1] $ \i -> do
        if i > 0
            then replaceNbrM mg (par!!i) (path!!i) (nv+i-1)
            else return ()
        replaceNbrM mg (par!!i) (path!!(i+1)) (nv+i)

    -- Update nv and deg5
    let deg5' = sort $ nv : (nv+pathlength-1) :
            filter (\v -> v /= path!!0 && v /= par!!pathlength) oldDeg5
    writeSTRef (mgNv mg) (nv + pathlength)
    writeSTRef (mgDeg5 mg) deg5'

    return (mkUndo saved pathlength nv oldDeg5)

-- | Apply B_{0,0} expansion (adds 3 vertices).
applyBentZeroM :: MutGraph s -> PathInfo -> Dir -> ST s UndoInfo
applyBentZeroM mg (PathInfo path par) dir = do
    (nv, oldDeg5) <- beginUndo mg
    let a = nv; b = nv + 1; c = nv + 2
        p0 = path!!0; p1 = path!!1; p2 = path!!2; p3 = path!!3; p4 = path!!4
        q0 = par!!0;  q1 = par!!1;  q2 = par!!2

    -- Save modified vertices
    saved <- mapM (saveVertex mg) [p0, p1, p2, p3, p4, q0, q1, q2]

    -- Set new vertex neighbor lists
    let aNbrs = case dir of
            DRight -> [p0, q0, q1, b, p1]
            DLeft  -> [p0, p1, b, q1, q0]
        bNbrs = case dir of
            DRight -> [a, q1, c, p3, p2, p1]
            DLeft  -> [a, p1, p2, p3, c, q1]
        cNbrs = case dir of
            DRight -> [b, q1, q2, p4, p3]
            DLeft  -> [b, p3, p4, q2, q1]
    setVertexNbrs mg a aNbrs
    setVertexNbrs mg b bNbrs
    setVertexNbrs mg c cNbrs

    -- insertAfter at endpoints
    case dir of
        DRight -> insertAfterM mg p0 q0 a
        DLeft  -> insertAfterM mg p0 p1 a
    case dir of
        DRight -> insertAfterM mg p4 p3 c
        DLeft  -> do
            prev <- prevCWM mg p4 p3
            insertAfterM mg p4 prev c

    -- Replace neighbors on path vertices
    replaceNbrM mg p1 q0 a
    replaceNbrM mg p1 q1 b
    replaceNbrM mg p2 q1 b
    replaceNbrM mg p3 q1 b
    replaceNbrM mg p3 q2 c

    -- Replace neighbors on parallel path vertices
    replaceNbrM mg q0 p1 a
    replaceNbrM mg q1 p1 a
    replaceNbrM mg q1 p2 b
    replaceNbrM mg q1 p3 c
    replaceNbrM mg q2 p3 c

    -- Update nv and deg5
    let deg5' = sort $ a : c : filter (\v -> v /= p0 && v /= p4) oldDeg5
    writeSTRef (mgNv mg) (nv + 3)
    writeSTRef (mgDeg5 mg) deg5'

    return (mkUndo saved 3 nv oldDeg5)

-- | Apply B_{i,j} expansion (i+j > 0, adds bentLen+3 vertices).
applyBentM :: MutGraph s -> PathInfo -> Dir -> Int -> Int -> ST s UndoInfo
applyBentM mg (PathInfo path par) dir bentPos bentLen = do
    (nv, oldDeg5) <- beginUndo mg
    let numNew = bentLen + 3
        bendI = bentPos + 2

    -- Vertices to save: endpoints + all path/par vertices touched
    let endpointVerts = [path !! 0, path !! (bentLen + 4)]
        beforePath = [path !! i | i <- [1..bentPos+1]]
        beforePar  = [par !! i  | i <- [0..bentPos]]
        bendPath   = [path !! bendI]
        bendPar    = [par !! (bendI - 1)]
        afterPath  = [path !! i | i <- [bendI+1..bentLen+3]]
        afterPar   = [par !! i  | i <- [bentPos+2..bentLen+2]]
    saved <- mapM (saveVertex mg) (endpointVerts ++ beforePath ++ beforePar
                                   ++ bendPath ++ bendPar ++ afterPath ++ afterPar)

    -- Set new vertex neighbor lists (following applyBent's newVertexNbrsDirect)
    forM_ [0..numNew-1] $ \i -> do
        let nbrsI
              -- Before the bend (i = 0..bentPos)
              | i <= bentPos && i == 0 = case dir of
                  DRight -> [path!!0, par!!0, par!!1, nv+1, path!!1]
                  DLeft  -> [path!!0, path!!1, nv+1, par!!1, par!!0]
              | i <= bentPos = case dir of
                  DRight -> [nv+i-1, par!!i, par!!(i+1), nv+i+1, path!!(i+1), path!!i]
                  DLeft  -> [nv+i-1, path!!i, path!!(i+1), nv+i+1, par!!(i+1), par!!i]
              -- Bend vertex (i = bentPos+1)
              | i == bentPos + 1 = case dir of
                  DRight -> [nv+i-1, par!!i, nv+i+1, path!!(i+2), path!!(i+1), path!!i]
                  DLeft  -> [nv+i-1, path!!i, path!!(i+1), path!!(i+2), nv+i+1, par!!i]
              -- After bend, last vertex
              | i == numNew - 1 = case dir of
                  DRight -> [nv+i-1, par!!(i-1), par!!i, path!!(i+2), path!!(i+1)]
                  DLeft  -> [nv+i-1, path!!(i+1), path!!(i+2), par!!i, par!!(i-1)]
              -- After bend, middle
              | otherwise = case dir of
                  DRight -> [nv+i-1, par!!(i-1), par!!i, nv+i+1, path!!(i+2), path!!(i+1)]
                  DLeft  -> [nv+i-1, path!!(i+1), path!!(i+2), nv+i+1, par!!i, par!!(i-1)]
        setVertexNbrs mg (nv + i) nbrsI

    -- insertAfter at endpoints
    case dir of
        DRight -> insertAfterM mg (path!!0) (par!!0) nv
        DLeft  -> insertAfterM mg (path!!0) (path!!1) nv
    case dir of
        DRight -> insertAfterM mg (path!!(bentLen+4)) (path!!(bentLen+3)) (nv+numNew-1)
        DLeft  -> do
            prev <- prevCWM mg (path!!(bentLen+4)) (path!!(bentLen+3))
            insertAfterM mg (path!!(bentLen+4)) prev (nv+numNew-1)

    -- Replace neighbors: before the bend
    forM_ [1..bentPos+1] $ \i -> do
        replaceNbrM mg (path!!i) (par!!(i-1)) (nv+i-1)
        replaceNbrM mg (path!!i) (par!!i) (nv+i)

    forM_ [0..bentPos] $ \i -> do
        if i > 0
            then replaceNbrM mg (par!!i) (path!!i) (nv+i-1)
            else return ()
        replaceNbrM mg (par!!i) (path!!(i+1)) (nv+i)

    -- Bend area
    replaceNbrM mg (path!!bendI) (par!!(bendI-1)) (nv+bendI-1)
    replaceNbrM mg (par!!(bendI-1)) (path!!(bendI-1)) (nv+bendI-2)
    replaceNbrM mg (par!!(bendI-1)) (path!!bendI) (nv+bendI-1)
    replaceNbrM mg (par!!(bendI-1)) (path!!(bendI+1)) (nv+bendI)

    -- After the bend: path[i] replace par[i-2] → nv+i-2, par[i-1] → nv+i-1
    forM_ [bendI+1..bentLen+3] $ \i -> do
        replaceNbrM mg (path!!i) (par!!(i-2)) (nv+i-2)
        replaceNbrM mg (path!!i) (par!!(i-1)) (nv+i-1)

    -- par[i] after bend: replace path[i+1] → nv+i, path[i+2] → nv+i+1
    forM_ [bentPos+2..bentLen+2] $ \i -> do
        replaceNbrM mg (par!!i) (path!!(i+1)) (nv+i)
        if i < bentLen + 2
            then replaceNbrM mg (par!!i) (path!!(i+2)) (nv+i+1)
            else return ()

    -- Update nv and deg5
    let deg5' = sort $ nv : (nv+numNew-1) :
            filter (\v -> v /= path!!0 && v /= path!!(bentLen+4)) oldDeg5
    writeSTRef (mgNv mg) (nv + numNew)
    writeSTRef (mgDeg5 mg) deg5'

    return (mkUndo saved numNew nv oldDeg5)

-- | Apply F (nanotube ring) expansion (adds 5 degree-6 vertices).
applyRingM :: MutGraph s -> [Vertex] -> [Vertex] -> ST s UndoInfo
applyRingM mg ring outer = do
    (nv, oldDeg5) <- beginUndo mg

    -- Save ring and outer vertices
    saved <- mapM (saveVertex mg) (ring ++ outer)

    -- Set new vertex neighbor lists
    forM_ [0..4] $ \i -> do
        let nbrsI = [ ring !! i
                     , nv + ((i - 1 + 5) `mod` 5)
                     , outer !! i
                     , outer !! ((i + 1) `mod` 5)
                     , nv + ((i + 1) `mod` 5)
                     , ring !! ((i + 1) `mod` 5)
                     ]
        setVertexNbrs mg (nv + i) nbrsI

    -- Update ring vertices
    forM_ [0..4] $ \i -> do
        let prev = (i - 1 + 5) `mod` 5
        replaceNbrM mg (ring !! i) (outer !! prev) (nv + prev)
        replaceNbrM mg (ring !! i) (outer !! i)    (nv + i)

    -- Update outer (cap) vertices
    forM_ [0..4] $ \i -> do
        let prev = (i - 1 + 5) `mod` 5
        replaceNbrM mg (outer !! i) (ring !! i)              (nv + prev)
        replaceNbrM mg (outer !! i) (ring !! ((i + 1) `mod` 5)) (nv + i)

    -- Update nv (deg5 unchanged — all new verts are degree 6)
    writeSTRef (mgNv mg) (nv + 5)

    return (mkUndo saved 5 nv oldDeg5)

---------------------------------------------------------------------
-- Dispatch
---------------------------------------------------------------------

-- | Apply an expansion to the mutable graph.
-- Computes PathInfo from the frozen parent graph (pure), then dispatches.
applyExpansionM :: MutGraph s -> Expansion -> DualGraph -> ST s (UndoInfo, PathInfo)
applyExpansionM mg (Exp (L 0) edge dir) parentGraph = do
    let pi = computeStraightPath parentGraph edge dir 3
    undo <- applyL0M mg pi dir
    return (undo, pi)
applyExpansionM mg (Exp (L i) edge dir) parentGraph = do
    let pi = computeStraightPath parentGraph edge dir (i + 3)
    undo <- applyStraightM mg pi dir (i + 2)
    return (undo, pi)
applyExpansionM mg (Exp (B 0 0) edge dir) parentGraph = do
    let pi = computeBentZeroPath parentGraph edge dir
    undo <- applyBentZeroM mg pi dir
    return (undo, pi)
applyExpansionM mg (Exp (B i j) edge dir) parentGraph = do
    let pi = computeBentPath parentGraph edge dir i (i + j)
    undo <- applyBentM mg pi dir i (i + j)
    return (undo, pi)
applyExpansionM _ (Exp F _ _) _ =
    error "applyExpansionM: F expansion uses applyRingM directly"
