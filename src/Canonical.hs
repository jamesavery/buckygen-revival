{-# LANGUAGE StrictData #-}
-- | Canonical construction path test for fullerene enumeration.
--
-- Implements McKay's canonical construction path: each fullerene has a unique
-- "canonical reduction" — the reduction with the lexicographically smallest
-- 5-tuple (x0, x1, x2, x3, x4). An expansion is accepted iff its inverse
-- is equivalent to the canonical reduction of the child graph.
--
-- The 5-tuple cascade uses lazy evaluation: x0-x3 are O(1) combinatorial
-- invariants that resolve >99% of cases; x4 (full BFS canonical form) is
-- O(n) but rarely needed.

module Canonical
    ( -- * Types
      BFSCode(..)
    , bfsCodeLength
    , CanOrd
    , Orientation(..)
    , Automorphism(..)
    , CanonResult(..)
      -- * Canonicity test
    , canonOrd
    , allReductions
    , allReductionsUpTo
    , isCanonical
      -- * BFS canonical form
    , bfsCanonicalForm
    , bfsWithNumbering
      -- * Automorphism group
    , canonicalBFSAndGroup
      -- * Rule 2: expansion orbit filtering
    , applyAutToExp
    , computeOtherTriple
    , filterByRule2
      -- * Bounding lemmas
    , buildMaxStraightLengths
    , maxExpansionLength
    , reductionFiveVertices
      -- * Instrumented canonicity test
    , CanonVerdict(..)
    , isCanonicalV
    ) where

import Seeds (DualGraph(..))
import Expansion
import Data.IntSet (IntSet)
import qualified Data.IntSet as IS
import Data.List (nub, foldl')
import Data.Bits ((.|.), shiftL)
import qualified Data.Set as Set
import Data.Array.Unboxed (UArray, (!), bounds, array, elems)
import Data.Array.ST (STUArray, runSTUArray, newArray, readArray, writeArray, freeze)
import Control.Monad (forM_)
import Control.Monad.ST (ST, runST)

---------------------------------------------------------------------
-- Types
---------------------------------------------------------------------

-- | BFS canonical form: a sequence of integers encoding the graph
-- structure from a given starting edge and direction.
-- Stored as an unboxed array for O(1) access and zero list allocation.
newtype BFSCode = BFS (UArray Int Int)

instance Eq BFSCode where
    BFS a == BFS b = a == b

instance Ord BFSCode where
    compare (BFS a) (BFS b) =
        let (_, hi) = bounds a
        in go 0 hi
      where
        go i hi'
            | i > hi'   = EQ
            | otherwise = case compare (a ! i) (b ! i) of
                            EQ -> go (i + 1) hi'
                            c  -> c

instance Show BFSCode where
    show (BFS arr) = "BFS " ++ show (elems arr)

-- | Length of BFS code array.
bfsCodeLength :: BFSCode -> Int
bfsCodeLength (BFS arr) = let (lo, hi) = bounds arr in hi - lo + 1

-- | Canonical ordering tuple. Smaller = more canonical.
-- (reductionLength, negate longestStraight, colourPair, pathColour, BFS code)
type CanOrd = (Int, Int, (Int, Int), Int, BFSCode)

-- | Orientation of an automorphism relative to the canonical orientation.
data Orientation = Preserving | Reversing
    deriving (Eq, Ord, Show)

-- | An automorphism: vertex permutation + orientation.
data Automorphism = Aut
    { autPerm        :: !(UArray Int Int)  -- vertex → vertex (O(1) lookup)
    , autOrientation :: !Orientation
    } deriving (Show)

-- | Result of canonical BFS computation with automorphism group.
data CanonResult = CanonResult
    { canonCode  :: !BFSCode
    , canonAuts  :: ![Automorphism]   -- includes identity
    , canonNbop  :: !Int              -- count of orientation-preserving auts
    , canonEdgesTotal    :: !Int      -- starting edges before colour filter
    , canonEdgesFiltered :: !Int      -- starting edges after colour filter
    , canonBfsGT :: !Int              -- BFS comparisons: GT (rejected)
    , canonBfsEQ :: !Int              -- BFS comparisons: EQ (automorphism)
    , canonBfsLT :: !Int              -- BFS comparisons: LT (new minimum)
    } deriving (Show)

---------------------------------------------------------------------
-- Colour functions
--
-- Direction convention (verified against buckygen.c):
--   DRight (Haskell) = use_prev (C) = CW walk
--   DLeft  (Haskell) = use_next (C) = CCW walk
--   C's ->next = CCW = Haskell's prevCW
--   C's ->prev = CW  = Haskell's nextCW
---------------------------------------------------------------------

-- | Walk 5 CW steps from start at vertex w, compute 5-bit colour.
-- Matches C code's get_colour_prev_5 (which walks via ->prev = CW).
colourCW5 :: DualGraph -> Vertex -> Vertex -> Int
colourCW5 g w start = go 4 start 0
  where
    go (-1) _ acc = acc
    go i nbr acc =
        let acc' = if deg g nbr == 5 then acc .|. (1 `shiftL` i) else acc
        in go (i - 1) (nextCW g w nbr) acc'

-- | Walk 5 CCW steps from start at vertex w, compute 5-bit colour.
-- Matches C code's get_colour_next_5 (which walks via ->next = CCW).
colourCCW5 :: DualGraph -> Vertex -> Vertex -> Int
colourCCW5 g w start = go 4 start 0
  where
    go (-1) _ acc = acc
    go i nbr acc =
        let acc' = if deg g nbr == 5 then acc .|. (1 `shiftL` i) else acc
        in go (i - 1) (prevCW g w nbr) acc'

-- | Path colour: follow straight-ahead path for 7 steps in direction dir,
-- recording degree-5 pattern as a 7-bit integer.
-- Matches C code's get_path_colour_prev (dir=DRight) /
-- get_path_colour_next (dir=DLeft).
pathColour :: DualGraph -> Dir -> Vertex -> Vertex -> Int
pathColour g dir u v = go 6 v u 0
  where
    go (-1) _ _ acc = acc
    go i atVtx fromVtx acc =
        let next = straightAhead g dir atVtx fromVtx
            acc' = if deg g next == 5
                   then acc .|. (1 `shiftL` i) else acc
        in go (i - 1) next atVtx acc'

---------------------------------------------------------------------
-- Second colour (x2): vertex colour pair for straight reductions
---------------------------------------------------------------------

-- | Compute colour pair at both endpoints of a straight reduction.
-- Returns (negate max, negate min) so higher colour = smaller tuple.
--
-- For DRight (use_prev in C):
--   C: get_colour_prev_5(e->prev->prev->invers->prev)
--   Haskell: w = advanceCW(u, v, 2), colourCW5 w (nextCW w u)
--
-- For DLeft (use_next in C):
--   C: get_colour_next_5(e->next->next->invers->next)
--   Haskell: w = advanceCW(u, v, deg-2), colourCCW5 w (prevCW w u)
straightSecondColour :: DualGraph -> Edge -> Edge -> Dir -> (Int, Int)
straightSecondColour g (u1, v1) (u2, v2) DRight =
    let w1 = advanceCW g u1 v1 2
        c1 = colourCW5 g w1 (nextCW g w1 u1)
        w2 = advanceCW g u2 v2 2
        c2 = colourCW5 g w2 (nextCW g w2 u2)
    in (negate (max c1 c2), negate (min c1 c2))
straightSecondColour g (u1, v1) (u2, v2) DLeft =
    let w1 = advanceCW g u1 v1 (deg g u1 - 2)
        c1 = colourCCW5 g w1 (prevCW g w1 u1)
        w2 = advanceCW g u2 v2 (deg g u2 - 2)
        c2 = colourCCW5 g w2 (prevCW g w2 u2)
    in (negate (max c1 c2), negate (min c1 c2))

-- | Second colour for a reduction.
--
-- L0: (negate max, negate min) of colour pair — C code has explicit
--     second colour in is_best_L0_reduction (lines 10315-10320).
--
-- L_i (i>=1): (negate max, 0) — C code's is_best_straight_reduction
--     has NO second colour; goes directly from first colour to path colour.
--     The min component must be 0 so it doesn't create a spurious tiebreaker.
--
-- B00: (negate max, negate min) with ASYMMETRIC colour functions.
--     C code's is_best_bent_zero_reduction (lines 10751-10779):
--     For use_next: endpoint1 uses get_colour_next_5, endpoint2 uses get_colour_prev_5.
--     Starting edges: test_edge->invers->next/prev (one hop, not two).
--
-- B_{i,j} (i+j>0): (negate max, 0) — same as L_i, no second colour.
secondColour :: DualGraph -> ExpKind -> Edge -> Dir -> (Int, Int)
secondColour g (L 0) (u, v) dir =
    straightSecondColour g (u, v) (v, u) dir
secondColour g (L i) (u, v) dir =
    let pi' = computeStraightPath g (u, v) dir (i + 3)
        mp = mainPath pi'
        u2 = last mp
        v2 = mp !! (length mp - 2)
        (negMax, _) = straightSecondColour g (u, v) (u2, v2) dir
    in (negMax, 0)  -- No second colour for L_i (i>=1)
secondColour g (B 0 0) (u, v) dir =
    bentZeroSecondColour g (u, v) dir
secondColour g (B i j) (u, v) dir =
    let pi' = computeBentPath g (u, v) dir i (i + j)
        mp = mainPath pi'
        u2 = last mp
        v2 = mp !! (length mp - 2)
        -- For bent types with i+j>0, use the same asymmetric B00 colour
        -- pattern but no second colour (just max).
        (negMax, _) = bentSecondColour g (u, v) (u2, v2) dir
    in (negMax, 0)
secondColour _ F _ _ = (0, 0)

---------------------------------------------------------------------
-- B00 second colour: asymmetric endpoint colours
---------------------------------------------------------------------

-- | Compute B00 colour pair. The B00 reduction has TWO endpoints
-- connected by a bent path. The colour functions are ASYMMETRIC:
-- one endpoint uses get_colour_next_5 (CCW), the other uses
-- get_colour_prev_5 (CW). Starting edges are just ->invers->next/prev
-- (one hop from the edge endpoint, not two like straight types).
--
-- C code (is_best_bent_zero_reduction, lines 10751-10754):
--   if(use_next) {
--       colour1 = get_colour_next_5(test_edge1->invers->next);
--       colour2 = get_colour_prev_5(test_edge2->invers->prev);
--   }
-- test_edge1 = edge from first 5-vertex to its path neighbor
-- test_edge2 = edge from last 5-vertex to its path neighbor
bentZeroSecondColour :: DualGraph -> Edge -> Dir -> (Int, Int)
bentZeroSecondColour g (u, v) dir =
    let pi' = computeBentZeroPath g (u, v) dir
        mp = mainPath pi'
        -- u = first 5-vertex, mp[1] = first path neighbor
        -- u2 = last 5-vertex (mp[4]), v2 = last path neighbor (mp[3])
        u2 = mp !! 4
        v2 = mp !! 3
    in bentSecondColour g (u, v) (u2, v2) dir

-- | Compute asymmetric bent colour pair for any bent reduction.
-- (u1, v1) = first endpoint edge, (u2, v2) = second endpoint edge.
--
-- For DLeft (use_next in C):
--   colour1 = get_colour_next_5(test_edge1->invers->next)
--     = colourCCW5 at v1, starting from prevCW(g, v1, u1)
--   colour2 = get_colour_prev_5(test_edge2->invers->prev)
--     = colourCW5 at v2, starting from nextCW(g, v2, u2)
--
-- For DRight (use_prev in C):
--   colour1 = get_colour_prev_5(test_edge1->invers->prev)
--     = colourCW5 at v1, starting from nextCW(g, v1, u1)
--   colour2 = get_colour_next_5(test_edge2->invers->next)
--     = colourCCW5 at v2, starting from prevCW(g, v2, u2)
bentSecondColour :: DualGraph -> Edge -> Edge -> Dir -> (Int, Int)
bentSecondColour g (u1, v1) (u2, v2) DLeft =
    let c1 = colourCCW5 g v1 (prevCW g v1 u1)
        c2 = colourCW5  g v2 (nextCW g v2 u2)
    in (negate (max c1 c2), negate (min c1 c2))
bentSecondColour g (u1, v1) (u2, v2) DRight =
    let c1 = colourCW5  g v1 (nextCW g v1 u1)
        c2 = colourCCW5 g v2 (prevCW g v2 u2)
    in (negate (max c1 c2), negate (min c1 c2))

---------------------------------------------------------------------
-- Third colour (x3): path colour for straight reductions
---------------------------------------------------------------------

-- | Path colour at both endpoints, return negate of max.
--
-- For DRight (use_prev in C):
--   C: get_path_colour_prev(e->next->next)
--   Haskell: w = advanceCW(u, v, deg-2), pathColour DRight u w
--
-- For DLeft (use_next in C):
--   C: get_path_colour_next(e->prev->prev)
--   Haskell: w = advanceCW(u, v, 2), pathColour DLeft u w
straightThirdColour :: DualGraph -> Edge -> Edge -> Dir -> Int
straightThirdColour g (u1, v1) (u2, v2) DRight =
    let w1 = advanceCW g u1 v1 (deg g u1 - 2)
        pc1 = pathColour g DRight u1 w1
        w2 = advanceCW g u2 v2 (deg g u2 - 2)
        pc2 = pathColour g DRight u2 w2
    in negate (max pc1 pc2)
straightThirdColour g (u1, v1) (u2, v2) DLeft =
    let w1 = advanceCW g u1 v1 2
        pc1 = pathColour g DLeft u1 w1
        w2 = advanceCW g u2 v2 2
        pc2 = pathColour g DLeft u2 w2
    in negate (max pc1 pc2)

-- | Third colour for a reduction.
-- Straight types use straightThirdColour (symmetric, both endpoints same dir).
-- Bent types use bentThirdColour (asymmetric, endpoints use opposite dirs).
thirdColour :: DualGraph -> ExpKind -> Edge -> Dir -> Int
thirdColour g (L 0) (u, v) dir =
    straightThirdColour g (u, v) (v, u) dir
thirdColour g (L i) (u, v) dir =
    let pi' = computeStraightPath g (u, v) dir (i + 3)
        mp = mainPath pi'
        u2 = last mp
        v2 = mp !! (length mp - 2)
    in straightThirdColour g (u, v) (u2, v2) dir
thirdColour g (B 0 0) (u, v) dir =
    let pi' = computeBentZeroPath g (u, v) dir
        mp = mainPath pi'
        u2 = mp !! 4
        v2 = mp !! 3
    in bentThirdColour g (u, v) (u2, v2) dir
thirdColour g (B i j) (u, v) dir =
    let pi' = computeBentPath g (u, v) dir i (i + j)
        mp = mainPath pi'
        u2 = last mp
        v2 = mp !! (length mp - 2)
    in bentThirdColour g (u, v) (u2, v2) dir
thirdColour _ F _ _ = 0

---------------------------------------------------------------------
-- Bent third colour (path colour for bent reductions)
---------------------------------------------------------------------

-- | Path colour for bent reduction endpoints.
--
-- C code is_best_third_colour_bent (lines 10214-10278):
-- The starting edges for bent path colour are different from straight:
--   test_edge = (u→v), ->invers = (v→u)
--   For use_next: get_path_colour_next(test_edge->invers->prev->prev)
--     ->prev->prev at v from u = 2 CW steps from u at v = advanceCW(g, v, u, 2)
--     Then path colour in CCW (DLeft) direction from u through that vertex.
--   For use_prev: get_path_colour_prev(test_edge->invers->next->next)
--     ->next->next at v from u = 2 CCW steps = advanceCW(g, v, u, deg-2)
--     Then path colour in CW (DRight) direction from u through that vertex.
bentThirdColour :: DualGraph -> Edge -> Edge -> Dir -> Int
bentThirdColour g (u1, v1) (u2, v2) DLeft =
    let w1 = advanceCW g v1 u1 (deg g v1 - 2)  -- 2 CCW from u1 at v1
        pc1 = pathColour g DLeft v1 w1           -- CCW path from v1
        w2 = advanceCW g v2 u2 2                 -- 2 CW from u2 at v2
        pc2 = pathColour g DRight v2 w2          -- CW path from v2
    in negate (max pc1 pc2)
bentThirdColour g (u1, v1) (u2, v2) DRight =
    let w1 = advanceCW g v1 u1 2                 -- 2 CW from u1 at v1
        pc1 = pathColour g DRight v1 w1          -- CW path from v1
        w2 = advanceCW g v2 u2 (deg g v2 - 2)   -- 2 CCW from u2 at v2
        pc2 = pathColour g DLeft v2 w2           -- CCW path from v2
    in negate (max pc1 pc2)

---------------------------------------------------------------------
-- Canonical ordering
---------------------------------------------------------------------

-- | Compute the canonical ordering key for a reduction.
-- The 5-tuple is compared lexicographically; smaller = more canonical.
-- Haskell's lazy evaluation ensures x4 (BFS) is only computed when
-- x0-x3 tie — matching the C code's cascade.
canonOrd :: Reduction -> DualGraph -> CanOrd
canonOrd (Red kind (u, v) dir) g = (x0, x1, x2, x3, x4)
  where
    x0 = reductionLength kind
    x1 = negate (longestStraight kind)
    x2 = secondColour g kind (u, v) dir
    x3 = thirdColour g kind (u, v) dir
    x4 = bfsCanonicalForm g u v dir

-- | Cheap canonical ordering: only (x0, x1, x2, x3), no BFS.
-- All four components are O(1) combinatorial invariants.
-- Used by the two-phase isCanonical to avoid BFS for >99% of reductions.
type CheapOrd = (Int, Int, (Int, Int), Int)

cheapCanonOrd :: Reduction -> DualGraph -> CheapOrd
cheapCanonOrd (Red kind (u, v) dir) g = (x0, x1, x2, x3)
  where
    x0 = reductionLength kind
    x1 = negate (longestStraight kind)
    x2 = secondColour g kind (u, v) dir
    x3 = thirdColour g kind (u, v) dir

-- | Extract just the BFS canonical form (x4) from a reduction.
redBFS :: Reduction -> DualGraph -> BFSCode
redBFS (Red _ (u, v) dir) g = bfsCanonicalForm g u v dir

---------------------------------------------------------------------
-- BFS Canonical Form (STUArray-based, matching C code's approach)
---------------------------------------------------------------------

-- | Find index of v in w's CW neighbor list (flat array).
-- Linear scan over 5-6 entries. Used instead of Expansion.indexOf
-- to avoid the module dependency and allow inlining.
indexOfFlat :: UArray Int Int -> Int -> Int -> Int -> Int
indexOfFlat flat w d v = go 0
  where
    base = w * 6
    go i | i >= d    = error $ "indexOfFlat: " ++ show v ++ " not in " ++ show w
         | flat ! (base + i) == v = i
         | otherwise = go (i + 1)
{-# INLINE indexOfFlat #-}

-- | Compute the "first-row colour" for starting edge (u, v, dir).
-- The first row of the BFS code is entirely determined by deg(u) and the
-- degree sequence of u's non-v neighbors in BFS traversal order. Since
-- all degrees are 5 or 6, we pack this into a single Int:
--   bits = deg(u) * 32 + binary encoding of neighbor degrees (0=deg5, 1=deg6)
-- This preserves lexicographic ordering of the first BFS row.
-- Only edges with the minimum edgeColour can produce the minimum full BFS code.
edgeColour :: UArray Int Int -> UArray Int Int -> Int -> Int -> Dir -> Int
edgeColour flat dfl u v dir = dU * 32 + go startIdx (dU - 1) 0
  where
    dU = dfl ! u
    step = case dir of { DRight -> 1; DLeft -> -1 }
    refIdx = indexOfFlat flat u dU v
    startIdx = (refIdx + step + dU) `mod` dU
    go _ 0 acc = acc
    go idx remaining acc =
        let nbr = flat ! (u * 6 + idx)
            bit = dfl ! nbr - 5   -- 0 for deg-5, 1 for deg-6
            nextIdx = (idx + step + dU) `mod` dU
        in go nextIdx (remaining - 1) (acc * 2 + bit)
{-# INLINE edgeColour #-}

-- | Compute BFS canonical form starting from directed edge (u -> v)
-- walking in direction dir.
--
-- DRight = CW walk (matches C's testcanon_mirror, used with use_prev)
-- DLeft  = CCW walk (matches C's testcanon, used with use_next)
--
-- Uses STUArray internally — no IntMap, no list allocation.
-- Matches C code's testcanon/testcanon_mirror (lines 2246-2522).
bfsCanonicalForm :: DualGraph -> Vertex -> Vertex -> Dir -> BFSCode
bfsCanonicalForm g u v dir = BFS $ runSTUArray $ do
    let nv = numVertices g
        codeLen = 6 * nv - 12
        maxN = nv + 1
        flat = adjFlat g
        dfl = degFlat g
        step = case dir of { DRight -> 1; DLeft -> -1 }

    numArr <- newArray (0, nv - 1) 0 :: ST s (STUArray s Int Int)
    queueV <- newArray (1, nv) 0     :: ST s (STUArray s Int Int)
    queueR <- newArray (1, nv) 0     :: ST s (STUArray s Int Int)
    codeArr <- newArray (0, codeLen - 1) 0 :: ST s (STUArray s Int Int)

    writeArray numArr u 1
    writeArray numArr v 2
    writeArray queueV 1 u
    writeArray queueR 1 v
    writeArray queueV 2 v
    writeArray queueR 2 u

    let process curNum nextNum codeIdx
          | curNum > nv = return codeArr
          | otherwise = do
              w <- readArray queueV curNum
              ref <- readArray queueR curNum
              let d = dfl ! w
                  refIdx = indexOfFlat flat w d ref
                  startIdx = (refIdx + step + d) `mod` d
              (nextNum', codeIdx') <-
                  walkNbrs step flat dfl maxN numArr queueV queueR codeArr
                           w d startIdx (d - 1) nextNum codeIdx
              writeArray codeArr codeIdx' 0
              process (curNum + 1) nextNum' (codeIdx' + 1)

    process 1 3 0

-- | BFS from a single starting edge, returning both code and vertex numbering.
-- Returns (BFSCode, numArray) where numArray[v] = BFS number of vertex v.
bfsWithNumbering :: DualGraph -> Vertex -> Vertex -> Dir
                 -> (BFSCode, UArray Int Int)
bfsWithNumbering g u v dir = runST $ do
    let nv = numVertices g
        codeLen = 6 * nv - 12
        maxN = nv + 1
        flat = adjFlat g
        dfl = degFlat g
        step = case dir of { DRight -> 1; DLeft -> -1 }

    numArr <- newArray (0, nv - 1) 0 :: ST s (STUArray s Int Int)
    queueV <- newArray (1, nv) 0     :: ST s (STUArray s Int Int)
    queueR <- newArray (1, nv) 0     :: ST s (STUArray s Int Int)
    codeArr <- newArray (0, codeLen - 1) 0 :: ST s (STUArray s Int Int)

    writeArray numArr u 1
    writeArray numArr v 2
    writeArray queueV 1 u
    writeArray queueR 1 v
    writeArray queueV 2 v
    writeArray queueR 2 u

    let process curNum nextNum codeIdx
          | curNum > nv = return ()
          | otherwise = do
              w <- readArray queueV curNum
              ref <- readArray queueR curNum
              let d = dfl ! w
                  refIdx = indexOfFlat flat w d ref
                  startIdx = (refIdx + step + d) `mod` d
              (nextNum', codeIdx') <-
                  walkNbrs step flat dfl maxN numArr queueV queueR codeArr
                           w d startIdx (d - 1) nextNum codeIdx
              writeArray codeArr codeIdx' 0
              process (curNum + 1) nextNum' (codeIdx' + 1)

    process 1 3 0
    frozenCode <- freeze codeArr
    frozenNum  <- freeze numArr
    return (BFS frozenCode, frozenNum)

-- | Walk neighbors of vertex w in BFS order, writing codes to codeArr.
-- Shared between bfsCanonicalForm and bfsWithNumbering.
-- Uses flat array index stepping (no indexOf per neighbor).
walkNbrs :: Int -> UArray Int Int -> UArray Int Int -> Int
         -> STUArray s Int Int -> STUArray s Int Int -> STUArray s Int Int
         -> STUArray s Int Int
         -> Int -> Int -> Int -> Int -> Int -> Int
         -> ST s (Int, Int)
walkNbrs step flat dfl maxN numArr queueV queueR codeArr
         w d idx remaining nextNum codeIdx
    | remaining <= 0 = return (nextNum, codeIdx)
    | otherwise = do
        let nbr = flat ! (w * 6 + idx)
        num <- readArray numArr nbr
        if num > 0
            then do
                writeArray codeArr codeIdx num
                let nextIdx = (idx + step + d) `mod` d
                walkNbrs step flat dfl maxN numArr queueV queueR codeArr
                         w d nextIdx (remaining - 1) nextNum (codeIdx + 1)
            else do
                writeArray codeArr codeIdx (dfl ! nbr + maxN)
                writeArray numArr nbr nextNum
                writeArray queueV nextNum nbr
                writeArray queueR nextNum w
                let nextIdx = (idx + step + d) `mod` d
                walkNbrs step flat dfl maxN numArr queueV queueR codeArr
                         w d nextIdx (remaining - 1) (nextNum + 1) (codeIdx + 1)

---------------------------------------------------------------------
-- Automorphism group computation
---------------------------------------------------------------------

-- | Compute canonical BFS form and automorphism group.
-- Tries all (vertex, neighbor, direction) starting configurations.
-- Starting edges that produce the minimum BFS code define automorphisms.
--
-- Matches C code's canon() (lines 3267-3402): pre-allocates work arrays
-- in a single runST block, does interleaved BFS + comparison (like testcanon),
-- aborting early for GT cases (~99%). Only freezes numbering for automorphisms.
--
-- Optimization: edge colour pre-filter. The first BFS row is determined by
-- deg(u) and the degree pattern of u's neighbors. Only edges with the minimum
-- first-row colour can produce the minimum BFS code. This typically reduces
-- starting edges from ~390 to ~20 (at nv=35).
canonicalBFSAndGroup :: DualGraph -> CanonResult
canonicalBFSAndGroup g = runST $ do
    let nv = numVertices g
        codeLen = 6 * nv - 12
        maxN = nv + 1
        flat = adjFlat g
        dfl = degFlat g
        -- Edge colour pre-filter: compute first-row colour for each starting
        -- edge, keep only those with the minimum colour.
        allWithColour = [ (edgeColour flat dfl u v d, (u, v, d))
                        | u <- [0..nv-1]
                        , i <- [0 .. dfl ! u - 1]
                        , let v = flat ! (u * 6 + i)
                        , d <- [DLeft, DRight]
                        ]
        !minCol = minimum (map fst allWithColour)
        allStarts = [s | (c, s) <- allWithColour, c == minCol]
        !nTotal = length allWithColour
        !nFiltered = length allStarts

    -- Pre-allocate reusable work arrays (like C's stack-allocated arrays)
    numArr  <- newArray (0, nv - 1) 0     :: ST s (STUArray s Int Int)
    queueV  <- newArray (1, nv) 0         :: ST s (STUArray s Int Int)
    queueR  <- newArray (1, nv) 0         :: ST s (STUArray s Int Int)
    codeArr <- newArray (0, codeLen - 1) 0 :: ST s (STUArray s Int Int)

    -- Freeze STUArray to UArray (type-disambiguated helper)
    let freezeNum :: STUArray s Int Int -> ST s (UArray Int Int)
        freezeNum = freeze

    -- Clear only the numbered vertices in numArr (tracked via queueV)
    let clearNum lastNum =
            forM_ [1..lastNum] $ \k -> do
                vtx <- readArray queueV k
                writeArray numArr vtx 0

    -- Full BFS: compute code into codeArr, numbering into numArr.
    -- Returns the number of vertices numbered (= nv when complete).
    let fullBFS u0 v0 dir0 = do
            let step = case dir0 of { DRight -> 1; DLeft -> -1 }
            writeArray numArr u0 1
            writeArray numArr v0 2
            writeArray queueV 1 u0
            writeArray queueR 1 v0
            writeArray queueV 2 v0
            writeArray queueR 2 u0
            let process curNum nextNum codeIdx
                  | curNum > nv = return (nextNum - 1)
                  | otherwise = do
                      w <- readArray queueV curNum
                      ref <- readArray queueR curNum
                      let d = dfl ! w
                          refIdx = indexOfFlat flat w d ref
                          startIdx = (refIdx + step + d) `mod` d
                      (nextNum', codeIdx') <-
                          walkNbrs step flat dfl maxN numArr queueV queueR codeArr
                                   w d startIdx (d - 1) nextNum codeIdx
                      writeArray codeArr codeIdx' 0
                      process (curNum + 1) nextNum' (codeIdx' + 1)
            process 1 3 0

    -- Interleaved BFS + comparison against refCode (like C's testcanon).
    -- Returns (Ordering, lastNumbered). Aborts early on mismatch (GT/LT).
    let compareBFS u0 v0 dir0 refCode = do
            let step = case dir0 of { DRight -> 1; DLeft -> -1 }
            writeArray numArr u0 1
            writeArray numArr v0 2
            writeArray queueV 1 u0
            writeArray queueR 1 v0
            writeArray queueV 2 v0
            writeArray queueR 2 u0
            let process curNum nextNum codeIdx
                  | curNum > nv = return (EQ, nextNum - 1)
                  | otherwise = do
                      w <- readArray queueV curNum
                      ref <- readArray queueR curNum
                      let d = dfl ! w
                          refIdx = indexOfFlat flat w d ref
                          startIdx = (refIdx + step + d) `mod` d
                      result <- walkCmp step flat dfl maxN numArr queueV queueR
                                        refCode w d startIdx (d - 1) nextNum codeIdx
                      case result of
                          Left (ord, lastNum) -> return (ord, lastNum)
                          Right (nextNum', codeIdx') -> do
                              -- Check separator against refCode
                              let refVal = refCode ! codeIdx'
                              if 0 < refVal then return (LT, nextNum' - 1)
                              else if 0 > refVal then return (GT, nextNum' - 1)
                              else process (curNum + 1) nextNum' (codeIdx' + 1)
            process 1 3 0

    -- Process all starting edges
    case allStarts of
        [] -> error "canonicalBFSAndGroup: no starting edges"
        ((u0,v0,d0):rest) -> do
            _ <- fullBFS u0 v0 d0
            bestRef <- freezeNum codeArr
            nm0 <- freezeNum numArr
            clearNum nv

            -- Main loop: compare each starting edge against best.
            -- Track GT/EQ/LT counts for instrumentation.
            let go bestCode matches gtC eqC ltC [] = return (bestCode, matches, gtC, eqC, ltC)
                go bestCode matches gtC eqC ltC ((u,v,d):rs) = do
                    (ord, lastNum) <- compareBFS u v d bestCode
                    case ord of
                        GT -> do
                            clearNum lastNum
                            go bestCode matches (gtC+1) eqC ltC rs
                        EQ -> do
                            nm <- freezeNum numArr
                            clearNum lastNum
                            go bestCode ((nm, d) : matches) gtC (eqC+1) ltC rs
                        LT -> do
                            -- New minimum found. Redo full BFS (rare).
                            clearNum lastNum
                            _ <- fullBFS u v d
                            newBest <- freezeNum codeArr
                            nm <- freezeNum numArr
                            clearNum nv
                            go newBest [(nm, d)] gtC eqC (ltC+1) rs

            (bestCode, matches, gtCount, eqCount, ltCount) <- go bestRef [(nm0, d0)] 0 0 0 rest

            -- Build automorphisms from UArray numberings
            let (refNumMap, bestDir) = head matches
                refInv = array (1, nv) [(refNumMap ! vv, vv) | vv <- [0..nv-1]]
                         :: UArray Int Int
                auts = [ Aut perm ori
                       | (numMap, dd) <- matches
                       , let perm = array (0, nv-1)
                               [(vv, refInv ! (numMap ! vv)) | vv <- [0..nv-1]]
                               :: UArray Int Int
                             ori = if dd == bestDir then Preserving else Reversing
                       ]
                nbop = length [a | a <- auts, autOrientation a == Preserving]

            return $ CanonResult (BFS bestCode) auts nbop
                                 nTotal nFiltered gtCount eqCount ltCount

-- | Walk neighbors comparing each BFS code element against a reference.
-- Returns Left (Ordering, lastNumbered) on mismatch (early abort),
-- Right (nextNum, codeIdx) when all neighbors matched.
walkCmp :: Int -> UArray Int Int -> UArray Int Int -> Int
        -> STUArray s Int Int -> STUArray s Int Int -> STUArray s Int Int
        -> UArray Int Int
        -> Int -> Int -> Int -> Int -> Int -> Int
        -> ST s (Either (Ordering, Int) (Int, Int))
walkCmp step flat dfl maxN numArr queueV queueR refCode
        w d idx remaining nextNum codeIdx
    | remaining <= 0 = return $ Right (nextNum, codeIdx)
    | otherwise = do
        let nbr = flat ! (w * 6 + idx)
        num <- readArray numArr nbr
        let val = if num > 0 then num else dfl ! nbr + maxN
            refVal = refCode ! codeIdx
        if val < refVal then return $ Left (LT, nextNum - 1)
        else if val > refVal then return $ Left (GT, nextNum - 1)
        else do
            if num == 0
                then do
                    writeArray numArr nbr nextNum
                    writeArray queueV nextNum nbr
                    writeArray queueR nextNum w
                    let nextIdx = (idx + step + d) `mod` d
                    walkCmp step flat dfl maxN numArr queueV queueR refCode
                            w d nextIdx (remaining - 1) (nextNum + 1) (codeIdx + 1)
                else do
                    let nextIdx = (idx + step + d) `mod` d
                    walkCmp step flat dfl maxN numArr queueV queueR refCode
                            w d nextIdx (remaining - 1) nextNum (codeIdx + 1)

---------------------------------------------------------------------
-- Rule 2: Expansion orbit filtering
---------------------------------------------------------------------

-- | Apply an automorphism to an expansion site.
-- Orientation-preserving: keep direction. Orientation-reversing: flip direction.
applyAutToExp :: Automorphism -> Expansion -> Expansion
applyAutToExp (Aut sigma Preserving) (Exp kind (u, v) dir) =
    Exp kind (sigma ! u, sigma ! v) dir
applyAutToExp (Aut sigma Reversing) (Exp kind (u, v) dir) =
    Exp kind (sigma ! u, sigma ! v) (flipDir dir)

-- | Compute the "other triple" for an expansion — the expansion from the
-- other endpoint of the central path that produces the same child graph.
--
-- L types: same direction (path is symmetric)
-- B types: flip direction and swap (i,j) (bend reverses direction)
computeOtherTriple :: DualGraph -> Expansion -> Maybe Expansion
computeOtherTriple g (Exp (L 0) (u, v) d) =
    let pi' = computeStraightPath g (u, v) d 3
        par = parallelPath pi'
        mp = mainPath pi'
        -- The other degree-5 vertex is par[2] (= q2)
        otherU = par !! 2
        -- The starting neighbor from otherU: from C code's edge traversal chain
        -- e->invers->next->next->invers->next, the neighbor is one step in the
        -- flipDir direction from mp[2] (= p2) at otherU's neighbor list.
        otherV = sideNbr g (flipDir d) otherU (mp !! 2)
    in Just (Exp (L 0) (otherU, otherV) d)  -- same direction for L

computeOtherTriple g (Exp (L i) (u, v) d) =
    let pi' = computeStraightPath g (u, v) d (i + 3)
        mp = mainPath pi'
        par = parallelPath pi'
        -- Other degree-5 vertex is par[last] (= q_{i+2})
        otherU = last par
        -- Starting neighbor: one step in flipDir direction from mp[last] at otherU
        otherV = sideNbr g (flipDir d) otherU (last mp)
    in Just (Exp (L i) (otherU, otherV) d)  -- same direction for L

computeOtherTriple g (Exp (B 0 0) (u, v) d) =
    let pi' = computeBentZeroPath g (u, v) d
        mp = mainPath pi'
        -- Other degree-5 vertex is path[4]
        otherU = mp !! 4
        otherV = mp !! 3  -- penultimate on path
    in Just (Exp (B 0 0) (otherU, otherV) (flipDir d))  -- flip for B

computeOtherTriple g (Exp (B i j) (u, v) d) =
    case computeBentPathSafe g (u, v) d i (i + j) of
        Nothing -> Nothing
        Just pi' ->
            let mp = mainPath pi'
                otherU = last mp
                otherV = mp !! (length mp - 2)
            in Just (Exp (B j i) (otherU, otherV) (flipDir d))  -- swap i,j + flip

computeOtherTriple _ (Exp F _ _) = Nothing  -- F has no other triple

-- | Filter expansion sites by Rule 2: keep one representative per orbit.
-- An orbit is generated by: (1) automorphism group images, and
-- (2) other-triple equivalence.
filterByRule2 :: DualGraph -> [Automorphism] -> [Expansion] -> [Expansion]
filterByRule2 g auts exps = go Set.empty exps
  where
    go _ [] = []
    go seen (e:rest)
        | expKey e `Set.member` seen = go seen rest
        | otherwise =
            let orbit = fullOrbit e
                seen' = foldl' (\s x -> Set.insert x s) seen orbit
            in e : go seen' rest

    -- Key for Set membership: directed edges are distinct, so use full key.
    expKey (Exp kind (u, v) dir) = (kind, u, v, dir)

    fullOrbit e =
        let autImages = map (\a -> expKey (applyAutToExp a e)) auts
            otherTriple = computeOtherTriple g e
            otherImages = case otherTriple of
                Nothing -> []
                Just e' -> map (\a -> expKey (applyAutToExp a e')) auts
        in autImages ++ otherImages

---------------------------------------------------------------------
-- All reductions of a graph
---------------------------------------------------------------------

-- | Enumerate all valid reduction sites of graph g.
-- Each reduction is a (kind, edge, direction) triple representing
-- a site where the graph could be reduced to a smaller fullerene.
allReductions :: DualGraph -> [Reduction]
allReductions = allReductionsUpTo 5

-- | Enumerate all reductions of length up to maxLen.
-- Using a smaller maxLen skips expensive bent reduction enumeration
-- when shorter reductions are guaranteed to dominate in canonOrd.
allReductionsUpTo :: Int -> DualGraph -> [Reduction]
allReductionsUpTo maxLen g = l0Reds ++ lReds ++ b00Reds ++ bigBReds
  where
    l0Reds =
        [ Red (L 0) (u, v) d
        | u <- degree5 g
        , i <- [0 .. deg g u - 1]
        , let v = nbrAt g u i
        , deg g v == 5
        , d <- [DLeft, DRight]
        , isValidL0Direction g (u, v) d
        ]

    lReds = if maxLen >= 2 then straightReductions maxLen g else []

    b00Reds =
        if maxLen >= 2
        then [ Red (B 0 0) (u, v) d
             | u <- degree5 g
             , i <- [0 .. deg g u - 1]
             , let v = nbrAt g u i
             , d <- [DLeft, DRight]
             , isValidB00Reduction g (u, v) d
             ]
        else []

    bigBReds = if maxLen >= 3 then bentReductions maxLen g else []

-- | Fast check: does the graph have any valid L0 reduction?
-- Used for early exit in isCanonical: when the inverse has reductionLength > 1
-- but the child has an L0 reduction (x0=1), the inverse can never be canonical.
hasAnyL0Reduction :: DualGraph -> Bool
hasAnyL0Reduction g = any hasL0Edge (degree5 g)
  where
    hasL0Edge u = any (\i -> let v = nbrAt g u i
                             in deg g v == 5 && anyDir u v)
                      [0 .. deg g u - 1]
    anyDir u v = isValidL0Direction g (u, v) DLeft
              || isValidL0Direction g (u, v) DRight

-- | Check if L0 direction is valid.
-- Matches C's has_L0_reduction: flanking vertices 2 positions away
-- from the edge endpoints (in the walk direction) must both be degree 6.
--
-- DRight (use_prev): 2 CW from opposite endpoint at each vertex.
-- DLeft  (use_next): 2 CCW from opposite endpoint at each vertex.
isValidL0Direction :: DualGraph -> Edge -> Dir -> Bool
isValidL0Direction g (u, v) DRight =
    deg g (advanceCW g u v 2) == 6 && deg g (advanceCW g v u 2) == 6
isValidL0Direction g (u, v) DLeft =
    deg g (advanceCW g u v (deg g u - 2)) == 6 &&
    deg g (advanceCW g v u (deg g v - 2)) == 6

-- | Check if B_{0,0} site is valid for REDUCTION.
-- Uses a mix of structural and reduction-site criteria.
-- The C code's has_bent_zero_reduction checks specific vertex degrees
-- along the path. We approximate with:
--   1. No self-intersection in path
--   2. Far endpoint is degree-5
--   3. First and last path vertices differ
--   4. Path and parallel path are disjoint (conservative: prevents
--      spurious B00 reductions at tight configurations)
isValidB00Reduction :: DualGraph -> Edge -> Dir -> Bool
isValidB00Reduction g (u, v) dir =
    let pi' = computeBentZeroPath g (u, v) dir
        path = mainPath pi'
        par  = parallelPath pi'
        pathSet = IS.fromList path
    in IS.size pathSet == length path    -- no duplicates (was O(L²) nub)
    && deg g (last path) == 5
    && head path /= last path
    && IS.null (pathSet `IS.intersection` IS.fromList par)  -- was O(L²) intersect

---------------------------------------------------------------------
-- Straight reduction enumeration (L_i, i >= 1)
---------------------------------------------------------------------

-- | Enumerate all L_i (i >= 1) reduction sites.
-- Matches C code's has_short_straight_reduction / has_short_straight_reduction_L1:
-- for each degree-5 vertex u and each neighbor v, follow straight-ahead
-- to find the next degree-5 vertex. If found, check flanking degree-6
-- at both endpoints. The path direction doesn't matter since all
-- intermediate vertices must be degree > 5 (otherwise a shorter
-- reduction would exist, and we stop at the first degree-5 vertex).
straightReductions :: Int -> DualGraph -> [Reduction]
straightReductions maxRedLen g =
    [ Red (L (dist - 1)) (u, v) d
    | u <- degree5 g
    , i <- [0 .. deg g u - 1]
    , let v = nbrAt g u i
    , deg g v /= 5  -- distance 1 = L0, handled separately
    , Just (w, prevW, dist) <- [followStraightToFive g u v maxRedLen]
    , d <- [DLeft, DRight]
    , isValidStraightReduction g (u, v) (w, prevW) d
    ]

-- | Follow straight-ahead from u through v until a degree-5 vertex is found.
-- Returns Just (endpoint, prevVertex, distance) or Nothing if not found
-- within maxDist steps. Stops at first degree-5 vertex (shorter reductions
-- take priority). Includes cycle detection for longer paths.
followStraightToFive :: DualGraph -> Vertex -> Vertex -> Int
                     -> Maybe (Vertex, Vertex, Int)
followStraightToFive g u v maxDist = go u v 1 (IS.singleton u)
  where
    go from to dist visited
        | dist > maxDist = Nothing
        | IS.member to visited = Nothing  -- cycle (was O(dist) elem)
        | deg g to == 5 = Just (to, from, dist)
        | otherwise =
            let next = straightAhead g DLeft to from  -- dir irrelevant at degree-6
            in go to next (dist + 1) (IS.insert to visited)

-- | Check if a straight reduction direction is valid (flanking degree-6).
-- Matches C's has_L0_reduction / has_short_straight_reduction: the vertex
-- 2 CW (DRight) or 2 CCW (DLeft) from the path direction at BOTH endpoints
-- must be degree-6.
isValidStraightReduction :: DualGraph -> Edge -> Edge -> Dir -> Bool
isValidStraightReduction g (u, v) (w, prevW) DRight =
    deg g (advanceCW g u v 2) == 6 && deg g (advanceCW g w prevW 2) == 6
isValidStraightReduction g (u, v) (w, prevW) DLeft =
    deg g (advanceCW g u v (deg g u - 2)) == 6 &&
    deg g (advanceCW g w prevW (deg g w - 2)) == 6

---------------------------------------------------------------------
-- Bent reduction enumeration (B_{i,j}, i+j > 0)
---------------------------------------------------------------------

-- | Enumerate all B_{i,j} (i+j > 0) reduction sites using reduction-site
-- criteria. Matches C code's has_short_bent_reduction:
-- - Path non-self-intersection (via computeBentPathSafe)
-- - Far endpoint degree-5
-- - Flanking vertex at far endpoint degree-6 (OPPOSITE side from straight
--   reductions, because the bend flips which side is the "colour" side)
-- - Path/parallel path disjoint
bentReductions :: Int -> DualGraph -> [Reduction]
bentReductions maxRedLen g =
    [ Red (B i j) (u, v) d
    | u <- degree5 g
    , k <- [0 .. deg g u - 1]
    , let v = nbrAt g u k
    , d <- [DLeft, DRight]
    , bentLen <- [1 .. maxRedLen - 2]  -- i+j >= 1, i+j+2 <= maxRedLen
    , bentPos <- [0 .. bentLen]
    , let i = bentPos
          j = bentLen - bentPos
    , isValidBentReduction g (u, v) d i j
    ]

-- | Check if a B_{i,j} reduction site is valid using reduction-site criteria.
-- This is DIFFERENT from expansion-site criteria (canBentPath) which checks
-- parallel path degrees and specific adjacencies for the replacement ops.
-- Reduction-site criteria (matching C code's has_short_bent_reduction):
--   1. Path computable without self-intersection
--   2. All intermediate vertices are degree-6 (degree-5 means shorter reduction)
--   3. Far endpoint is degree-5
--   4. Flanking vertex at far endpoint is degree-6 and not on the path
-- NOTE: path/parallel disjointness is NOT checked — that's expansion-specific.
isValidBentReduction :: DualGraph -> Edge -> Dir -> Int -> Int -> Bool
isValidBentReduction g (u0, v0) dir i j =
    case computeBentPathSafe g (u0, v0) dir i (i + j) of
        Nothing -> False
        Just (PathInfo path _par) ->
            let lastV = last path
                prevW = path !! (length path - 2)
                intermediates = init (tail path)  -- all except first and last
                flankV = bentEndpointFlank g lastV prevW dir
            in deg g lastV == 5
            && head path /= lastV
            && all (\v -> deg g v == 6) intermediates
            && deg g flankV == 6
            && flankV `notElem` path

-- | Flanking vertex at the far endpoint of a bent reduction.
-- Bent reductions use the OPPOSITE side from straight reductions because
-- the bend flips the "colour" side of the path.
-- C code (has_short_bent_reduction):
--   use_next (DLeft):  e->invers->prev->prev->end = 2 CW from incoming
--   use_prev (DRight): e->invers->next->next->end  = 2 CCW from incoming
bentEndpointFlank :: DualGraph -> Vertex -> Vertex -> Dir -> Vertex
bentEndpointFlank g w prevW DLeft  = advanceCW g w prevW 2
bentEndpointFlank g w prevW DRight = advanceCW g w prevW (deg g w - 2)

---------------------------------------------------------------------
-- Bounding lemmas
---------------------------------------------------------------------

-- | Extract the pair of degree-5 vertices involved in a reduction.
-- For L0: both vertices of the edge are degree-5.
-- For L_i (i>=1): start vertex u and far endpoint (found by straight walk).
-- For B_{0,0}: start vertex u and far endpoint of bent path.
-- For B_{i,j} (i+j>0): start vertex u and far endpoint of bent path.
reductionFiveVertices :: DualGraph -> Reduction -> (Vertex, Vertex)
reductionFiveVertices g (Red (L 0) (u, v) _) = (u, v)
reductionFiveVertices g (Red (L i) (u, v) _) =
    case followStraightToFive g u v (i + 1) of
        Just (w, _, _) -> (u, w)
        Nothing        -> (u, u)  -- shouldn't happen for valid reduction
reductionFiveVertices g (Red (B 0 0) (u, v) d) =
    let pi' = computeBentZeroPath g (u, v) d
    in (u, last (mainPath pi'))
reductionFiveVertices g (Red (B i j) (u, v) d) =
    case computeBentPathSafe g (u, v) d i (i + j) of
        Just (PathInfo path _) -> (u, last path)
        Nothing                -> (u, u)  -- shouldn't happen for valid reduction
reductionFiveVertices _ (Red F _ _) = (-1, -1)  -- F has no 5-vertex pair

-- | Precomputed table: max_straight_lengths[nv] = maximum possible minimum
-- pentagon distance minus 1 for a fullerene dual with nv vertices.
-- Matches C code's initialize_fuller_arrays() (lines 4732-4774).
-- See doc/MAX_STRAIGHT_LENGTHS.md for the geometric argument.
buildMaxStraightLengths :: Int -> UArray Int Int
buildMaxStraightLengths maxNV = runSTUArray $ do
    arr <- newArray (0, maxNV) 0
    let go dist minReqPrev
          | minReqPrev > maxNV = return arr
          | otherwise =
              let r = (dist - 1) `div` 2
                  base = 12 * (1 + 5 * (r + 1) * r `div` 2)
                  minReq = if even dist && dist > 2 then base + 20 else base
              in do forM_ [minReqPrev .. min (minReq - 1) maxNV] $ \i ->
                      writeArray arr i (dist - 1)
                    go (dist + 1) minReq
    go 1 0

-- | Compute the maximum expansion length for a graph, applying all
-- bounding lemmas. The result is the maximum reduction length of any
-- expansion worth trying. Compose via minimum — soundness is preserved.
--
-- C code: determine_max_pathlength_straight() (lines 4868-4966)
maxExpansionLength :: Int -> UArray Int Int -> DualGraph -> [Reduction] -> Int
maxExpansionLength maxDV mslTable g reds =
    let nv = numVertices g
        -- Base: can't expand beyond target size.
        -- L_i adds i+2 vertices, B_{i,j} adds i+j+3, so the worst case
        -- (fewest vertices per reduction length) is L_i: redLen = i+1, adds i+2.
        -- So maxRedLen such that nv + (maxRedLen+1) <= maxDV gives:
        baseBound = maxDV - nv - 1

        -- Table-based geometric bound (C code lines 4873-4882).
        -- For each candidate pathlength p, check if a reduction of length p-1
        -- is geometrically possible in a child with nv+p vertices.
        -- Subsumes both the old geoMaxRedLen and Lemma 6 (hasTwoDistantL0s).
        -- See doc/MAX_STRAIGHT_LENGTHS.md for the full analysis.
        tableBound = last $ 1 : [ p | p <- [2..baseBound]
                                    , nv + p <= maxDV
                                    , p - 1 <= mslTable ! (nv + p) ]

        hasL0  = any (\(Red k _ _) -> k == L 0) reds

        -- Length-2 reductions: L1 and B00
        len2Reds = [ r | r <- reds, reductionLength (redKind r) == 2 ]
        -- Length-2 five-vertex pairs, deduplicated by sorted pair
        len2Pairs = nub [ normalize (reductionFiveVertices g r) | r <- len2Reds ]
          where normalize (a, b) = (min a b, max a b)

        hasLen2 = not (null len2Pairs)

        -- Lemma 3: L0 exists → max length 2
        lemma3 = if hasL0 then Just 2 else Nothing

        -- Lemma 5: ≥3 independent length-2 reductions (pairwise disjoint
        -- 5-vertex sets) → max length 2
        lemma5 = if hasThreeIndependent len2Pairs then Just 2 else Nothing

        -- Lemma 4: ≥2 length-2 reductions with different 5-vertex sets
        -- → max length 3
        lemma4 = if length len2Pairs >= 2 then Just 3 else Nothing

        -- Lemma 2: parent has reduction of length d ≤ 2 → children
        -- have length ≤ d+2. With length-2 (L1/B00): max = 4
        -- With length-1 (L0): max = 3 (subsumed by Lemma 3's tighter bound)
        lemma2 = if hasLen2 then Just 4
                 else if hasL0 then Just 3
                 else Nothing

        allBounds = [baseBound, tableBound]
                 ++ catMaybes [lemma3, lemma5, lemma4, lemma2]

    in maximum [1, minimum allBounds]  -- at least 1

-- | Check if there exist 3 pairs with pairwise disjoint vertex sets.
-- C code: contains_three_indep_L1s / contains_three_indep_B00s
-- Uses backtracking search: try taking each pair, then recurse.
hasThreeIndependent :: [(Vertex, Vertex)] -> Bool
hasThreeIndependent pairs
    | length pairs < 3 = False
    | otherwise = search 0 IS.empty 0
  where
    n = length pairs
    search count _ _
        | count >= 3 = True
    search _ _ start
        | start >= n = False
    search count used start =
        let (a, b) = pairs !! start
        in if IS.member a used || IS.member b used
           then search count used (start + 1)
           else -- Try taking this pair
                search (count + 1) (IS.insert a (IS.insert b used)) (start + 1)
                -- If that doesn't find 3, try skipping this pair
                || search count used (start + 1)

catMaybes :: [Maybe a] -> [a]
catMaybes [] = []
catMaybes (Just x : xs) = x : catMaybes xs
catMaybes (Nothing : xs) = catMaybes xs

---------------------------------------------------------------------
-- Canonical test
---------------------------------------------------------------------

-- | Verdict from canonical test with rejection reason for instrumentation.
data CanonVerdict = CanonAccept | RejectL0Early | RejectNoInverse
                  | RejectPhase1 | RejectPhase2
    deriving (Eq, Show)

-- | Check if an expansion produces a canonical child.
isCanonical :: Expansion -> Int -> DualGraph -> Bool
isCanonical e parentNV g' = isCanonicalV e parentNV g' == CanonAccept

-- | Like isCanonical but returns the specific verdict for instrumentation.
--
-- After applying expansion e to graph g with parentNV vertices
-- (producing child graph g'), checks whether e's inverse is the
-- canonical reduction of g'.
--
-- Uses lazy short-circuiting to match C code's early-exit strategy:
--   Guard 1: Longer inverse + L0 exists → reject immediately.
--   Guard 2: No inverse reduction found → reject.
--   Guard 3: any non-inverse has strictly better (x0,x1,x2,x3) → reject.
--            `any` stops at the first better reduction (like C's `return 0`).
--   Guard 4: No non-inverse ties the inverse at (x0,x1,x2,x3) → accept.
--   Guard 5: BFS tiebreaker for the few tied reductions.
--
-- allReds yields L0 first, so for an L1/B00 inverse, guard 3 finds
-- a better L0 reduction after examining just 1-3 reductions.
isCanonicalV :: Expansion -> Int -> DualGraph -> CanonVerdict
isCanonicalV (Exp kind _ dir) parentNV g'
    -- Early exit: longer inverse but L0 exists → can never be canonical
    | reductionLength kind > 1 && hasAnyL0Reduction g' = RejectL0Early
    -- No inverse reduction in the child → reject
    | null inverseReds = RejectNoInverse
    -- Phase 1: any reduction strictly better than inverse at (x0,x1,x2,x3)?
    -- `any` short-circuits: stops at the first better reduction.
    | any (\r -> cheapCanonOrd r g' < invCheap) allReds = RejectPhase1
    -- No non-inverse ties the inverse → inverse is the unique minimum
    | null tiedNonInv = CanonAccept
    -- Phase 2: BFS tiebreaker among reductions tied at invCheap
    | minInvBFS <= minAllBFS = CanonAccept
    | otherwise = RejectPhase2
  where
    -- Only enumerate reductions up to the expansion's length: longer
    -- reductions can never beat the inverse in lexicographic canonOrd
    -- (since x0 = reductionLength is the first component).
    allReds = allReductionsUpTo (reductionLength kind) g'
    numNew = numNewVertices kind
    newV1 = parentNV
    newV2 = parentNV + numNew - 1

    -- Bug #14: Both vertices of the reduction's starting edge must be
    -- within the newly-added vertex range (prevents spurious matches
    -- with existing degree-5 vertices).
    -- Bug #15: Direction must match the expansion's direction — the C
    -- code tests the specific directed edge, not all directions.
    isInverseRed (Red k (a, b) d) =
        k == kind && a >= newV1 && a <= newV2
                  && b >= newV1 && b <= newV2 && d == dir
    inverseReds = filter isInverseRed allReds

    -- Cheapest canonical ordering among inverse reductions (1-4 items)
    invCheap = minimum (map (`cheapCanonOrd` g') inverseReds)

    -- Phase 2 (reached ~3% of the time): find tied non-inverse reductions
    tiedNonInv = [ r | r <- allReds, not (isInverseRed r)
                     , cheapCanonOrd r g' == invCheap ]
    tiedInv    = [ r | r <- inverseReds, cheapCanonOrd r g' == invCheap ]
    minInvBFS  = minimum (map (`redBFS` g') tiedInv)
    minAllBFS  = minimum (map (`redBFS` g') (tiedInv ++ tiedNonInv))
