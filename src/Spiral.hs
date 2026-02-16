-- | Spiral algorithm for oriented triangulations (regular and generalized).
--
-- Given an oriented triangulation (a planar graph where every face is a
-- triangle, with neighbors stored in consistent cyclic order), this module
-- computes canonical spirals — face-degree sequences that uniquely identify
-- the graph up to isomorphism.
--
-- For fullerene duals (12 degree-5 vertices, rest degree-6), the canonical
-- spiral computation is O(N).

module Spiral
  ( Graph(..), Node, mkGraph
  , deg, nbrs, next, prev, edge
  , Orientation(..), orient, apex
  , regularSpiral, canonicalSpiral
  , GeneralSpiral(..), generalSpiral, canonicalGeneralSpiral
  , startingTriples
  , rotateFwd
  ) where

import           Control.Monad       (guard, foldM)
import qualified Data.IntSet         as IS
import           Data.Maybe          (mapMaybe)
import           Data.List           (group, minimumBy, sort)
import           Data.Ord            (comparing)
import           Data.Sequence       (Seq, ViewL(..), ViewR(..), (<|), (|>))
import qualified Data.Sequence       as Seq
import qualified Data.Vector.Unboxed as VU

-- Unchecked array index (bounds guaranteed by algorithm invariants).
(!) :: VU.Unbox a => VU.Vector a -> Int -> a
(!) = VU.unsafeIndex
{-# INLINE (!) #-}
infixl 9 !

-- ════════════════════════════════════════════════════════════════
-- Oriented planar graph (CSR format)
-- ════════════════════════════════════════════════════════════════

type Node = Int

-- | Oriented triangulation in CSR (Compressed Sparse Row) format.
--   All neighbor lists concatenated into one flat unboxed array.
data Graph = Graph
  { gAdj    :: !(VU.Vector Node)  -- ^ all neighbors, concatenated
  , gOffset :: !(VU.Vector Int)   -- ^ start offset per vertex; length N+1
  }

mkGraph :: [[Node]] -> Graph
mkGraph lists = Graph (VU.fromList (concat lists))
                      (VU.fromList (scanl (+) 0 (map length lists)))

order :: Graph -> Int
order g = VU.length (gOffset g) - 1

deg :: Graph -> Node -> Int
deg g v = gOffset g ! (v+1) - gOffset g ! v

-- | Neighbor slice for vertex v (O(1), no copy).
nbrSlice :: Graph -> Node -> VU.Vector Node
nbrSlice g v = VU.unsafeSlice (gOffset g ! v) (deg g v) (gAdj g)

nbrs :: Graph -> Node -> [Node]
nbrs g v = VU.toList (nbrSlice g v)

-- | Cyclic successor of w in v's neighbor list.
next :: Graph -> Node -> Node -> Maybe Node
next g v w = do
  i <- VU.elemIndex w nv
  pure $ nv ! mod (i + 1) (VU.length nv)
  where nv = nbrSlice g v

-- | Cyclic predecessor of w in v's neighbor list.
prev :: Graph -> Node -> Node -> Maybe Node
prev g v w = do
  i <- VU.elemIndex w nv
  pure $ nv ! mod (i - 1 + n) n
  where nv = nbrSlice g v; n = VU.length nv

-- | Are u and v adjacent?
edge :: Graph -> Node -> Node -> Bool
edge g u v = VU.elem v (nbrSlice g u)

-- ════════════════════════════════════════════════════════════════
-- Orientation and triangle navigation
-- ════════════════════════════════════════════════════════════════

data Orientation = CCW | CW

-- | Detect orientation of a starting triangle.
orient :: Graph -> Node -> Node -> Node -> Maybe Orientation
orient g f1 f2 f3
  | next g f1 f2 == Just f3 = Just CCW
  | prev g f1 f2 == Just f3 = Just CW
  | otherwise                = Nothing

-- | The vertex across boundary edge (u -> w): the third vertex of the
--   adjacent triangle, found by looking backward (CCW) or forward (CW).
--   On a GPU this is a scan of <= D_max entries with no branching.
apex :: Graph -> Orientation -> Node -> Node -> Maybe Node
apex g CCW = prev g
apex g CW  = next g

-- ════════════════════════════════════════════════════════════════
-- Drain: the core boundary operation
-- ════════════════════════════════════════════════════════════════

-- | Drain from front: consume one connection, cascade-pop saturated nodes.
--   Returns (updated boundary, total connections consumed).
--   Cascade depth <= D_max - 3 (at most 3 for fullerenes).
drainFront :: Seq (Node, Int) -> (Seq (Node, Int), Int)
drainFront s = case Seq.viewl s of
  EmptyL          -> (Seq.empty, 0)
  (v, r) :< rest
    | r > 1     -> ((v, r - 1) <| rest, 1)
    | otherwise -> let (s', k) = drainFront rest in (s', k + 1)

-- | Drain from back (symmetric).
drainBack :: Seq (Node, Int) -> (Seq (Node, Int), Int)
drainBack s = case Seq.viewr s of
  EmptyR          -> (Seq.empty, 0)
  rest :> (v, r)
    | r > 1     -> (rest |> (v, r - 1), 1)
    | otherwise -> let (s', k) = drainBack rest in (s', k + 1)

-- ════════════════════════════════════════════════════════════════
-- Regular spiral extraction (fast path — no cut-vertex checks)
-- ════════════════════════════════════════════════════════════════

-- | Boundary + removed-set, threaded through the fold.
data St = St !(Seq (Node, Int)) !IS.IntSet

-- | One peel step.  Given boundary [... u] [w ...], find the apex v
--   of the triangle across edge (u,w), drain both ends, push v.
peel :: Graph -> Orientation -> St -> Maybe (St, Int)
peel g ori (St bnd removed) = do
  (w, _) <- front bnd
  (u, _) <- back  bnd
  v      <- apex g ori u w
  guard (v `IS.notMember` removed)
  let (b1, n1) = drainFront bnd
      (b2, n2) = drainBack b1
      open     = deg g v - n1 - n2
  guard (open >= 1)
  pure (St (b2 |> (v, open)) (IS.insert v removed), deg g v)

-- | Extract a regular (jump-free) spiral from a starting triple.
regularSpiral :: Graph -> (Node, Node, Node) -> Maybe [Int]
regularSpiral g (f1, f2, f3) = do
  ori <- orient g f1 f2 f3
  let n  = order g
      s0 = St (Seq.fromList [(f1, deg g f1 - 2), (f2, deg g f2 - 2), (f3, deg g f3 - 2)])
              (IS.fromList [f1, f2, f3])

  (St bnd removed, revDs) <-
    foldM (\(st, acc) _ -> do (st', dv) <- peel g ori st
                              pure (st', dv : acc))
          (s0, []) [3 .. n - 2]

  let lastV = IS.findMin (IS.fromDistinctAscList [0..n-1] `IS.difference` removed)
  guard (Seq.length bnd == deg g lastV)
  guard (all ((== 1) . snd) bnd)
  pure $ map (deg g) [f1, f2, f3] ++ reverse revDs ++ [deg g lastV]

-- ════════════════════════════════════════════════════════════════
-- Generalized spiral extraction (handles cut vertices via jumps)
-- ════════════════════════════════════════════════════════════════

data GeneralSpiral = GeneralSpiral
  { gsJumps  :: ![(Int, Int)]   -- ^ (position, jump_length) pairs
  , gsSpiral :: ![Int]          -- ^ face-degree sequence
  } deriving (Eq, Show)

-- | Canonical ordering: fewest jumps > smallest jump list > smallest spiral.
instance Ord GeneralSpiral where
  compare (GeneralSpiral j1 s1) (GeneralSpiral j2 s2) =
    compare (length j1) (length j2) <>
    compare j1 j2 <>
    compare s1 s2

data GSt = GSt
  { gstBnd     :: !(Seq (Node, Int))
  , gstRemoved :: !IS.IntSet
  , gstJumpSt  :: !Int             -- current consecutive rotation count
  , gstJumps   :: ![(Int, Int)]    -- accumulated jumps (reversed)
  }

-- | Is v a cut vertex of the remaining (un-removed) graph?
--   Count edges between consecutive remaining neighbors in cyclic order;
--   v is a cut vertex iff this count < (number of remaining neighbors) - 1.
isCutVertex :: Graph -> IS.IntSet -> Node -> Bool
isCutVertex g removed v =
  let rn = [w | w <- nbrs g v, w `IS.notMember` removed]
      n  = length rn
      ce = length [() | (a, b) <- zip rn (tail rn ++ [head rn]), edge g a b]
  in n >= 2 && ce < n - 1

-- | One peel step with jump handling.  If the apex is a cut vertex,
--   rotate the boundary and retry instead of peeling.
generalPeelStep :: Graph -> Orientation
                -> Int        -- step index (for recording jumps)
                -> GSt -> Maybe (GSt, Int)
generalPeelStep g ori stepIdx = tryPeel 0
  where
    tryPeel !rotations (GSt bnd removed jumpSt jumpLog) = do
      guard (rotations < Seq.length bnd)  -- stop after full boundary rotation
      (w, _) <- front bnd
      (u, _) <- back  bnd
      v      <- apex g ori u w
      guard (v `IS.notMember` removed)

      if isCutVertex g removed v
        then tryPeel (rotations + 1) (GSt (rotateFwd bnd) removed (jumpSt + 1) jumpLog)
        else do
          let jumpLog' = if jumpSt > 0
                then (stepIdx, jumpSt) : jumpLog
                else jumpLog
              (b1, n1) = drainFront bnd
              (b2, n2) = drainBack b1
              open     = deg g v - n1 - n2
          guard (open >= 1)
          pure (GSt (b2 |> (v, open)) (IS.insert v removed) 0 jumpLog', deg g v)

-- | Rotate a deque: move front to back.
rotateFwd :: Seq a -> Seq a
rotateFwd s = case Seq.viewl s of
  EmptyL    -> Seq.empty
  a :< rest -> rest |> a

-- | Extract a generalized spiral from a starting triple.
generalSpiral :: Graph -> (Node, Node, Node) -> Maybe GeneralSpiral
generalSpiral g (f1, f2, f3) = do
  ori <- orient g f1 f2 f3
  let n  = order g
      s0 = GSt (Seq.fromList [(f1, deg g f1 - 2), (f2, deg g f2 - 2), (f3, deg g f3 - 2)])
               (IS.fromList [f1, f2, f3]) 0 []

  (GSt bnd removed _ jumpLog, revDs) <-
    foldM (\(st, acc) i -> do (st', dv) <- generalPeelStep g ori i st
                              pure (st', dv : acc))
          (s0, []) [3 .. n - 2]

  let lastV = IS.findMin (IS.fromDistinctAscList [0..n-1] `IS.difference` removed)
  guard (Seq.length bnd == deg g lastV)
  guard (all ((== 1) . snd) bnd)
  pure $ GeneralSpiral (reverse jumpLog)
           (map (deg g) [f1, f2, f3] ++ reverse revDs ++ [deg g lastV])

-- ════════════════════════════════════════════════════════════════
-- Canonical spirals
-- ════════════════════════════════════════════════════════════════

-- | Canonical regular spiral, or Nothing if every starting triple needs jumps.
canonicalSpiral :: Graph -> Maybe [Int]
canonicalSpiral g = case mapMaybe (regularSpiral g) (startingTriples g) of
  []      -> Nothing
  spirals -> Just (minimum spirals)

-- | Canonical generalized spiral.  Tries regular spirals first (fast path);
--   falls back to generalized only if all regular attempts fail.
canonicalGeneralSpiral :: Graph -> Maybe GeneralSpiral
canonicalGeneralSpiral g =
  case mapMaybe (regularSpiral g) triples of
    spirals@(_:_) -> Just (minimum [GeneralSpiral [] sp | sp <- spirals])
    [] -> case mapMaybe (generalSpiral g) triples of
            []      -> Nothing
            spirals -> Just (minimum spirals)
  where triples = startingTriples g

-- | All valid starting triples, restricted to the rarest special face type.
startingTriples :: Graph -> [(Node, Node, Node)]
startingTriples g = do
  f1 <- rarestSpecial g
  f2 <- nbrs g f1
  f3 <- mapMaybe id [prev g f2 f1, next g f2 f1]
  pure (f1, f2, f3)

-- | Vertices with the rarest non-hexagonal degree (the 12 pentagons
--   in a fullerene dual).  Falls back to all vertices if every vertex
--   has degree 6.
rarestSpecial :: Graph -> [Node]
rarestSpecial g
  | null specials = [0 .. order g - 1]
  | otherwise     = [v | v <- [0 .. order g - 1], deg g v == target]
  where
    specials = [deg g v | v <- [0 .. order g - 1], deg g v /= 6]
    target   = head . minimumBy (comparing length) . group . sort $ specials

-- ════════════════════════════════════════════════════════════════
-- Helpers
-- ════════════════════════════════════════════════════════════════

front :: Seq a -> Maybe a
front s = case Seq.viewl s of { a :< _ -> Just a; _ -> Nothing }

back :: Seq a -> Maybe a
back s = case Seq.viewr s of { _ :> a -> Just a; _ -> Nothing }
