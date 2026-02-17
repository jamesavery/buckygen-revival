-- | Spiral windup: construct an oriented triangulation from a face-degree sequence.
--
-- This is the inverse of spiral extraction: given a sequence of face degrees
-- (and optionally a jump list), build the oriented adjacency lists of the
-- corresponding triangulation.

module Windup
  ( windupSpiral
  , windupGeneralSpiral
  , fromRSPI
  ) where

import           Data.List             (foldl')
import qualified Data.IntMap.Strict    as IM
import qualified Data.IntSet           as IS
import           Data.Sequence         (Seq, ViewL(..), ViewR(..), (|>))
import qualified Data.Sequence         as Seq
import qualified Data.Vector.Unboxed   as VU
import           Spiral                (Graph, Node, mkGraph, rotateFwd)

-- Unchecked array index (bounds guaranteed by algorithm invariants).
(!) :: VU.Unbox a => VU.Vector a -> Int -> a
(!) = VU.unsafeIndex
{-# INLINE (!) #-}
infixl 9 !

type Nbrs = IM.IntMap [Node]
type OV   = Seq (Node, Int)

-- | Construct an oriented triangulation from a face-degree sequence (no jumps).
windupSpiral :: [Int] -> Graph
windupSpiral spiral = windupGeneralSpiral spiral []

-- | Construct an oriented triangulation from a face-degree sequence with jumps.
windupGeneralSpiral :: [Int] -> [(Int, Int)] -> Graph
windupGeneralSpiral spiral jumps =
  let n    = length spiral
      sp   = VU.fromList spiral
      jm   = IM.fromList jumps
      nbs0 = IM.fromList [(i, []) | i <- [0..n-1]]
      nbs1 = insEdge 0 1 (-1) (-1) nbs0
      ov0  = Seq.fromList [(0, sp ! 0 - 1), (1, sp ! 1 - 1)]
      (nbsF, ovF) = foldl' (\(nb, ov) k ->
                       let ov' = applyJump jm k ov
                       in stepK sp k nb ov')
                     (nbs1, ov0) [2..n-2]
      nbsL = closeLast (n-1) nbsF ovF
  in mkGraph [IM.findWithDefault [] i nbsL | i <- [0..n-1]]

-- | Construct a fullerene dual triangulation from 1-indexed RSPI and atom count.
--   E.g. @fromRSPI 20 [1..12]@ builds the icosahedron (C20 dual).
fromRSPI :: Int -> [Int] -> Graph
fromRSPI nAtoms rspi1 = windupSpiral spiral
  where
    nFaces  = nAtoms `div` 2 + 2
    pentSet = IS.fromList (map (subtract 1) rspi1)
    spiral  = [if IS.member i pentSet then 5 else 6 | i <- [0..nFaces-1]]

-- ════════════════════════════════════════════════════════════════
-- Windup internals
-- ════════════════════════════════════════════════════════════════

-- | Apply a jump at step k if one exists: rotate boundary by jump_length.
applyJump :: IM.IntMap Int -> Int -> OV -> OV
applyJump jm k ov = case IM.lookup k jm of
  Nothing -> ov
  Just n  -> iterate rotateFwd ov !! n

-- | Insert edge (u,v) into the planar embedding: place v in u's neighbor
--   list before sucUV, and u in v's list before sucVU.
insEdge :: Node -> Node -> Node -> Node -> Nbrs -> Nbrs
insEdge u v sucUV sucVU =
  IM.adjust (insBefore v sucUV) u .
  IM.adjust (insBefore u sucVU) v

insBefore :: Eq a => a -> a -> [a] -> [a]
insBefore new _      []     = [new]
insBefore new before (x:xs)
  | x == before = new : x : xs
  | otherwise   = x : insBefore new before xs

-- | One step of the windup: add face k to the partially-built triangulation.
--   Connect k to front and back of boundary, cascade-drain saturated nodes,
--   push k with remaining open valencies.
stepK :: VU.Vector Int -> Node -> Nbrs -> OV -> (Nbrs, OV)
stepK sp k nbs0 ov0 =
  let fv = fst (seqHead ov0)
      bv = fst (seqLast ov0)

      -- Connect k to front (insert before back in both lists)
      nbs1 = insEdge k fv bv bv nbs0
      ov1  = decFront ov0

      -- Connect k to back (insert before front / second-to-last)
      bv'  = fst (seqLast ov1)
      s2l  = fst (sec2last ov1)
      nbs2 = insEdge k bv' fv s2l nbs1
      ov2  = decBack ov1

      -- Cascade drain (2 connections already made)
      (nbs3, ov3, pu3) = cascFwd k nbs2 ov2 2
      (nbs4, ov4, pu4) = cascBwd k nbs3 ov3 pu3

  in (nbs4, ov4 |> (k, sp ! k - pu4))

-- | Cascade forward: while front is saturated, pop it and connect k to new front.
--   Stops when boundary has <= 1 element, or when both ends are saturated
--   with only 2 elements left (to avoid cascading past the back end).
cascFwd :: Node -> Nbrs -> OV -> Int -> (Nbrs, OV, Int)
cascFwd k nbs ov pu
  | Seq.length ov <= 1       = (nbs, ov, pu)
  | snd (seqHead ov) /= 0   = (nbs, ov, pu)
  | Seq.length ov == 2
    && snd (seqLast ov) == 0 = (nbs, ov, pu)  -- both ends saturated, don't cross
  | otherwise =
      let old  = fst (seqHead ov)
          ov'  = Seq.drop 1 ov
          nbs' = insEdge k (fst (seqHead ov')) (fst (seqLast ov')) old nbs
      in cascFwd k nbs' (decFront ov') (pu + 1)

-- | Cascade backward: while back is saturated, pop it and connect k to new back.
--   Symmetric stopping condition to cascFwd.
cascBwd :: Node -> Nbrs -> OV -> Int -> (Nbrs, OV, Int)
cascBwd k nbs ov pu
  | Seq.length ov <= 1       = (nbs, ov, pu)
  | snd (seqLast ov) /= 0   = (nbs, ov, pu)
  | Seq.length ov == 2
    && snd (seqHead ov) == 0 = (nbs, ov, pu)  -- both ends saturated, don't cross
  | otherwise =
      let old  = fst (seqLast ov)
          ov'  = Seq.take (Seq.length ov - 1) ov
          nbs' = insEdge k (fst (seqLast ov')) old (fst (sec2last ov')) nbs
      in cascBwd k nbs' (decBack ov') (pu + 1)

-- | Close the last face: connect it to all remaining boundary nodes.
closeLast :: Node -> Nbrs -> OV -> Nbrs
closeLast lastN nbs ov =
  let nodes = fmap fst ov
      len   = Seq.length nodes
  in foldl' (\nb i ->
       let ni   = Seq.index nodes i
           succ = Seq.index nodes ((i + 1) `mod` len)
           pred = Seq.index nodes ((i - 1 + len) `mod` len)
       in insEdge lastN ni succ pred nb
     ) nbs [0..len-1]

-- ════════════════════════════════════════════════════════════════
-- Seq helpers
-- ════════════════════════════════════════════════════════════════

seqHead :: Seq a -> a
seqHead s = case Seq.viewl s of { a Seq.:< _ -> a; _ -> error "seqHead: empty" }

seqLast :: Seq a -> a
seqLast s = case Seq.viewr s of { _ Seq.:> a -> a; _ -> error "seqLast: empty" }

sec2last :: Seq a -> a
sec2last s = Seq.index s (Seq.length s - 2)

decFront :: Seq (Node, Int) -> Seq (Node, Int)
decFront = Seq.adjust (\(v, r) -> (v, r - 1)) 0

decBack :: Seq (Node, Int) -> Seq (Node, Int)
decBack ov = Seq.adjust (\(v, r) -> (v, r - 1)) (Seq.length ov - 1) ov
