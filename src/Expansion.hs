{-# LANGUAGE StrictData #-}
module Expansion
    ( -- * Types
      Vertex, Edge, Dir(..), ExpKind(..), Expansion(..), Reduction(..)
    , PathInfo(..)
    , inverse
      -- * Graph primitives
    , nbrAt, nbrs, deg
    , nextCW, prevCW, advanceCW, sideNbr, straightAhead, turnAhead
    , insertAfter, replaceNbr, removeNbr
      -- * Path computation
    , computeStraightPath, computeBentZeroPath, computeBentPath
    , computeBentPathSafe, canBentPath, isValidStraightSite
      -- * Expansion / reduction
    , applyExpansion, applyReduction
      -- * Enumeration
    , expansions, expansionsL0, expansionsL, expansionsB
      -- * Nanotube ring (F) expansion / reduction
    , findNanotubeRing, applyRing, reduceRing
      -- * Reduction helpers
    , reductionLength, longestStraight, numNewVertices
      -- * Direction
    , flipDir
    ) where

import Seeds (DualGraph(..), EdgeList, initEdgeList, mkDualGraph, mkDualGraphLite)
import qualified Data.Array as A
import Data.Array.Unboxed (UArray)
import qualified Data.Array.Unboxed as UA
import Data.IntMap.Strict (IntMap)
import qualified Data.IntMap.Strict as IM
import qualified Data.IntSet as IS
import Data.List (sort, elemIndex, foldl', intersect, nub)
import Data.Maybe (fromJust)

---------------------------------------------------------------------
-- Types
---------------------------------------------------------------------

type Vertex = Int
type Edge   = (Vertex, Vertex)          -- directed: (start, end)

data Dir = DLeft | DRight deriving (Eq, Ord, Show)

data ExpKind
    = L Int           -- L_i: straight path, i >= 0
    | B Int Int       -- B_{i,j}: bent path, i+j >= 0
    | F               -- Nanotube ring expansion (adds 5 degree-6 vertices)
    deriving (Eq, Ord, Show)

data Expansion = Exp
    { expKind   :: ExpKind
    , startEdge :: Edge
    , direction :: Dir
    } deriving (Eq, Show)

data Reduction = Red
    { redKind :: ExpKind
    , redEdge :: Edge
    , redDir  :: Dir
    } deriving (Eq, Show)

-- | Inverse relationship: same triple, different semantic role
inverse :: Expansion -> Reduction
inverse (Exp k e d) = Red k e d

-- | Path and parallel path arrays computed from the graph.
-- For straight L_i: mainPath has i+3 entries, parallelPath has i+3 entries.
-- For B_{0,0}: mainPath has 5 entries, parallelPath has 3 entries.
-- For B_{i,j} (i+j>0): mainPath has i+j+5 entries, parallelPath has i+j+3 entries.
data PathInfo = PathInfo
    { mainPath     :: [Vertex]
    , parallelPath :: [Vertex]
    } deriving (Eq, Show)

---------------------------------------------------------------------
-- Graph primitives: cyclic neighbor list operations
---------------------------------------------------------------------

-- | O(1) unboxed access to the i-th CW neighbor of vertex u.
nbrAt :: DualGraph -> Vertex -> Int -> Vertex
nbrAt g u i = adjFlat g UA.! (u * 6 + i)
{-# INLINE nbrAt #-}

-- | Get the neighbor list of vertex u (O(1) via pre-built boxed Array).
nbrs :: DualGraph -> Vertex -> [Vertex]
nbrs g u = adjArray g A.! u

-- | Degree of vertex u (O(1) unboxed lookup).
deg :: DualGraph -> Vertex -> Int
deg g u = degFlat g UA.! u
{-# INLINE deg #-}

-- | Find index of v in u's CW neighbor list.
-- Linear scan over 5-6 contiguous unboxed Ints (faster than IntMap).
indexOf :: DualGraph -> Vertex -> Vertex -> Int
indexOf g u v = go 0
  where
    base = u * 6
    d = deg g u
    flat = adjFlat g
    go i | i >= d    = error $ "indexOf: " ++ show v ++ " not in " ++ show u
                             ++ "'s neighbors " ++ show (nbrs g u)
         | flat UA.! (base + i) == v = i
         | otherwise = go (i + 1)
{-# INLINE indexOf #-}

-- | Next neighbor clockwise from v in u's neighbor list.
nextCW :: DualGraph -> Vertex -> Vertex -> Vertex
nextCW g u v = nbrAt g u ((indexOf g u v + 1) `mod` deg g u)
{-# INLINE nextCW #-}

-- | Previous neighbor (counter-clockwise) from v in u's neighbor list.
prevCW :: DualGraph -> Vertex -> Vertex -> Vertex
prevCW g u v = let d = deg g u in nbrAt g u ((indexOf g u v - 1 + d) `mod` d)
{-# INLINE prevCW #-}

-- | Advance n positions clockwise from v in u's neighbor list.
advanceCW :: DualGraph -> Vertex -> Vertex -> Int -> Vertex
advanceCW g u v n' = nbrAt g u ((indexOf g u v + n') `mod` deg g u)
{-# INLINE advanceCW #-}

-- | Side neighbor: for DRight use prevCW, for DLeft use nextCW
-- This gives the parallel-path vertex when looking from 'from' toward 'to'.
sideNbr :: DualGraph -> Dir -> Vertex -> Vertex -> Vertex
sideNbr g DRight from to = prevCW g from to
sideNbr g DLeft  from to = nextCW g from to

-- | Straight-ahead advance: 3 CW for DRight (next^3), (deg-3) CW for DLeft (prev^3).
-- At degree-6 these are equivalent (3 = 6-3). At degree-5 they differ (3 vs 2).
straightAhead :: DualGraph -> Dir -> Vertex -> Vertex -> Vertex
straightAhead g dir u from = advanceCW g u from n
  where n = case dir of { DRight -> 3; DLeft -> deg g u - 3 }

-- | Turn advance for bent paths: 2 CW for DRight (next^2), (deg-2) CW for DLeft (prev^2).
turnAhead :: DualGraph -> Dir -> Vertex -> Vertex -> Vertex
turnAhead g dir u from = advanceCW g u from n
  where n = case dir of { DRight -> 2; DLeft -> deg g u - 2 }

---------------------------------------------------------------------
-- IntMap-level neighbor list modification
---------------------------------------------------------------------

-- | Insert 'new' clockwise after 'after' in vertex u's neighbor list
insertAfter :: IntMap [Int] -> Vertex -> Vertex -> Vertex -> IntMap [Int]
insertAfter adj u after new = IM.adjust go u adj
  where
    go ns = let (before, rest) = span (/= after) ns
            in case rest of
                 (a:as) -> before ++ [a, new] ++ as
                 []     -> error $ "insertAfter: " ++ show after
                                ++ " not in " ++ show u ++ "'s list"

-- | Replace 'old' with 'new' in vertex u's neighbor list
replaceNbr :: IntMap [Int] -> Vertex -> Vertex -> Vertex -> IntMap [Int]
replaceNbr adj u old new = IM.adjust (map (\x -> if x == old then new else x)) u adj

-- | Remove v from vertex u's neighbor list
removeNbr :: IntMap [Int] -> Vertex -> Vertex -> IntMap [Int]
removeNbr adj u v = IM.adjust (filter (/= v)) u adj

-- | Flip direction
flipDir :: Dir -> Dir
flipDir DLeft  = DRight
flipDir DRight = DLeft

-- | CW-position-based lookup on raw IntMap (for use during graph surgery)
nextCW' :: IntMap [Int] -> Vertex -> Vertex -> Vertex
nextCW' adj u v =
    let ns = adj IM.! u
        n  = length ns
        i  = case elemIndex v ns of
               Just idx -> idx
               Nothing  -> error $ "nextCW': " ++ show v ++ " not in "
                                ++ show u ++ "'s list " ++ show ns
    in ns !! ((i + 1) `mod` n)

prevCW' :: IntMap [Int] -> Vertex -> Vertex -> Vertex
prevCW' adj u v =
    let ns = adj IM.! u
        n  = length ns
        i  = case elemIndex v ns of
               Just idx -> idx
               Nothing  -> error $ "prevCW': " ++ show v ++ " not in "
                                ++ show u ++ "'s list " ++ show ns
    in ns !! ((i - 1 + n) `mod` n)

sideNbr' :: IntMap [Int] -> Dir -> Vertex -> Vertex -> Vertex
sideNbr' adj DRight from to = prevCW' adj from to
sideNbr' adj DLeft  from to = nextCW' adj from to

---------------------------------------------------------------------
-- Edge list: slot-based replacement (matching C code's edge_list)
-- EdgeList and initEdgeList are imported from Seeds.
---------------------------------------------------------------------

-- | Replace by slot position, not by value search.
-- Looks up the position of orgNbr in vertex's list (via edge_list),
-- replaces whatever is currently at that position with newNbr,
-- and registers newNbr at that position in the edge_list.
replaceNbrEL :: IntMap [Int] -> EdgeList -> Vertex -> Vertex -> Vertex
             -> (IntMap [Int], EdgeList)
replaceNbrEL adj el vertex orgNbr newNbr =
    let elV = el IM.! vertex
        pos = case IM.lookup orgNbr elV of
                Just p  -> p
                Nothing -> error $ "replaceNbrEL: " ++ show orgNbr
                    ++ " not in edge_list for vertex " ++ show vertex
                    ++ " (edge_list keys: " ++ show (IM.keys elV) ++ ")"
        ns  = adj IM.! vertex
        ns' = take pos ns ++ [newNbr] ++ drop (pos + 1) ns
        elV' = IM.insert newNbr pos elV
    in (IM.insert vertex ns' adj, IM.insert vertex elV' el)

-- | Register a new vertex and its neighbors in the edge list.
registerVertex :: EdgeList -> Vertex -> [Vertex] -> EdgeList
registerVertex el v ns = IM.insert v (IM.fromList (zip ns [0..])) el

-- | Insert 'new' clockwise after 'after' in vertex u's neighbor list,
-- and update the edge list accordingly (shift positions, register new).
-- Uses edge_list for position lookup (slot-based, not value-based).
insertAfterEL :: IntMap [Int] -> EdgeList -> Vertex -> Vertex -> Vertex
              -> (IntMap [Int], EdgeList)
insertAfterEL adj el u after new =
    let ns = adj IM.! u
        elU = case IM.lookup u el of
                Just m  -> m
                Nothing -> IM.fromList (zip ns [0..])
        afterPos = case IM.lookup after elU of
                Just p  -> p
                Nothing -> error $ "insertAfterEL: " ++ show after
                               ++ " not in " ++ show u ++ "'s edge_list"
                               ++ " (keys: " ++ show (IM.keys elU) ++ ")"
        pos = afterPos + 1  -- insert position (after 'after')
        ns' = take pos ns ++ [new] ++ drop pos ns
        elU' = IM.map (\p -> if p >= pos then p + 1 else p) elU
        elU'' = IM.insert new pos elU'
    in (IM.insert u ns' adj, IM.insert u elU'' el)


---------------------------------------------------------------------
-- Path computation
---------------------------------------------------------------------

-- | Compute path and parallel path for a straight expansion (L_i).
-- numEntries = i+3 (the number of vertices in path and parallelPath).
-- For L_0: numEntries = 3
-- For L_1: numEntries = 4
-- For L_i: numEntries = i+3
computeStraightPath :: DualGraph -> Edge -> Dir -> Int -> PathInfo
computeStraightPath g (u0, v0) dir numEntries = PathInfo path par
  where
    -- Generate directed edges along the strip.
    -- edge[k] = (from, to): from = path[k], to = direction toward path[k+1].
    edges = take numEntries $ iterate advance (u0, v0)

    advance (from, to) = (to, straightAhead g dir to from)

    path = map fst edges
    par  = map (\(f, t) -> sideNbr g dir f t) edges

-- | Compute path and parallel path for a B_{0,0} expansion.
-- path has 5 entries, parallelPath has 3 entries.
-- The bent path goes: straight step, turn (advance 2), straight step.
computeBentZeroPath :: DualGraph -> Edge -> Dir -> PathInfo
computeBentZeroPath g (u0, v0) dir = PathInfo path par
  where
    -- Phase 1: one straight step (producing 2 path entries)
    -- edge0 = (u0, v0), edge1 = (v0, straightAhead)
    e0 = (u0, v0)
    (_, v1) = e0
    e1 = (v1, straightAhead g dir v1 u0)

    -- Phase 2: turn at e1's target
    (from2, turnV) = e1
    afterTurn = turnAhead g dir turnV from2
    e2 = (turnV, afterTurn)

    -- Phase 3: one straight step after turn
    (from3, v3) = e2
    e3 = (v3, straightAhead g dir v3 from3)

    -- Phase 4: final endpoint
    (from4, v4) = e3

    -- path: [u0, v1, turnV, v3, v4]
    path = [u0, v1, turnV, v3, v4]

    -- Parallel path: 3 entries.  The turn edge is SKIPPED for parallel path
    -- (matching the C code's extend_bent_zero convention).
    -- par[0] = sideNbr at edge path[0]->path[1] (before turn)
    -- par[1] = sideNbr at edge path[2]->path[3] (after turn)
    -- par[2] = sideNbr at edge path[3]->path[4] (after turn)
    par = [ sideNbr g dir (path !! 0) (path !! 1)
          , sideNbr g dir (path !! 2) (path !! 3)
          , sideNbr g dir (path !! 3) (path !! 4)
          ]

-- | Compute path and parallel path for a B_{i,j} expansion with i+j > 0.
-- bentPos = i, bentLen = i+j.
-- path has bentLen+5 entries, parallelPath has bentLen+3 entries.
computeBentPath :: DualGraph -> Edge -> Dir -> Int -> Int -> PathInfo
computeBentPath g (u0, v0) dir bentPos bentLen = PathInfo pathVerts parVerts
  where
    -- Straight step advance
    adv (from, to) = (to, straightAhead g dir to from)

    -- Phase 1: straight steps producing bentPos+2 path entries
    preBendEdges = take (bentPos + 2) $ iterate adv (u0, v0)

    -- Path[0..bentPos+1]
    prePath = map fst preBendEdges
    -- Par[0..bentPos]: side of edges path[k]->path[k+1] for k=0..bentPos
    prePar = map (\(f, t) -> sideNbr g dir f t) (take (bentPos + 1) preBendEdges)

    -- At the end of pre-bend: the edge pointing into the turn vertex
    (preLast, turnV) = last preBendEdges

    -- Phase 2: the turn
    afterTurn = turnAhead g dir turnV preLast
    -- turnV = path[bentPos+1] is already in prePath
    -- path[bentPos+2] = turnV (wait -- turnV was the TARGET of the last edge,
    --   but we need the next vertex AFTER the turn)

    -- Actually: the last edge in preBendEdges is (path[bentPos], path[bentPos+1]).
    -- wait, no: preBendEdges has bentPos+2 entries. Let me trace carefully.
    --   preBendEdges[0] = (u0, v0)          → path[0] = u0
    --   preBendEdges[1] = (v0, next)         → path[1] = v0
    --   ...
    --   preBendEdges[bentPos+1] = (prev, X)  → path[bentPos+1] = prev
    -- So the TARGET of the last edge = X = the direction beyond path[bentPos+1]
    -- But path[bentPos+1] = fst of... hmm, no:
    -- preBendEdges[k] = (path[k], direction_to_path[k+1])
    -- So fst of preBendEdges[bentPos+1] = path[bentPos+1]
    -- And snd of preBendEdges[bentPos+1] = direction toward path[bentPos+2]
    -- This "direction" vertex IS the next path vertex in a straight walk.
    -- But we're about to turn, so path[bentPos+2] is NOT that vertex.

    -- Let me re-examine. The turn happens AFTER bentPos+1 straight steps.
    -- After phase 1, work_edge in C points from path[bentPos+1] toward
    -- the next straight-ahead vertex. Then:
    --   path[bentPos+2] = work_edge->end (the straight-ahead vertex)
    --   work_edge = work_edge->invers
    --   work_edge = work_edge->next->next (turn)
    -- This means path[bentPos+2] is the vertex BEFORE the turn,
    -- and the turn happens at that vertex.

    -- Let me re-trace the C code for the turn phase:
    --   After phase 1 loop, work_edge = (path[bentPos+1], straightAhead)
    --   where straightAhead was set by the last advance in the loop.
    --
    --   Turn:
    --     path_bent[length++] = work_edge->end;  → path[bentPos+2] = straightAhead
    --     work_edge = work_edge->invers;          → (straightAhead, path[bentPos+1])
    --     work_edge = work_edge->next->next;      → advance 2 CW from path[bentPos+1]
    --                                               at vertex straightAhead

    -- So the turn vertex is straightAhead = snd(last preBendEdges) = turnV
    -- And path[bentPos+2] = turnV

    -- After the turn, work_edge points from turnV toward the turned direction.

    turnEdge = (turnV, afterTurn)

    -- Phase 3: straight steps after the turn
    -- (bentLen - bentPos + 1) more path entries
    postBendEdges = take (bentLen - bentPos + 1) $ iterate adv turnEdge

    postPath = map fst postBendEdges
    postPar  = map (\(f, t) -> sideNbr g dir f t) postBendEdges

    -- Phase 4: final endpoint
    (lastFrom, lastTo) = last postBendEdges
    finalVertex = straightAhead g dir lastTo lastFrom
    -- Final parallel: sideNbr at (lastTo -> finalVertex)
    -- But C code uses work_edge->prev->end (for use_next) where work_edge = (lastTo, finalVertex)
    -- This is prevCW(lastTo, finalVertex) for DRight = sideNbr DRight lastTo finalVertex
    -- But wait: the C code for the last parallel uses the OPPOSITE side function!
    -- For use_next: work_edge->prev->end (prev, not next)
    -- For !use_next: work_edge->next->end (next, not prev)
    -- This is the SAME as sideNbr but at a different vertex pairing.
    -- Actually: sideNbr DRight lastTo finalVertex = prevCW lastTo finalVertex
    -- Which IS what the C code computes. Let me verify:
    -- C code: work_edge = (lastTo, finalVertex), work_edge->prev->end
    --   = prevCW of finalVertex in lastTo's list? No:
    --   work_edge->prev->end: work_edge starts at lastTo pointing to finalVertex.
    --   prev of that edge is the edge just before in CW order at lastTo.
    --   prev->end = the vertex just before finalVertex in lastTo's CW list
    --   = prevCW(lastTo, finalVertex)
    -- And sideNbr DRight lastTo finalVertex = prevCW g lastTo finalVertex. ✓

    finalPar = sideNbr g dir lastTo finalVertex

    -- Assemble path
    pathVerts = prePath ++ postPath ++ [lastTo, finalVertex]

    -- Parallel path: the turn edge is SKIPPED (matching C code convention).
    -- par[0..bentPos] = sideNbr at pre-bend edges
    -- par[bentPos+1..bentLen+2] = sideNbr at post-bend edges (shifted +1 in path)
    parVerts = prePar ++ postPar ++ [finalPar]

---------------------------------------------------------------------
-- Expansion operations
---------------------------------------------------------------------

-- | Apply an expansion, returning the expanded graph and the path info used.
applyExpansion :: Expansion -> DualGraph -> (DualGraph, PathInfo)
applyExpansion (Exp (L 0) edge dir) g =
    let pi = computeStraightPath g edge dir 3
    in (applyL0 g pi dir, pi)
applyExpansion (Exp (L i) edge dir) g =
    let pi = computeStraightPath g edge dir (i + 3)
    in (applyStraight g pi dir (i + 2), pi)
applyExpansion (Exp (B 0 0) edge dir) g =
    let pi = computeBentZeroPath g edge dir
    in (applyBentZero g pi dir, pi)
applyExpansion (Exp (B i j) edge dir) g =
    let pi = computeBentPath g edge dir i (i + j)
    in (applyBent g pi dir i (i + j), pi)

-- | Apply a reduction (inverse of expansion) using pre-computed path info.
applyReduction :: Expansion -> PathInfo -> DualGraph -> DualGraph
applyReduction (Exp (L 0)   _ _) pi g  = reduceL0 g pi
applyReduction (Exp (L _)   _ _) pi g = reduceStraight g pi
applyReduction (Exp (B 0 0) _ _) pi g = reduceBentZero g pi
applyReduction (Exp (B i j) _ _) pi g = reduceBent g pi i (i + j)

---------------------------------------------------------------------
-- L0 expansion: adds 2 degree-5 vertices
---------------------------------------------------------------------

applyL0 :: DualGraph -> PathInfo -> Dir -> DualGraph
applyL0 g (PathInfo path par) dir = mkDualGraph nv' adj' deg5' el'
  where
    nv = numVertices g
    a  = nv          -- first new vertex (degree-5)
    b  = nv + 1      -- second new vertex (degree-5)
    nv' = nv + 2

    p0 = path !! 0;  p1 = path !! 1;  p2 = path !! 2
    q0 = par  !! 0;  q1 = par  !! 1;  q2 = par  !! 2

    aNbrs = case dir of
        DRight -> [p0, q0, q1, b, p1]
        DLeft  -> [p0, p1, b, q1, q0]
    bNbrs = case dir of
        DRight -> [a, q1, q2, p2, p1]
        DLeft  -> [a, p1, p2, q2, q1]

    -- Start from persistent edge list
    adj0 = neighbours g
    el_  = edgeList g

    -- Add new vertices and register them
    adj1 = IM.insert a aNbrs $ IM.insert b bNbrs adj0
    el0  = registerVertex (registerVertex el_ a aNbrs) b bNbrs

    -- p0: insert a after the appropriate neighbor (becomes degree-6)
    (adj2, el1) = case dir of
        DRight -> insertAfterEL adj1 el0 p0 q0 a
        DLeft  -> insertAfterEL adj1 el0 p0 p1 a

    -- q2: insert b after the appropriate neighbor (becomes degree-6)
    (adj3, el2) = case dir of
        DRight -> insertAfterEL adj2 el1 q2 p2 b
        DLeft  -> insertAfterEL adj2 el1 q2 (prevCW g q2 p2) b

    -- p1: replace q0 with a, replace q1 with b
    (adj4, el3) = let (a1, e1) = replaceNbrEL adj3 el2 p1 q0 a
                  in  replaceNbrEL a1 e1 p1 q1 b

    -- p2: replace q1 with b
    (adj5, el4) = replaceNbrEL adj4 el3 p2 q1 b

    -- q0: replace p1 with a
    (adj6, el5) = replaceNbrEL adj5 el4 q0 p1 a

    -- q1: replace p1 with a, replace p2 with b
    (adj', el') = let (a1, e1) = replaceNbrEL adj6 el5 q1 p1 a
                  in  replaceNbrEL a1 e1 q1 p2 b

    -- Degree-5 list: remove p0 and q2 (now degree-6), add a and b
    deg5' = sort $ a : b : filter (\v -> v /= p0 && v /= q2) (degree5 g)

---------------------------------------------------------------------
-- L0 reduction: removes 2 vertices
---------------------------------------------------------------------

reduceL0 :: DualGraph -> PathInfo -> DualGraph
reduceL0 g (PathInfo path par) = mkDualGraphLite nv' adj' deg5'
  where
    nv = numVertices g
    a  = nv - 2       -- first added vertex
    b  = nv - 1       -- second added vertex
    nv' = nv - 2

    p0 = path !! 0;  p1 = path !! 1;  p2 = path !! 2
    q0 = par  !! 0;  q1 = par  !! 1;  q2 = par  !! 2

    adj0 = neighbours g

    -- Remove new vertices
    adj1 = IM.delete a $ IM.delete b adj0

    -- p0: remove a (restore to degree-5)
    adj2 = removeNbr adj1 p0 a

    -- p1: replace a with q0, replace b with q1
    adj3 = replaceNbr (replaceNbr adj2 p1 a q0) p1 b q1

    -- p2: replace b with q1
    adj4 = replaceNbr adj3 p2 b q1

    -- q0: replace a with p1
    adj5 = replaceNbr adj4 q0 a p1

    -- q1: replace a with p1, replace b with p2
    adj6 = replaceNbr (replaceNbr adj5 q1 a p1) q1 b p2

    -- q2: remove b (restore to degree-5)
    adj' = removeNbr adj6 q2 b

    -- Degree-5 list: remove a and b, add back p0 and q2
    deg5' = sort $ p0 : q2 : filter (\v -> v /= a && v /= b) (degree5 g)

---------------------------------------------------------------------
-- Straight expansion (L_i, i >= 1): adds pathlength new vertices
---------------------------------------------------------------------

applyStraight :: DualGraph -> PathInfo -> Dir -> Int -> DualGraph
applyStraight g (PathInfo path par) dir pathlength = mkDualGraph nv' adj' deg5' el'
  where
    nv  = numVertices g
    nv' = nv + pathlength

    newVertexNbrs i
        | i == 0 && i == pathlength - 1 = error "applyStraight: pathlength must be >= 2"
        | i == 0 = case dir of
            DRight -> [path!!0, par!!0, par!!1, nv+1, path!!1]
            DLeft  -> [path!!0, path!!1, nv+1, par!!1, par!!0]
        | i == pathlength - 1 = case dir of
            DRight -> [nv+i-1, par!!i, par!!(i+1), path!!(i+1), path!!i]
            DLeft  -> [nv+i-1, path!!i, path!!(i+1), par!!(i+1), par!!i]
        | otherwise = case dir of
            DRight -> [nv+i-1, par!!i, par!!(i+1), nv+i+1, path!!(i+1), path!!i]
            DLeft  -> [nv+i-1, path!!i, path!!(i+1), nv+i+1, par!!(i+1), par!!i]

    -- Start from persistent edge list
    adj0 = neighbours g
    el_  = edgeList g

    -- Add new vertices and register them
    adj1 = foldl' (\acc i -> IM.insert (nv + i) (newVertexNbrs i) acc)
                  adj0 [0..pathlength-1]
    el0  = foldl' (\acc i -> registerVertex acc (nv+i) (newVertexNbrs i))
                  el_ [0..pathlength-1]

    -- Insert edges at endpoints (degree-5 → degree-6)
    (adj2, el1) = case dir of
        DRight -> insertAfterEL adj1 el0 (path!!0) (par!!0) nv
        DLeft  -> insertAfterEL adj1 el0 (path!!0) (path!!1) nv
    (adj3, el2) = case dir of
        DRight -> insertAfterEL adj2 el1 (par!!pathlength) (path!!pathlength) (nv+pathlength-1)
        DLeft  -> insertAfterEL adj2 el1 (par!!pathlength) (prevCW g (par!!pathlength) (path!!pathlength)) (nv+pathlength-1)

    -- Replace neighbors on path vertices
    (adj4, el3) = foldl' (\(acc, el) i ->
        let (acc', el') = replaceNbrEL acc el (path!!i) (par!!(i-1)) (nv+i-1)
        in if i < pathlength
           then replaceNbrEL acc' el' (path!!i) (par!!i) (nv+i)
           else (acc', el'))
        (adj3, el2) [1..pathlength]

    -- Replace neighbors on parallel path vertices
    (adj', el') = foldl' (\(acc, el) i ->
        let (acc', el') = if i > 0
                          then replaceNbrEL acc el (par!!i) (path!!i) (nv+i-1)
                          else (acc, el)
        in replaceNbrEL acc' el' (par!!i) (path!!(i+1)) (nv+i))
        (adj4, el3) [0..pathlength-1]

    -- Degree-5 list
    deg5' = sort $ (nv) : (nv+pathlength-1) :
        filter (\v -> v /= path!!0 && v /= par!!pathlength) (degree5 g)

---------------------------------------------------------------------
-- Straight reduction (L_i, i >= 1)
---------------------------------------------------------------------

reduceStraight :: DualGraph -> PathInfo -> DualGraph
reduceStraight g (PathInfo path par) = mkDualGraphLite nv' adj' deg5'
  where
    -- Determine pathlength from the path array
    -- path has pathlength+1 entries, par has pathlength+1 entries
    pathlength = length path - 1
    nv  = numVertices g
    nv' = nv - pathlength

    -- The new vertices that were added are nv'..nv'-1+pathlength-1 = nv'..nv-1
    adj0 = neighbours g

    -- Remove new vertices
    adj1 = foldl' (\acc i -> IM.delete (nv' + i) acc) adj0 [0..pathlength-1]

    -- Remove inserted edges at endpoints
    adj2 = removeNbr adj1 (path!!0) nv'
    adj3 = removeNbr adj2 (par!!pathlength) (nv' + pathlength - 1)

    -- Restore path vertex neighbors
    adj4 = foldl' (\acc i ->
        let acc' = replaceNbr acc (path!!i) (nv'+i-1) (par!!(i-1))
        in if i < pathlength
           then replaceNbr acc' (path!!i) (nv'+i) (par!!i)
           else acc')
        adj3 [1..pathlength]

    -- Restore parallel path vertex neighbors
    adj' = foldl' (\acc i ->
        let acc' = if i > 0
                   then replaceNbr acc (par!!i) (nv'+i-1) (path!!i)
                   else acc
        in replaceNbr acc' (par!!i) (nv'+i) (path!!(i+1)))
        adj4 [0..pathlength-1]

    -- Degree-5 list
    deg5' = sort $ (path!!0) : (par!!pathlength) :
        filter (\v -> v /= nv' && v /= (nv'+pathlength-1)) (degree5 g)

---------------------------------------------------------------------
-- B_{0,0} expansion: adds 3 vertices
---------------------------------------------------------------------

applyBentZero :: DualGraph -> PathInfo -> Dir -> DualGraph
applyBentZero g (PathInfo path par) dir = mkDualGraph nv' adj' deg5' el'
  where
    nv  = numVertices g
    a   = nv           -- degree-5
    b   = nv + 1       -- degree-6 (bend vertex)
    c   = nv + 2       -- degree-5
    nv' = nv + 3

    p0 = path!!0; p1 = path!!1; p2 = path!!2; p3 = path!!3; p4 = path!!4
    q0 = par!!0;  q1 = par!!1;  q2 = par!!2

    aNbrs = case dir of
        DRight -> [p0, q0, q1, b, p1]
        DLeft  -> [p0, p1, b, q1, q0]
    bNbrs = case dir of
        DRight -> [a, q1, c, p3, p2, p1]
        DLeft  -> [a, p1, p2, p3, c, q1]
    cNbrs = case dir of
        DRight -> [b, q1, q2, p4, p3]
        DLeft  -> [b, p3, p4, q2, q1]

    -- Start from persistent edge list
    adj0 = neighbours g
    el_  = edgeList g

    -- Add new vertices and register them
    adj1 = IM.insert a aNbrs $ IM.insert b bNbrs $ IM.insert c cNbrs adj0
    el0  = registerVertex (registerVertex (registerVertex el_ a aNbrs) b bNbrs) c cNbrs

    -- Insert edges at endpoints
    (adj2, el1) = case dir of
        DRight -> insertAfterEL adj1 el0 p0 q0 a
        DLeft  -> insertAfterEL adj1 el0 p0 p1 a
    (adj3, el2) = case dir of
        DRight -> insertAfterEL adj2 el1 p4 p3 c
        DLeft  -> insertAfterEL adj2 el1 p4 (prevCW g p4 p3) c

    -- path[1]: replace par[0] → a, par[1] → b
    (adj4, el3) = let (a1, e1) = replaceNbrEL adj3 el2 p1 q0 a
                  in  replaceNbrEL a1 e1 p1 q1 b
    -- path[2]: replace par[1] → b
    (adj5, el4) = replaceNbrEL adj4 el3 p2 q1 b
    -- path[3]: replace par[1] → b, par[2] → c
    (adj6, el5) = let (a1, e1) = replaceNbrEL adj5 el4 p3 q1 b
                  in  replaceNbrEL a1 e1 p3 q2 c

    -- par[0]: replace path[1] → a
    (adj7, el6) = replaceNbrEL adj6 el5 q0 p1 a
    -- par[1]: replace p1 → a, p2 → b, p3 → c
    (adj8, el7) = let (a1, e1) = replaceNbrEL adj7 el6 q1 p1 a
                      (a2, e2) = replaceNbrEL a1 e1 q1 p2 b
                  in  replaceNbrEL a2 e2 q1 p3 c
    -- par[2]: replace path[3] → c
    (adj', el') = replaceNbrEL adj8 el7 q2 p3 c

    deg5' = sort $ a : c : filter (\v -> v /= p0 && v /= p4) (degree5 g)

---------------------------------------------------------------------
-- B_{0,0} reduction: removes 3 vertices
---------------------------------------------------------------------

reduceBentZero :: DualGraph -> PathInfo -> DualGraph
reduceBentZero g (PathInfo path par) = mkDualGraphLite nv' adj' deg5'
  where
    nv  = numVertices g
    a   = nv - 3;  b = nv - 2;  c = nv - 1
    nv' = nv - 3

    p0 = path!!0; p1 = path!!1; p2 = path!!2; p3 = path!!3; p4 = path!!4
    q0 = par!!0;  q1 = par!!1;  q2 = par!!2

    adj0 = neighbours g
    adj1 = IM.delete a $ IM.delete b $ IM.delete c adj0

    -- Remove inserted edges
    adj2 = removeNbr (removeNbr adj1 p0 a) p4 c

    -- Restore path neighbors
    adj3 = replaceNbr (replaceNbr adj2 p1 a q0) p1 b q1
    adj4 = replaceNbr adj3 p2 b q1
    adj5 = replaceNbr (replaceNbr adj4 p3 b q1) p3 c q2

    -- Restore parallel neighbors
    adj6 = replaceNbr adj5 q0 a p1
    adj7 = replaceNbr (replaceNbr (replaceNbr adj6 q1 a p1) q1 b p2) q1 c p3
    adj' = replaceNbr adj7 q2 c p3

    deg5' = sort $ p0 : p4 : filter (\v -> v /= a && v /= c) (degree5 g)

---------------------------------------------------------------------
-- B_{i,j} expansion (i+j > 0): adds bentLen+3 vertices
---------------------------------------------------------------------

applyBent :: DualGraph -> PathInfo -> Dir -> Int -> Int -> DualGraph
applyBent g (PathInfo path par) dir bentPos bentLen = mkDualGraph nv' adj' deg5' el'
  where
    nv  = numVertices g
    numNew = bentLen + 3
    nv' = nv + numNew

    -- Index into new vertices: nv+0 .. nv+numNew-1
    -- nv+0 and nv+numNew-1 are degree-5, rest degree-6

    -- The bend vertex is at index bentPos+1 in the new vertices.
    -- Before bend: new vertices 0..bentPos (standard straight pattern)
    -- Bend vertex: bentPos+1 (connects to 3 path vertices)
    -- After bend: bentPos+2..numNew-1 (shifted straight pattern)

    -- Build neighbor lists for new vertices.
    -- Before bend (i=0..bentPos): same as straight
    -- Bend vertex (i=bentPos+1): 6 neighbors, 3 path vertices
    -- After bend (i=bentPos+2..numNew-1): shifted pattern

    newVertexNbrs i
        -- Before the bend
        | i <= bentPos =
            let pred = if i > 0 then [nv+i-1] else []
                succ' = [nv+i+1]
                core = [path!!i, path!!(i+1)]
                sides = [par!!(i+1), par!!i]
            in case dir of
                DRight -> pred ++ reverse sides ++ succ' ++ reverse core
                    -- Actually, let me derive this more carefully from the C code.
                    -- This is getting complex, let me use a different approach.
                DLeft  -> pred ++ core ++ succ' ++ sides

        -- Bend vertex
        | i == bentPos + 1 =
            let pred = [nv+i-1]
                succ' = [nv+i+1]
                -- Bend vertex connects to 3 path vertices:
                -- path[i], path[i+1], path[i+2]
                -- And 1 parallel vertex: par[i] (= par[bentPos+1])
                cores = [path!!i, path!!(i+1), path!!(i+2)]
                side = [par!!i]
            in case dir of
                DRight -> pred ++ reverse side ++ succ' ++ reverse cores
                DLeft  -> pred ++ cores ++ succ' ++ side

        -- After the bend (shifted indices: path[i+1] and path[i+2])
        | otherwise =
            let pred = [nv+i-1]
                succ' = if i < numNew - 1 then [nv+i+1] else []
                core = [path!!(i+1), path!!(i+2)]
                sides = [par!!i, par!!(i-1)]
            in case dir of
                DRight -> pred ++ reverse sides ++ succ' ++ reverse core
                DLeft  -> pred ++ core ++ succ' ++ sides
    -- OK, deriving CW orders from the C code is error-prone with this approach.
    -- Let me directly implement the C code's wiring pattern instead.

    -- Direct implementation following C code extend_bent
    -- For use_next (DRight), walking via ->prev from firstedge:

    newVertexNbrsDirect i
        -- Before the bend (i = 0..bentPos)
        | i <= bentPos && i == 0 = case dir of
            -- degree-5, no predecessor
            -- DRight prev traversal: path[0], path[1], nv+1, par[1], par[0]
            -- CW: [path[0], par[0], par[1], nv+1, path[1]]
            DRight -> [path!!0, par!!0, par!!1, nv+1, path!!1]
            -- DLeft next traversal: path[0], path[1], nv+1, par[1], par[0]
            -- CW: same order
            DLeft  -> [path!!0, path!!1, nv+1, par!!1, par!!0]

        | i <= bentPos = case dir of
            -- degree-6, has predecessor and successor
            -- DRight prev traversal: nv+i-1, path[i], path[i+1], nv+i+1, par[i+1], par[i]
            -- CW: [nv+i-1, par[i], par[i+1], nv+i+1, path[i+1], path[i]]
            DRight -> [nv+i-1, par!!i, par!!(i+1), nv+i+1, path!!(i+1), path!!i]
            DLeft  -> [nv+i-1, path!!i, path!!(i+1), nv+i+1, par!!(i+1), par!!i]

        -- Bend vertex (i = bentPos+1), always degree-6
        | i == bentPos + 1 = case dir of
            -- DRight prev traversal: nv+i-1, path[i], path[i+1], path[i+2], nv+i+1, par[i]
            -- CW: [nv+i-1, par[i], nv+i+1, path[i+2], path[i+1], path[i]]
            DRight -> [nv+i-1, par!!i, nv+i+1, path!!(i+2), path!!(i+1), path!!i]
            DLeft  -> [nv+i-1, path!!i, path!!(i+1), path!!(i+2), nv+i+1, par!!i]

        -- After the bend (i = bentPos+2..numNew-1)
        | i == numNew - 1 = case dir of
            -- Last vertex, degree-5, no successor
            -- DRight prev traversal: nv+i-1, path[i+1], path[i+2], par[i], par[i-1]
            -- CW: [nv+i-1, par[i-1], par[i], path[i+2], path[i+1]]
            DRight -> [nv+i-1, par!!(i-1), par!!i, path!!(i+2), path!!(i+1)]
            DLeft  -> [nv+i-1, path!!(i+1), path!!(i+2), par!!i, par!!(i-1)]

        | otherwise = case dir of
            -- degree-6, has predecessor and successor
            -- DRight prev traversal: nv+i-1, path[i+1], path[i+2], nv+i+1, par[i], par[i-1]
            -- CW: [nv+i-1, par[i-1], par[i], nv+i+1, path[i+2], path[i+1]]
            DRight -> [nv+i-1, par!!(i-1), par!!i, nv+i+1, path!!(i+2), path!!(i+1)]
            DLeft  -> [nv+i-1, path!!(i+1), path!!(i+2), nv+i+1, par!!i, par!!(i-1)]

    -- Start from persistent edge list
    adj0 = neighbours g
    el_  = edgeList g

    -- Add new vertices and register them
    adj1 = foldl' (\acc i -> IM.insert (nv+i) (newVertexNbrsDirect i) acc)
                  adj0 [0..numNew-1]
    el0  = foldl' (\acc i -> registerVertex acc (nv+i) (newVertexNbrsDirect i))
                  el_ [0..numNew-1]

    -- Insert edges at endpoints
    (adj2, el1) = case dir of
        DRight -> insertAfterEL adj1 el0 (path!!0) (par!!0) nv
        DLeft  -> insertAfterEL adj1 el0 (path!!0) (path!!1) nv
    (adj3, el2) = case dir of
        DRight -> insertAfterEL adj2 el1 (path!!(bentLen+4)) (path!!(bentLen+3)) (nv+numNew-1)
        DLeft  -> insertAfterEL adj2 el1 (path!!(bentLen+4)) (prevCW g (path!!(bentLen+4)) (path!!(bentLen+3))) (nv+numNew-1)

    -- Replace neighbors: before the bend
    (adj4, el3) = foldl' (\(acc, el) i ->
        let (acc', el') = replaceNbrEL acc el (path!!i) (par!!(i-1)) (nv+i-1)
        in replaceNbrEL acc' el' (path!!i) (par!!i) (nv+i))
        (adj3, el2) [1..bentPos+1]

    -- par[i]: replace path[i] → nv+i-1 (for i>0), path[i+1] → nv+i
    (adj5, el4) = foldl' (\(acc, el) i ->
        let (acc', el') = if i > 0
                          then replaceNbrEL acc el (par!!i) (path!!i) (nv+i-1)
                          else (acc, el)
        in replaceNbrEL acc' el' (par!!i) (path!!(i+1)) (nv+i))
        (adj4, el3) [0..bentPos]

    -- Bend area (i = bentPos + 2 in path indexing)
    bendI = bentPos + 2
    (adj6, el5) = replaceNbrEL adj5 el4 (path!!bendI) (par!!(bendI-1)) (nv+bendI-1)
    (adj7, el6) = replaceNbrEL adj6 el5 (par!!(bendI-1)) (path!!(bendI-1)) (nv+bendI-2)
    (adj8, el7) = replaceNbrEL adj7 el6 (par!!(bendI-1)) (path!!bendI) (nv+bendI-1)
    (adj9, el8) = replaceNbrEL adj8 el7 (par!!(bendI-1)) (path!!(bendI+1)) (nv+bendI)

    -- After the bend: path[i] replace par[i-2] → nv+i-2, par[i-1] → nv+i-1
    (adj10, el9) = foldl' (\(acc, el) i ->
        let (acc', el') = replaceNbrEL acc el (path!!i) (par!!(i-2)) (nv+i-2)
        in replaceNbrEL acc' el' (path!!i) (par!!(i-1)) (nv+i-1))
        (adj9, el8) [bendI+1..bentLen+3]

    -- par[i] (after bend): replace path[i+1] → nv+i, path[i+2] → nv+i+1
    (adj', el') = foldl' (\(acc, el) i ->
        let (acc', el') = replaceNbrEL acc el (par!!i) (path!!(i+1)) (nv+i)
        in if i < bentLen + 2
           then replaceNbrEL acc' el' (par!!i) (path!!(i+2)) (nv+i+1)
           else (acc', el'))
        (adj10, el9) [bentPos+2..bentLen+2]

    deg5' = sort $ nv : (nv+numNew-1) :
        filter (\v -> v /= path!!0 && v /= path!!(bentLen+4)) (degree5 g)

---------------------------------------------------------------------
-- B_{i,j} reduction (i+j > 0)
---------------------------------------------------------------------

reduceBent :: DualGraph -> PathInfo -> Int -> Int -> DualGraph
reduceBent g (PathInfo path par) bentPos bentLen = mkDualGraphLite nv' adj' deg5'
  where
    numNew = bentLen + 3
    nv  = numVertices g
    nv' = nv - numNew

    adj0 = neighbours g

    -- Remove new vertices
    adj1 = foldl' (\acc i -> IM.delete (nv'+i) acc) adj0 [0..numNew-1]

    -- Remove inserted edges
    adj2 = removeNbr adj1 (path!!0) nv'
    adj3 = removeNbr adj2 (path!!(bentLen+4)) (nv'+numNew-1)

    -- Restore before-bend path neighbors
    adj4 = foldl' (\acc i ->
        let acc' = replaceNbr acc (path!!i) (nv'+i-1) (par!!(i-1))
        in replaceNbr acc' (path!!i) (nv'+i) (par!!i))
        adj3 [1..bentPos+1]

    adj5 = foldl' (\acc i ->
        let acc' = if i > 0
                   then replaceNbr acc (par!!i) (nv'+i-1) (path!!i)
                   else acc
        in replaceNbr acc' (par!!i) (nv'+i) (path!!(i+1)))
        adj4 [0..bentPos]

    -- Restore bend area
    bendI = bentPos + 2
    adj6 = replaceNbr adj5 (path!!bendI) (nv'+bendI-1) (par!!(bendI-1))
    adj7 = replaceNbr adj6 (par!!(bendI-1)) (nv'+bendI-2) (path!!(bendI-1))
    adj8 = replaceNbr adj7 (par!!(bendI-1)) (nv'+bendI-1) (path!!bendI)
    adj9 = replaceNbr adj8 (par!!(bendI-1)) (nv'+bendI) (path!!(bendI+1))

    -- Restore after-bend path neighbors
    adj10 = foldl' (\acc i ->
        let acc' = replaceNbr acc (path!!i) (nv'+i-2) (par!!(i-2))
        in replaceNbr acc' (path!!i) (nv'+i-1) (par!!(i-1)))
        adj9 [bendI+1..bentLen+3]

    adj' = foldl' (\acc i ->
        let acc' = replaceNbr acc (par!!i) (nv'+i) (path!!(i+1))
        in if i < bentLen + 2
           then replaceNbr acc' (par!!i) (nv'+i+1) (path!!(i+2))
           else acc')
        adj10 [bentPos+2..bentLen+2]

    deg5' = sort $ (path!!0) : (path!!(bentLen+4)) :
        filter (\v -> v /= nv' && v /= (nv'+numNew-1)) (degree5 g)

---------------------------------------------------------------------
-- Expansion enumeration
---------------------------------------------------------------------

-- | Enumerate all valid expansions of a graph up to the given max reduction length.
expansions :: Int -> DualGraph -> [Expansion]
expansions maxLen g =
       expansionsL0 g
    ++ (if maxLen >= 2 then expansionsL maxLen g else [])
    ++ (if maxLen >= 2 then expansionsB maxLen g else [])

-- | Enumerate all L0 expansions.
-- Matches C code's find_L0_extensions_next/prev: for each degree-5 vertex u,
-- iterate over ALL neighbors v (v can be any degree), compute the straight
-- path, and check that par[2] (= q2) is degree-5.
--
-- Note: v does NOT need to be degree-5. The L0 expansion only promotes
-- path[0] (= u, degree-5) and par[2] (= q2, degree-5) to degree-6.
-- Vertex v (= path[1]) keeps its original degree.
--
-- We enumerate BOTH (u, v) and (v', u') orderings because the expansion
-- path extends past v, so different orderings produce different children.
expansionsL0 :: DualGraph -> [Expansion]
expansionsL0 g =
    [ Exp (L 0) (u, v) d
    | u <- degree5 g
    , v <- nbrs g u         -- v can be any degree (C code doesn't filter)
    , d <- [DLeft, DRight]
    , isValidStraightSite g (u, v) d 3
    ]

-- | Enumerate all L_i (i >= 1) expansions up to the given max reduction length.
-- path[0] and par[pathlength] must be degree-5 (they gain an edge via insertAfter).
-- Path and parallel path must be disjoint.
expansionsL :: Int -> DualGraph -> [Expansion]
expansionsL maxLen g =
    [ Exp (L i) (u, v) d
    | u <- degree5 g
    , v <- nbrs g u
    , d <- [DLeft, DRight]
    , i <- [1 .. maxLen - 1]
    , let numE = i + 3
    , canStraightPath g (u, v) d numE
    , isValidStraightSite g (u, v) d numE
    ]

-- | Check if a straight expansion site is valid: par[last] is degree-5,
-- path/parallel path are disjoint, and all adjacencies needed by
-- replaceNbrEL hold (essential for non-triangulated graphs).
isValidStraightSite :: DualGraph -> Edge -> Dir -> Int -> Bool
isValidStraightSite g edge dir numEntries =
    let pi = computeStraightPath g edge dir numEntries
        p = mainPath pi
        q = parallelPath pi
        adj = neighbours g
        isAdj u v = v `elem` (adj IM.! u)
        pathlength = numEntries - 1
        qLast = last q
        -- path[k] must be adj to par[k-1] and par[k] (for replaceNbrEL)
        pathAdj = all (\k -> isAdj (p!!k) (q!!(k-1))
                          && isAdj (p!!k) (q!!k))
                      [1..pathlength]
        -- par[k] must be adj to path[k+1] (for replaceNbrEL)
        parAdj = all (\k -> isAdj (q!!k) (p!!(k+1))) [0..pathlength-1]
    in deg g qLast == 5
    && null (p `intersect` q)
    && pathAdj && parAdj

-- | Check if a straight path of the given length exists.
canStraightPath :: DualGraph -> Edge -> Dir -> Int -> Bool
canStraightPath g (u0, v0) dir numEntries = go (u0, v0) numEntries []
  where
    advance (from, to) = (to, straightAhead g dir to from)
    go _ 0 _ = True
    go (from, to) n visited
        | from `elem` visited = False  -- self-intersecting
        | otherwise = go (advance (from, to)) (n - 1) (from : visited)

-- | Enumerate all B_{i,j} expansions up to the given max reduction length.
expansionsB :: Int -> DualGraph -> [Expansion]
expansionsB maxLen g =
    [ Exp (B i j) (u, v) d
    | u <- degree5 g
    , v <- nbrs g u
    , d <- [DLeft, DRight]
    , let maxBentLen = maxLen - 2  -- i+j+2 <= maxLen, so i+j <= maxLen-2
    , bentLen <- [0 .. maxBentLen]
    , bentPos <- [0 .. bentLen]
    , let i = bentPos
          j = bentLen - bentPos
    -- B_{0,0} (bent path of length 2) is a valid expansion; must be included
    , canBentPath g (u, v) d bentPos bentLen
    ]

-- | Check if a bent path of the given parameters exists and produces a valid expansion site.
-- Checks: no self-intersection, last path vertex is degree-5,
-- path and parallel path are disjoint, and all required adjacencies exist for the
-- replacement operations (important on non-triangulated graphs from previous expansions).
--
-- NOTE: Unlike straight expansions, bent expansions do NOT require par[last] to be
-- degree-5. The C code (find_bent_zero_extensions_next line 6036, find_bent_extensions_next
-- line 6351) only checks that the main path endpoint is degree-5. The parallel path
-- vertices serve as the "opposite side" of the local patch and have no degree constraint.
canBentPath :: DualGraph -> Edge -> Dir -> Int -> Int -> Bool
canBentPath g (u0, v0) dir bentPos bentLen =
    case computeBentPathSafe g (u0, v0) dir bentPos bentLen of
        Nothing -> False
        Just (PathInfo path par) ->
            let lastV = last path
                adj = neighbours g
                isAdj u v = v `elem` (adj IM.! u)
                bendI = bentPos + 2
                -- Check all adjacencies required by applyBent replacements:
                -- Before bend: path[k] adj par[k-1] and par[k]; par[k] adj path[k+1]
                beforeOk = all (\k -> isAdj (path!!k) (par!!(k-1))
                                   && isAdj (path!!k) (par!!k))
                               [1..bentPos+1]
                parBeforeOk = all (\k -> (k == 0 || isAdj (par!!k) (path!!k))
                                      && isAdj (par!!k) (path!!(k+1)))
                                  [0..bentPos]
                -- Bend area: path[bendI] adj par[bendI-1];
                -- par[bendI-1] adj path[bendI-1], path[bendI], path[bendI+1]
                bendOk = isAdj (path!!bendI) (par!!(bendI-1))
                      && isAdj (par!!(bendI-1)) (path!!(bendI-1))
                      && isAdj (par!!(bendI-1)) (path!!bendI)
                      && isAdj (par!!(bendI-1)) (path!!(bendI+1))
                -- After bend: path[k] adj par[k-2] and par[k-1]
                afterOk = all (\k -> isAdj (path!!k) (par!!(k-2))
                                  && isAdj (path!!k) (par!!(k-1)))
                              [bendI+1..bentLen+3]
                -- par after bend: par[k] adj path[k+1] and path[k+2] (if k < bentLen+2)
                parAfterOk = all (\k -> isAdj (par!!k) (path!!(k+1))
                                     && (k >= bentLen+2 || isAdj (par!!k) (path!!(k+2))))
                                 [bentPos+2..bentLen+2]
                -- DLeft insertion at endpoint: par[bentLen+2] must be in path[bentLen+4]'s list
                endpointOk = case dir of
                    DRight -> True  -- uses path[bentLen+3] which is always a neighbor
                    DLeft  -> isAdj (path!!(bentLen+4)) (par!!(bentLen+2))
            in deg g lastV == 5
            && head path /= lastV
            && null (path `intersect` par)
            && beforeOk && parBeforeOk && bendOk && afterOk && parAfterOk
            && endpointOk

-- | Safe version of computeBentPath that returns Nothing on self-intersection.
computeBentPathSafe :: DualGraph -> Edge -> Dir -> Int -> Int -> Maybe PathInfo
computeBentPathSafe g (u0, v0) dir bentPos bentLen = result
  where
    adv (from, to) = (to, straightAhead g dir to from)

    result = do
        -- Phase 1
        let preBendEdges = take (bentPos + 2) $ iterate adv (u0, v0)
            prePath = map fst preBendEdges
        -- Check no self-intersection in pre-bend path
        if hasDups prePath then Nothing else Just ()

        let (preLast, turnV) = last preBendEdges
            afterTurn = turnAhead g dir turnV preLast
            postBendEdges = take (bentLen - bentPos + 1) $
                            iterate adv (turnV, afterTurn)
            postPath = map fst postBendEdges
            (lastFrom, lastTo) = last postBendEdges
            finalV = straightAhead g dir lastTo lastFrom

        -- Check no self-intersection overall
        let allPath = prePath ++ postPath ++ [lastTo, finalV]
        if hasDups allPath then Nothing else do
            -- Compute parallel path uniformly as sideNbr(path[k], path[k+1])
            -- for k = 0..bentLen+2. This matches computeBentPath.
            let parVerts = [ sideNbr g dir (allPath !! k) (allPath !! (k+1))
                           | k <- [0 .. bentLen + 2]
                           ]
            Just (PathInfo allPath parVerts)

    hasDups xs = go IS.empty xs
      where go _ []     = False
            go s (x:rest) = IS.member x s || go (IS.insert x s) rest

---------------------------------------------------------------------
-- F expansion: nanotube ring (adds 5 degree-6 vertices)
---------------------------------------------------------------------

-- | Find the equatorial ring and outer path of a (5,0) nanotube.
-- Returns (ring, outer) where ring is a 5-cycle of degree-6 vertices
-- and outer[i] is the vertex between ring[i] and ring[(i+1)%5] on one cap.
-- Returns Nothing if no valid ring is found.
--
-- Algorithm (matching buckygen's has_lm_path):
-- 1. Try every directed edge between degree-6 vertices
-- 2. Follow "straightAhead" for 5 steps (3 positions = opposite in hexagon)
-- 3. If the path closes back to the starting edge, it's an equatorial ring
-- 4. Compute outers as common neighbors of consecutive ring vertices on a
--    consistent cap side; outers need NOT be degree-5 (for deeper nanotubes
--    they are degree-6 vertices from a previous Ring expansion)
findNanotubeRing :: DualGraph -> Maybe ([Vertex], [Vertex])
findNanotubeRing g
    | length deg6 < 5 = Nothing
    | otherwise = firstJust
        [ traceFromEdge u v
        | u <- deg6
        , v <- nbrs g u
        , deg g v == 6
        ]
  where
    nv = numVertices g
    deg6 = [v | v <- [0..nv-1], deg g v == 6]

    firstJust [] = Nothing
    firstJust (Nothing : rest) = firstJust rest
    firstJust (Just x : _) = Just x

    -- Trace a ring using straightAhead from edge (u0, u1).
    traceFromEdge u0 u1 =
        case traceRingSA u0 u1 of
            Nothing   -> Nothing
            Just ring -> computeOuter ring

    traceRingSA u0 u1 = step u0 u1 [u0, u1] 2
      where
        step prev cur ring n
            | n == 5 =
                let next = straightAhead g DRight cur prev
                in if next == u0 then Just ring else Nothing
            | otherwise =
                let next = straightAhead g DRight cur prev
                in if deg g next /= 6 || next `elem` ring
                   then Nothing
                   else step cur next (ring ++ [next]) (n + 1)

    -- Compute outers: for each ring edge, find common neighbors not in
    -- the ring (two per edge, one per cap side).  Pick a consistent cap
    -- by requiring each outer to be adjacent to the previous one.
    -- Validate with outersAdjacent.
    computeOuter ring =
        let commonNbrs i =
                let u = ring !! i
                    v = ring !! ((i + 1) `mod` 5)
                in [w | w <- nbrs g u, w `elem` nbrs g v, w `notElem` ring]
            allCommon = map commonNbrs [0..4]
        in if any null allCommon
           then Nothing
           else let results = [ ov
                              | first <- head allCommon
                              , let ov = resolveOuters allCommon first
                              , length ov == 5
                              , length (nub ov) == 5
                              , outersAdjacent ring ov
                              ]
                in case results of
                    (ov:_) -> Just (ring, ov)
                    []     -> Nothing

    resolveOuters allCommon first = go [first] 1
      where
        go outers 5 = outers
        go outers i =
            let prev = last outers
                candidates = allCommon !! i
                pick = case filter (`elem` nbrs g prev) candidates of
                           (c:_) -> c
                           []    -> case candidates of
                                     [c] -> c
                                     _   -> head candidates
            in go (outers ++ [pick]) (i + 1)

    -- Verify that for each ring vertex, its two outer vertices are
    -- CW-adjacent (no extra non-ring vertices between them).
    outersAdjacent ring' ov = all check [0..4]
      where
        check i =
            let rv    = ring' !! i
                oPrev = ov !! ((i - 1 + 5) `mod` 5)
                oCurr = ov !! i
                nbs   = nbrs g rv
            in case (elemIndex oPrev nbs, elemIndex oCurr nbs) of
                (Just pi, Just ci) ->
                    let n  = length nbs
                        between = [ nbs !! (j `mod` n)
                                  | j <- [pi + 1 .. (if ci > pi then ci else ci + n) - 1] ]
                        extras  = [v | v <- between, v `notElem` ring']
                    in null extras
                _ -> False

-- | Apply the F (nanotube ring) expansion.
-- Adds 5 new degree-6 vertices between the equatorial ring and one cap.
-- ring = [r0, r1, r2, r3, r4] (the 5-cycle of degree-6 vertices)
-- outer = [o0, o1, o2, o3, o4] (the cap vertices, where outer[i]
--         is between ring[i] and ring[(i+1)%5])
applyRing :: DualGraph -> [Vertex] -> [Vertex] -> DualGraph
applyRing g ring outer = mkDualGraph nv' adj' deg5' el'
  where
    nv  = numVertices g
    nv' = nv + 5

    -- New vertex indices
    new i = nv + i

    -- CW neighbor list for each new vertex
    -- nv+i connects to: ring[i], then CW around the vertex:
    --   new[(i-1+5)%5], outer[i], outer[(i+1)%5], new[(i+1)%5], ring[(i+1)%5]
    newNbrs i = [ ring !! i
                , new ((i - 1 + 5) `mod` 5)
                , outer !! i
                , outer !! ((i + 1) `mod` 5)
                , new ((i + 1) `mod` 5)
                , ring !! ((i + 1) `mod` 5)
                ]

    adj0 = neighbours g
    el_  = edgeList g

    -- Add new vertices
    adj1 = foldl' (\acc i -> IM.insert (new i) (newNbrs i) acc) adj0 [0..4]
    el0  = foldl' (\acc i -> registerVertex acc (new i) (newNbrs i)) el_ [0..4]

    -- Update ring vertices: replace cap-side degree-5 neighbors with new vertices
    -- ring[i]: replace outer[(i-1+5)%5] with nv+((i-1+5)%5),
    --          replace outer[i] with nv+i
    (adj2, el1) = foldl' (\(acc, el) i ->
        let prev = (i - 1 + 5) `mod` 5
            (a1, e1) = replaceNbrEL acc el (ring !! i) (outer !! prev) (new prev)
        in  replaceNbrEL a1 e1 (ring !! i) (outer !! i) (new i))
        (adj1, el0) [0..4]

    -- Update outer (cap) vertices: replace ring neighbors with new vertices
    -- outer[i]: replace ring[i] with nv+((i-1+5)%5),
    --           replace ring[(i+1)%5] with nv+i
    (adj', el') = foldl' (\(acc, el) i ->
        let prev = (i - 1 + 5) `mod` 5
            (a1, e1) = replaceNbrEL acc el (outer !! i) (ring !! i) (new prev)
        in  replaceNbrEL a1 e1 (outer !! i) (ring !! ((i + 1) `mod` 5)) (new i))
        (adj2, el1) [0..4]

    -- degree-5 list unchanged (all new vertices are degree 6)
    deg5' = degree5 g

-- | Reduce a nanotube ring (F) expansion — inverse of applyRing.
-- Takes the ring and outer arrays from the original expansion.
-- The 5 added vertices are assumed to be (nv-5)..(nv-1).
--
-- Reversal of applyRing:
--   applyRing replaced in ring[i]: outer[(i-1+5)%5] -> newV[(i-1+5)%5],
--                                  outer[i] -> newV[i]
--   applyRing replaced in outer[i]: ring[i] -> newV[(i-1+5)%5],
--                                   ring[(i+1)%5] -> newV[i]
-- reduceRing undoes these replacements and deletes newV[0..4].
reduceRing :: DualGraph -> [Vertex] -> [Vertex] -> DualGraph
reduceRing g ring outer = mkDualGraphLite nv' adj' deg5'
  where
    nv  = numVertices g
    nv' = nv - 5

    newV i = nv' + i

    adj0 = neighbours g

    -- Delete the 5 new ring vertices
    adj1 = foldl' (\acc i -> IM.delete (newV i) acc) adj0 [0..4]

    -- Restore ring vertices: reverse the two replacements per ring vertex
    adj2 = foldl' (\acc i ->
        let prev = (i - 1 + 5) `mod` 5
            acc' = replaceNbr acc (ring !! i) (newV prev) (outer !! prev)
        in  replaceNbr acc' (ring !! i) (newV i) (outer !! i))
        adj1 [0..4]

    -- Restore outer (cap) vertices: reverse the two replacements per outer vertex
    adj' = foldl' (\acc i ->
        let prev = (i - 1 + 5) `mod` 5
            acc' = replaceNbr acc (outer !! i) (newV prev) (ring !! i)
        in  replaceNbr acc' (outer !! i) (newV i) (ring !! ((i + 1) `mod` 5)))
        adj2 [0..4]

    -- degree-5 list unchanged (removed vertices were all degree 6)
    deg5' = degree5 g

---------------------------------------------------------------------
-- Reduction helpers
---------------------------------------------------------------------

-- | Number of vertices added by an expansion (= removed by its inverse).
numNewVertices :: ExpKind -> Int
numNewVertices (L 0)   = 2
numNewVertices (L i)   = i + 2
numNewVertices (B 0 0) = 3
numNewVertices (B i j) = i + j + 3
numNewVertices F       = 5

-- | Reduction length (distance between the two degree-5 endpoints).
reductionLength :: ExpKind -> Int
reductionLength (L i)   = i + 1
reductionLength (B i j) = i + j + 2
reductionLength F       = maxBound

-- | Length of the longest straight segment.
longestStraight :: ExpKind -> Int
longestStraight (L i)   = i + 1
longestStraight (B i j) = max (i + 1) (j + 1)
longestStraight F       = 0
