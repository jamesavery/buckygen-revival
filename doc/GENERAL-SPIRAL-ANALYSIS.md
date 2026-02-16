# Generalized Spiral: Analysis for Haskell Implementation

## What Changes from Regular to General

The regular spiral fails when the candidate face `v = findV(back, front)` would, if
removed, disconnect the remaining graph. The generalized spiral resolves this by
**rotating the boundary deque** — moving front to back — until the candidate face is
safe to remove.

The only differences from the regular spiral are:

1. Before peeling v, check if v is a **cut vertex** of the remaining graph.
2. If yes: rotate the boundary (cyclic shift by 1) and retry. Track the number of
   consecutive rotations as `jump_state`.
3. When a non-cut-vertex v is found and `jump_state > 0`: record `(step_index,
   jump_state)` in the jump list, reset `jump_state = 0`.
4. Then peel v normally (drain front, drain back, push).

Everything else — the drain operations, the boundary deque, the orientation logic — is
identical.

## The Cut Vertex Test

**C++ reference:** `PlanarGraph::is_cut_vertex` in `planargraph.cc:161`

For vertex v with neighbors n_0, n_1, ..., n_{k-1} in cyclic order (in the remaining
graph), count edges between consecutive pairs (n_i, n_{i+1 mod k}). Vertex v is a cut
vertex iff this count < k-1.

**Why this works:** In a triangulation, the neighbors of v form a cycle (the "link" of
v). Each edge (n_i, n_{i+1}) corresponds to a triangle (v, n_i, n_{i+1}). If all k
triangles exist, the link is a complete cycle — removing v leaves the neighbors
connected. If fewer than k-1 edges exist, the link has a gap — removing v disconnects
the graph.

**Key constraint:** This test requires at most one non-triangular face incident to v.
Since the remaining graph is a triangulation with one "hole" (the peeled region), this
constraint is satisfied: there's exactly one non-triangular face adjacent to the boundary
vertices.

### Haskell Implementation Strategy

The C++ code maintains a mutable `remaining_graph` and physically removes nodes. For
Haskell, I avoid this by using the **original graph** plus the **removed set**:

1. **v's remaining neighbors:** Filter v's original neighbors, removing any in the
   removed set. The cyclic order is preserved since we just skip elements.
2. **Edge existence:** Check in the original graph. Edges between non-removed vertices
   are unchanged (only edges incident to removed vertices are gone).

This is correct because `remove_node(u)` only removes edges incident to u. Edges
between surviving vertices are untouched.

```haskell
isCutVertex :: Graph -> IS.IntSet -> Node -> Bool
isCutVertex g removed v =
  let remNbrs = [w | w <- nbrs g v, IS.notMember w removed]
      n = length remNbrs
  in n >= 2 && countConsecutiveEdges g remNbrs < n - 1

countConsecutiveEdges :: Graph -> [Node] -> Int
countConsecutiveEdges g ns =
  length [() | (a, b) <- zip ns (tail ns ++ [head ns]), isEdge g a b]

isEdge :: Graph -> Node -> Node -> Bool
isEdge g u v = V.elem v (adj g ! u)
```

## Jump Mechanics

From `triangulation.cc:601-656`:

```
for i = 3 to N-2:
  u = back(boundary), w = front(boundary)
  v = findV(u, w)

  if v is cut vertex of remaining graph:
    rotate boundary: move front to back
    do NOT advance i
    jump_state++
    continue

  if jump_state > 0:
    record jump (i, jump_state)
    jump_state = 0

  peel v normally
```

The rotation `[w, a, b, ..., u] → [a, b, ..., u, w]` changes which edge we look across.
After rotation, the new front (a) and new back (w) are consecutive in the boundary, so
they share an edge in the original graph — guaranteeing `findV` will find a valid v.

**Maximum rotations per step:** At most `|boundary| - 1` rotations before returning to
the original configuration. The Brinkmann & Dress theorem guarantees that a non-cut
vertex will be found before completing a full cycle (for 3-connected planar graphs).

## Canonical Ordering

From `spiral.hh:18-22`:

```
(jumps1, spiral1) < (jumps2, spiral2) iff:
  |jumps1| < |jumps2|, OR
  |jumps1| == |jumps2| AND jumps1 < jumps2  (lexicographic on pairs), OR
  jumps1 == jumps2 AND spiral1 < spiral2    (lexicographic)
```

Regular spirals (0 jumps) are always preferred over general spirals.

## Canonical Search Strategy

From `triangulation.cc:735-840`:

1. **Phase 1:** Try all starting triples with regular spiral only (general=false).
   If any succeed, take the lexicographic minimum. Done.
2. **Phase 2:** Only if ALL regular spirals failed, try all starting triples with
   general spiral (general=true). Take the canonical minimum.

This two-phase approach ensures regular spirals are always preferred when available.

## Windup with Jumps

The windup constructor (`triangulation.cc:275`) handles jumps:

```
if jumps is non-empty and k == jumps[0].position:
  rotate boundary by jumps[0].length positions
  remove this jump from the list
```

The rotation in the windup is identical to the rotation in the unwinding: move front to
back, repeated `jump_length` times.

## Data Types

```haskell
data GeneralSpiral = GeneralSpiral
  { gsJumps  :: [(Int, Int)]   -- (position, jump_length) pairs
  , gsSpiral :: [Int]          -- face-degree sequence
  } deriving (Eq, Show)

-- Canonical ordering matching C++ general_spiral::operator<
instance Ord GeneralSpiral where
  compare (GeneralSpiral j1 s1) (GeneralSpiral j2 s2) =
    compare (length j1) (length j2) <>
    compare j1 j2 <>
    compare s1 s2
```

## New API Functions

```haskell
-- Spiral.hs additions:
generalSpiral        :: Graph -> (Node, Node, Node) -> Maybe GeneralSpiral
canonicalGeneralSpiral :: Graph -> Maybe GeneralSpiral

-- Windup.hs addition:
windupGeneralSpiral  :: [Int] -> [(Int, Int)] -> Graph
```

## Test Plan

1. **Td-C100:** `canonicalGeneralSpiral` returns `Just` with non-empty jumps (currently
   `canonicalSpiral` returns `Nothing`).
2. **Round-trip:** Wind up the Td-C100 general spiral, extract it again, get the same
   result.
3. **Existing tests still pass:** Regular spirals on C20, C24, C28 are unchanged.
4. **Consistency:** For graphs where regular spirals exist, `canonicalGeneralSpiral`
   returns the same spiral as `canonicalSpiral` (with empty jump list).
