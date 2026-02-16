# The Spiral Algorithm for Triangulations

This document describes the spiral algorithm for oriented triangulations, as implemented
in `Spiral.hs` (unwinding) and `Windup.hs` (construction). The algorithm computes a
canonical face-degree sequence that uniquely identifies a triangulation up to isomorphism.

## Context and Purpose

**Goal:** Given two triangulations (or their dual cubic graphs), determine whether they
are isomorphic. Compute a *canonical spiral* for each; the graphs are isomorphic iff
their canonical spirals are equal.

**Complexity:** For fullerene duals (12 degree-5 vertices, rest degree-6), the canonical
spiral computation is O(N). There are at most 12 x 5 x 2 = 120 starting configurations,
each requiring one O(N) pass.

## Duality: Triangulations and Cubic Graphs

The spiral operates on **triangulations** (planar graphs where every face is a triangle),
not directly on cubic graphs.

- A **fullerene graph** is a cubic planar graph (every vertex has degree 3) with exactly
  12 pentagonal faces and the rest hexagonal.
- Its **dual** is a triangulation where:
  - Pentagon faces become degree-5 vertices
  - Hexagon faces become degree-6 vertices

The spiral records the degree of each vertex in this dual triangulation. For fullerenes,
each entry is 5 or 6.

## The Oriented Adjacency List

A triangulation is stored as oriented adjacency lists:

```haskell
newtype Graph = Graph { adj :: Vector (Vector Node) }
```

For each vertex `v`, `adj ! v` lists v's neighbors in a consistent cyclic order (all CCW
or all CW when viewed from outside the surface). The key navigation operations:

```haskell
cNext :: Graph -> Node -> Node -> Maybe Node  -- cyclic successor of w in v's list
cPrev :: Graph -> Node -> Node -> Maybe Node  -- cyclic predecessor of w in v's list
```

**Orientation invariant:** For a CCW-oriented triangulation, triangle (a, b, c) traversed
counter-clockwise satisfies: the predecessor of c in adj[a] is b. Equivalently, adj[a]
contains `...b, c, ...` in order.

## The Spiral (Output Format)

A **regular spiral** is a face-degree sequence of length N (number of vertices in the
triangulation):

```
[deg(f1), deg(f2), ..., deg(fN)]
```

where f1, f2, ..., fN is a specific ordering of all vertices determined by peeling faces
from the triangulation one at a time.

For **fullerenes**, the compact **RSPI** (Ring-Spiral Pentagon Indices) representation
stores only the 1-indexed positions where degree = 5:

```
rspi = [i | (i, d) <- zip [1..] spiral, d == 5]   -- 12 entries for fullerenes
```

## Unwinding: Extracting a Spiral from a Graph

**Module:** `Spiral.hs`

### Surface Orientation

The starting triple (f1, f2, f3) determines the surface orientation. The `orient`
function detects this and returns the appropriate face-finder:

```haskell
orient g f1 f2 f3
  | cNext g f1 f2 == Just f3 = Just (cPrev g)   -- CCW surface → use prev to find faces
  | cPrev g f1 f2 == Just f3 = Just (cNext g)   -- CW surface  → use next to find faces
  | otherwise                 = Nothing          -- f1,f2,f3 not a valid triangle
```

The intuition: if f3 follows f2 in f1's CCW neighbor list, the triple traverses their
shared face CCW. To find the *next* face to peel (on the opposite side of the boundary
edge), we look in the *reverse* rotational direction, hence `cPrev`.

### The Boundary

The algorithm maintains a deque (sequence) of `(node, remainingValencies)` pairs:

```haskell
data St = St !(Seq (Node, Int)) !IS.IntSet
```

Each entry represents a vertex already added to the spiral that still has unsaturated
connections. The `IntSet` tracks which vertices have been removed.

### Drain: The Core Primitive

The **drain** operation connects the newly peeled face to one end of the boundary,
cascading through any faces that become fully saturated:

```haskell
drainFront ⟨⟩           = (⟨⟩, 0)
drainFront ⟨(v,r), B'⟩  = if r > 1  then (⟨(v, r-1), B'⟩, 1)
                           else       let (B'', k) = drainFront B' in (B'', k+1)
```

**Cascade bound:** At each step, the total number of cascading pops is bounded by
`deg(v) - 3`, because the new face must retain at least 1 open valency:

```
pre_used = 1 (front drain) + front_cascades + 1 (back drain) + back_cascades
remaining = deg(v) - pre_used ≥ 1
∴ total_cascades ≤ deg(v) - 3
```

For fullerene duals (max degree 6): at most 3 cascading pops per step.

### Peel: One Spiral Step

```haskell
peel g findV (St bnd removed) = do
  (w, _) <- viewFront bnd              -- front of boundary
  (u, _) <- viewBack  bnd              -- back of boundary
  v      <- findV u w                  -- next face: prev(back, front) or next(back, front)
  guard (IS.notMember v removed)       -- must not already be removed
  let (b1, n1) = drainFront bnd       -- connect to front, cascade
      (b2, n2) = drainBack b1         -- connect to back, cascade
      remaining = deg g v - n1 - n2
  guard (remaining >= 1)               -- must have open valencies left
  pure (St (b2 |> (v, remaining)) (IS.insert v removed), deg g v)
```

### Full Regular Spiral

```haskell
regularSpiral :: Graph -> (Node, Node, Node) -> Maybe [Int]
```

1. **Orient** from the starting triple to get the face-finder function.
2. **Initialize** boundary with `[(f1, deg-2), (f2, deg-2), (f3, deg-2)]` (each face
   uses 2 connections for the starting triangle).
3. **Fold** `peel` over steps 3 through N-2, collecting face degrees.
4. **Finalize:** identify the last unremoved vertex, validate that the boundary has
   exactly `deg(lastV)` entries each with 1 remaining valency.
5. Return `[deg(f1), deg(f2), deg(f3)] ++ middle ++ [deg(lastV)]`.

Returns `Nothing` if any step fails (the regular spiral doesn't exist from this start).

### Canonical Spiral

```haskell
canonicalSpiral :: Graph -> Maybe [Int]
canonicalSpiral g = case mapMaybe (regularSpiral g) (startingTriples g) of
  []      -> Nothing
  spirals -> Just (minimum spirals)
```

Starting triples are restricted to the **rarest non-hexagonal degree** (for fullerenes:
the 12 pentagons). For each such vertex f1, try every neighbor f2 and both orientations
(f3 = cPrev or cNext of f2 w.r.t. f1). Take the lexicographic minimum of all successful
spirals.

Returns `Nothing` when no regular spiral exists from any starting triple. This first
occurs for fullerenes at Td-C100. Among ~2.7 x 10^12 isomers up to C400, the longest
needed jump sequence (in the generalized spiral, not implemented here) is 4.

## Windup: Constructing a Graph from a Spiral

**Module:** `Windup.hs`

The windup is the inverse of unwinding: given a face-degree sequence, construct the
oriented adjacency lists of the corresponding triangulation.

### Core Operation: `insEdge`

```haskell
insEdge :: Node -> Node -> Node -> Node -> Nbrs -> Nbrs
insEdge u v sucUV sucVU nbrs
```

Insert edge (u, v): place v in u's neighbor list immediately before `sucUV`, and place u
in v's neighbor list immediately before `sucVU`. If the "before" node isn't found, append
at end. This maintains the planar embedding's orientation.

### Algorithm

```haskell
windupSpiral :: [Int] -> Graph
```

1. **Initialize:** Create edge {0, 1}. Boundary = `[(0, spiral[0]-1), (1, spiral[1]-1)]`.

2. **For each face k = 2 to N-2:**
   - **connect_forward:** Insert edge {k, front}. In front's list, k goes before back.
     In k's list, front goes before back. Decrement front's open valency.
   - **connect_backward:** Insert edge {k, back}. In back's list, k goes before
     second-to-last. In k's list, back goes before front. Decrement back's valency.
   - **Cascade forward:** While front is saturated (valency = 0), pop it and connect k
     to the new front. The popped node becomes the "insert before" hint for k's list.
   - **Cascade backward:** Symmetric, from the back end.
   - **Push** k onto the boundary with `spiral[k] - connectionsUsed` remaining valencies.

3. **Close last face (N-1):** Connect it to all remaining boundary nodes (which should
   each have exactly 1 open valency), forming the final polygon.

### RSPI to Graph

```haskell
fromRSPI :: Int -> [Int] -> Graph
fromRSPI nAtoms rspi1 = windupSpiral spiral
  where
    nFaces  = nAtoms `div` 2 + 2
    pentSet = IS.fromList (map (subtract 1) rspi1)   -- convert 1-indexed to 0-indexed
    spiral  = [if IS.member i pentSet then 5 else 6 | i <- [0..nFaces-1]]
```

Example: `fromRSPI 20 [1..12]` builds the icosahedron (C20 dual, 12 vertices all degree 5).

## GPU Lockstep Analysis

Since jumps are rare (first at Td-C100), regular spirals can be computed for an entire
isomer space in lockstep with zero branch divergence:

1. **Fixed iteration count:** Main loop runs exactly N-3 times for every thread.
2. **Unrolled cascades:** The while-loops become nested ifs with depth ≤ 3 (for max
   degree 6). These compile to predicated instructions — no warp divergence.
3. **Failure masking:** Failed threads set `ok = false` and continue with frozen state.
4. **Bounded scans:** `prev(u, w)` scans at most D_max = 6 entries — same for all threads.
5. **No graph mutation:** Original adjacency lists are read-only. Mutable state is just
   the boundary buffer and a removed-node bitset.

### Canonical Spiral via Parallel Reduction

```
Phase 1: For each of 120 starting configs, run regular_spiral for all isomers in lockstep
Phase 2: Parallel reduction to find lexicographic minimum across successful spirals
Phase 3: Handle rare failures (Td-C100 onward) with generalized spiral or CPU fallback
```

## Generalized Spiral

The generalized spiral handles cases where the regular spiral fails (candidate face
is a cut vertex). It adds one mechanism: **boundary rotation** (jumps).

### Cut Vertex Test

Before peeling face v, check if removing v would disconnect the remaining graph.
For vertex v with remaining neighbors n_0, ..., n_{k-1} in cyclic order, count edges
between consecutive pairs (n_i, n_{i+1 mod k}). v is a cut vertex iff count < k-1.

The Haskell implementation avoids maintaining a mutable graph copy. Instead, it
filters v's original neighbors to exclude removed nodes, then checks edges in the
original graph (edges between non-removed vertices are unchanged by node removal).

### Jump Mechanics

When v is a cut vertex:
1. Rotate boundary deque: move front to back
2. Increment `jump_state` counter
3. Retry (do not advance the step counter)

When v is NOT a cut vertex and `jump_state > 0`:
1. Record `(step_index, jump_state)` in the jump list
2. Reset `jump_state = 0`
3. Peel v normally

### Canonical Ordering

```
(jumps1, spiral1) < (jumps2, spiral2) iff:
  |jumps1| < |jumps2|, OR
  |jumps1| == |jumps2| AND jumps1 < jumps2, OR
  jumps1 == jumps2 AND spiral1 < spiral2
```

Regular spirals (0 jumps) are always preferred. The canonical search first tries all
regular spirals; only if all fail does it try generalized spirals.

### Windup with Jumps

`windupGeneralSpiral spiral jumps` applies jump rotations during construction: at
step k, if there's a jump `(k, n)`, rotate the boundary by n positions before adding
node k.

## Test Results

The implementation passes all tests (`TestSpiral.hs`):

| Test | Vertices | Expected | Result |
|------|----------|----------|--------|
| Tetrahedron | 4 (all deg 3) | `[3,3,3,3]` | PASS |
| Octahedron | 6 (all deg 4) | `[4,4,4,4,4,4]` | PASS |
| Icosahedron | 12 (all deg 5) | `[5,5,...,5]` | PASS |
| C20 from RSPI | 12 (all deg 5) | `[5,5,...,5]` | PASS |
| D6d-C24 from RSPI | 14 (mixed 5/6) | Just _ | PASS |
| Td-C28 from RSPI | 16 (mixed 5/6) | Just _ | PASS |
| Td-C100 regular | 52 (mixed 5/6) | Nothing | PASS |
| Td-C100 general | 52 (mixed 5/6) | `jumps=[(42,2)]` | PASS |
| Td-C100 round-trip | 52 | wind up + re-extract | PASS |
| C20 general == regular | 12 | no jumps, same spiral | PASS |
| D6d-C24 general == regular | 14 | no jumps, same spiral | PASS |

All 120/120 starting triples succeed for C20, D6d-C24, and Td-C100 (general).

## API Summary

### Spiral.hs (unwinding)

| Function | Type | Purpose |
|----------|------|---------|
| `mkGraph` | `[[Node]] -> Graph` | Construct from adjacency lists |
| `regularSpiral` | `Graph -> (Node,Node,Node) -> Maybe [Int]` | Extract regular spiral from one start |
| `canonicalSpiral` | `Graph -> Maybe [Int]` | Canonical regular spiral |
| `generalSpiral` | `Graph -> (Node,Node,Node) -> Maybe GeneralSpiral` | Extract generalized spiral from one start |
| `canonicalGeneralSpiral` | `Graph -> Maybe GeneralSpiral` | Canonical spiral (regular preferred, general fallback) |
| `startingTriples` | `Graph -> [(Node,Node,Node)]` | All valid starting configurations |

### Windup.hs (construction)

| Function | Type | Purpose |
|----------|------|---------|
| `windupSpiral` | `[Int] -> Graph` | Build graph from face-degree sequence |
| `windupGeneralSpiral` | `[Int] -> [(Int,Int)] -> Graph` | Build graph from spiral + jumps |
| `fromRSPI` | `Int -> [Int] -> Graph` | Build fullerene dual from RSPI |

## Dependencies

Only the `containers` and `vector` packages (beyond `base`). Compiles with GHC 8.8+:

```bash
cabal install --lib vector
ghc -O2 Spiral.hs Windup.hs TestSpiral.hs -o test-spiral
./test-spiral
```
