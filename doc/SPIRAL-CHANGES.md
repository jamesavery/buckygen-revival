# Spiral Haskell: Changes Since Last Release

## Bug Fix: Windup Cascade Collision

**Problem**: `windupSpiral` produced invalid graphs for most spirals with 16+
vertices. The cascade functions (`cascFwd`/`cascBwd`) could eat past the
opposite end of the boundary deque when both ends reached remaining=0
simultaneously, creating duplicate edges and corrupted adjacency lists.

Smallest failing case: C28 D2, spiral `[5,5,5,5,5,5,6,5,6,6,6,5,5,5,5,5]`.
At step k=14, boundary `[9(0), 10(1), 11(1), 13(0)]` — `cascFwd` would pop
through 9, 10, 11 and connect k to 13 a second time (13 was already connected
by the initial connect-back), producing a duplicate edge (14,13).

**Fix**: Added two stopping guards to both `cascFwd` and `cascBwd`, matching
the C++ reference in `triangulation.cc`:

```haskell
cascFwd k nbs ov pu
  | Seq.length ov <= 1       = (nbs, ov, pu)  -- only back node left
  | snd (seqHead ov) /= 0   = (nbs, ov, pu)  -- front not saturated
  | Seq.length ov == 2
    && snd (seqLast ov) == 0 = (nbs, ov, pu)  -- both saturated, don't cross
  | otherwise = ...
```

The symmetric guard in `cascBwd` checks `snd (seqHead ov) == 0`.

**Verified**: 17 expansion-generated graphs (C20-C40) all pass
spiral-extract -> windup -> re-extract round-trips.

## Efficiency: CSR Graph Representation

Replaced the old `Vector (Vector Node)` graph representation with CSR
(Compressed Sparse Row) format:

```haskell
-- Old: one heap-allocated boxed vector per vertex
data Graph = Graph { adj :: Vector (Vector Node) }

-- New: two flat unboxed arrays, zero per-vertex allocation
data Graph = Graph
  { gAdj    :: !(VU.Vector Node)   -- all neighbors, concatenated
  , gOffset :: !(VU.Vector Int)    -- start offset per vertex; length N+1
  }
```

Benefits:
- Eliminates N pointer-chasing indirections and N small heap objects.
- `nbrSlice g v` returns an O(1) zero-copy slice into `gAdj`.
- Same data layout as GPU-side CSR, so CPU-GPU transfer is a memcpy.
- All navigation (`next`, `prev`, `edge`, `deg`) works on flat unboxed arrays.

The `mkGraph` constructor is unchanged: `mkGraph [[1,2,3], [0,3,2], ...]`.

## Readability: Cleaner Abstractions

### `(!)` operator
```haskell
(!) :: VU.Unbox a => VU.Vector a -> Int -> a
(!) = VU.unsafeIndex
```
Replaces verbose `VU.unsafeIndex arr i` with `arr ! i` throughout.

### Orientation and apex
```haskell
data Orientation = CCW | CW

orient :: Graph -> Node -> Node -> Node -> Maybe Orientation
apex   :: Graph -> Orientation -> Node -> Node -> Maybe Node
```
The starting triple's winding determines `CCW` or `CW`. Then `apex g ori u w`
finds the triangle vertex across edge (u,w) — it's `prev g u w` for CCW or
`next g u w` for CW. This replaces an anonymous closure that was hard to read.

### Renamed navigation
- `cNext` -> `next`, `cPrev` -> `prev`, `isEdge` -> `edge`
- `viewFront`/`viewBack` -> `front`/`back`

### Windup helpers
- `decFront`/`decBack` replace `decAt 0`/`decAt (len-1)` — clearer intent.
- Windup uses `VU.Vector Int` (unboxed) for the spiral array instead of boxed.

## Architecture: Separate Fast and Slow Paths

Regular and generalized spiral extraction are deliberately kept separate:

- **`peel`** (fast path): No cut-vertex check. Used by `regularSpiral`. This is
  the hot path — sufficient for all fullerenes up to C98.
- **`generalPeelStep`** (slow path): Calls `isCutVertex` on every step (O(deg^2)).
  Used by `generalSpiral`. Only needed starting at Td-C100 (exceedingly rare).
- **`canonicalGeneralSpiral`**: Two-phase — tries all 120 regular spirals first;
  only if ALL fail does it try generalized. This means the common case never
  pays the cut-vertex cost.

## File Inventory

| File | Purpose |
|------|---------|
| `Spiral.hs` | Core module: CSR graph, regular + generalized spiral extraction |
| `Windup.hs` | Inverse: construct graph from spiral (with cascade bug fix) |
| `TestSpiral.hs` | Test suite: hand-crafted graphs, RSPI round-trips, WindupTestdata round-trips |
| `WindupTestdata.hs` | 17 known-correct graphs (C20-C40) with canonical spirals |
| `SPIRAL-CHANGES.md` | This file |

## API Quick Reference

```haskell
-- Graph construction
mkGraph :: [[Node]] -> Graph

-- Navigation
deg  :: Graph -> Node -> Int
nbrs :: Graph -> Node -> [Node]
next :: Graph -> Node -> Node -> Maybe Node   -- cyclic successor
prev :: Graph -> Node -> Node -> Maybe Node   -- cyclic predecessor
edge :: Graph -> Node -> Node -> Bool

-- Orientation
orient :: Graph -> Node -> Node -> Node -> Maybe Orientation
apex   :: Graph -> Orientation -> Node -> Node -> Maybe Node

-- Spiral extraction
regularSpiral          :: Graph -> (Node, Node, Node) -> Maybe [Int]
canonicalSpiral        :: Graph -> Maybe [Int]
generalSpiral          :: Graph -> (Node, Node, Node) -> Maybe GeneralSpiral
canonicalGeneralSpiral :: Graph -> Maybe GeneralSpiral
startingTriples        :: Graph -> [(Node, Node, Node)]

-- Windup (inverse)
windupSpiral        :: [Int] -> Graph
windupGeneralSpiral :: [Int] -> [(Int, Int)] -> Graph
fromRSPI            :: Int -> [Int] -> Graph   -- from atom count + 1-indexed RSPI
```
