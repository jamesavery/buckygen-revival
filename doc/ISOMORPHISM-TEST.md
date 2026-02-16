# Fullerene Isomorphism Testing via Canonical Spirals

This document explains how to use the Haskell spiral implementation (`Spiral.hs` +
`Windup.hs`) for fast isomorphism equivalence testing of fullerene graphs. It is written
for a Claude session implementing buckygen in Haskell.

## Quick Start

### Prerequisites

```bash
# GHC and the vector package
apt install ghc cabal-install
cabal update && cabal install --lib vector
```

### Build and test

```bash
ghc -O2 Spiral.hs Windup.hs TestSpiral.hs -o test-spiral
./test-spiral
```

Expected output: all PASS, including C20, D6d-C24, Td-C28, and the Td-C100 negative test.

## How Isomorphism Testing Works

Two fullerene graphs are isomorphic **iff** their canonical spirals are equal.

The canonical spiral is the lexicographically smallest face-degree sequence obtainable
by peeling faces from the dual triangulation, over all valid starting configurations.

### The Algorithm in One Paragraph

Take the dual of a cubic fullerene graph to get a triangulation (pentagons become
degree-5 vertices, hexagons become degree-6). Pick a starting triple of three mutually
adjacent vertices. Peel faces one at a time from the boundary, recording each vertex's
degree. The result is a sequence like `[5,5,5,5,5,6,5,5,6,5,5,5,5,5]`. Try all valid
starting triples (up to 120 for fullerenes), take the lexicographic minimum. This is the
canonical spiral — a complete isomorphism invariant.

## API for Isomorphism Testing

### Core types

```haskell
import Spiral (Graph, Node, mkGraph, canonicalSpiral, regularSpiral, startingTriples)
import Windup (windupSpiral, fromRSPI)

type Node = Int
newtype Graph = Graph { adj :: Vector (Vector Node) }  -- oriented adjacency lists
```

### Testing two graphs for isomorphism

```haskell
isIsomorphic :: Graph -> Graph -> Bool
isIsomorphic g1 g2 = canonicalSpiral g1 == canonicalSpiral g2
```

`canonicalSpiral` returns `Maybe [Int]`. Two graphs are isomorphic iff both return the
same `Just spiral`. (If both return `Nothing`, regular spirals failed for both — this
means the generalized spiral algorithm is needed, which is not implemented here. First
occurrence: Td-C100.)

### Constructing a graph from RSPI

If you have a fullerene specified by its RSPI (Ring-Spiral Pentagon Indices, 1-indexed):

```haskell
-- C60 Buckminsterfullerene (Ih symmetry)
-- RSPI: positions of the 12 pentagons in the canonical spiral
let g = fromRSPI 60 [1,7,9,11,13,20,22,24,26,32,34,36]
```

Arguments: `fromRSPI nAtoms rspi` where:
- `nAtoms` = number of carbon atoms (20, 24, 28, ..., 60, ...)
- `rspi` = list of 12 integers, 1-indexed positions where the spiral has a pentagon

The function builds the full face-degree sequence (`nAtoms/2 + 2` entries of 5 or 6)
and constructs the oriented triangulation via the windup algorithm.

### Constructing a graph from oriented adjacency lists

If your buckygen implementation produces oriented adjacency lists directly:

```haskell
-- Each inner list is the neighbors of vertex i in consistent cyclic order
let g = mkGraph
      [ [1, 2, 3, 4, 5]       -- vertex 0's neighbors, CCW
      , [2, 0, 5, 10, 6]      -- vertex 1's neighbors, CCW
      , ...
      ]
```

**Important orientation rule:** For a CCW-oriented triangulation with triangle (a, b, c)
traversed counter-clockwise, `adj[a]` must contain `...b, c, ...` in that order (b is
the predecessor of c in a's cyclic neighbor list).

### Getting just the spiral (without canonicalization)

```haskell
-- Try one specific starting triple
regularSpiral g (f1, f2, f3) :: Maybe [Int]

-- Get all valid starting triples for a graph
startingTriples g :: [(Node, Node, Node)]
```

## Integrating with Buckygen

### If buckygen produces triangulations directly

Pass them directly to `canonicalGeneralSpiral`, which always succeeds (tries regular
first, falls back to generalized):

```haskell
let g = mkGraph adjacencyLists
case canonicalGeneralSpiral g of
  Just gs -> ...  -- gs :: GeneralSpiral, with gsJumps and gsSpiral fields
  Nothing -> ...  -- should never happen for valid 3-connected planar graphs
```

Or if you only need regular spirals (sufficient for fullerenes up to C98):

```haskell
case canonicalSpiral g of
  Just spiral -> ...  -- spiral :: [Int], the face-degree sequence
  Nothing     -> ...  -- need generalized spiral (rare: first at C100)
```

### Comparing isomers during generation

A typical buckygen workflow:
1. Generate candidate triangulations
2. For each, compute `canonicalGeneralSpiral`
3. Use the spiral as a hash key to detect and reject duplicates

```haskell
import qualified Data.Set as Set

deduplicate :: [Graph] -> [Graph]
deduplicate = go Set.empty
  where
    go _    []     = []
    go seen (g:gs) =
      case canonicalGeneralSpiral g of
        Nothing -> g : go seen gs  -- shouldn't happen for valid graphs
        Just gs'
          | Set.member gs' seen -> go seen gs       -- duplicate, skip
          | otherwise           -> g : go (Set.insert gs' seen) gs
```

## Performance Characteristics

| Fullerene | Dual vertices | Starting triples | Time per canonicalization |
|-----------|---------------|------------------|--------------------------|
| C20       | 12            | 120              | ~microseconds            |
| C60       | 32            | 120              | ~microseconds            |
| C100      | 52            | 120              | ~microseconds            |
| C200      | 102           | 120              | ~tens of microseconds     |

The number of starting triples is always ≤ 120 for fullerenes (12 pentagons x 5
neighbors x 2 orientations). Each `regularSpiral` call is O(N) where N = number of
dual vertices. Total: O(N) per isomorphism test.

## Limitations

1. **No cubic-to-dual conversion.** You must provide the dual triangulation. If your
   buckygen works with cubic graphs, you need to dualize first.

2. **Assumes valid oriented triangulation.** The input must have consistently oriented
   cyclic neighbor lists forming a valid triangulation. Invalid input may produce wrong
   results or crash.

## File Inventory

| File | Purpose |
|------|---------|
| `Spiral.hs` | Core module: unwinding, canonical spiral computation |
| `Windup.hs` | Inverse: construct triangulation from spiral or RSPI |
| `TestSpiral.hs` | Test suite: tetrahedron, octahedron, icosahedron, C20-C100 |
| `SPIRAL-ALGORITHM.md` | Detailed algorithm description with proofs |
| `ISOMORPHISM-TEST.md` | This file |

## Worked Example

```haskell
import Spiral
import Windup

-- Two isomers of C28: Td symmetry and D2 symmetry
let td_c28 = fromRSPI 28 [1,2,3,5,7,9,10,11,12,13,14,15]
let d2_c28 = fromRSPI 28 [1,2,3,4,5,6,8,12,13,14,15,16]

-- Compute canonical spirals
canonicalSpiral td_c28  -- Just [5,5,5,6,5,6,5,6,5,5,5,5,5,5,5,6]
canonicalSpiral d2_c28  -- Just [5,5,5,5,5,5,6,6,6,6,5,5,5,5,5,5]

-- Different spirals → not isomorphic
canonicalSpiral td_c28 == canonicalSpiral d2_c28  -- False

-- Same RSPI twice → isomorphic (trivially)
let g1 = fromRSPI 24 [1,2,3,4,5,7,8,10,11,12,13,14]
let g2 = fromRSPI 24 [1,2,3,4,5,7,8,10,11,12,13,14]
canonicalSpiral g1 == canonicalSpiral g2  -- True
```
