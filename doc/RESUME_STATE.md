# Resume State — Step 2 COMPLETE, Test Data Generated

## What's done
- All expansion/reduction operations implemented and tested (L0, Straight, BentZero, Bent, Ring)
- 461 tests pass across 8 test sections (TestExpansion.hs), 0 failures
- Generation verified through C60: all 5770 fullerene isomers match known counts, 0 spiral failures
- Cross-check (spiral dedup vs backtracking isomorphism) passes through C40
- Test data generated: TestGraphs60.hs with all 5770 graphs + canonical generalized spirals
- 10 bugs found and fixed across 10 sessions

## Files
| File | Lines | Purpose |
|------|-------|---------|
| `Seeds.hs` | 920 | Step 1: 49 seed graphs + EdgeList type, initEdgeList, 4-field DualGraph |
| `Expansion.hs` | 1320 | All expansion/reduction operations including F (nanotube ring) |
| `Spiral.hs` | 297 | Canonical generalized spiral computation |
| `TestExpansion.hs` | 509 | Comprehensive test suite (461 tests) |
| `Generate.hs` | 252 | Generation sanity check with isomorphism deduplication |
| `GenTestData60.hs` | 139 | Test data generator (all graphs up to C60) |
| `TestGraphs60.hs` | 5790 | Generated test data (5770 graphs with canonical spirals, ~3.4 MB) |

## Key implementation details

### findNanotubeRing (Expansion.hs)
Uses `straightAhead` tracing (matching buckygen's `has_lm_path`) to find the (5,0)
nanotube equatorial ring. Tries all degree-6 vertex pairs as starting edges.
`computeOuter` finds common neighbors of consecutive ring vertices (no degree filter),
resolves a consistent cap side via adjacency (`resolveOuters`), and validates with
`outersAdjacent` to prevent false positives on non-nanotube graphs.

### applyRing (Expansion.hs)
Adds 5 new degree-6 vertices around the equatorial ring. Each new vertex i has
neighbors: `[ring[i], new[(i-1)%5], outer[i], outer[(i+1)%5], new[(i+1)%5], ring[(i+1)%5]]`.
Ring and outer vertices get their neighbor lists patched via `replaceNbrEL`.

### Expansion enumeration
- L0: bidirectional (both directions for each degree-5 pair, 120 on C20)
- L_i (i>=1): path of length i+1 between two degree-5 vertices
- B_{0,0}: included (was previously excluded by mistake)
- B_{i,j} (i+j>0): bent path with a turn
- F (Ring): detected via `findNanotubeRing`, only for (5,0) nanotube family from C30

### Spiral deduplication
`Spiral.canonicalGeneralSpiral` computes a canonical generalized spiral (with jumps
for non-Hamiltonian cases). Used as the isomorphism invariant for deduplication.
All 5770 graphs through C60 have successful spiral computations (0 failures).

## Generation counts (all match known isomer numbers)
```
C20=1, C24=1, C26=1, C28=2, C30=3, C32=6, C34=6, C36=15, C38=17, C40=40,
C42=45, C44=89, C46=116, C48=199, C50=271, C52=437, C54=580, C56=924,
C58=1205, C60=1812
Total: 5770
```

## Next step: Step 3 — Canonicity test

The 6-tuple `(x0, x1, x2, x3, x4, x5)` cascade with lazy evaluation and early exit.
- x0-x4 are O(1) combinatorial invariants resolving >99.9% of cases
- x5 (full BFS canonical form) is O(n), needed <0.1% of the time
- Automorphism group falls out as a byproduct of x5

Reference: `bucky_pseudocode.tex` for the formal spec, `buckygen.c` for the C implementation.

## How to resume
1. Read this file and CLAUDE.md for project context
2. Read PROGRESS.md for bug history and architecture details
3. For Step 3: read `bucky_pseudocode.tex` for canonicity test spec
4. Implement canonicity test in a new module (e.g. `Canonical.hs`)
5. Test against TestGraphs60.hs data
