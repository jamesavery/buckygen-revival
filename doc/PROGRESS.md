# Progress — Buckygen Rewrite

## Status: C60 generation CORRECT — Parallel evaluation WORKING — Search module is sole API

## Latest: EdgeList removed from DualGraph (2026-02-18)

### EdgeList removal
- **Removed `EdgeList` type** (`IntMap (IntMap Int)`) from `DualGraph` — it mapped
  `(vertex, neighbor) → position` in the CW adjacency list, used only during expansion
  operations. The C code's `edge_list[MAXN][MAXN]` is a flat O(1) 2D array; our IntMap
  equivalent was O(log n) with pointer chasing — slower than the existing `indexOf`
  which scans 5-6 contiguous UArray elements.
- **Seeds.hs**: Removed `EdgeList` type alias, `initEdgeList`, `edgeList` field from
  `DualGraph` record, merged `mkDualGraph`/`mkDualGraphLite` into single `mkDualGraph`.
  DualGraph is now 5 fields: `nv`, `neighbours`, `degree5`, `adjFlat`, `degFlat`.
- **Expansion.hs**: Deleted 3 EL helper functions (`replaceNbrEL`, `registerVertex`,
  `insertAfterEL`). All 5 apply functions (applyL0, applyStraight, applyBentZero,
  applyBent, applyRing) now use the existing `replaceNbr`/`insertAfter` helpers.
  Eliminated all `(adj, el)` tuple threading — just plain `adj` threading.
  All 4 reduce functions changed from `mkDualGraphLite` to `mkDualGraph`.
- **MutGraph.hs**: Removed `EdgeList` import, removed extra `IM.empty` from DG constructors.
- **Generate.hs**: Removed `initEdgeList` import.
- **Search.hs**: Updated NFData instance (5 fields instead of 6).
- **All tests pass**: expansion round-trips, canonical generation through C60 (5770 isomers),
  all 6 demo cross-checks, F expansion tests.

### Benchmark: before/after EdgeList removal (C80, 42 dv, 131,200 isomers, M2 Ultra, 10 cores)

```
Scheme                        Before    After     Change
Sequential DFS                 59.8s    52.8s     -12%
Sequential BFS                 59.7s    52.6s     -12%
Parallel Pure (d=9)              7.2s     6.4s     -11%
Parallel MutGraph (d=9)          5.2s     5.2s       0%
Work-queue Pure (10w, d=9)       8.5s     7.7s      -9%
Work-queue MutGraph (10w, d=9)   6.6s     6.4s      -3%
Hierarchical MutGraph (d1=4)     5.0s     5.2s      +4%
BFS+DFS MutGraph (d=7)           5.8s     6.7s     +16%
```

Pure backends show a consistent ~11% speedup from removing EdgeList allocation/copy
overhead on every DualGraph. MutGraph backends are flat (within noise) because they
bypass `mkDualGraph` for the 97% of candidates that fail canonical testing via
zero-copy `unsafeFreeze` — they never paid the EdgeList cost.

## GenForest retired, Search is the single module (2026-02-18)

### GenForest → Search migration
- **`GenForest.hs` retired** to `obsolete/GenForest.hs` via `git mv`
- **`Search.hs` is now self-contained** (~680 lines): absorbs all core definitions
  (GenConfig, mkConfig, GenTree, genForest, seedsWithAuts, generateChildren,
  NFData instances, depthProfile, subtreeSizes, allGraphs, graphsOfSize, bfsByLevel)
  plus the Fold/Schedule/search abstraction and all 8 schedule backends
- **All consumers rewritten** to import only Search:
  - `DemoForest.hs` — all 12 items use Search API (was 10 GenForest + 6 Search)
  - `BenchPar.hs` — all modes use `search countBySize <Schedule> cfg`
  - `BenchSearch.hs` — pure Search benchmark (removed old/new comparison)
  - `TestCanonical.hs` — forest cross-check uses `Search.searchPure`
- **Zero functional changes**: all 6 cross-checks pass at C60 (5770 isomers, 1812 at C60)
- **Slightly faster** than old GenForest+Search: SeqDFS 59.8s vs 67.4s at C80, 132.5s vs 143.2s at C86

### Benchmark: self-contained Search module (M2 Ultra, GHC 9.6.7 -O2, 10 cores)

**C80 (131,200 isomers)**
```
Scheme                              Time        Speedup
Sequential DFS                      59.8s         1.0x
Sequential BFS                      59.7s         1.0x
Parallel Pure (d=9)                  7.2s         8.3x
Parallel MutGraph (d=9)              5.2s        11.5x
Work-queue Pure (10w, d=9)           8.5s         7.0x
Work-queue MutGraph (10w, d=9)       6.6s         9.1x
Hierarchical MutGraph (d1=4,d2=5)    5.0s        12.0x
BFS+DFS MutGraph (d=7)               5.8s        10.3x
```

**C86 (286,272 isomers)**
```
Scheme                              Time        Speedup
Sequential DFS                     132.5s         1.0x
Sequential BFS                     132.0s         1.0x
Parallel Pure (d=9)                 16.5s         8.0x
Parallel MutGraph (d=9)             12.1s        11.0x
Work-queue Pure (10w, d=9)          17.2s         7.7x
Work-queue MutGraph (10w, d=9)      12.7s        10.4x
Hierarchical MutGraph (d1=4,d2=5)   11.4s        11.6x
BFS+DFS MutGraph (d=7)              13.9s         9.5x
```

HierMut is fastest at both sizes (12.0x and 11.6x on 10 cores).

### Files changed
- **Moved**: `src/GenForest.hs` → `obsolete/GenForest.hs`
- **Rewritten**: `src/Search.hs` (self-contained, ~680 lines)
- **Rewritten**: `src/DemoForest.hs` (Search-only)
- **Rewritten**: `src/BenchPar.hs` (Search-only)
- **Rewritten**: `src/BenchSearch.hs` (Search-only benchmark)
- **Modified**: `src/TestCanonical.hs` (GenForest → Search)
- **Modified**: `buckygen-revival.cabal` (removed GenForest from all module lists)

### Previous: Search module created (2026-02-18)

#### SearchMonad abstraction: `Search.hs`
- **`Fold r`** record: `foldVisit`, `foldMerge`, `foldEmpty` — captures what to collect
- **`Schedule`** sum type: `SeqDFS | SeqBFS | ParFlat | ParMut | WorkQ | WorkQMut | HierMut | BfsDfs`
- **`search`**: single entry point dispatching to per-schedule backends
- **`searchPure`**: pure sequential DFS (no IO)
- **Predefined folds**: `countBySize`, `collectAll`, `collectOfSize`, `mapNodes`
- **`SearchM` monad**: `yield` + `searchForest` + `runSearchM` for custom pure DFS
- **`mutFoldSubtree`**: generalized `mutCountSubtree` parameterized by `Fold r`
- All 8 schedule backends implemented, all cross-checked against DFS at C60
- Custom fold (max |Aut| per size) matches `foldForest` result
- `SearchM` monad produces identical counts
- Performance regression test at C80 and C86: all results match, zero overhead

## Previous: 10 traversal schemes + work-queue parallel (2026-02-17)

### Work-queue parallel generation
- `workQueueCountBySize`: TQueue + forkIO workers + per-worker IORef accumulators
- Pre-splits tree at fork depth (like flat-fork), puts chunks in TQueue
- Workers pull chunks and DFS each one — automatic load balancing
- Matches pure flat-fork performance: C80 N=8: 9.5s (vs 9.6s flat-fork)
- Initial per-graph STM implementation deadlocked (bounded TBQueue) and was 10x slow;
  redesigned to pre-split + coarse-grained DFS per chunk

### 10 forest traversal schemes (all cross-checked at C60)
1. **Sequential DFS** (`dfsCountBySize`) — baseline
2. **BFS level-by-level** (`bfsByLevel`, `bfsCountBySize`)
3. **Generic fold** (`foldForest`) — arbitrary accumulator
4. **Lazy DFS stream** (`allGraphs`) — on-demand [(graph, auts)]
5. **Parallel flat-fork pure** (`parCountBySize`)
6. **Parallel flat-fork MutGraph** (`parMutCountBySize`) — fastest
7. **Target-size generation** (`graphsOfSize`) — pruned for specific C_n
8. **Parallel map** (`parMapForest`) — apply function to every graph
9. **Work-queue** (`workQueueCountBySize`) — dynamic load balance
10. **DemoForest.hs** exercises all schemes and cross-checks: ALL PASS

### Hybrid parallel: pure forest at top, MutGraph DFS within each subtree
- `parMutCountBySize`: uses pure `generateChildren` to force tree to fork depth,
  collects subtree roots, spawns independent `runST` blocks with fresh MutGraph
  per subtree, evaluates via `parList rdeepseq`
- Combines best of both worlds: pure forest for parallelism (zero shared state),
  MutGraph for sequential speed within each subtree (in-place surgery + unsafeFreeze)
- 1.45x faster than pure forest at any core count; 25% less allocation

**Benchmark results — C80 (131,200 isomers), all 3 parallel backends, N=8:**
```
Backend                 Elapsed     CPU util
Pure flat-fork d=6       9.5s        612%
MutGraph hybrid d=6      6.7s        601%      (fastest)
Work-queue N=8           9.5s        607%
```

**Pure forest vs MutGraph hybrid (both N=8 d=6):**
```
                C60         C70         C80
Pure forest     0.55s       —           9.5s
MutGraph        0.45s       —           6.7s       (1.45x faster)
```

### Flat-fork parallel strategy
- `parCountBySize` / `parMutCountBySize` use flat-fork: force tree to fork depth,
  collect all subtrees into a single flat list, spark each subtree's DFS count
- Eliminates recursive spark fizzling (229/230 sparks converted vs 83/416 with recursive approach)
- Optimal fork depth: 6 (230 subtrees) for up to ~10 cores
- Scaling improves with problem size (more work per subtree → lower overhead ratio)
- N>8 gives diminishing returns (bottleneck: work imbalance, not GC)

### GenForest: explicit lazy search forest for parallel evaluation
- Module `GenForest.hs` (~530 lines): decouples tree structure from traversal strategy
- `generateChildren` — pure function computing all canonical children
- `genForest` — lazy unfold from 3 seeds
- 10 traversal schemes (see above)
- Each subtree is independent — purity guarantees no cross-tree communication

### Selective unsafeFreezeGraph for rejected children
- 97% of children fail `isCanonical` but previously each paid the O(7n) `freeze` cost
- Added `unsafeFreezeGraph`: zero-copy alias of MutGraph's mutable arrays as immutable UArrays
- Safe because `isCanonicalV` returns strict `CanonVerdict` (nullary constructors) — all graph
  reads complete before `undoMutation` modifies the shared arrays
- Only the 3% of children that pass `isCanonical` use safe `freezeGraph` (for recursion)
- C60: 1.45s → 1.35s (7% speedup), 7.7GB → 7.5GB. C70: 7.73s → 7.06s (9% speedup), 43.0 → 41.6GB.

## Previous: Comprehensive instrumentation + edge colour pre-filter (2026-02-17)

### Edge colour pre-filter for canonicalBFSAndGroup
- Computes "first-row colour" for each BFS starting edge: `deg(u)*32 + degree_pattern`
- Only edges with the minimum colour can produce the minimum full BFS code
- Filter effectiveness: 58x reduction at C60 (335→5.8 avg), 69x at C70 (390→5.6 avg)
- Wall-clock improvement: ~5.6% (modest because interleaved BFS comparison already aborted cheaply)

### Performance instrumentation
- `CanonResult` now includes BFS stats: edges total/filtered, GT/EQ/LT comparison counts
- `CanonVerdict` type for `isCanonicalV`: tracks rejection stage (L0Early/NoInverse/Phase1/Phase2)
- `SizeStats` expanded with 11 new counters for complete generation pipeline analysis
- All stats displayed in `printStats` with canonical test breakdown + BFS group analysis

### Critical fix: `freezeGraph` switched from `unsafeFreeze` to safe `freeze`
- `unsafeFreeze` caused data corruption when frozen graph arrays (aliasing MutGraph) were
  read after `undoMutation` had restored parent state. The BFS stats fields in `CanonResult`
  could be lazily evaluated after mutations, reading corrupted data.
- Safe `freeze` copies the active array portion (~200 bytes per graph), adding only 3% overhead.
- C60: 1.45s / 7.7GB (was 1.4s / 7.5GB). C70: 7.73s / 43.0GB (was 7.5s / 41.4GB).

## Previous: BUCKYGEN-ALGORITHM.md written (2026-02-16)
- Self-contained 954-line algorithm specification in `doc/BUCKYGEN-ALGORITHM.md`
- Fills all 15 gaps identified in `doc/BUCKYPAPER.md` (paper analysis)
- Designed to be implementable from scratch without any other reference
- Covers: graph representation, seeds, expansions, reductions, canonical test (5-tuple),
  Rule 2 orbit filtering, generation loop, bounding lemmas, worked example, validation data

## Reference Documents
- `BUCKYGEN-ALGORITHM.md` — **COMPLETE SPECIFICATION**: self-contained algorithm description
- `C-CODE-ALGORITHM.md` — algorithm as reverse-engineered from C code
- `BUCKYGEN-MAP.md` — mapping between paper, C code, and Haskell (function index, line numbers)
- `BUCKYPAPER.md` — analysis of 15 gaps/errors in the original buckygen paper
- `buckygen_paper_source/` — original paper (has significant gaps; see BUCKYPAPER.md)

## Files
- `Seeds.hs` (954 lines) — DualGraph type, planar code decoder, 49 seed graphs
- `Expansion.hs` (1354 lines) — Pure graph surgery: 5 expansion types × apply + reduce + enumerate
- `Canonical.hs` (1206 lines) — 5-tuple canonical test, BFS, automorphisms, Rule 2, bounding lemmas
- `MutGraph.hs` (551 lines) — Mutable graph for sequential DFS optimization (apply/freeze/undo)
- `Search.hs` (~680 lines) — **Unified search**: Fold + Schedule + search, generation forest, all 8 backends, SearchM monad, convenience traversals, analysis
- `DemoForest.hs` (~120 lines) — Exercises all traversal schemes with cross-checks (Search-only)
- `BenchPar.hs` (~75 lines) — Parallel benchmark: `pure`, `mut`, `wq`, `wqmut`, `hier`, `bfsdfs`, `analyze` modes (Search-only)
- `BenchSearch.hs` (~100 lines) — Benchmark exercising all 8 schedule backends with timing
- `Spiral.hs` (297 lines) — Canonical generalized spiral computation
- `TestCanonical.hs` (555 lines) — Generation tests with MutGraph + pure forest cross-check
- `TestExpansion.hs` (509 lines) — Expansion round-trip tests
- `Generate.hs` (252 lines) — Brute-force generation with isomorphism dedup
- `GenTestData60.hs` (139 lines) — Test data generator (all graphs up to C60)
- `TestGraphs60.hs` (5790 lines) — Generated test data (5770 graphs with canonical spirals)

## Generation Counts (ALL 21 sizes pass through C60)

```
Dual V  C_n     Found       Expected    Status
----------------------------------------------------
12      C20     1           1           PASS
14      C24     1           1           PASS
15      C26     1           1           PASS
16      C28     2           2           PASS
17      C30     3           3           PASS
18      C32     6           6           PASS
19      C34     6           6           PASS
20      C36     15          15          PASS
21      C38     17          17          PASS
22      C40     40          40          PASS
23      C42     45          45          PASS
24      C44     89          89          PASS
25      C46     116         116         PASS
26      C48     199         199         PASS
27      C50     271         271         PASS
28      C52     437         437         PASS
29      C54     580         580         PASS
30      C56     924         924         PASS
31      C58     1205        1205        PASS
32      C60     1812        1812        PASS
```

## Current Architecture

### Canonical test (Rule 1): 5-tuple cascade
```
(x0, x1, x2, x3, x4) = (reductionLength, -longestStraight, colourPair, pathColour, BFS)
```

### Automorphism group: `canonicalBFSAndGroup`
- Tries all (vertex, neighbor, direction) starting edges
- Edges producing the minimum BFS code define automorphisms
- Returns `CanonResult` with code, automorphism list, and nbop count

### Rule 2: `filterByRule2`
- Partitions expansion sites into equivalence classes
- Each class generated by: automorphism images + other-triple equivalence
- Eliminates ~45% of expansion sites

### Rule 1 (isCanonical): direction-specific inverse test
- Identifies the inverse reduction by: same kind, both vertices new, same direction
- The direction-specific test (Bug #15) is the key insight: the C code tests the
  expansion's edge in its SPECIFIC direction, rejecting if the other direction is better
- This single check eliminates all duplicates — no additional Rule 2 orbit extensions needed

### Generation loop (TestCanonical.hs)
1. Compute automorphism group for seeds
2. For each parent: enumerate expansions, filter by Rule 2, apply canonical test
3. Compute child's automorphism group for recursion
4. Spiral dedup as safety net (currently catches 0 duplicates)

## Bugs Fixed

### Bug #15: Direction-specific canonical test (session 18)
The C code's `is_best_L0_reduction` tests the new edge in the SPECIFIC direction used
by the expansion. If the other direction has a better canonOrd, the expansion is rejected.
Our `isInverseRed` was matching ALL directions at new vertices.

**Discovery**: C28's `(0,1) DLeft` and `(0,1) DRight` both passed isCanonical. C code
instrumentation showed `L0 colour_reject: nv=18 from edge 4-0 use_next=1`.

**Fix**: Added `d == dir` to `isInverseRed`, extracting `dir` from the expansion pattern.

**Effect**: Dedup catches 4,268 → 0. This single fix subsumes the L0 reverse equivalence
issue (Bug #16 attempted) — the direction-specific test already rejects the duplicate direction.

### Bug #14: isInverseRed too loose — required both vertices new (session 17)
Previously only checked that the start vertex was new. A new degree-5 vertex adjacent to
an existing degree-5 vertex created spurious L0 reduction matches.

**Fix**: Require both `a` and `b` of the reduction's edge to be within the new vertex range.

**Effect**: Dedup catches 10,109 → 4,268.

### Bug #13: canBentPath incorrectly required par[last] degree-5 (session 16)
**Fix**: Removed `deg g qLast == 5` from `canBentPath`. Now 1812/1812 C60 isomers found.

### Bug #12: L_i reduction enumeration used expansion-site criteria (session 12)
### Bug #11: L0 expansion enumeration too restrictive (session 11)

## Performance (C60 run)

### Current (with selective unsafeFreezeGraph, 2026-02-17)
```
Time: 1.35s  |  Allocation: 7.5 GB  |  Peak mem: 75 MB  |  Zero duplicates
Expansion sites (bounded):     380,060
After Rule 2:                  180,958 (Rule 2 eliminated 52%)
Canon acceptance rate:           5,767 / 180,958 (3%)
C70: 7.06s / 41.6 GB / 30,579 isomers through C70
```

### Optimization history
| Optimization | Time | Alloc | Speedup |
|---|---|---|---|
| Baseline (no bounds) | 310s | ~1 TB | 1.0x |
| + Bounding lemmas 1-5 | 170s | 467 GB | 1.8x |
| + allReductionsUpTo | 25.5s | 66.6 GB | **12.2x** |
| + boxed Array caches | 23.5s | 66.6 GB | 13.2x |
| + flat UArray (adjFlat/degFlat) | 12.1s | 71 GB | **25.6x** |
| + incremental BFS + IntSet fixes | 15.4s | 64 GB | **20.1x** |
| + two-phase isCanonical | 13.8s | 57 GB | 22.5x |
| + MutGraph (mutable surgery + unsafeFreeze) | 12.4s | 55 GB | 25.0x |
| + profiling-guided fixes (3 below) | 4.2s | 35.2 GB | 73.8x |
| + geometric table (remove BFS distance) | 2.5s | 17.3 GB | 124x |
| + STUArray BFS + interleaved comparison | 2.0s | 12.1 GB | 155x |
| + lazy isCanonical + nbrs elimination | 1.4s | 7.5 GB | 221x |
| + edge colour pre-filter + instrumentation | 1.45s | 7.7 GB | 214x |
| + selective unsafeFreezeGraph | **1.35s** | **7.5 GB** | **230x** |
| C reference (buckygen.c) | 0.02s | — | ~15000x |

### Profiling-guided fixes (2026-02-17, 13.8s → 4.2s, 3.3x speedup)

GHC profiling (`-fprof-auto`) revealed three major bottlenecks:

1. **Parent allReductions used maxLen=5, only ≤2 needed** (27.5% time, 27.8% alloc)
   `expandChildrenST` called `allReductions g` (= `allReductionsUpTo 5 g`), which
   enumerates expensive bent reductions of length 3-5. But `maxExpansionLength` only
   ever inspects L0 (length 1) and L1/B00 (length 2) reductions.
   **Fix**: Changed to `allReductionsUpTo 2 g`.

2. **isCanonical enumerated B00/L1 when L0 dominates** (19.5% time, 16.8% alloc)
   When the inverse has reductionLength > 1 (L1, B00, or longer) but the child has
   an L0 reduction (x0=1), the inverse can never be canonical since x0=1 always beats
   x0≥2 lexicographically. The old code still enumerated all B00 reductions.
   **Fix**: Added `hasAnyL0Reduction` early exit in `isCanonical`.

3. **bfsDistance traversed entire graph, only needed depth ≤ 5** (10.7% time, 20.9% alloc)
   `hasTwoDistantL0s` calls `bfsDistance` to check if distance > 4, but the old BFS
   visited every reachable vertex. Graphs at C60 have 32 dual vertices.
   **Fix**: Added depth bound to `bfsDistanceBounded g src tgt 5`.

Key wins (earlier):
- `allReductionsUpTo`: Skip bent reduction enumeration when expansion length ≤ 2 (6.7x)
- Flat `UArray Int Int` for navigation: unboxed contiguous memory, linear scan for
  indexOf (5-6 entries), `nbrAt` = single array index. Eliminated all navigation
  primitives from the profile (were 25%+ of time). (1.9x over boxed Array)
- Incremental BFS in `canonicalBFSAndGroup`: instead of computing all ~360 BFS codes
  then taking `minimum`, compare each candidate against best-so-far. Lazy list
  comparison aborts most candidates after 1-3 elements. Peak memory dropped from
  796 MB to 31 MB (25x). Key bug: `bestDir` must come from same match as `refNumMap`.
- Two-phase `isCanonical`: Phase 1 compares cheap (x0-x3) tuples for all reductions.
  Only reductions that tie the minimum cheap tuple proceed to Phase 2 (BFS). Reduces
  x4 evaluations from 512K to 28K (18x). Overall speedup: 11% (BFS was only ~15% of
  total work due to lazy list comparison aborting early in the old code).
- IntSet for `isValidB00Reduction` (nub/intersect → O(L log L)),
  `followStraightToFive` (elem → O(log n)), `bfsDistance` (Seq queue → O(n)).
- STUArray BFS with interleaved comparison (matching C code's `testcanon` approach):
  pre-allocate work arrays in a single `runST` block, compare each BFS element
  against the reference during construction. GT cases (~99%) abort after 1-3 elements
  with zero allocation. Only EQ (automorphisms) and LT (rare new minimum) freeze arrays.
  Eliminated IntMap from BFS entirely. Also changed `Automorphism.autPerm` from
  `IntMap Int` to `UArray Int Int` for O(1) permutation lookup.
  C60: 2.5s → 2.0s (20%), 17.3GB → 12.1GB (30%). C70: 13.7s → 10.4s (24%).
- Geometric table (doc/MAX_STRAIGHT_LENGTHS.md): precomputed `max_straight_lengths[nv]`
  replaces BFS-based `hasTwoDistantL0s` (Lemma 6). O(1) table lookup vs O(k²·n) BFS.
  C60: 4.2s → 2.5s (40%), 35.2GB → 17.3GB (51%).
- Lazy `isCanonical` with `any` short-circuiting (matching C's `return 0` pattern):
  reformulated `isCanonical` from `minimum (map cheapCanonOrd allReds)` (forces all) to
  `any (\r -> cheapCanonOrd r < invCheap) allReds` (stops at first better). Guards cascade:
  L0 early exit → null check → `any` phase 1 → tied check → BFS phase 2.
  97% of calls reject after examining 1-5 reductions instead of all 20-60.
  C60: 2.0s → 1.4s (30%), 12.1GB → 7.5GB (38%). C70: 10.4s → 7.5s (28%).
- `nbrs` list elimination: replaced `v <- nbrs g u` (allocates 5-6 cons cells per call)
  with indexed iteration `i <- [0..deg g u - 1], let v = nbrAt g u i` in all reduction
  and expansion enumerations. Added `isAdj g u v` (direct flat-array scan) to replace
  `v `elem` nbrs g u`. Eliminated `nbrs` from entire Canonical.hs hot path.
- Remaining gap to C: ~70x

### Profile distribution (C40, post flat UArray)
```
straightAhead 7.3%  |  isValidB00 5.4%  |  bfsDistance 5.0%  |  bfsCF.procNbrs 4.6%
bentPath.adv 4.2%   |  hasDups 3.8%     |  bentZeroPath 3.1% |  bentPath 2.7%
```
Navigation primitives (indexOf, advanceCW, deg, nextCW, prevCW) have disappeared
from the profile. Remaining time is distributed across algorithm logic and BFS.

## Key Decisions
- Direction-specific canonical test (Bug #15) eliminates all duplicates without
  needing L0 reverse equivalence or or_same_edge_found in Rule 2
- Spiral dedup remains as safety net but currently catches 0 duplicates

## Next Steps (in priority order)

### 1. DONE: Lazy incremental rejection in `isCanonical` + `nbrs` elimination

Implemented 2026-02-17. C60: 2.0s → 1.4s (30%), 12.1GB → 7.5GB (38%).
C70: 10.4s → 7.5s (28%), 65.2GB → 41.4GB (36%). See optimization history above.

**What was done:**
- Reformulated `isCanonical` from `minimum (map cheapCanonOrd allReds)` to guard-based
  `any (\r -> cheapCanonOrd r < invCheap) allReds`. Five guards cascade: L0 early exit →
  null check → `any` short-circuit (phase 1) → unique minimum check → BFS tiebreaker (phase 2).
- Replaced `a `elem` newRange` with O(1) bounds check `a >= newV1 && a <= newV2`.
- Replaced all `v <- nbrs g u` with indexed iteration `i <- [0..deg g u - 1], let v = nbrAt g u i`
  in Canonical.hs (reduction enumerations) and Expansion.hs (expansion enumerations).
- Added `isAdj :: DualGraph -> Vertex -> Vertex -> Bool` (direct flat-array scan) to
  replace all `v `elem` nbrs g u` patterns. Eliminated `nbrs` from entire Canonical.hs.

**Remaining from the analysis (not yet done):**
- Path sharing via `AnnotatedRed` (eliminates 3x redundant path computation for B00)
- Cross-type priority guards (`hasAnyL1Reduction` for B00 rejection)
- Parent-to-child information carrying

### 2. DONE: Colour-pair pre-filter for BFS

Implemented 2026-02-17. 58x edge reduction at C60 (334→5.8 avg). See edge colour
pre-filter section above.

### 3. DONE: Avoid `freezeGraph` for rejected children

Implemented 2026-02-17. Added `unsafeFreezeGraph` (zero-copy alias) for the 97% of
children that fail `isCanonical`. Safe `freezeGraph` only for the 3% that pass.
C60: 1.45s → 1.35s (7%), C70: 7.73s → 7.06s (9%).

**Remaining optimization opportunities:**
- Path sharing via `AnnotatedRed` (eliminates 3x redundant path computation for B00)
- Cross-type priority guards (`hasAnyL1Reduction` for B00 rejection)
- Parent-to-child information carrying
- Remaining gap to C: ~67x

### 4. DONE: Parallel forest traversal

Implemented 2026-02-17. Flat-fork `parList rdeepseq` at configurable fork depth.
C80: 47.9s → 9.7s (4.95x on 8 cores). See parallel benchmark results above.

### 5. Future optimization opportunities
- Path sharing via `AnnotatedRed` (eliminates 3x redundant path computation for B00)
- Cross-type priority guards (`hasAnyL1Reduction` for B00 rejection)
- Parent-to-child information carrying
- Work-stealing thread pool (for better load balance at high core counts)
- Remaining gap to C (sequential): ~70x
