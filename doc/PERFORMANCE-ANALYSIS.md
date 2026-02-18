# Performance Analysis: Closing the 600x Gap

## Current state

| | Haskell | C (buckygen.c) | Gap |
|---|---|---|---|
| C60 (1812 isomers) | 12.1s | 0.02s | 605x |
| Allocation | 71 GB | ~0 (arena) | ∞ |

## Root cause: two algorithmic violations

### Problem 1: `canonicalBFSAndGroup` evaluates ALL full BFS codes

**Where:** `Canonical.hs:445-476`, called once per accepted isomer.

**What it does:** For every starting edge (u, v, dir) — roughly `nv × avgDeg × 2 ≈ 360`
candidates at C60 — it computes the **full** BFS canonical form (a ~180-element
list), then finds `minimum` over all 360 codes.

**The cost:**
- 360 full BFS traversals × O(n log n) each = O(n² log n) per isomer
- `minimum` compares 360 lists of length ~6n pairwise = O(n² · numStarts)
- For 5770 isomers at C60: ~1.93M BFS calls, each producing ~180 boxed Ints

**What C does instead (`testcanon` / `testcanon_mirror`):**
1. **Pre-filter by colour pair:** Only test starting edges whose
   `(colour[u], colour[v])` matches the minimum colour pair found so far.
   This reduces ~360 candidates to ~10-20.
2. **Interleaved BFS + compare:** Compute BFS elements one at a time,
   comparing against the current best. Abort as soon as the current code
   exceeds the best. Most candidates are rejected after examining 1-3 elements.
3. **Result:** O(n) amortized per isomer, not O(n² log n).

**Impact of fixing:** Estimated **10-50x speedup** from this alone.

### Problem 2: O(n) work in operations that should be O(1)

Several inner-loop operations do O(n) or O(n²) work where O(1) is expected:

| Operation | Location | Actual | Expected | Called |
|---|---|---|---|---|
| `mkDualGraph` array rebuild | Seeds.hs:52 | O(n) | O(1) amortized | Every expansion |
| `intersect` on path lists | Expansion.hs:1005 | O(L²) | O(L) | ~7500/graph |
| `nub` on paths | Canonical.hs:627 | O(L²) | O(L log L) | ~200/graph |
| `bfsDistance` queue append | Canonical.hs:871 | O(n²) total | O(n) | ~10/graph |
| `elem` on growing visited | Expansion.hs:1015 | O(k) per step | O(1) | ~1500/graph |
| `sort` on degree5 list | Expansion.hs:483+ | O(12 log 12) | O(1) insert | Every expansion |

The `mkDualGraph` issue is particularly bad: every candidate expansion pays O(n) to
build 3 arrays, but ~97% of candidates fail the canonicity test. The arrays are wasted.

## Why x4 (BFS tiebreaker) fires 8.6% of the time, not <1%

Profiling at C40 shows the 5-tuple cascade resolution:

```
x0 evaluated: 50,664 (100%)
x1 evaluated: 17,858 (35.2%)     — x0 resolves 65%
x2 evaluated: 16,994 (33.5%)     — x1 resolves 1.7%
x3 evaluated:  4,382 (8.65%)     — x2 resolves 25%
x4 evaluated:  4,340 (8.57%)     — x3 resolves 0.08% ← almost nothing!
```

**x3 (pathColour) has near-zero discriminating power.** 99% of cases reaching x3
fall through to x4, which forces a full BFS computation. This means BFS is not
a rare tiebreaker — it's computed for ~1 in 12 reductions.

## Fix plan (priority order)

### Fix A: Interleaved BFS with early abort in `canonicalBFSAndGroup`

Replace the current "compute all, then minimum" pattern:

```haskell
-- CURRENT (O(n² log n) per isomer):
results = [(bfsWithNumbering g u v d, d) | (u,v,d) <- allStarts]
bestCode = minimum [c | (c,_) <- results]
```

with C's testcanon strategy:

```haskell
-- NEW (O(n) amortized per isomer):
-- For each starting edge, compute BFS one element at a time.
-- Compare against current best. Abort early if worse.
-- Only complete BFS for codes that tie the best through all elements.
```

This requires a new `bfsCompare` function that interleaves BFS traversal with
lexicographic comparison, returning `LT | EQ | GT` without materializing the
full code.

### Fix B: Colour-pair pre-filter for starting edges

Before trying all starting edges, compute the colour pair for each and only
test edges with the minimum colour pair. This is what the C code calls
`canon_edge_colour` / `better_colour` — it reduces candidates from ~360 to ~20.

### Fix C: Delay `mkDualGraph` until after canonicity test

Introduce a lightweight graph that skips array construction:

```haskell
-- Test canonicity on IntMap-only graph (no arrays):
isCanonical e nv partialGraph

-- Only build full DualGraph for the ~3% that pass:
fullGraph = mkDualGraph nv adj deg5 edgeList
```

### Fix D: Replace list operations with IntSet in validation

- `intersect` → `IS.intersection`
- `nub` → `IS.fromList`
- Growing `visited` list → `IntSet`

### Fix E: Fix `bfsDistance` O(n²) queue

Replace `queue ++ next` with `Data.Sequence` for O(1) amortized enqueue.

## Expected impact

| Fix | Expected speedup | Cumulative |
|---|---|---|
| A: Interleaved BFS + early abort | 5-10x | 1.2-2.4s |
| B: Colour-pair pre-filter | 2-5x | 0.2-1.2s |
| C: Delay mkDualGraph | 1.3-2x | 0.1-0.9s |
| D: IntSet in validation | 1.1-1.2x | 0.1-0.8s |
| E: Fix bfsDistance queue | 1.01x | 0.1-0.8s |
| **Combined** | **15-100x** | **0.1-0.8s** |

The C reference at 0.02s would still be ~5-40x faster due to:
- O(1) mutable array access vs O(log n) IntMap
- Zero allocation (global arena) vs per-graph construction
- Cache-friendly pointer-based edge traversal
- C compiler optimizations (inlining, register allocation)

These constant-factor gaps are addressable later via the do/undo pattern
(mutable graph with backtracking) but are not the priority.

## Profile snapshot (C40, post flat-UArray optimization)

```
straightAhead          7.3%  |  isValidB00Reduction  5.4%
bfsDistance.go         5.0%  |  bfsCF.processNeighbors 4.6%
computeBentPathSafe    4.2%  |  hasDups.go           3.8%
computeBentZeroPath    3.1%  |  turnAhead            2.3%
```

No single function > 8%. Navigation primitives (indexOf, advanceCW, deg,
nextCW, prevCW) have been eliminated from the profile by the flat UArray
optimization. Remaining time is algorithm logic and BFS.
