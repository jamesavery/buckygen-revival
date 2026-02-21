# The Generation Forest — Search Module Edition

## The forest structure

Buckygen enumerates fullerene isomers by growing them from small seed graphs.
Each growth step applies one of five expansion operations (L0, L1, ..., B00, B01,
..., F) to a parent graph, adding 2--7 new dual vertices. The canonical
construction path principle guarantees that every isomer is generated exactly
once: a child is accepted only if its *inverse reduction* is the
lexicographically smallest reduction of that child.

This defines a forest --- a collection of rooted trees with no cross-edges:

```
C20 -- C24 -- C28a -- ...
  |      \---- C28b -- ...
  \---- C26 -- C30a -- ...
               \---- ...
C28 -- C32a -- ...
  \---- C32b -- ...

C30 -- C34a -- ...
  \---- C34b -- ...
         \---- C38a -- ...
```

Three roots: the dodecahedron C20, the D2 fullerene C28, and the (5,0)
nanotube cap C30. Every other fullerene isomer appears as exactly one node
in exactly one of these three trees.

### Key structural properties

**Subtree independence.** A node's children depend only on the node itself (its
graph and its automorphism group). There is no communication between subtrees,
no shared mutable state, no global counters. This is not a design choice ---
it is forced by the mathematics. The canonical test for a child graph examines
only the child's own reductions.

**Monotone growth.** Every child has strictly more dual vertices than its parent.
Children of a node with `nv` dual vertices have `nv+2` to `nv+7` vertices.
This means the tree has bounded depth (at most `(maxDV - 12)/2` levels from the
smallest seed), and every path from root to leaf is a strictly increasing
sequence of graph sizes.

**Uneven branching.** The branching factor varies wildly. Some nodes (especially
small, highly symmetric graphs) have dozens of children; others have few or
none. The three root subtrees are also very different sizes: at C60, C20's
subtree contains 97% of all isomers. This unevenness is the central challenge
for parallel evaluation.

### Node computation

Each node in the forest stores three things:

```haskell
data GenTree = GenTree
    { gnGraph    :: !DualGraph          -- the fullerene dual graph
    , gnAuts     :: ![Automorphism]     -- its automorphism group
    , gnChildren :: [GenTree]           -- lazily computed children
    }
```

The children are computed by `generateChildren`, a pure function that
implements the five-stage pipeline:

1. **Bound.** Enumerate reductions of the parent (length <= 2). Apply
   bounding lemmas to determine the maximum useful expansion length.
2. **Enumerate.** List all expansion sites up to that length.
3. **Filter (Rule 2).** Partition expansions into orbits under the parent's
   automorphism group. Keep one representative per orbit.
4. **Apply.** Apply each surviving expansion via pure graph surgery,
   producing a candidate child graph.
5. **Accept (Rule 1).** Test whether the expansion's inverse reduction is
   the canonical (lexicographically smallest) reduction of the child.
   If so, compute the child's automorphism group and accept it.

Steps 1--3 depend on the parent graph and its automorphisms.
Steps 4--5 depend on the parent graph and the expansion.
Nothing depends on any other node in the forest.

### Laziness

The `gnChildren` field is a lazy list. `genForest` constructs the entire
infinite-depth tree structure without evaluating any children:

```haskell
genForest cfg = map mkTree seedsWithAuts
  where
    mkTree (g, auts) = GenTree g auts
        [mkTree (c, ca) | (c, ca) <- generateChildren cfg g auts]
```

No `generateChildren` call executes until something demands a particular
node's children. This is what allows the same forest to support DFS, BFS,
parallel, streaming, and filtered traversals without changing the generation
logic.


## The Search abstraction

The `Search` module separates **what to collect** from **how to traverse** the
forest. You plug in a `Fold` (the collection algebra) and a `Schedule` (the
traversal strategy), and `search` does the rest.

### `Fold r` --- what to collect

```haskell
data Fold r = Fold
    { foldVisit :: !(DualGraph -> [Automorphism] -> r)   -- node -> result
    , foldMerge :: !(r -> r -> r)                        -- combine results
    , foldEmpty :: !r                                     -- identity element
    }
```

The three fields form a monoid homomorphism: `foldVisit` maps each graph to a
value, `foldMerge` combines values (must be associative), and `foldEmpty` is
the identity. Because the forest's subtrees are independent, results can merge
in any order --- this is what makes parallelism safe.

Predefined folds:

| Fold | Type | Description |
|------|------|-------------|
| `countBySize` | `Fold (Map Int Int)` | Count isomers per dual vertex count |
| `collectAll` | `Fold [(DualGraph, [Automorphism])]` | Collect every graph with its automorphism group |
| `collectOfSize n` | `Fold [DualGraph]` | Collect only graphs with `n` dual vertices |
| `mapNodes f` | `Fold [b]` | Apply `f` to every node, collect in a list |

Writing a custom fold is straightforward:

```haskell
-- Maximum automorphism group size per vertex count:
maxAutFold = Fold
    { foldVisit = \g auts -> Map.singleton (numVertices g) (length auts)
    , foldMerge = Map.unionWith max
    , foldEmpty = Map.empty
    }

-- Just count everything:
totalCount = Fold { foldVisit = \_ _ -> 1, foldMerge = (+), foldEmpty = 0 }
```

### `Schedule` --- how to traverse

```haskell
data Schedule
    = SeqDFS                  -- sequential depth-first
    | SeqBFS                  -- sequential breadth-first
    | ParFlat  depth          -- parallel flat-fork, pure DFS per subtree
    | ParMut   depth          -- parallel flat-fork, MutGraph DFS per subtree
    | WorkQ    nWorkers depth -- TQueue work-stealing, pure DFS backend
    | WorkQMut nWorkers depth -- TQueue work-stealing, MutGraph DFS backend
    | HierMut  nWorkers d1 d2 -- 3-phase hierarchical parallel
    | BfsDfs   nBfs nDfs depth -- parallel BFS feeds parallel DFS
```

All schedules produce **identical results** for a given fold and config. They
differ only in evaluation order, memory profile, and parallelism.

### `search` --- the entry point

```haskell
search     :: NFData r => Fold r -> Schedule -> GenConfig -> IO r
searchPure :: Fold r -> GenConfig -> r
```

`search` dispatches to the appropriate backend. `searchPure` is a convenience
for sequential DFS without IO.


## Traversal schemes

### 1. Sequential DFS (`SeqDFS`)

The baseline. A strict left fold over the forest in pre-order. Each node is
visited, folded, then its children are visited left-to-right. Memory usage is
proportional to tree depth, not total size.

```haskell
searchPure countBySize cfg
-- or: search countBySize SeqDFS cfg
```

### 2. Sequential BFS (`SeqBFS`)

Processes all nodes at generation depth *d* before any at depth *d+1*. The
entire frontier at the widest level must be held in memory simultaneously.
Useful for algorithms that need level-order processing.

```haskell
search countBySize SeqBFS cfg
```

### 3. Parallel flat-fork, pure (`ParFlat depth`)

Two phases:

**Phase 1 --- Sequential split.** Walk the lazy tree from roots to fork depth
*d*, folding nodes above and collecting subtree roots at depth *d*.

**Phase 2 --- Parallel evaluation.** Each subtree is folded via pure DFS,
sparked in parallel via `parList rdeepseq`.

```haskell
search countBySize (ParFlat 9) cfg   -- +RTS -N10
```

### 4. Parallel flat-fork, MutGraph (`ParMut depth`)

Same split-and-spark structure, but each subtree uses a mutable graph
(`MutGraph`) for in-place surgery with apply/unsafeFreeze/undo. Each subtree
gets its own `MutGraph` in a separate `runST` block --- zero shared state.

The 97% of children that fail the canonical test use zero-copy `unsafeFreeze`.
Only the 3% that pass allocate a safe copy.

```haskell
search countBySize (ParMut 9) cfg    -- fastest for most sizes
```

### 5. Work-queue, pure (`WorkQ nWorkers depth`)

Pre-splits the tree at a fork depth, puts subtree roots into a `TQueue`.
Workers pull chunks and pure-DFS each one. Dynamic load balancing: fast workers
that finish small subtrees immediately grab the next chunk.

```haskell
search countBySize (WorkQ 10 9) cfg
```

### 6. Work-queue, MutGraph (`WorkQMut nWorkers depth`)

Same as WorkQ but each chunk is processed using MutGraph DFS. Combines
dynamic load balancing with MutGraph's low-allocation inner loop.

```haskell
search countBySize (WorkQMut 10 9) cfg
```

### 7. Hierarchical parallel (`HierMut nWorkers d1 d2`)

Three phases: (1) sequential split to small depth *d1*, (2) parallel split
of each level-1 subtree to depth *d2* via `parList rdeepseq`, (3) work-queue
distributes all resulting chunks to MutGraph DFS workers.

Parallelizes the split phase itself. At C400+ the sequential split would
dominate; this scheme avoids that bottleneck.

```haskell
search countBySize (HierMut 10 4 5) cfg
```

### 8. Parallel BFS+DFS (`BfsDfs nBfs nDfs depth`)

Self-bootstrapping: BFS workers generate children concurrently, routing them
to either the BFS queue (if above fork depth) or the DFS queue (at fork depth).
DFS workers pull chunks and count them with MutGraph. Even the initial 3 seeds
are processed concurrently --- no sequential bottleneck at any scale.

```haskell
search countBySize (BfsDfs 10 10 7) cfg
```

### `SearchM` monad

For custom traversal logic within pure DFS:

```haskell
-- Count only IPR fullerenes
let iprCount = runSearchM countBySize cfg $
      searchForest $ \g auts ->
        when (isIPR g) (yield g auts)
```

`searchForest` walks every node; your action decides whether to `yield`. The
fold accumulates yielded results. Runs pure sequential DFS only.


## Cross-validation

`DemoForest.hs` exercises all traversal schemes at C60 and cross-checks
identical results:

```
=== All cross-checks ===
   PASS: BFS == DFS
   PASS: fold == DFS
   PASS: stream == DFS
   PASS: par == DFS
   PASS: mut == DFS
   PASS: searchPure == DFS
   PASS: search ParFlat == DFS
   PASS: search ParMut == DFS
   PASS: SearchM == DFS
   PASS: collectOfSize == DFS

ALL SCHEMES AGREE
```


## Performance

All benchmarks on Apple M2 Ultra (arm64), GHC 9.6.7 -O2, 10 cores
(`+RTS -N10`), 3 runs per test (median reported).

### C70 (37 dv, 30,579 isomers, 10 cores, fork depth 9)

```
Scheme                              Time        Speedup vs SeqDFS
Sequential DFS                       8.6s         1.0x
Sequential BFS                       8.4s         1.0x
Parallel Pure (d=9)                  1.2s         7.2x
Parallel MutGraph (d=9)              1.2s         7.2x
Work-queue Pure (10w, d=9)           2.0s         4.3x
Work-queue MutGraph (10w, d=9)       1.9s         4.5x
Hierarchical MutGraph (d1=4,d2=5)    950ms        9.1x
BFS+DFS MutGraph (d=7)               1.1s         7.8x
```

All results match (30,579 isomers).

### C80 (42 dv, 131,200 isomers, 10 cores, fork depth 9)

```
Scheme                              Time        Speedup vs SeqDFS
Sequential DFS                      39.5s         1.0x
Sequential BFS                      37.7s         1.0x
Parallel Pure (d=9)                  4.3s         9.2x
Parallel MutGraph (d=9)              3.9s        10.1x
Work-queue Pure (10w, d=9)           5.3s         7.5x
Work-queue MutGraph (10w, d=9)       5.1s         7.7x
Hierarchical MutGraph (d1=4,d2=5)    3.9s        10.1x
BFS+DFS MutGraph (d=7)               5.2s         7.6x
```

All results match (131,200 isomers).

### Analysis

HierMut and ParMut are the fastest schedules, achieving near-linear scaling
on 10 cores (10.1x at C80).  The three-phase design of HierMut ---
sequential bootstrap, parallel split, work-queue distribution --- eliminates
both the sequential bottleneck (that limits ParMut at very large sizes) and
the load imbalance (that limits flat-fork schemes).

The `Fold r` parameterization introduces no measurable overhead --- GHC
specializes the fold record's function fields at compile time, producing the
same inner loops as hardcoded versions would.


## Choosing a schedule

| Situation | Schedule | Why |
|-----------|----------|-----|
| Small runs, debugging | `SeqDFS` | Deterministic, no threading overhead |
| Need level-order output | `SeqBFS` | Processes all depth-d nodes before depth d+1 |
| Multicore, moderate size | `ParMut depth` | Best wall-clock: MutGraph avoids allocation |
| Multicore, unbalanced tree | `WorkQMut nw depth` | Dynamic load balancing via TQueue |
| Very large (C400+) | `HierMut nw d1 d2` | Parallelizes the expensive split phase too |
| Experimental | `BfsDfs nb nd depth` | Self-generating work queue, no sequential bottleneck |

The `depth` parameter controls where to split the tree. Higher depth = more
chunks = better load balance, but more splitting overhead. Typical values:
4--8 for C60--C200, 9+ for larger sizes.


## Module interface

```haskell
module Search
    ( -- Configuration
      GenConfig(..), mkConfig
      -- Seeds
    , seedsWithAuts
      -- Pure child generation
    , generateChildren
      -- Lazy forest
    , GenTree(..), genForest
      -- Collection algebra
    , Fold(..)                         -- Fold { foldVisit, foldMerge, foldEmpty }
    , countBySize                      -- Fold (Map Int Int)
    , collectAll                       -- Fold [(DualGraph, [Automorphism])]
    , collectOfSize                    -- Int -> Fold [DualGraph]
    , mapNodes                         -- (DualGraph -> [Automorphism] -> b) -> Fold [b]
      -- Scheduling strategy
    , Schedule(..)                     -- SeqDFS | SeqBFS | ParFlat | ... | BfsDfs
      -- Unified search
    , search                           -- NFData r => Fold r -> Schedule -> GenConfig -> IO r
    , searchPure                       -- Fold r -> GenConfig -> r
      -- Monadic interface (pure DFS only)
    , SearchM                          -- opaque
    , yield                            -- DualGraph -> [Automorphism] -> SearchM r ()
    , searchForest                     -- (...) -> SearchM r ()
    , runSearchM                       -- Fold r -> GenConfig -> SearchM r a -> r
      -- Convenience traversals
    , allGraphs                        -- GenConfig -> [(DualGraph, [Automorphism])]
    , graphsOfSize                     -- GenConfig -> Int -> [DualGraph]
    , bfsByLevel                       -- GenConfig -> [[(DualGraph, [Automorphism])]]
      -- Analysis
    , depthProfile                     -- GenConfig -> [(Int, Int)]
    , subtreeSizes                     -- GenConfig -> Int -> [Int]
    ) where
```


## File layout

```
src/Search.hs          Self-contained: forest + Fold + Schedule + search (~680 lines)
src/DemoForest.hs      Demo: exercises all schemes, cross-checks at C60
src/BenchSearch.hs     Benchmark: all 8 schedule backends with timing
src/BenchPar.hs        CLI benchmark: per-mode parallel execution
obsolete/GenForest.hs  Retired monolithic module (replaced by Search.hs)
```
