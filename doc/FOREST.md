# The Generation Forest

## The forest structure

Buckygen enumerates fullerene isomers by growing them from small seed graphs.
Each growth step applies one of five expansion operations (L0, L1, ..., B00, B01,
..., F) to a parent graph, adding 2-7 new dual vertices. The canonical
construction path principle guarantees that every isomer is generated exactly
once: a child is accepted only if its *inverse reduction* is the
lexicographically smallest reduction of that child.

This defines a forest --- a collection of rooted trees with no cross-edges:

```
C20 ── C24 ── C28a ── ...
  │      └──── C28b ── ...
  └──── C26 ── C30a ── ...
               └──── ...
C28 ── C32a ── ...
  └──── C32b ── ...

C30 ── C34a ── ...
  └──── C34b ── ...
         └──── C38a ── ...
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


## Traversal schemes

All traversal schemes consume the same lazy `genForest`. The generation
logic is written once; the traversal strategy is chosen independently.

### 1. Sequential DFS (`dfsCountBySize`)

The baseline. A strict left fold over the forest in pre-order, accumulating
a `Map Int Int` of counts by dual vertex size:

```haskell
dfsCountBySize cfg = foldl' countTree Map.empty (genForest cfg)
```

Each node is visited once, counted, then its children are visited
left-to-right. Memory usage is proportional to the tree depth (the
stack of pending sibling lists), not the total tree size. At C60 this
processes 5,770 nodes in 1.9s.

### 2. BFS by level (`bfsByLevel`, `bfsCountBySize`)

Returns graphs grouped by generation depth: level 0 = seeds, level 1 =
children of seeds, and so on.

```haskell
bfsByLevel cfg = go (genForest cfg)
  where
    go [] = []
    go level =
        let valid   = [... | GenTree g auts children <- level, nv <= maxDV]
            graphs  = [(g, auts) | (g, auts, _) <- valid]
            nextLevel = concatMap (\(_, _, cs) -> cs) valid
        in graphs : go nextLevel
```

All nodes at depth *d* are materialized before any at depth *d+1*. This
means the entire frontier at the widest level must be held in memory
simultaneously. For C60 the widest level has ~1,400 nodes --- manageable,
but for larger targets BFS uses significantly more memory than DFS.

BFS is useful for understanding the tree's shape (how many nodes at each
generation depth) and for algorithms that need to process all graphs of
a given generation before proceeding.

### 3. Lazy DFS stream (`allGraphs`)

Flattens the forest into a lazy list of `(DualGraph, [Automorphism])` pairs
in DFS pre-order:

```haskell
allGraphs cfg = concatMap go (genForest cfg)
  where
    go (GenTree g auts children)
        | numVertices g > maxDV = []
        | otherwise = (g, auts) : concatMap go children
```

Graphs are produced one at a time on demand. You can take the first 10
graphs, filter for a property, or pipe into a downstream consumer without
ever materializing the full tree. This is the most flexible interface for
composition with other lazy operations.

### 4. Generic fold (`foldForest`)

A strict left fold with an arbitrary accumulator, generalizing
`dfsCountBySize`:

```haskell
foldForest :: (a -> DualGraph -> [Automorphism] -> a) -> a -> GenConfig -> a
```

The function receives the accumulator, the graph, and its automorphism group
at each node. Traversal is DFS pre-order. Examples:

```haskell
-- Count by size (reimplements dfsCountBySize):
foldForest (\m g _ -> Map.insertWith (+) (numVertices g) 1 m) Map.empty cfg

-- Max |Aut| per size:
foldForest (\m g auts -> Map.insertWith max (numVertices g) (length auts) m)
           Map.empty cfg
```

### 5. Target-size generation (`graphsOfSize`)

Returns only graphs with exactly the target number of dual vertices, pruning
subtrees whose roots already exceed the target:

```haskell
graphsOfSize cfg target = concatMap go (genForest cfg)
  where
    go (GenTree g _ children)
        | nv > target  = []     -- prune: all descendants are larger
        | nv == target = [g]    -- emit, don't recurse
        | otherwise    = concatMap go children
```

Because growth is monotone (children are always larger), any node exceeding
the target can be pruned along with its entire subtree. This is substantially
faster than generating all isomers up to `target` and filtering --- it avoids
expanding the large population of intermediate-size graphs that aren't needed.

### 6. Parallel flat-fork, pure (`parCountBySize`)

The first parallel scheme. Two phases:

**Phase 1 --- Sequential split.** Walk the lazy tree from the roots to a
configurable fork depth *d*, counting nodes above the fork and collecting
all subtree roots at depth *d* into a flat list:

```
         root
        / | \
       /  |  \          depth 0..d: sequential count
      .   .   .
     /|  /|\  |\
    s1 s2 s3 ... sK     depth d: collect subtree roots
```

At C60 with *d* = 6, this yields ~230 subtrees.

**Phase 2 --- Parallel evaluation.** Each subtree is DFS-counted
sequentially, and all subtree counts are sparked in parallel:

```haskell
map countSeqTree subtrees `using` parList rdeepseq
```

`parList` creates one GHC spark per subtree. The RTS distributes sparks
across capabilities (OS threads). Each spark runs a sequential DFS of its
subtree with zero shared state.

**Phase 3 --- Merge.** Union the per-subtree maps with the top-level count.

**Why flat-fork?** An earlier recursive approach --- spark children at every
level --- suffered from spark fizzling: the parent thread's fold evaluated
sparks before worker threads could pick them up (80%+ fizzle rate). The
flat-fork creates all sparks at once in a single flat list, achieving
229/230 spark conversion.

**Limitation.** The assignment of subtrees to capabilities is essentially
static. If one subtree has 10x more nodes than average, the capability that
gets it runs 10x longer while others idle.

### 7. Parallel flat-fork, MutGraph (`parMutCountBySize`)

Same split-and-spark structure as scheme 6, but the sequential DFS within
each subtree uses a mutable graph (`MutGraph`) instead of pure graph surgery.

The pure `generateChildren` creates a fresh immutable `DualGraph` for every
expansion via IntMap surgery and array construction. The MutGraph backend
instead:

1. Pre-allocates a single mutable flat array (Nf x 6 neighbor slots).
2. Applies each expansion by mutating ~10 array entries in place.
3. Freezes the array to an immutable snapshot for the canonical test.
4. Undoes the mutation to restore the parent state.

This eliminates most allocation in the inner loop. The 97% of children
that fail the canonical test use a zero-copy `unsafeFreeze` (safe because
the strict `CanonVerdict` result is fully evaluated before the undo).

Each subtree gets its own `MutGraph` in a separate `runST` block, so
there is still zero shared state between sparks:

```haskell
mutCountSubtree (g, auts) = runST $ do
    mg <- newMutGraph maxDV
    loadGraph mg g
    mutDFS mg g auts
```

This is the fastest backend: 1.45x faster than the pure forest at any
core count, with 25% less allocation.

### 8. Parallel map (`parMapForest`)

Applies a function to every graph in the forest and collects the results,
using flat-fork parallelism:

```haskell
parMapForest :: NFData b => GenConfig -> Int
             -> (DualGraph -> [Automorphism] -> b) -> [b]
```

Uses the same split-at-depth pattern: results above the fork are collected
sequentially; subtrees below the fork are mapped in parallel via
`parList rdeepseq`.

Examples:

```haskell
-- Automorphism group size for every graph:
parMapForest cfg 4 (\_ auts -> length auts) :: [Int]

-- Minimum degree-5 distance per graph:
parMapForest cfg 4 (\g _ -> minDeg5Distance g) :: [Int]
```

### 9. Work-queue (`workQueueCountBySize`)

An explicit thread pool with dynamic load balancing via STM.

**Phase 1 --- Split.** Same as flat-fork: walk the tree to a fork depth,
collect subtree roots. The fork depth is chosen automatically as
`ceil(log2(nWorkers * 8))`, clamped to at least 4.

**Phase 2 --- Enqueue.** Put all subtree roots into a `TQueue`. Set
`inFlight` = number of chunks.

**Phase 3 --- Worker loop.** Each of `nWorkers` threads (via `forkIO`)
runs:

```
loop:
    atomically:
        try read from queue
        if got item    -> return it
        if empty, inFlight > 0 -> retry (block until queue changes)
        if empty, inFlight = 0 -> return Nothing (all done)

    DFS the subtree locally (pure generateChildren)
    merge count into per-worker IORef
    atomically: decrement inFlight
    goto loop
```

**Phase 4 --- Merge.** Read all per-worker `IORef` accumulators and union
them with the top-level counts.

**Termination.** The `inFlight` TVar counts chunks that have been dequeued
but not yet fully processed. When `inFlight` reaches 0 and the queue is
empty, every worker's `retry` resolves to `Nothing` and they all exit.

**Comparison with flat-fork.** Both schemes split the tree at the same
depth and produce the same chunks. The difference is scheduling:

| | Flat-fork (`parList`) | Work-queue (`TQueue`) |
|---|---|---|
| Scheduler | GHC spark pool | Explicit forkIO + TQueue |
| Assignment | Static | Dynamic: fast workers take more |
| Per-chunk overhead | ~0 (spark = pointer) | 1 STM read + 1 STM write |
| Tail latency | Waits for slowest spark | Less likely to strand work |
| Purity | Pure (`Map Int Int`) | IO (returns `IO (Map Int Int)`) |

In practice at C80 with 8 cores and ~230 chunks, both give identical
wall-clock times (~9.5s) because the chunks are numerous enough that
statistical imbalance is small. The work-queue's advantage would emerge at
higher core counts or with fewer, more uneven chunks.

**Design history.** The first implementation enqueued individual graphs
(not subtree chunks), with workers re-enqueueing children. This had two
problems: (1) a bounded `TBQueue` deadlocked when the frontier exceeded
capacity and a single worker blocked on enqueue while being the only
consumer; (2) even with an unbounded `TQueue`, per-graph STM transactions
made it 10x slower than flat-fork. The final design pre-splits into chunks
(eliminating STM from the inner loop) and uses `TQueue` (eliminating the
deadlock).

### 10. Analysis tools (`depthProfile`, `subtreeSizes`)

Not traversal schemes per se, but tools for understanding the tree's shape
to inform the choice of fork depth.

`depthProfile` counts nodes at each generation depth:

```
depth 0:     3 nodes     (seeds)
depth 1:    11 nodes
depth 2:    43 nodes
  ...
depth 9:  1406 nodes     (widest level at C60)
  ...
depth 14:    2 nodes     (deepest)
```

`subtreeSizes` computes the node count of each subtree at a given fork
depth, sorted descending. This reveals work imbalance:

```
subtreeSizes cfg 6 = [182, 156, 143, ..., 1, 1, 1]
-- 230 subtrees, biggest has 182 nodes, smallest has 1
-- ratio biggest/average: ~7.3x
```


## Cross-validation

`DemoForest.hs` exercises schemes 1--8 at C60 and cross-checks that all
produce identical isomer counts:

```
=== All cross-checks ===
   PASS: BFS == DFS
   PASS: fold == DFS
   PASS: stream == DFS
   PASS: par == DFS
   PASS: mut == DFS

ALL SCHEMES AGREE
```

The work-queue (scheme 9) is validated separately via `BenchPar.hs` since
it returns `IO`:

```
Total isomers (work-queue, 8 workers): 5770    -- matches DFS
```


## Performance summary (C80, 131,200 isomers, 8 cores)

```
Scheme                          Wall clock    CPU utilization
Sequential DFS (1 core)           48s            100%
Parallel flat-fork, pure          9.5s           612%
Parallel flat-fork, MutGraph      6.7s           601%
Work-queue, 8 workers             9.5s           607%
```

The MutGraph backend is the fastest because it eliminates allocation in the
inner loop (in-place surgery + unsafeFreeze). The pure flat-fork and
work-queue are equivalent because the coordination overhead (sparks vs STM)
is negligible relative to the per-subtree DFS work. All three parallel
schemes achieve near-linear scaling from 1 to 8 cores.
