# Search: Unified Fold + Schedule Abstraction

`Search.hs` separates **what to collect** from **how to traverse** the generation
forest.  The core recursion -- seed expansion, bounding, canonical testing -- is
written once; you plug in a `Fold` (the collection algebra) and a `Schedule`
(the traversal strategy), and `search` does the rest.

## The problem it solves

`GenForest.hs` grew 10+ traversal functions -- `dfsCountBySize`,
`parMutCountBySize`, `workQueueCountBySize`, etc. -- all duplicating the same
fold logic over the search tree with the collection pattern (count by size)
baked in.  Want a new collection pattern (collect graphs, compute a histogram,
find the maximum)?  Write another traversal function.  Want it parallel?  Copy
the parallel plumbing too.

The Search module eliminates this M-collection-patterns x N-strategies
duplication.

## Design

Three concepts, one entry point:

```
Fold r          what to collect (monoid homomorphism from nodes to results)
Schedule        how to traverse (DFS, BFS, parallel, work-queue, ...)
search          the single entry point that combines them
```

### `Fold r` -- the collection algebra

```haskell
data Fold r = Fold
    { foldVisit :: !(DualGraph -> [Automorphism] -> r)   -- node -> result
    , foldMerge :: !(r -> r -> r)                        -- combine results
    , foldEmpty :: !r                                     -- identity element
    }
```

The three fields form a monoid homomorphism: `foldVisit` maps each graph in the
generation forest to a value of type `r`, and `foldMerge`/`foldEmpty` must form
a monoid (associative merge, identity element).  Because the generation forest's
subtrees are independent, results can merge in any order -- this is what makes
parallelism safe.

#### Predefined folds

| Fold | Type | Description |
|------|------|-------------|
| `countBySize` | `Fold (Map Int Int)` | Count isomers per dual vertex count |
| `collectAll` | `Fold [(DualGraph, [Automorphism])]` | Collect every graph with its automorphism group |
| `collectOfSize n` | `Fold [DualGraph]` | Collect only graphs with `n` dual vertices |
| `mapNodes f` | `Fold [b]` | Apply `f` to every node, collect in a list |

#### Writing a custom fold

Any monoid works.  Example -- maximum automorphism group size per vertex count:

```haskell
maxAutFold :: Fold (Map Int Int)
maxAutFold = Fold
    { foldVisit = \g auts -> Map.singleton (numVertices g) (length auts)
    , foldMerge = Map.unionWith max
    , foldEmpty = Map.empty
    }
```

Example -- total count (just an integer):

```haskell
totalCount :: Fold Int
totalCount = Fold
    { foldVisit = \_ _ -> 1
    , foldMerge = (+)
    , foldEmpty = 0
    }
```

The monoid laws must hold: `foldMerge` is associative, `foldEmpty` is its
identity.  If they don't, results will differ between sequential and parallel
schedules (where merge order varies).

### `Schedule` -- traversal strategy

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

All schedules produce **identical results** for a given fold and config.  They
differ only in evaluation order, memory profile, and parallelism.

**How to choose:**

| Situation | Schedule | Why |
|-----------|----------|-----|
| Small runs, debugging | `SeqDFS` | Deterministic, no threading overhead |
| Need level-order output | `SeqBFS` | Processes all depth-d nodes before depth d+1 |
| Multicore, moderate size | `ParMut depth` | Best wall-clock: MutGraph avoids allocation |
| Multicore, unbalanced tree | `WorkQMut nw depth` | Dynamic load balancing via TQueue |
| Very large (C400+) | `HierMut nw d1 d2` | Parallelizes the expensive split phase too |
| Experimental | `BfsDfs nb nd depth` | Self-generating work queue, no sequential bottleneck |

The `depth` parameter controls where to split the tree: nodes above the depth
are processed to generate work units; nodes at and below the depth become
independent chunks.  Higher depth = more chunks = better load balance, but more
splitting overhead.  Typical values: 4--8 for C60--C200.

### `search` -- the entry point

```haskell
search     :: NFData r => Fold r -> Schedule -> GenConfig -> IO r
searchPure :: Fold r -> GenConfig -> r
```

`search` dispatches to the appropriate backend.  Pure schedules (`SeqDFS`,
`SeqBFS`, `ParFlat`, `ParMut`) are wrapped in `return $!` for a uniform IO
interface.  `searchPure` is a convenience for when you know you want sequential
DFS and don't want IO.

## Usage examples

### Count C60 isomers (sequential)

```haskell
let cfg = mkConfig 32        -- 32 dual vertices = C60
    counts = searchPure countBySize cfg
-- counts = fromList [(12,1),(14,1),(16,2),(18,5),...,(32,1812)]
```

### Count C60 isomers (parallel, 4-core)

```haskell
counts <- search countBySize (ParMut 6) (mkConfig 32)
-- Same result, uses all cores.  Run with: +RTS -N4
```

### Collect all C60 graphs

```haskell
graphs <- search (collectOfSize 32) (ParMut 6) (mkConfig 32)
-- graphs :: [DualGraph], length 1812
```

### Compute a custom statistic in parallel

```haskell
-- Histogram of min degree-5 vertex distance
let distFold = Fold
      { foldVisit = \g _ -> Map.singleton (minDeg5Dist g) 1
      , foldMerge = Map.unionWith (+)
      , foldEmpty = Map.empty
      }
hist <- search distFold (WorkQMut 8 6) (mkConfig 42)
```

### Use the SearchM monad for custom traversal logic

```haskell
-- Count only IPR fullerenes
let iprCount = runSearchM countBySize (mkConfig 42) $
      searchForest $ \g auts ->
        when (isIPR g) (yield g auts)
```

`SearchM` provides a monadic interface over pure DFS.  `searchForest` walks
every node; your action decides whether to `yield`.  The fold accumulates
yielded results.

## Architecture

### How `Fold` works with parallel schedules

Every parallel schedule follows the same three-phase pattern:

1. **Split**: walk the tree to a fork depth, applying the fold to nodes above
2. **Fan out**: distribute subtree roots to workers (sparks, TQueue, etc.)
3. **Fold per worker**: each worker folds its subtree independently
4. **Merge**: combine all worker results with `foldMerge`

Because `foldMerge` is associative and `foldEmpty` is its identity, the merge
order doesn't matter.  Workers never communicate during the fold phase -- no
locks, no shared mutable state.

### `mutFoldSubtree` -- the MutGraph backend

The key performance function.  Generalizes `GenForest.mutCountSubtree` to
accept any `Fold r`:

```haskell
mutFoldSubtree :: Fold r -> GenConfig -> (DualGraph, [Automorphism]) -> r
```

Runs in a fresh `runST` block with a pre-allocated `MutGraph`.  The DFS loop:

1. Apply expansion (mutate graph in place)
2. `unsafeFreezeGraph` for the canonical test (zero-copy; safe because the test
   is pure and fully evaluates before any further mutation)
3. If canonical: `freezeGraph` (safe copy), call `foldVisit` on the child, recurse
4. `undoMutation` to restore parent state
5. `foldMerge` to accumulate child results

The fold's `foldVisit` and `foldMerge` are pure functions called from within ST.
This is safe: ST enforces that no mutable references escape, and the fold
functions don't use any.  The `Fold r` type has no `ST` in it at all -- it is
the same fold used by the pure DFS backend.

### `SearchM` -- the monadic interface

```haskell
newtype SearchM r a = SearchM (Fold r -> GenConfig -> (a, r))
```

A reader-like monad carrying the fold and config.  Bind (`>>=`) merges the
accumulated `r` from both sides via `foldMerge`.  `yield` calls `foldVisit` to
produce a result fragment.  `searchForest` performs DFS internally, calling your
action at each node.

`SearchM` runs as pure sequential DFS only -- it does not support parallel
schedules.  For parallel execution of custom logic, build a `Fold` and pass it
to `search` with a parallel `Schedule`.

### Why not a type class?

The generation forest uses three fundamentally different evaluation strategies:

- **Pure lazy trees** (`GenTree`): nodes are thunks; DFS = pattern matching
- **ST-based MutGraph**: in-place mutation with apply/freeze/undo
- **IO-based work-queues**: TQueue + forkIO for dynamic load balancing

A type class `SearchMonad m` would need to unify `Identity`, `ST s`, and `IO`
under one interface.  The result type `r` would need to thread through the class
as a parameter.  The MutGraph backend's apply/undo protocol doesn't fit
naturally into a monadic bind.

The `Fold r` + `Schedule` approach avoids all of this.  The fold is a plain
record of pure functions.  The schedule selects the traversal machinery.  The
two are orthogonal -- every fold works with every schedule, with no type-level
gymnastics.

## Relationship to GenForest

`Search` is a layer on top of `GenForest`.  It imports:

- `GenConfig`, `mkConfig`, `GenTree`, `genForest` -- the lazy tree
- `generateChildren` -- the pure generation pipeline
- `seedsWithAuts` -- the three roots (C20, C28, C30)

The existing `GenForest` traversal functions (`dfsCountBySize`,
`parMutCountBySize`, etc.) remain and continue to work.  `Search` provides the
same functionality through a unified interface.  The equivalences are verified
by cross-checks in `DemoForest.hs`:

```
searchPure countBySize cfg             ==  dfsCountBySize cfg
search countBySize (ParFlat 4) cfg     ==  parCountBySize cfg 4
search countBySize (ParMut 4) cfg      ==  parMutCountBySize cfg 4
runSearchM countBySize cfg (searchForest $ \g a -> yield g a)
                                       ==  dfsCountBySize cfg
```

All pass at C60 (5770 total isomers across C20--C60, 1812 at C60).

## Module interface

```haskell
module Search
    ( -- Collection algebra
      Fold(..)                         -- Fold { foldVisit, foldMerge, foldEmpty }
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
    , searchForest                     -- (DualGraph -> [Automorphism] -> SearchM r ()) -> SearchM r ()
    , runSearchM                       -- Fold r -> GenConfig -> SearchM r a -> r

      -- Re-exports
    , GenConfig(..), mkConfig
    ) where
```

## File layout

```
src/Search.hs          The module described here (~640 lines)
src/GenForest.hs       Generation forest (one new export: collectSubtreeRoots)
src/DemoForest.hs      Demo program (items 11-16 exercise the Search API)
```
