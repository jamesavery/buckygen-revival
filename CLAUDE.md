# Buckygen Rewrite — Mission & Architecture

## Repository layout

```
buckygen-revival/
  CLAUDE.md              — This file (project instructions)
  src/                   — Haskell source files
    Seeds.hs             — Seed graphs (C20, C28, C30) + DualGraph type
    Expansion.hs         — Expansion/reduction operations (L, B, F types)
    Canonical.hs         — Canonical test (Rule 1), automorphisms, Rule 2
    Spiral.hs            — Canonical generalized spiral computation
    TestCanonical.hs     — Main test: generation with canonical test + Rule 2
    TestExpansion.hs     — Expansion/reduction unit tests
    Generate.hs          — Brute-force generation with isomorphism dedup
    GenTestData60.hs     — Test data generator (all graphs up to C60)
    TestGraphs60.hs      — Generated test data (5770 graphs with spirals)
  doc/                   — Algorithm documentation
    C-CODE-ALGORITHM.md  — AUTHORITATIVE: complete algorithm from C code
    BUCKYGEN-MAP.md      — Mapping between paper, C code, and Haskell
    PROGRESS.md          — Current status and bug history
    (other .md files)    — Algorithm notes (BFS, spiral, automorphisms, etc.)
```

The original C code (`buckygen.c`) and paper source (`buckygen_paper_source/`)
remain in the parent directory (`../`).

## Goal

Rewrite buckygen — a 15,000-line C program that enumerates fullerene isomers via
canonical construction path — as a clean, high-level Haskell implementation (and
eventually C++23), where the algebraic structure is explicit rather than buried in
global mutable state.

## Why rewrite

The C code is fast but architecturally hostile:

- **All state is global** — 11 independent mark systems, a static edge arena,
  megabytes of static arrays. You can't even run two instances in one address space.
- **Parallelism is a hack** — `res/mod` partitioning at a fixed split level,
  coordinated by System V message queues between separate processes.
- **Pruning is hardwired** — bounding lemmas are interleaved with generation logic;
  adding a new bound (e.g. pentagon-distance, HOMO-LUMO gap) means surgery on the
  core loop.
- **Scaling hits a wall** — the `numbering[2*MAXE][MAXE]` array is O(n^2); at 500
  dual vertices it's ~137 MB alone.

## The key insight: the algorithm is already algebraically clean

The pseudocode spec (`bucky_pseudocode.tex`) reveals that underneath the C hairball,
the algorithm decomposes into:

1. **A canonical search forest** — three roots (C20, C28, C30), plus IPR seeds. Each
   subtree is completely independent — no cross-tree communication.

2. **A per-node step** that's a pure 5-function pipeline:
   ```
   bound -> expansions -> apply -> isCanonical -> children
   ```

3. **Monoidal result collection** — fullerenes are gathered by list concatenation.
   Evaluation order doesn't matter.

4. **Bounds as a composable list** — each bounding lemma is a function
   `DualGraph -> Maybe Int`. They compose via `minimum` (meet in the lattice
   N union {infinity}). Adding a new bound = appending to a list.

5. **Parallelism as a strategy annotation** — `using parList` is semantically the
   identity. Depth-limited sparking: above a threshold, fork; below, run
   sequentially. The `SearchMonad` abstraction lets the same core code run as DFS,
   BFS, or parallel.

6. **Branch-and-bound as the one impure extension** — a shared monotonically-tightened
   bound (`MVar`/`TVar`). Safe because a stale read only causes wasted work, never
   missed results.

## The canonicity test is the performance heart

The 6-tuple `(x0, x1, x2, x3, x4, x5)` is evaluated lazily:

- x0-x4 are O(1) combinatorial invariants that resolve >99.9% of cases at 152 dual
  vertices
- x5 (full BFS canonical form) is O(n) but needed <0.1% of the time
- The automorphism group falls out as a byproduct of x5

## Deliverable path

1. **Seeds** -- DONE. `src/Seeds.hs` with all 49 seed graphs decoded and validated.
2. **Core types + expansion/reduction operations** — DONE. `src/Expansion.hs`.
3. **Canonicity test** — DONE. `src/Canonical.hs`. 5-tuple cascade, Rule 2 orbit
   filtering, direction-specific inverse test.
4. **Forest traversal** — `unfold` with the `SearchMonad` abstraction,
   strategy-parameterized.
5. **Bounding lemmas** — the 5+ lemmas as a `[Bound]` list, trivially extensible.
6. **Tests against the C implementation** — DONE through C60. All 1812 isomers correct.

## Reference documents

| File | Purpose |
|------|---------|
| `doc/C-CODE-ALGORITHM.md` | **AUTHORITATIVE.** Complete algorithm as implemented in the C code, reverse-engineered over multiple days. Covers the canonical test pipeline (4-stage colour cascade), Rule 2 marking systems, bounding lemmas, direction conventions, and function mappings. This is the primary reference for understanding what the code actually does. |
| `doc/BUCKYGEN-MAP.md` | **READ THIS REGULARLY.** Comprehensive mapping between paper, C code, and Haskell. Use this to locate specific C functions/line numbers and their paper references. Covers direction conventions, colour functions, BFS, bounding lemmas, edge traversal chains, and debugging notes. |
| `../bucky_pseudocode.tex` | Formal functional spec in Haskell notation |
| `../bucky_overview.tex` | Structural analysis mapping C code to math |
| `../buckygen_paper_source/` | Original Brinkmann-Goedgebeur-McKay 2012 paper (has significant gaps — see below) |
| `../buckygen.c` | The 15,000-line C implementation (source of truth for seed data) |

## Working with the algorithm

**For understanding what the algorithm does**: consult `doc/C-CODE-ALGORITHM.md`. It documents
the complete algorithm structure including details not in the paper:
- Section 4: Rule 1 canonical test pipeline (4-stage colour cascade with exact starting edges)
- Section 5: Rule 2 orbit filtering (marking systems, or_same_edge_found, direction-swapping)
- Section 6: Bounding lemmas (Lemmas 1-6 with C code line numbers)
- Section 7: `find_extensions_fuller()` full pseudocode
- Section 8: `is_best_L0_reduction()` detailed algorithm

**For locating specific C functions/line numbers and paper references**: consult `doc/BUCKYGEN-MAP.md`.
It contains:
- Section 0: Direction convention (DRight/DLeft ↔ use_prev/use_next, CW/CCW mapping)
- Section 3: Isomorphism rejection (the 6-tuple cascade, Rules 1 and 2)
- Section 4: Colour functions (exact C→Haskell edge traversal chains for each reduction type)
- Section 5: BFS canonical form (testcanon/testcanon_mirror ↔ bfsCanonicalForm)
- Section 6: Bounding lemmas (Lemmas 1-6 with paper/C/Haskell references)
- Section 8: Reduction validation (allReductions, each `is_best_*` function)
- Section 11: File/function index (line numbers for every key function)
- Section 12: Known issues and debugging notes

**Warning about the paper**: The buckygen paper (Brinkmann-Goedgebeur-McKay 2012) has
significant gaps that make it insufficient for a correct reimplementation on its own.
Key omissions include: direction-specific canonical testing, L0 edge-reversal handling,
colour function specifications, and the or_same_edge_found mechanism. The C code and
`doc/C-CODE-ALGORITHM.md` are the authoritative sources.

## Key algorithm concepts

### Expansion operations

- **L_i (straight)**: path of length i+1 between two 5-vertices, adds i+2 vertices.
  L0 always produces adjacent 5-vertices (never IPR).
- **B_{i,j} (bent)**: path with a turn, adds i+j+3 vertices, length i+j+2.
- **F (nanotube ring)**: adds a ring of degree-6 vertices around a (5,0)-nanotube
  equator. Only for the nanotube family from C30.

### Bounding lemmas

1. Every reducible non-IPR fullerene has a reduction of length <= 2.
2. If G has a reduction of length d <= 2, every child has length at most d+2.
3. If G has an L0 reduction, all canonical children have length at most 2.
4. Three independent L1 reductions with disjoint 5-vertex sets => length at most 2.
5. Geometric bound: relates minimum 5-vertex distance to minimum vertex count.

Bounds compose via `minimum` — soundness is preserved.

### Dual representation

Buckygen works entirely in the dual. A fullerene with n primal vertices has
n/2 + 2 dual vertices. The dual is a triangulation where every vertex has degree
5 or 6, with exactly 12 degree-5 vertices. Adjacency lists store neighbours in
cyclic planar order (essential for patch-replacement operations).

## Context Management

After every significant step or decision, update doc/PROGRESS.md with:
- What was just completed
- Current approach/strategy
- Key decisions made and why
- What's next / remaining TODOs
- Any important file paths or details

Do this proactively without being asked. This file is our safety net if context is lost.
