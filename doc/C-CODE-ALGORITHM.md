# Buckygen C Code: Algorithm Structure and Canonical Construction

This document captures the complete algorithm as implemented in `buckygen.c`
(~15,000 lines), based on multi-day reverse-engineering. The purpose is to
document what we learned from the C code that was NOT obvious from the paper,
enabling future paper-vs-code comparison.

---

## 1. Overall Architecture

The C code generates all fullerene isomers with a given number of vertices using
**McKay's canonical construction path** method. The core idea:

- Each fullerene has a unique "canonical parent" obtained by applying its
  minimum reduction.
- Only expand from the canonical parent → each isomorphism class generated
  exactly once.
- No hash tables, no isomorphism checking, no deduplication. Correctness
  relies entirely on the canonical test.

### Main entry point: `scansimple_fuller(int nbtot, int nbop)`

```
if nv == maxnv:
    got_one(nbtot, nbop)      // Output/count the fullerene
    return

find_extensions_fuller()      // Find all expansion sites, apply Rule 2 marking

for each unmarked expansion site:
    extend_*()                // Apply expansion (adds 2-7 new vertices)
    if is_best_*_reduction():  // Rule 1: canonical test (3-colour cascade)
        if test_canon:
            if canon_edge_oriented():  // Full BFS canonical form test
                scansimple_fuller()     // Recurse into child
        else:
            scansimple_fuller()   // Single best colour → skip BFS (fast path)
    reduce_*()                // Undo expansion
```

### Parameters `nbtot` and `nbop`

- `nbtot`: Total number of canonical numberings (= |Aut(G)|)
- `nbop`: Number of orientation-preserving canonical numberings (= |Aut+(G)|)
- These are used to correctly count graphs modulo symmetry
- Propagated from parent to child through the canonical test

---

## 2. The Three Seeds

| Seed | Dual vertices | Primal vertices | Symmetry | Role |
|------|---------------|-----------------|----------|------|
| C20 (dodecahedron) | 12 | 20 | Ih | Root for most fullerenes |
| C28 (Td) | 16 | 28 | Td | Irreducible, not reachable from C20 |
| C30 (D5h) | 17 | 30 | D5h | Root for (5,0) nanotube family |

C20 and C28 are expanded via L and B operations.
C30 is expanded via F (nanotube ring) only.

---

## 3. Expansion Operations

| Type | New vertices | Path length | C extend function | C reduce function |
|------|-------------|-------------|--------------------|--------------------|
| L0 | 2 | 1 | `extend_L0` (line 1184) | `reduce_L0` (line 1338) |
| L_i (i>=1) | i+2 | i+1 | `extend_straight` (line 1395) | `reduce_straight` (line 1541) |
| B_{0,0} | 3 | 2 | `extend_bent_zero` (line 1602) | `reduce_bent_zero` (line 1826) |
| B_{i,j} (i+j>0) | i+j+3 | i+j+2 | `extend_bent` (line 1886) | `reduce_bent` (line 2160) |
| F (ring) | 5 | - | `extend_F` | `reduce_F` |

### Expansion priority (shorter reductions block longer ones)

The C code enforces a strict priority ordering:
- L0 (length 2) blocks everything longer
- L1 (length 3) and B00 (length 3) block length >= 4
- L2 (length 4) and B10/B01 (length 4) block length >= 5
- And so on

This is checked inside `is_best_straight_reduction` and `is_best_bent_reduction`:
```c
if((pathlength > 3 && has_B00_reductions())
    || (pathlength > 4 && has_B10_reductions())
    || (pathlength > 5 && has_bent_reductions(pathlength - 4)))
    return 0;  // Shorter reduction exists → prune
```

In our Haskell code, this is handled by `x0 = reductionLength kind` in `canonOrd` —
shorter reductions naturally sort first.

---

## 4. Rule 1: The Canonical Test Pipeline

The canonical test checks: "Is this expansion the canonical parent of the child?"
It's a multi-stage cascade where each stage filters more candidates.

### Stage 0: Reduction type priority (x0 in our code)

L0 has the shortest reduction (length 2). If the child has an L0 reduction,
no longer expansion (L1, B00, etc.) can be canonical. This is checked by:
- `is_best_straight_reduction`: checks `has_B00_reductions()` etc. for shorter types
- `is_best_bent_reduction`: same

In our code: `x0 = reductionLength kind` ensures shorter types always win.

### Stage 1: First colour — 5-bit degree pattern (x2 in our code)

The C function `get_colour_next_5(edge)` / `get_colour_prev_5(edge)` computes
a 5-bit integer by examining 5 consecutive edges and checking if each endpoint
has degree 5 (bit=1) or degree 6 (bit=0).

```c
static int get_colour_next_5(EDGE *edge) {
    EDGE *e = edge;
    int i, colour = 0;
    for(i = 4; i >= 0; i--) {
        colour |= ((degree[e->end] == 5) << i);
        e = e->next;  // next = CCW direction
    }
    return colour;
}
```

For L0, both endpoints of the reduction are checked (test_edge1 and test_edge2).
The colour is `MAX(colour1, colour2)` as first comparison, then `MIN(colour1, colour2)`
as tiebreaker.

For each L0 reduction in the child, the same colours are computed. If ANY other
L0 reduction has strictly better (higher) colour → reject. If tied, continue to
Stage 2.

**In our code**: `secondColour` computes `(negate (max c1 c2), negate (min c1 c2))`.
The negation is because we want SMALLER = MORE canonical, while C uses MAX = BEST.

**Starting edges for L0 colour computation:**
```
use_next (DLeft):  test_edge->next->next->invers->next
  = advanceCW(g, u, v, deg-2), then prevCW(w, u)  [2 CCW from v at u, cross, 1 CCW]
  → colourCCW5

use_prev (DRight): test_edge->prev->prev->invers->prev
  = advanceCW(g, u, v, 2), then nextCW(w, u)  [2 CW from v at u, cross, 1 CW]
  → colourCW5
```

### Stage 2: Third colour — path colour (x3 in our code)

The C function `get_path_colour_next(edge)` / `get_path_colour_prev(edge)` computes
a longer (7-bit) invariant by following a straight-ahead path and recording the
degree pattern.

Called from `is_best_third_colour()` (line 10138). Same MAX logic: if any other
reduction has better path colour → reject.

**Starting edges for L0 path colour:**
```
use_next (DLeft):  test_edge->prev->prev = 2 CW from v at u
  = advanceCW(g, u, v, 2), then pathColour DLeft

use_prev (DRight): test_edge->next->next = 2 CCW from v at u
  = advanceCW(g, u, v, deg-2), then pathColour DRight
```

**In our code**: `straightThirdColour` computes `negate (max pc1 pc2)`.

### Stage 3: Full BFS canonical form (x4 in our code)

If colours 1-3 are tied among multiple candidates, the full BFS canonical numbering
is computed via `canon_edge_oriented()` (line 4089) or the optimized
`canon_edge_oriented_short()` for smaller graphs.

This computes the lexicographically minimal BFS numbering over all candidate
starting edges. Returns 1 if the test edge is among the canonical edges (i.e.,
gives a BFS form that's minimal).

**Optimization**: `test_canon_edge_short()` (line 10099):
- If only 1 candidate edge remains → skip BFS entirely (trivially canonical)
- If nv < NV_CANON_SHORT + 20 → skip BFS (graph too small for it to matter)
- Otherwise → full BFS canonical form

**In our code**: `bfsCanonicalForm g u v dir` computes the BFS code starting
from edge (u→v) walking in direction dir.

### Key difference: cascade early-exit vs 5-tuple comparison

The C code exits early at each stage — if a strictly better colour is found
at Stage 1, Stages 2-3 are never computed. Our code relies on Haskell's
lazy evaluation to achieve the same effect: `canonOrd` returns a 5-tuple,
and comparing tuples lexicographically means later components are only
evaluated if earlier ones tie.

---

## 5. Rule 2: Orbit Filtering (Expansion Equivalence Classes)

Rule 2 prevents trying multiple expansions from the same parent that produce
isomorphic children. The C code uses an intricate marking system.

### 5.1 The numbering array

```c
static EDGE *numbering[2*MAXE][MAXE];
```

This stores all canonical numberings of the current graph:
- `numbering[0..numb_pres-1]`: **Orientation-Preserving (OP)** automorphisms
- `numbering[numb_pres..numb_total-1]`: **Orientation-Reversing (OR)** automorphisms

These are computed by `canon_edge_oriented()` as a byproduct of the BFS
canonical form test.

### 5.2 Main marking loop in `find_extensions_fuller()` (line ~6784)

For each edge at position `index` in the canonical numbering:

```c
// Mark all OP-equivalent edges (same direction)
for(j = 1; j < npres; j++)
    MARKLO(numbering[j][index]);

// Mark all OR-equivalent edges (note: direction will be swapped)
for(; j < numb_total; j++) {
    MARKLO(numbering[j][index]);
    if(numbering[j][index] == e)
        or_same_edge_found = 1;  // Edge is fixed by an OR automorphism!
}
```

The MARKLO marks prevent the main loop from processing the same orbit twice.
Edges at position `index` under all automorphisms map to the same edge under
the automorphism. Already-marked edges are skipped.

### 5.3 Direction filtering via `or_same_edge_found`

```c
if(!or_same_edge_found)
    find_L0_extensions_next(e, ...);   // Only if NOT fixed by OR aut
find_L0_extensions_prev(e, ...);       // Always tried
```

**Key insight**: If an edge is fixed by an orientation-reversing automorphism,
then expanding in direction `next` and `prev` from that edge produces isomorphic
children. The code only tries `prev` in this case.

**In our code**: `applyAutToExp` with a `Reversing` automorphism flips the
direction. If an OR automorphism maps (u,v) to itself, then:
- `applyAutToExp(Reversing, Exp L0 (u,v) DRight) = Exp L0 (u,v) DLeft`
- These are in the same orbit → `filterByRule2` keeps only one

### 5.4 Expansion-specific marking: `mark_edges_L0()` (line 4976)

When an L0 expansion is accepted, its orbit equivalents are marked:

```c
static void mark_edges_L0(EDGE *e, int use_next, int numb_total, int npres) {
    if(use_next) MARK_L0_NEXT(e->label);
    else         MARK_L0_PREV(e->label);

    // Mark OP automorphism images (SAME direction)
    for(i = 1; i < npres; i++) {
        if(use_next) MARK_L0_NEXT(numbering[i][index]->label);
        else         MARK_L0_PREV(numbering[i][index]->label);
    }

    // Mark OR automorphism images (SWAPPED direction!)
    for(; i < numb_total; i++) {
        if(use_next) MARK_L0_PREV(numbering[i][index]->label);  // ← PREV, not NEXT
        else         MARK_L0_NEXT(numbering[i][index]->label);  // ← NEXT, not PREV
    }
}
```

**Critical detail**: For OR automorphisms, the direction is **swapped**.
An orientation-reversing automorphism maps CW walks to CCW walks and vice versa.
So an expansion in the `next` direction, under an OR automorphism, becomes an
expansion in the `prev` direction.

**In our code**: This direction swap is handled by `applyAutToExp`:
```haskell
applyAutToExp (Aut sigma Reversing) (Exp kind (u, v) dir) =
    Exp kind (sigma IM.! u, sigma IM.! v) (flipDir dir)
```

### 5.5 The "other triple" (Section 2.2 of the paper)

For each expansion type, there's a second expansion that produces the same child:

| Type | "Other triple" relationship |
|------|---------------------------|
| L0 | Edge (u,v): other = edge (u',v') at other end, same direction |
| L_i | Same: other endpoint of the straight path, same direction |
| B_{0,0} | Other endpoint of bent path, **flipped** direction |
| B_{i,j} | Other endpoint, **swapped** (i,j) and **flipped** direction |

**C code implementation**: The `find_L0_extensions_next/prev` functions mark
the "other" edge when finding an extension:
```c
mark_edges_L0(other_startedge, 1, numb_total, npres);  // line 5506
```

**In our code**: `computeOtherTriple` computes the other expansion explicitly.
`filterByRule2.fullOrbit` includes both automorphism images AND other-triple
images.

### 5.6 Independent mark systems

The C code uses SEPARATE mark arrays for each expansion type and direction:

| Mark system | Purpose |
|-------------|---------|
| `MARKLO` / `MARKHI` | Main edge marks (prevent revisiting orbits) |
| `MARK_L0_NEXT` / `MARK_L0_PREV` | L0-specific expansion marks |
| `MARK_DOUBLE_NEXT` / `MARK_DOUBLE_PREV` | L_i (straight) expansion marks |
| `MARK_B00_NEXT` / `MARK_B00_PREV` | B_{0,0} expansion marks |
| `MARK_BENT_NEXT` / `MARK_BENT_PREV` | B_{i,j} expansion marks |

Each uses epoch-based clearing (`markvalue += 2`) to avoid resetting arrays.

---

## 6. Bounding Lemmas (Pruning the Search Tree)

The bounding lemmas restrict the maximum expansion length needed at each node.
They're checked before generating expansion sites.

### Lemma 1 (line 12139)
Every reducible non-IPR fullerene has a reduction of length ≤ 2.
(Implemented as: if graph has no IPR violation, skip long expansions.)

### Lemma 2 (line 12183)
If G has a reduction of length d ≤ 2, every child has length ≤ d+2.
(The C code tracks `straight_length` and uses it to bound expansion length.)

### Lemma 3 (line 12192)
If G has an L0 reduction, all canonical children have length ≤ 2.
(Saves ~26% of expansion sites in practice.)

### Lemma 4
Three independent L1 reductions with disjoint 5-vertex sets ⟹ length ≤ 2.

### Lemma 5 (geometric bound, line 12142)
`maxlength = (nv - 16) / 2 + 3` for non-IPR. Relates pentagon distance to vertex count.

### Lemma 6 (IPR bound)
For IPR fullerenes: `maxlength = (nv - 26) / 2 + 3`.

These compose via `min()` — the tightest bound wins. In our Haskell design,
bounds are `[DualGraph -> Maybe Int]` that compose via `minimum`.

---

## 7. `find_extensions_fuller()` — The Full Extension Finding Loop

This is the most complex function (~260 lines, starting at line 6784).
It combines orbit marking and extension finding.

### Pseudocode (simplified):

```
RESETMARKS
compute degree_5_vertices[]  // exactly 12 vertices with degree 5

for each edge e at each degree-5 vertex:
    if ISMARKEDLO(e): continue  // already in a processed orbit

    // Mark all OP-equivalent edges
    for j = 1..numb_pres-1:
        MARKLO(numbering[j][e.index])

    // Mark all OR-equivalent edges, detect self-mapping
    or_same_edge_found = false
    for j = numb_pres..numb_total-1:
        MARKLO(numbering[j][e.index])
        if numbering[j][e.index] == e:
            or_same_edge_found = true

    // Try L0 extensions from this edge
    if can_perform_L0(e):
        find_L0_extensions_prev(e)
        if !or_same_edge_found:
            find_L0_extensions_next(e)

    // Try L1 extensions from this edge
    if can_perform_L1(e):
        find_L1_extensions_prev(e)
        if !or_same_edge_found:
            find_L1_extensions_next(e)

    // Try B00 extensions, straight extensions, bent extensions...
    // (same or_same_edge_found pattern for each type)
```

### Key observations:

1. **Edge-centric iteration**: The outer loop iterates over edges at degree-5
   vertices, not over vertex pairs. This naturally handles the "two directions
   per edge" issue.

2. **Marking before trying**: Automorphic edges are marked BEFORE extensions
   are tried. This prevents the same orbit from being processed twice.

3. **Direction filtering is per-edge**: `or_same_edge_found` is set once per
   edge and applies to ALL expansion types from that edge.

4. **Within find_L0_extensions_next/prev**: These functions iterate over
   neighbors to find valid L0 reduction sites. When one is found, they mark
   equivalent sites via `mark_edges_L0()` before continuing the search.

---

## 8. Graph Representation in C

### Edge structure (doubly-linked circular lists)

```c
typedef struct e {
    int start, end;         // endpoints
    struct e *prev, *next;  // CW/CCW around start vertex
    struct e *invers;       // reverse edge (end→start)
    int index;              // position in canonical numbering
    int label;              // for marking systems
    int mark;               // epoch-based mark value
} EDGE;

static EDGE *firstedge[MAXN];  // first edge for each vertex
static int degree[MAXN];       // vertex degrees
```

### Direction convention (CRITICAL)

| C code | Haskell | Meaning |
|--------|---------|---------|
| `e->next` | `prevCW g v w` | CCW successor of w around v |
| `e->prev` | `nextCW g v w` | CW successor of w around v |
| `use_next` | `DLeft` | CCW walk direction |
| `use_prev` | `DRight` | CW walk direction |

C's `->next` is CCW; Haskell's neighbor lists are CW-ordered. So C's `next`
corresponds to going backwards in Haskell's neighbor list.

### Edge traversal chains (most common patterns)

```
C: e->next->next->end           = advanceCW(g, u, v, deg(u)-2)  [2 CCW]
C: e->prev->prev->end           = advanceCW(g, u, v, 2)         [2 CW]
C: e->invers->next->next->next  = straight-ahead CCW through v
C: e->invers->prev->prev->prev  = straight-ahead CW through v
```

---

## 9. `is_best_L0_reduction()` — Detailed Algorithm (line 10292)

This is the canonical test for L0 expansions. Parameters:
- `test_edge1`: edge nv→nv+1 (between the two NEW vertices)
- `test_edge2`: edge nv+1→nv (inverse of test_edge1)
- `use_next`: direction flag (1 = CCW/DLeft, 0 = CW/DRight)

### Step 1: Compute colour pair for test edges

```c
if(use_next) {
    colour_testedge1 = get_colour_next_5(test_edge1->next->next->invers->next);
    colour_testedge2 = get_colour_next_5(test_edge2->next->next->invers->next);
}
best_colour = MAX(colour_testedge1, colour_testedge2);
best_colour_two = MIN(colour_testedge1, colour_testedge2);
```

### Step 2: Compare against ALL other L0 reductions

Loop through all 12 degree-5 vertices and their edges:
```c
for(i = 0; i < 12; i++) {
    e = firstedge[degree_5_vertices[i]];
    do {
        if(e->start < e->end && has_L0_reduction(e, &direction_bitvector)) {
            // Check next direction
            if(bitvector & NEXT) {
                colour_tmp = get_colour_next_5(e->next->next->invers->next);
                colour_tmp_other = get_colour_next_5(e->invers->next->next->invers->next);
                if(MAX(colour_tmp, colour_tmp_other) > best_colour)
                    return 0;  // REJECT: strictly better L0 reduction exists
                // Tiebreaker: compare second-best colours
                // ... collect tied candidates into good_next_tmp[]
            }
            // Same for prev direction
        }
        e = e->next;
    } while(e != firstedge[degree_5_vertices[i]]);
}
```

### Step 3: Third colour test and BFS

```c
return is_best_third_colour(test_edge1, test_edge2, ...,
    good_next_tmp, num_good_next_tmp, ...);
```

This computes path colours via `get_path_colour_next/prev()` for the
remaining tied candidates, then delegates to `test_canon_edge_short()`
which may call `canon_edge_oriented_short()` for full BFS.

### `has_L0_reduction()` (line 8990)

Checks if edge e connects two degree-5 vertices AND the 4 surrounding
vertices (2 CCW from e at u, 2 CW from e at u, same at v) all have degree 6.

```c
degree[e->end] == 5
  && degree[e->prev->prev->end] == 6   // 2 CW from v at u
  && degree[e->invers->prev->prev->end] == 6  // 2 CW from u at v
  && degree[e->next->next->end] == 6   // 2 CCW from v at u
  && degree[e->invers->next->next->end] == 6  // 2 CCW from u at v
```

---

## 10. Processing Flow for L0 in `scansimple_fuller()`

Lines 12213-12270:

```c
for(i = 0; i < num_ext_L0; i++) {
    startedge = ext_L0[i].startedge;     // edge at expansion site
    use_next = ext_L0[i].use_next;        // direction

    extend_L0(startedge, use_next);       // Add 2 new vertices
    // After extend: last_edge_L0 = edge nv→nv+1

    if(is_best_L0_reduction(last_edge_L0, last_edge_L0->invers, use_next,
                            good_next, &num_good_next, &can_edges_next,
                            good_prev, &num_good_prev, &can_edges_prev)) {
        test_canon = 1;
        if(num_good_next + num_good_prev == 1) {
            test_canon = 0;      // Only 1 candidate → skip BFS
            xnbtot = xnbop = 1;
        }
        if(!test_canon ||
           canon_edge_oriented(good_next, num_good_next, can_edges_next,
                              good_prev, num_good_prev, can_edges_prev,
                              degree, numbering, &xnbtot, &xnbop)) {
            if(nv < maxnv - 1)
                add_straight_extension_to_list(...);
            scansimple_fuller(xnbtot, xnbop);  // RECURSE
        }
    }

    reduce_L0();  // Undo expansion
}
```

### Important: `last_edge_L0`

After `extend_L0()`, `last_edge_L0` is set to the edge between the two
NEW vertices (nv → nv+1). This edge is passed to `is_best_L0_reduction`
as the test edge.

In our Haskell code, `isCanonical` checks if any reduction at the new
vertices has the minimum canonOrd. The "new vertices" for L0 are
`[parentNV, parentNV+1]`.

---

## 11. Mapping Between C Functions and Haskell Functions

### Canonical test

| C function | Haskell function | Purpose |
|-----------|-----------------|---------|
| `is_best_L0_reduction` | `isCanonical` + `canonOrd` | L0 canonical test |
| `is_best_straight_reduction` | `isCanonical` + `canonOrd` | L_i canonical test |
| `is_best_bent_zero_reduction` | `isCanonical` + `canonOrd` | B00 canonical test |
| `is_best_bent_reduction` | `isCanonical` + `canonOrd` | B_{i,j} canonical test |
| `get_colour_next_5` | `colourCCW5` | 5-bit CCW degree pattern |
| `get_colour_prev_5` | `colourCW5` | 5-bit CW degree pattern |
| `get_path_colour_next` | `pathColour DLeft` | CCW path colour |
| `get_path_colour_prev` | `pathColour DRight` | CW path colour |
| `canon_edge_oriented` | `bfsCanonicalForm` | Full BFS canonical numbering |
| `testcanon` | `bfsCanonicalForm _ _ _ DLeft` | CCW BFS |
| `testcanon_mirror` | `bfsCanonicalForm _ _ _ DRight` | CW BFS |

### Rule 2

| C mechanism | Haskell function | Purpose |
|-------------|-----------------|---------|
| `MARKLO` loop in `find_extensions_fuller` | `filterByRule2` | Orbit representatives |
| `mark_edges_L0` (OP: same dir) | `applyAutToExp Preserving` | OP aut images |
| `mark_edges_L0` (OR: swap dir) | `applyAutToExp Reversing` | OR aut images (dir flip) |
| `or_same_edge_found` | Implicit in `applyAutToExp` | Edge fixed by OR aut |
| `mark_edges_L0(other_startedge)` | `computeOtherTriple` | Other-triple marking |

### Expansion/reduction

| C function | Haskell function |
|-----------|-----------------|
| `extend_L0` | `applyL0` |
| `reduce_L0` | `reduceL0` |
| `extend_straight` | `applyStraight` |
| `reduce_straight` | `reduceStraight` |
| `extend_bent_zero` | `applyBentZero` |
| `reduce_bent_zero` | `reduceBentZero` |
| `extend_bent` | `applyBent` |
| `reduce_bent` | `reduceBent` |

---

## 12. Known Discrepancies and Open Questions

### 12.1 Dedup catches (4,268 at C60)

Our code generates 1812/1812 correct C60 isomers but has 4,268 "dedup catches"
where the same child is generated via two different expansions from the same
parent. Both expansions survive Rule 2 and the canonical test.

**Investigation findings:**
- All 4,268 cases involve L0 expansions
- The parents have trivial automorphism groups (|Aut|=1)
- The child has multiple L0 reductions tied at minimum canonOrd
- The expansion's "inverse" reduction (at new vertices nv, nv+1) has the same
  canonOrd as some other L0 reduction at old vertices
- The reduction at (nv, nv+1) does NOT reconstruct the parent graph — it's a
  coincidental L0 site with the same canonical ordering

**Open question**: Does the C code also generate these children twice, or does
some mechanism we haven't found prevent it? The C code has no hash-based dedup,
so either the canonical test prevents duplicates or it doesn't.

### 12.2 Precise isInverseRed filter

We changed `isInverseRed` to require BOTH vertices of the reduction to be in
the new vertex range `[parentNV .. parentNV + numNew - 1]`. This matches the
C code's approach of using edge `nv→nv+1` as the test edge.

The original (loose) filter required only one vertex to be new, which created
false positive matches where a new degree-5 vertex happened to be adjacent to
an old degree-5 vertex (creating an L0 site that's not the expansion's inverse).
This change reduced dedup catches from 10,109 to 4,268.

### 12.3 The secondColour tiebreaker

For L0, the C code uses a two-level colour comparison:
1. `best_colour = MAX(colour1, colour2)` — primary comparison
2. `best_colour_two` — secondary tiebreaker (second-best colour)

Both colours are computed from the two endpoints of the L0 pair. Our code
computes `(negate (max c1 c2), negate (min c1 c2))`, which should be equivalent.

However, the C code's tiebreaker logic is more nuanced: when comparing against
other L0 reductions, it checks each endpoint separately. An edge qualifies for
the "good" list only if one of its colours matches `best_colour` AND the other
matches or exceeds `best_colour_two`. This could produce different filtering
than our simple (max, min) comparison.

### 12.4 Extension recording for lookahead

The C code records accepted extensions for later use:
- `add_straight_extension_to_list()` in `scansimple_fuller()`
- `add_straight_colour_2_L0_extension_to_list()` in `is_best_L0_reduction()`
- These feed into the bounding lemma system for subsequent levels

Our code doesn't yet implement this optimization. It shouldn't affect correctness,
only performance.

---

## 13. The Automorphism Group Computation

The automorphism group is computed as a byproduct of BFS canonical form.
`canon_edge_oriented()` tries all starting edges and orientations:

1. For each starting edge (u,v) and orientation (next/prev):
   - Compute BFS numbering
   - Compare with current best
   - If equal: found an automorphism
   - If better: new best; reset automorphism count

2. Output:
   - `numbering[0..numb_pres-1]`: OP automorphisms
   - `numbering[numb_pres..numb_total-1]`: OR automorphisms
   - `nbtot = numb_total`: total number of canonical numberings
   - `nbop`: number of OP canonical numberings

In our code: `canonicalBFSAndGroup` performs the same computation, returning
a `CanonResult` with automorphism list and orientation flags.

---

## 14. Output: `got_one(nbtot, nbop)` (line 8847)

When `nv == maxnv`, the completed fullerene is output:
- `nbtot`: number of automorphisms (= |Aut(G)|)
- `nbop`: number of orientation-preserving automorphisms (= |Aut+(G)|)
- These are used for counting by symmetry group
- An optional FILTER macro can add user-defined output filtering
- Spiral computation can be optionally performed for validation

The function increments `numb_total_out` (total graphs found).
There is NO deduplication at this stage — the canonical construction path
guarantees each isomorphism class arrives exactly once.
