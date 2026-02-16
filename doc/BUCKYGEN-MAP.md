# Buckygen: Paper-to-Code Structural Map

Machine-readable mapping between the buckygen paper (fullerene-generation.tex),
the C implementation (buckygen.c), and the Haskell reimplementation.

---

## 0. Direction Convention (Critical — Read First)

The C code uses doubly-linked edge records with `->next` (CCW) and `->prev` (CW).
The Haskell code stores neighbor lists in **CW order**.

| C code | Haskell | Meaning |
|--------|---------|---------|
| `e->next` | `prevCW g v w` | CCW successor of w around v |
| `e->prev` | `nextCW g v w` | CW successor of w around v |
| `e->invers` | swap endpoints | reverse directed edge |
| `use_next` | `DLeft` | CCW walk direction |
| `use_prev` | `DRight` | CW walk direction |
| `testcanon` | `bfsCanonicalForm _ _ _ DLeft` | CCW BFS (->next) |
| `testcanon_mirror` | `bfsCanonicalForm _ _ _ DRight` | CW BFS (->prev) |

**Mnemonic**: C's "next" is CCW, Haskell's neighbor lists are CW. So C's `->next`
from vertex v past neighbor w = Haskell's `prevCW g v w` (going backwards in
the CW list = going forwards in the CCW direction).

**Edge traversal translation**: Given C code `startedge` from u to v:
- `startedge->next` = the edge from u whose end is `prevCW g u v` (one CCW step)
- `startedge->prev` = the edge from u whose end is `nextCW g u v` (one CW step)
- `startedge->next->next->end` = `advanceCW g u v (deg g u - 2)` (2 CCW = deg-2 CW)
- `startedge->prev->prev->end` = `advanceCW g u v 2` (2 CW steps past v)
- `startedge->invers` = directed edge from v to u
- `startedge->invers->next->next->next` = straight-ahead from u through v in CCW mode

**Straight-ahead step**: For a degree-d vertex entered from one neighbor,
straight-ahead = skip ceil(d/2)-1 neighbors in the walk direction after the
entry edge. For degree 6: `->invers->next->next->next` (CCW) or
`->invers->prev->prev->prev` (CW). In Haskell: `straightAhead g dir atVtx fromVtx`.

---

## 1. Construction Algorithm

### Paper: Section 2.1 (lines 166-273)

| Topic | Paper lines | Description |
|-------|------------|-------------|
| Expansion/reduction defs | 166-170 | Patch replacement: expansion replaces fragment with larger fragment. G' = child, G = parent. |
| Infinite families needed | 172-177 | No finite set of patch replacements suffices for all fullerenes (Brinkmann et al. 2006). |
| Three expansion families | 178-188 | L_i (straight), B_{i,j} (bent), F (nanotube ring). From Hasheminezhad et al. 2008. Mirror image of L must also be considered. All faces drawn completely or labeled f_k/g_k must be distinct. |
| Dual representation | 197-216 | Figure 2. 5-vertices = solid white, 6-vertices = solid black, shaded = either. |
| Irreducible fullerenes | 218-237 | C20 (dodecahedron), C28(Td), C30(D5h). Type-(5,0) nanotubes from C30 via F. All others are "reducible". |
| **Theorem 1.1** | 229-233 | Every fullerene except C28(Td) and (5,0)-nanotubes is constructible from C20 via L and B. Seeds: {C20, C28, C30}. |
| Algorithm overview | 235-237 | Apply L,B from C20 and C28(Td); separately generate (5,0)-nanotubes from C30 via F. |

### C code: Expansion operations

| Operation | C function | Lines | New verts | Path length | Distance between 5-verts |
|-----------|-----------|-------|-----------|-------------|------------------------|
| L0 expand | `extend_L0` | 1184-1337 | 2 | 1 | 1 (adjacent) |
| L0 reduce | `reduce_L0` | 1338-1393 | -2 | | |
| L_i (i≥1) expand | `extend_straight` | 1395-1540 | i+2 | i+1 | i+1 |
| L_i reduce | `reduce_straight` | 1541-1600 | -(i+2) | | |
| B_{0,0} expand | `extend_bent_zero` | 1602-1825 | 3 | 2 | 2 |
| B_{0,0} reduce | `reduce_bent_zero` | 1826-1884 | -3 | | |
| B_{i,j} expand | `extend_bent` | 1886-2159 | i+j+3 | i+j+2 | i+j+2 |
| B_{i,j} reduce | `reduce_bent` | 2160-2242 | -(i+j+3) | | |

### Haskell: Expansion.hs

| Operation | Function | Key helpers |
|-----------|----------|-------------|
| L0 expand | `applyL0` | `computeStraightPath g (u,v) dir 3` → PathInfo with mainPath, parallelPath |
| L0 reduce | `reduceL0` | Inverse of applyL0 |
| L_i expand | `applyStraight` | `computeStraightPath g (u,v) dir (i+3)` |
| L_i reduce | `reduceStraight` | |
| B_{0,0} expand | `applyBentZero` | `computeBentZeroPath g (u,v) dir` → specialized 4-vertex path |
| B_{0,0} reduce | `reduceBentZero` | |
| B_{i,j} expand | `applyBent` | `computeBentPath g (u,v) dir i j` |
| B_{i,j} reduce | `reduceBent` | |
| F expand | `applyRing` | `findNanotubeRing g` → (ring, outer) |
| F reduce | `reduceRing` | |

### Expansion geometry — L0 in detail

An L0 expansion at edge (u→v) with direction `use_next` (C code, line 12213-12226):

```
Before:  u (deg-5) ——— v (deg-5)
         |             |
    parallel_path[0]   parallel_path[2]

Path computation (C code line 12214-12224):
  path[0] = startedge->start = u
  path[1] = first straight-ahead = intermediate deg-6 vertex
  path[2] = second straight-ahead = v (the other deg-5 vertex)
  parallel_path[j] = startedge->prev->end (for use_next)
                   = startedge->next->end (for use_prev)

After:   u (deg-6) — nv (deg-5) — nv+1 (deg-5) — v' (deg-6)
```

New vertices nv and nv+1 are the new 5-vertices. The old 5-vertices u and
parallel_path[2] become 6-vertices (C line 1332-1333: `replace_degree_5_vertex`).

The edge between nv and nv+1 is stored as `last_edge_L0` (C line 1235,1251).
This is the edge passed to `is_best_L0_reduction` (line 12233):
```c
is_best_L0_reduction(last_edge_L0, last_edge_L0->invers, ext_L0_use_next[i], ...)
```
So `test_edge1` = nv→nv+1, `test_edge2` = nv+1→nv.

### Expansion geometry — B00 in detail

A B00 expansion at edge (u→v) with `use_next`:
- u is degree-5, v is degree-6 (NOT degree-5 — this is a key difference from L)
- Path goes: u → v → turn → w (degree-5)
- Total path length: 2 (distance between the two 5-vertices)

B00 path traversal (`has_bent_zero_reduction`, C line 10468-10487):
```c
e = startedge;               // from u (deg-5)
if(degree[e->end] == 5) return 0;   // v must be deg-6
e = e->next->next;            // 2 CCW steps around u from v
if(degree[e->end] == 5) return 0;   // must be deg-6
e = e->invers->prev->prev;   // flip, 2 CW steps = turn
if(degree[e->end] == 6) return 0;   // w must be deg-5
e = e->invers->next->next;   // flip, 2 CCW steps
if(degree[e->end] == 5) return 0;   // must be deg-6
*enarc = e;                   // this edge is the "other end"
return 1;
```

The two test edges passed to `is_best_bent_zero_reduction`:
- `test_edge1` = edge from new vertex at u's end of the bent path
- `test_edge2` = edge from new vertex at w's end of the bent path
These are at DIFFERENT ends and use DIFFERENT colour functions.

### Expansion geometry — Straight (L_i, i≥1) in detail

Path computation for L_i follows straight-ahead steps:
```c
// scansimple_fuller, line 12273-12285
work_edge = ext_L1[i];
for(j = 0; j < 3 + 1; j++) {   // 3 = pathlength for L1
    path[j] = work_edge->start;
    if(use_next) {
        parallel_path[j] = work_edge->prev->end;  // CW neighbor
        work_edge = work_edge->invers->next->next->next;  // straight-ahead (CCW)
    } else {
        parallel_path[j] = work_edge->next->end;  // CCW neighbor
        work_edge = work_edge->invers->prev->prev->prev;  // straight-ahead (CW)
    }
}
```

For L1 (pathlength=3): path = [u, mid1, mid2, v] where u,v are degree-5.
New vertices: nv, nv+1, nv+2 (3 new vertices).

Test edges for `is_best_straight_reduction` (line 12289-12290):
```c
is_best_straight_reduction(edge_list[nv - 3][nv - 2],
    edge_list[nv - 1][nv - 2], 3, ext_L1_use_next[i], ...)
```
- `test_edge1` = edge from first new vertex (nv-3) to second new vertex (nv-2)
- `test_edge2` = edge from last new vertex (nv-1) to second-to-last (nv-2)
These are the directed edges at both ends of the path between the two new 5-vertices.

---

## 2. Isomorphism Rejection (Canonical Construction Path)

### Paper: Section 2.2 (lines 275-423)

| Topic | Paper lines | Key details |
|-------|------------|-------------|
| Goal | 277-282 | Output one per isomorphism class. McKay's canonical construction path. |
| Canonical reduction | 284-290 | Unique (up to automorphisms) reduction for each reducible dual fullerene G. The "canonical parent" is obtained by applying this reduction. |
| Equivalence of expansions | 292-304 | Two expansions equivalent iff automorphism maps patches. For L-type: same-direction ↔ orientation-preserving auto; opposite-direction ↔ orientation-reversing auto. |
| **Rule 1** | 310-311 | Only accept G' if the last step was a canonical expansion (= inverse of canonical reduction of G'). |
| **Rule 2** | 313-314 | For each G, only apply one expansion per equivalence class. |
| Representing triples | 318-331 | Reduction = (e, x, d): e = directed start edge, x = parameter set, d = direction. Each reduction has exactly **two** representing triples (one from each end of the central path). Figure 3 shows example. |
| Equivalence on triples | 333-340 | Two relations generate equivalence: (1) two triples representing same reduction; (2) automorphism mapping edges (same direction ↔ orientation-preserving, different direction ↔ orientation-reversing). |
| **6-tuple** | 351-359 | Each triple gets (x0,...,x5). Canonical reduction = smallest 6-tuple. |
| x0: length | 361-368 | Distance between two 5-vertices before reduction. L_x: x+1. B_{x,y}: x+y+2. Shorter = more canonical. |
| x1: longest straight | 370-377 | Negative of longest straight path in reduction. L_x: x1 = -x0 (trivial). B_{x,y}: x1 = -(max(x,y)+1). Always distinguishes L from B with same x0. |
| x2, x3, x4: degree strings | 379-382 | Degrees in well-defined neighborhoods of **increasing constant size**. C code implements as: (colour1_max, colour1_min) pair, then path_colour. |
| Lazy cascade | 384-393 | x_i computed ONLY for triples that minimize (x0,...,x_{i-1}). Stop when unique minimum found or two triples represent same reduction. |
| x5: BFS canonical form | 394-406 | Full graph encoding relative to edge + direction. BFS-numbering, neighbors in rotational order. Same x5 ↔ orientation-preserving isomorphism (same dir) or orientation-reversing (diff dir). |
| Automorphism group | 415-418 | Byproduct of x5 computation. |
| Complexity | 419-423 | x0-x4: O(1) each. x5: O(n). Total: O(n). |
| Performance | 426-433 | At 152 dual verts (C300): x0-x4 resolve >99.9% of cases. |

### C code: 6-tuple cascade — how it actually works

The C code does NOT compute a 6-tuple and compare. Instead it uses an
**incremental early-exit rejection** approach:

**x0 and x1 are implicit in the type dispatch.** The generation loop
(scansimple_fuller) processes expansion types in strict priority order:
```
L0 (length=1, longest_straight=1)  →  most canonical
L1 (length=2, longest_straight=2)  →  next
B00 (length=2, longest_straight=1) →  same x0 as L1, but x1 worse
L2+ (length≥3, longest_straight=length)
B_{i,j} with i+j>0 (length≥3, longest_straight<length)
```

Each `is_best_*` function for a given type checks whether any reduction of
the same or shorter type has a better colour. If so, the expansion is rejected.

**Key cross-type priority checks:**

| In function | What it checks | C lines | Effect |
|------------|---------------|---------|--------|
| `is_best_straight_reduction` (L_i, i≥1) | `has_B00_reductions()` if pathlength>3 | 11915-11918 | If B00 reduction exists and pathlength of tested reduction is >3, reject (B00 has shorter length) |
| `is_best_straight_reduction` | `has_B10_reductions()` if pathlength>4 | 11916 | Same logic for B10 |
| `is_best_straight_reduction` | `has_bent_reductions(pathlength-4)` if pathlength>5 | 11917 | Same for longer bent |
| `is_best_bent_zero_reduction` (B00) | `is_L0_or_L1_reduction(e)` for all 5-verts | 10724-10731 | If ANY L0 or L1 reduction exists, reject B00 |

### C code: is_best_L0_reduction — complete cascade (10292-10433)

This is the simplest case because L0 is the shortest possible reduction (length=1).
No need to check for shorter reductions.

**Parameters**: `test_edge1` = nv→nv+1, `test_edge2` = nv+1→nv, `use_next` = direction.
These are the representing triples of the inverse reduction in the child graph.

**Step 1 — Compute colour pair of test edges (x2):**
```c
if(use_next) {  // DLeft
    colour_testedge1 = get_colour_next_5(test_edge1->next->next->invers->next);
    colour_testedge2 = get_colour_next_5(test_edge2->next->next->invers->next);
    best_colour = MAX(colour_testedge1, colour_testedge2);
} else {        // DRight
    colour_testedge1 = get_colour_prev_5(test_edge1->prev->prev->invers->prev);
    colour_testedge2 = get_colour_prev_5(test_edge2->prev->prev->invers->prev);
    best_colour = MAX(colour_testedge1, colour_testedge2);
}
best_colour_two = (both best) ? best_colour : MIN(colour_testedge1, colour_testedge2);
```

**Step 2 — Check all other L0 reductions:**
Loop over all 12 degree-5 vertices and their edges (10332-10422):
```c
for(i = 0; i < 12; i++) {
    e = firstedge[degree_5_vertices[i]];
    do {
        if(e->start < e->end) {  // avoid double-counting undirected edge
            if(has_L0_reduction(e, &direction_bitvector)) {
                // For each valid direction (NEXT and/or PREV):
                // Compute colour of both triples (e and e->invers)
                // If colour > best_colour → return 0 (better reduction exists)
                // If colour == best_colour, check second colour
                // If second colour > best_colour_two → return 0
                // If tied, add to good_*_tmp lists
            }
        }
        e = e->next;
    } while(e != ee);
}
```

**Step 3 — Path colour (x4) via is_best_third_colour (10137-10201):**
```c
if(use_next) {  // DLeft
    colour_testedge1 = get_path_colour_next(test_edge1->prev->prev);
    colour_testedge2 = get_path_colour_next(test_edge2->prev->prev);
} else {        // DRight
    colour_testedge1 = get_path_colour_prev(test_edge1->next->next);
    colour_testedge2 = get_path_colour_prev(test_edge2->next->next);
}
// Compare against all edges in good_next_tmp/good_prev_tmp lists
// Return 0 if any has better path colour
```

**Step 4 — BFS (x5) via test_canon_edge_short → canon_edge_oriented:**
If still tied after path colour, compute full BFS canonical form.
`canon_edge_oriented` (4088-4194) takes the surviving candidate edges
and runs `testcanon`/`testcanon_mirror` to find the true minimum.

### C code: is_best_straight_reduction — cascade for L_i, i≥1 (11812-11925)

**Extra step: check for shorter reductions of different types (lines 11851-11918).**

Before colour comparison, scan all degree-5 vertices for reductions of the
SAME length. Function `has_short_straight_reduction` (9017-9074) follows
straight-ahead paths from each degree-5 vertex:
```c
while(length < pathlength_cur_best - 1 && degree[e->end] != 5) {
    e = e->invers->next->next->next;  // straight-ahead step
    length++;
}
if(degree[e->end] == 5) {
    if(length < pathlength_cur_best - 1) return 2;  // SHORTER reduction found
    // else check direction validity, return 1 if same length valid reduction found
}
```

Return values: 0 = no reduction from this edge, 1 = same-length reduction,
2 = **shorter** reduction found → immediately reject.

After scanning all L_i reductions of same length, check for shorter B reductions:
```c
if(pathlength > 3 && has_B00_reductions()) return 0;
if(pathlength > 4 && has_B10_reductions()) return 0;
if(pathlength > 5 && has_bent_reductions(pathlength - 4)) return 0;
```

This implements the priority: for a given x0 (= pathlength), L has better x1
than B (because longest_straight(L_i) = i+1 = pathlength, while
longest_straight(B_{i,j}) = max(i,j)+1 < pathlength). But a SHORTER B
(smaller x0) beats a LONGER L.

Colour functions are the same as L0:
```c
// use_next: get_colour_next_5(e->next->next->invers->next)
// use_prev: get_colour_prev_5(e->prev->prev->invers->prev)
```

### C code: is_best_bent_zero_reduction — cascade for B00 (10716-10833)

**First step: reject if any L0 or L1 reduction exists (10724-10731).**
```c
for(i = 0; i < 12; i++) {
    e = firstedge[degree_5_vertices[i]];
    for(j = 0; j < 5; j++) {
        if(is_L0_or_L1_reduction(e)) return 0;  // shorter reduction exists
        e = e->next;
    }
}
```

`is_L0_or_L1_reduction` (10442-10452): check if startedge leads to L0 (end is
degree-5) or L1 (end is degree-6, one straight-ahead step reaches degree-5).

**Colour functions for B00 are DIFFERENT from straight types:**
```c
if(use_next) {
    colour_testedge1 = get_colour_next_5(test_edge1->invers->next);   // CCW, simpler chain
    colour_testedge2 = get_colour_prev_5(test_edge2->invers->prev);   // CW, different function!
} else {
    colour_testedge2 = get_colour_next_5(test_edge2->invers->next);
    colour_testedge1 = get_colour_prev_5(test_edge1->invers->prev);
}
```
Note: test_edge1 uses `get_colour_next_5` and test_edge2 uses `get_colour_prev_5`
(or vice versa depending on direction). This is because the bent path changes
direction at the bend point, so the two endpoints see the graph from different
rotational perspectives.

**Path colour for B00 also differs** (`is_best_third_colour_bent`, 10214-10278):
```c
if(use_next) {
    colour_testedge1 = get_path_colour_next(test_edge1->invers->next->next);
    colour_testedge2 = get_path_colour_prev(test_edge2->invers->prev->prev);
} else {
    colour_testedge1 = get_path_colour_prev(test_edge1->invers->prev->prev);
    colour_testedge2 = get_path_colour_next(test_edge2->invers->next->next);
}
```

### C code: BFS canonical form

| Function | Lines | Direction | Description |
|----------|-------|-----------|-------------|
| `testcanon_first_init` | ~2380 | CCW (->next) | Compute initial BFS representation into `representation[]` |
| `testcanon_init` | ~2430 | CCW (->next) | Compare new starting edge against `representation[]`, update if better |
| `testcanon` | 2246-2346 | CCW (->next) | Full BFS comparison. Returns 0=worse, 1=automorphism, 2=better |
| `testcanon_mirror_init` | ~2490 | CW (->prev) | Initial mirror BFS |
| `testcanon_mirror` | 2459-2522 | CW (->prev) | Mirror BFS comparison |
| `canon` | 3267-3402 | Both | Full canonicity: finds minimum BFS over all starting edges and both orientations. Returns automorphism group as byproduct. |
| `canon_edge_oriented` | 4088-4194 | Both | Canonicity restricted to specific candidate edges (from colour cascade). Takes `edgelist_or[]` (orientation-preserving, use testcanon) and `edgelist_inv[]` (orientation-reversing, use testcanon_mirror). Returns 1 iff one of the first `can_edges_or/inv` candidates is optimal. |
| `canon_edge_oriented_short` | 4211-... | Both | Like `canon_edge_oriented` but early-exit after scanning a limited prefix. |

**BFS encoding scheme** (`testcanon`, line 2278-2304):
```
For each vertex in BFS order (1..nv):
  For each neighbor in cyclic order (->next for CCW, ->prev for CW),
  skipping the edge we came from:
    If neighbor already numbered → emit its number
    If neighbor is new → emit colour[vertex] (= degree + MAXN), assign next number
  Emit 0 (separator)
```
- Vertex 1 = startedge->start, vertex 2 = startedge->end
- The first entry is omitted (colour of vertex 2, always the same)
- Colours = degree + MAXN ensures colours (≥ MAXN+5) can't be confused with
  vertex numbers (≤ nv)

### Haskell: Canonical.hs functions

| Function | Description | C equivalent |
|----------|-------------|-------------|
| `canonOrd :: Reduction -> DualGraph -> CanOrd` | 5-tuple (x0, x1, x2, x3, x4) | Combined logic of `is_best_*` functions |
| `bfsCanonicalForm :: DualGraph -> Vertex -> Vertex -> Dir -> BFSCode` | BFS encoding | `testcanon` / `testcanon_mirror` |
| `allReductions :: DualGraph -> [Reduction]` | Enumerate all valid reduction sites | The scanning loops inside each `is_best_*` |
| `isCanonical :: Expansion -> Int -> DualGraph -> Bool` | Full canonicity test | `is_best_*` + `canon_edge_oriented` |
| `colourCW5 :: DualGraph -> Vertex -> Vertex -> Int` | 5-bit CW colour walk | `get_colour_prev_5` |
| `colourCCW5 :: DualGraph -> Vertex -> Vertex -> Int` | 5-bit CCW colour walk | `get_colour_next_5` |
| `pathColour :: DualGraph -> Dir -> Vertex -> Vertex -> Int` | 7-bit straight-ahead | `get_path_colour_prev/next` |
| `straightSecondColour` | Colour pair for straight reductions | Inline in `is_best_L0_reduction` / `is_best_straight_reduction` |
| `straightThirdColour` | Path colour for straight reductions | `is_best_third_colour` |
| `secondColour` | Dispatch: straight uses `straightSecondColour`, bent returns (0,0) | Various `is_best_*` |
| `thirdColour` | Dispatch: straight uses `straightThirdColour`, bent returns 0 | Various |

---

## 3. Colour Functions — Detailed Edge Traversal Mapping

### get_colour_prev_5 / get_colour_next_5 (C lines 5170-5196)

These compute a 5-bit integer encoding the degree-5/degree-6 pattern of 5
consecutive neighbors reachable by walking in one direction.

```c
// get_colour_prev_5 (CW walk):
EDGE *e = edge;
for(i = 4; i >= 0; i--) {
    colour |= ((degree[e->end] == 5) << i);
    e = e->prev;   // CW step
}
// Haskell: colourCW5 g w start
// Walk 5 CW steps from 'start' around vertex w, recording degree-5 bits.

// get_colour_next_5 (CCW walk):
EDGE *e = edge;
for(i = 4; i >= 0; i--) {
    colour |= ((degree[e->end] == 5) << i);
    e = e->next;   // CCW step
}
// Haskell: colourCCW5 g w start
// Walk 5 CCW steps from 'start' around vertex w, recording degree-5 bits.
```

The starting edge is critical — it determines WHICH 5 neighbors are examined.

### Colour starting edges for STRAIGHT reductions (L0 and L_i)

Given directed edge e from u to v (the reduction starting edge):

**use_next (DLeft) — C code: `e->next->next->invers->next`**

Step-by-step C→Haskell translation:
1. `e->next` = one CCW step around u from v = edge to `prevCW g u v`
2. `e->next->next` = another CCW = edge to `advanceCW g u v (deg g u - 2)`
3. Let w = that vertex. `->invers` = flip to directed edge from w to u
4. `->next` = one CCW step around w from u = edge to `prevCW g w u`
5. Then `get_colour_next_5` walks 5 more CCW steps from there

**Haskell**:
```haskell
-- DLeft (use_next):
w = advanceCW g u v (deg g u - 2)     -- 2 CCW from v around u
startNbr = prevCW g w u                -- 1 CCW from u around w
colourCCW5 g w startNbr                -- walk 5 CCW steps
```

**use_prev (DRight) — C code: `e->prev->prev->invers->prev`**

1. `e->prev` = one CW step around u from v = edge to `nextCW g u v`
2. `e->prev->prev` = another CW = edge to `advanceCW g u v 2`
3. Let w = that vertex. `->invers` = flip to edge from w to u
4. `->prev` = one CW step around w from u = edge to `nextCW g w u`
5. Then `get_colour_prev_5` walks 5 more CW steps

**Haskell**:
```haskell
-- DRight (use_prev):
w = advanceCW g u v 2                  -- 2 CW from v around u
startNbr = nextCW g w u                -- 1 CW from u around w
colourCW5 g w startNbr                 -- walk 5 CW steps
```

**For the OTHER endpoint** (e->invers = edge from v to u):
Same computation but starting from v, with u as the "opposite" vertex:
```haskell
-- DLeft: w2 = advanceCW g v u (deg g v - 2), colourCCW5 g w2 (prevCW g w2 v)
-- DRight: w2 = advanceCW g v u 2, colourCW5 g w2 (nextCW g w2 v)
```

**Colour pair**: `(negate (max c1 c2), negate (min c1 c2))` where c1, c2 are
the colours at the two endpoints. Higher colour = more canonical in C code
(C uses MAX for best), so we negate for Haskell's min-is-canonical ordering.

### Colour starting edges for B00 reductions

**DIFFERENT from straight types!** (C comment at line 10846: "uses different
colour1 than is_best_bent_zero_reduction")

Given test_edge1 (from first endpoint) and test_edge2 (from second endpoint):

**use_next (DLeft):**
```c
colour_testedge1 = get_colour_next_5(test_edge1->invers->next);
colour_testedge2 = get_colour_prev_5(test_edge2->invers->prev);
```

Note the asymmetry: edge1 uses `get_colour_next_5` (CCW), edge2 uses
`get_colour_prev_5` (CW). This is because the bent path reverses direction
at the bend point.

`test_edge1->invers->next`:
1. `->invers`: flip edge1 (now pointing toward first endpoint)
2. `->next`: one CCW step
3. Then `get_colour_next_5` walks 5 CCW steps

`test_edge2->invers->prev`:
1. `->invers`: flip edge2
2. `->prev`: one CW step
3. Then `get_colour_prev_5` walks 5 CW steps

**use_prev (DRight):** — note the swap:
```c
colour_testedge2 = get_colour_next_5(test_edge2->invers->next);  // note: edge2 uses next
colour_testedge1 = get_colour_prev_5(test_edge1->invers->prev);  // edge1 uses prev
```

### Path colour starting edges for STRAIGHT reductions

From `is_best_third_colour` (10137-10201):

**use_next (DLeft):**
```c
get_path_colour_next(test_edge1->prev->prev);
```
- `test_edge1->prev->prev`: 2 CW steps around u from v = edge to `advanceCW g u v 2`
- Then `get_path_colour_next` does:
  ```c
  e = e->invers->next->next->next;  // first straight-ahead step (CCW direction)
  for(i = 6; i >= 0; i--) {
      colour |= ((degree[e->end] == 5) << i);
      e = e->invers->next->next->next;  // next straight-ahead
  }
  ```

**Haskell**:
```haskell
-- DLeft (use_next):
w = advanceCW g u v 2
pathColour g DLeft u w    -- straightAhead in DLeft (CCW) direction for 7 steps
```

Wait — this needs careful attention. C code uses `test_edge->prev->prev` with
`get_path_colour_next`. The `->prev->prev` goes 2 CW steps, but
`get_path_colour_next` then walks straight-ahead in the CCW direction.

**use_prev (DRight):**
```c
get_path_colour_prev(test_edge1->next->next);
```
- `->next->next`: 2 CCW steps around u from v = edge to `advanceCW g u v (deg-2)`
- Then `get_path_colour_prev` walks straight-ahead in CW direction

**Haskell**:
```haskell
-- DRight (use_prev):
w = advanceCW g u v (deg g u - 2)
pathColour g DRight u w    -- straightAhead in DRight (CW) direction for 7 steps
```

### Path colour starting edges for B00 reductions

From `is_best_third_colour_bent` (10214-10278):

**use_next (DLeft):**
```c
colour_testedge1 = get_path_colour_next(test_edge1->invers->next->next);
colour_testedge2 = get_path_colour_prev(test_edge2->invers->prev->prev);
```

Again asymmetric: edge1 uses `_next` (CCW straight-ahead), edge2 uses `_prev`
(CW straight-ahead).

**use_prev (DRight):**
```c
colour_testedge1 = get_path_colour_prev(test_edge1->invers->prev->prev);
colour_testedge2 = get_path_colour_next(test_edge2->invers->next->next);
```

### get_path_colour_next / get_path_colour_prev (C lines 5246-5279)

```c
// get_path_colour_next (CCW straight-ahead):
e = e->invers->next->next->next;       // first straight-ahead step
for(i = PATHLENGTH_COLOUR - 1; i >= 0; i--) {   // PATHLENGTH_COLOUR = 7
    colour |= ((degree[e->end] == 5) << i);
    e = e->invers->next->next->next;   // next straight-ahead step
}

// get_path_colour_prev (CW straight-ahead):
e = e->invers->prev->prev->prev;       // first straight-ahead step
for(i = PATHLENGTH_COLOUR - 1; i >= 0; i--) {
    colour |= ((degree[e->end] == 5) << i);
    e = e->invers->prev->prev->prev;   // next straight-ahead step
}
```

Straight-ahead step at a degree-6 vertex: from the entry edge, the
straight-ahead neighbor is exactly opposite (3 CW or 3 CCW steps from entry).
`->invers->next->next->next` = flip to entry vertex, advance 3 CCW = opposite.

Haskell `straightAhead g dir atVtx fromVtx` computes this:
- DRight (CW): `advanceCW g atVtx fromVtx (halfDeg)` where halfDeg = deg/2
- DLeft (CCW): `advanceCW g atVtx fromVtx (deg - halfDeg)` where halfDeg = deg/2

---

## 4. Bounding Lemmas

### Paper: Section 2.3 "Optimizations" (lines 512-698)

| Lemma | Paper lines | Statement | Bound on canonical expansion length |
|-------|------------|-----------|--------------------------------------|
| **Lemma 1** | 519-525 | Reducible non-IPR duals have L0, L1, or B00 reduction. | length ≤ 2 for non-IPR |
| **Geometric** | 528-542 | Length-d expansion not canonical if child has < 12·f(⌊(d-1)/2⌋) verts, f(x)=1+5x(x+1)/2. | Depends on target size |
| **Lemma 2** | 547-567 | G has reduction of length d≤2 ⟹ all children have reduction of length ≤ d+2. | If parent has L1/B00: child ≤ 4 |
| **Lemma 3** | 576-614 | G has L0 reduction ⟹ all canonical children have reduction of length ≤ 2. Proof uses 5-edge-connectivity. | If parent has L0: child ≤ 2 (i.e., max expansion = 3) |
| **Obs. 3.1** | 636-643 | ≥3 five-vertices in expansion patch ⟹ child has 5-verts at distance ≤ l/2+1. | Used in Lemma 4,5 proofs |
| **Lemma 4** | 646-660 | ≥2 length-2 reductions with different 5-vertex sets ⟹ canonical children have length ≤ 3. | max expansion = 4 |
| **Lemma 5** | 663-674 | ≥3 length-2 reductions with pairwise disjoint 5-vertex sets ⟹ canonical children have length ≤ 2. | max expansion = 3 |
| **Lemma 6** | 681-693 | Two L0 reductions R1, R2 with d(R1,R2)>4 ⟹ all canonical children have L0 reduction. | max expansion = 2 |
| Effectiveness | 695-698 | Lemmas 2-6 determine bound in 93.9% of cases at 152 dual verts. | |

### C code: determine_max_pathlength_straight (4868-4966)

Detailed line-by-line:
```
4869: max_pathlength = maxnv - nv              // can't exceed target size
4870-4871: if ≤ 3, return immediately
4873-4882: Refine via max_straight_lengths[] table (precomputed geometric bounds)
4885-4891: Adjust for parity (nv + max_pathlength != maxnv - 1 for non-start)
4892-4908: Non-IPR optimizations:
  4901: if straight_length==2 (= L0 exists) && num_straight_extensions>0
        → max_pathlength = 3  [Lemma 3]
  4902: if straight_length==3 && contains_three_indep_L1s()
        → max_pathlength = 3  [Lemma 5]
  4903: if contains_three_indep_B00s() → max_pathlength = 3
4910-4938: Further refinements (only when max_pathlength > 4):
  4919: if straight_length==3 && num_straight_extensions>2
        → max_pathlength = 4  [Lemma 2 with two L1s]
  4920: if num_bent_zero_extensions>1 → max_pathlength = 4  [two B00s]
  4932-4937: if straight_length==3 && num_straight_extensions>0 && max>5
        → max_pathlength = 5  [one L2 exists]
```

Key variables:
- `straight_length`: length of shortest straight reduction found in the PARENT
  (set during the parent's is_best_* call, propagated to children)
  - 2 means L0 exists (distance 1 between 5-verts, pathlength = distance+1 = 2)
  - 3 means L1 exists (distance 2, pathlength 3)
- `num_straight_extensions`: count of straight reductions at that length
- `num_bent_zero_extensions`: count of B00 reductions

### Haskell: Not yet implemented

Currently TestCanonical.hs uses `maxLen = maxDV - nv + 2` (simple geometric bound).

---

## 5. Generation Loop

### C code: scansimple_fuller (12100-12722) — detailed structure

```
12103-12108: if nv == maxnv → output graph, return
12113-12124: Split-level handling (parallel computation via mod/res)
12129-12131: Intermediate output (startswitch mode)
12138: max_pathlength_straight = determine_max_pathlength_straight()
12144-12147: if max_pathlength < 2 → return (no valid expansions possible)
12149-12197: FIND_EXTENSIONS_SIMPLE_FULLER (= find_extensions_fuller)
  Populates: ext_L0[], ext_L1[], ext_bent_zero[], ext_straight[], ext_bent[]
  Each with use_next[] flag and (for longer types) length[]/position[]

12213-12270: L0 expansions loop:
  for each L0:
    compute path[] and parallel_path[]
    extend_L0(...)
    if is_best_L0_reduction(...):
      if only 1 surviving candidate: skip BFS, accept with trivial group
      else: canon_edge_oriented(...) → get automorphism group
      if accepted: scansimple_fuller(xnbtot, xnbop)  // recurse
    reduce_L0(...)

12273-12325: L1 expansions loop:
  for each L1:
    compute path[0..3] and parallel_path[0..3]
    extend_straight(..., pathlength=3, ...)
    if is_best_straight_reduction(edge_list[nv-3][nv-2], edge_list[nv-1][nv-2], 3, ...):
      [...same pattern: optional BFS, recurse...]
    reduce_straight(...)

12327-...: B00 expansions loop:
  for each B00:
    extend_bent_zero(...)
    if is_best_bent_zero_reduction(...):
      [...same pattern...]
    reduce_bent_zero(...)

...: Longer straight (L2+) expansions loop:
  Same pattern with extend_straight/reduce_straight and pathlength > 3

...: Longer bent (B_{i,j}, i+j>0) expansions loop:
  Same pattern with extend_bent/reduce_bent
```

### Haskell: TestCanonical.hs generation

```haskell
expandChildren maxDV st g =
    let nv = numVertices g
        maxLen = maxDV - nv + 2
        exps = expansions maxLen g           -- all expansions (no Rule 2)
        children = [ g'
                   | e <- exps
                   , let (g', _) = applyExpansion e g
                   , numVertices g' <= maxDV
                   , isCanonical e nv g'     -- Rule 1
                   ]
        ringChildren = case findNanotubeRing g of  -- F expansion
            Just (ring, outer) | nv + 5 <= maxDV -> [applyRing g ring outer]
            _ -> []
    in foldl' (processTree maxDV) st (children ++ ringChildren)

processTree maxDV st g
    | numVertices g > maxDV = st
    | Just gs <- generalSpiralKey g
    , Set.member gs (gsSeen st) = st           -- spiral dedup (Rule 2 substitute)
    | Just gs <- generalSpiralKey g =
        let st' = st { gsSeen = Set.insert gs (gsSeen st)
                     , gsGraphs = Map.insertWith (++) nv [g] (gsGraphs st) }
        in expandChildren maxDV st' g
```

---

## 6. Extension Finding (Rule 2) — Detailed

### C code: find_extensions_fuller (6784-7047)

**Trivial group** (numb_total == 1, line 6844-6963):

Uses multiple independent mark arrays to avoid duplicate expansions:
- `MARK_L0_NEXT/PREV`: marks L0 edges already tried in next/prev direction
- `MARK_DOUBLE_NEXT/PREV`: marks L1+ straight edges
- `MARK_B00_NEXT/PREV`: marks B00 edges
- `MARK_BENT_NEXT/PREV`: marks longer bent edges

For L0 specifically (line 6847-6857):
```c
for(i = 0; i < 12; i++) {
    e = firstedge[degree_5_vertices[i]];
    for(k = 0; k < 5; k++) {
        find_L0_extensions_prev(e, ...);  // try CW direction
        find_L0_extensions_next(e, ...);  // try CCW direction
        e = e->next;
    }
}
```

`find_L0_extensions_next` (5432-5507):
```c
e = startedge->invers->next->next->next;
temp_vertex = e->invers->next->next->end;  // the other 5-vertex (3 straight-ahead steps away)
if(degree[temp_vertex] == 5 && !ISMARKED_L0_NEXT(startedge->label)) {
    // Lookahead: check if this L0 expansion can be canonical
    // by comparing its colour against the parent's best known colours
    if(can_be_canonical) {
        edge_ext_straight[*num_ext_straight] = startedge;
        ext_use_next[*num_ext_straight] = 1;
        (*num_ext_straight)++;
    }
    // Mark the other triple (e->invers->next->next->invers->next) and equivalents
    mark_edges_L0(other_startedge, 1, numb_total, npres);
}
```

The L0 marking ensures that for an undirected edge {u,v} between two degree-5
vertices, we only try the expansion from one directed edge in each direction
(not both u→v and v→u in the same direction, since they give the same expansion).

**Non-trivial group** (numb_total > 1, line 6968-7046):
- `numbering[j][index]` maps each edge to its image under the j-th automorphism
- `or_same_edge_found`: if an orientation-reversing automorphism fixes the edge,
  only one direction needs to be tried
- All images of each edge are marked before processing, preventing equivalent
  edges from being tried

### Haskell: Expansion.hs expansion enumeration

```haskell
expansions :: Int -> DualGraph -> [Expansion]
expansions maxLen g = expansionsL0 g ++ expansionsL maxLen g
                   ++ expansionsB maxLen g
```

This tries ALL possible expansions without Rule 2 filtering. The spiral
dedup in the generation loop catches duplicates.

---

## 7. Reduction Validation Functions

### L0: has_L0_reduction (C lines 8990-9002)

```c
has_L0_reduction(EDGE *startedge, unsigned char *direction_bitvector) {
    if(degree[startedge->end] == 5) {   // other end is also deg-5
        *direction_bitvector = 0;
        // CW flanking: startedge->prev->prev->end (2 CW from v around u)
        //              startedge->invers->prev->prev->end (2 CW from u around v)
        if(degree[startedge->prev->prev->end] == 6 && degree[startedge->invers->prev->prev->end] == 6)
            (*direction_bitvector) |= PREV_BITVECTOR;   // DRight valid
        // CCW flanking:
        if(degree[startedge->next->next->end] == 6 && degree[startedge->invers->next->next->end] == 6)
            (*direction_bitvector) |= NEXT_BITVECTOR;   // DLeft valid
        if(*direction_bitvector > 0)
            return 1;
    }
    return 0;
}
```

Precondition: `startedge->start` has degree 5 (only called from degree-5 vertex loop).

Why flanking must be degree-6: The L0 reduction removes 2 vertices and reconnects
their neighbors. The flanking vertices gain one neighbor each in the process.
If they're already degree-5, they'd become degree-6 (fine). If degree-6, they'd
become degree-7 (not a valid fullerene dual — only degrees 5 and 6 allowed).
Wait, actually it's the reverse: the flanking vertices LOSE connectivity. Let me
re-examine... Actually, the condition ensures the expansion is reversible: after
expanding, the L0 reduction back must produce a valid graph. The degree-6
flanking condition ensures the patch replacement doesn't create degree violations.

### L0/L1: is_L0_or_L1_reduction (C lines 10442-10452)

```c
is_L0_or_L1_reduction(EDGE *startedge) {
    EDGE *e = startedge->invers;
    if(degree[e->start] == 6) {           // endpoint is deg-6 → not L0, check L1
        e = e->invers->next->next->next;  // straight-ahead step
        if(degree[e->start] == 6)         // next vertex also deg-6 → no L1
            return 0;
    }
    // Now e points to a deg-5 vertex (either directly adjacent or 1 step away)
    // Check direction validity (flanking vertices must be deg-6):
    return (degree[startedge->prev->prev->end] == 6 && degree[e->prev->prev->end] == 6)
        || (degree[startedge->next->next->end] == 6 && degree[e->next->next->end] == 6);
}
```

Used by `is_best_bent_zero_reduction` to quickly check if any L0 or L1 reduction
exists (which would be shorter than B00, making B00 non-canonical).

### B00: has_bent_zero_reduction (C lines 10468-10487)

```c
has_bent_zero_reduction(EDGE *startedge, EDGE **enarc) {
    EDGE *e = startedge;       // from deg-5 vertex u
    if(degree[e->end] == 5) return 0;    // adjacent must be deg-6 (not L0)
    e = e->next->next;                    // 2 CCW around u
    if(degree[e->end] == 5) return 0;    // must be deg-6
    e = e->invers->prev->prev;           // flip + 2 CW = the "turn" step
    if(degree[e->end] == 6) return 0;    // must be deg-5 (other endpoint of bent path)
    e = e->invers->next->next;           // flip + 2 CCW
    if(degree[e->end] == 5) return 0;    // must be deg-6
    *enarc = e;                           // edge at the other end
    return 1;
}
```

Path traversal: u (deg-5) → v1 (deg-6, 1 step) → v2 (deg-6, 2 CCW around u)
→ TURN → v3 (deg-5 = other endpoint). The "turn" is implemented as
`->invers->prev->prev` instead of the straight-ahead `->invers->next->next->next`.

Note: "Searching in one direction is sufficient, since enarc is same but in
other direction" (C comment at line 10464-10465). This means each B00 reduction
is found once from each 5-vertex endpoint, and the `enarc` returned is the
edge at the OTHER endpoint. Only prev direction is searched (CW turn).

### Straight (L_i, i≥1): has_short_straight_reduction (C lines 9017-9074)

```c
has_short_straight_reduction(startedge, pathlength_cur_best, &enarc, &dir_bitvector) {
    e = startedge;
    length = 1;
    while(length < pathlength_cur_best - 1 && degree[e->end] != 5) {
        e = e->invers->next->next->next;  // straight-ahead (CCW)
        // Mark vertices to detect self-intersection
        if(ISMARKED_V(e->start)) return 0;
        MARK_V(e->start); MARK_V(e->prev->end); MARK_V(e->next->end);
        length++;
    }
    if(degree[e->end] == 5) {
        e = e->invers;
        // Check for self-intersection
        if(ISMARKED_V(e->start) || ISMARKED_V(e->prev->prev->end)
           || ISMARKED_V(e->next->next->end)) return 0;
        if(length < pathlength_cur_best - 1)
            return 2;  // SHORTER reduction found!
        // Same length: check direction validity
        *direction_bitvector = 0;
        if(flanking CW both deg-6) (*direction_bitvector) |= PREV_BITVECTOR;
        if(flanking CCW both deg-6) (*direction_bitvector) |= NEXT_BITVECTOR;
        if(*direction_bitvector > 0) { *enarc = e; return 1; }
        return 0;
    }
    return 0;
}
```

Returns 2 if a strictly shorter reduction is found (immediate rejection),
1 if same-length reduction found (continue with colour comparison), 0 if no
valid reduction from this edge.

### Haskell: allReductions (Canonical.hs)

```haskell
allReductions g = l0Reds ++ lReds ++ b00Reds ++ bigBReds
  where
    maxLen = numVertices g
    l0Reds = [ Red (L 0) (u, v) d
             | u <- degree5 g
             , v <- filter (\w -> deg g w == 5) (nbrs g u)
             , d <- [DLeft, DRight]
             , isValidL0Direction g (u, v) d ]
    lReds = map inverse (expansionsL maxLen g)      -- reuse expansion finder
    b00Reds = [ Red (B 0 0) (u, v) d
              | u <- degree5 g, v <- nbrs g u
              , d <- [DLeft, DRight]
              , isValidB00Site g (u, v) d ]
    bigBReds = map inverse
        [ e | e@(Exp (B i j) _ _) <- expansionsB maxLen g, i + j > 0 ]
```

---

## 8. IPR Generation

### Paper: Section 3 (lines 704+)

L0 expansion always creates adjacent 5-vertices → always non-IPR.
For IPR-only generation: skip L0, use modified colour functions and
bounding lemmas.

### C code: scansimple_fuller_ipr

Separate function with:
- `is_best_straight_reduction_ipr` (11939+): uses `get_colour_next_3` / `get_colour_prev_3`
  instead of the 5-variant (only 3-bit colour since IPR guarantees no adjacent
  5-vertices, so fewer informative bits needed)
- `is_best_bent_zero_reduction_ipr` (10848+): additional IPR validity checks
  (`has_bent_zero_reduction_ipr` checks that reduced graph is also IPR)
- Additional IPR-specific bounding
- 49 IPR seed graphs for large sizes

---

## 9. Key Data Structures

### C: EDGE record (doubly-linked, oriented)
```c
typedef struct e {
    int start;            // source vertex
    int end;              // target vertex
    struct e *prev;       // CW neighbor (same start vertex)
    struct e *next;       // CCW neighbor (same start vertex)
    struct e *invers;     // reverse edge (swap start/end)
    int label;            // unique edge label for marking
    int index;            // index in numbering[][] array
} EDGE;
```

### C: Global state variables

| Variable | Type | Purpose |
|----------|------|---------|
| `nv` | int | Current vertex count |
| `ne` | int | Current directed edge count (= 2 × undirected edges) |
| `maxnv` | int | Target vertex count |
| `degree[]` | int[MAXN] | Vertex degrees (5 or 6) |
| `firstedge[]` | EDGE*[MAXN] | First edge record per vertex |
| `degree_5_vertices[]` | int[12] | The 12 degree-5 vertices (always exactly 12) |
| `numbering[][]` | EDGE*[2·MAXE][MAXE] | Automorphism group as edge permutations |
| `straight_length` | int | Pathlength of shortest straight reduction (2=L0, 3=L1, etc.) |
| `straight_extensions[][]` | int[][] | Cached straight reduction vertex lists |
| `num_straight_extensions` | int | Count of cached straight reductions |
| `bent_zero_extensions[]` | int[][] | Cached B00 reduction vertex lists |
| `num_bent_zero_extensions` | int | Count of B00 reductions |
| `best_straight_colour_1` | int | Best colour1 value from parent (for lookahead) |
| `best_straight_colour_2` | int | Best colour2 value from parent |
| `last_edge_L0` | EDGE* | Edge between the two new vertices after L0 expansion |

### Haskell: DualGraph (Seeds.hs)
```haskell
data DualGraph = DG
    { numVertices :: !Int
    , neighbours  :: !(IntMap [Int])   -- CW-ordered neighbor lists
    , degree5     :: [Int]             -- 12 degree-5 vertices
    , edgeList    :: !EdgeList          -- O(1) edge lookup
    }
```

### Haskell: Key types (Expansion.hs)
```haskell
type Vertex = Int
type Edge   = (Vertex, Vertex)   -- directed: (start, end)
data Dir    = DLeft | DRight     -- DLeft = CCW (use_next), DRight = CW (use_prev)
data ExpKind = L Int             -- L_i: straight, parameter i
             | B Int Int         -- B_{i,j}: bent, parameters i,j
             | F                 -- nanotube ring
data Expansion = Exp ExpKind Edge Dir
data Reduction = Red ExpKind Edge Dir
data PathInfo = PI { mainPath :: [Vertex], parallelPath :: [Vertex] }
```

---

## 10. Lookahead Optimization (C-specific)

The C code caches reduction information from the PARENT graph to prune
children's expansion search early. This is not in the paper.

### Colour lookahead in find_L0_extensions_next (5432-5507)

When the parent has L0 reductions (`straight_length == 2` and
`best_straight_colour_1 > 0`), the extension finder pre-computes the colour
of each candidate L0 expansion and compares it against the parent's best
colours. If the candidate can't possibly beat the parent's canonical
reduction's colour, it's skipped.

```c
if(straight_length == 2 && best_straight_colour_1 > 0) {
    // Compute colour of this expansion's inverse:
    colour_one = get_colour_next_5(startedge);
    colour_one_enarc = get_colour_next_5(other_startedge);
    best_colour = MAX(colour_one, colour_one_enarc);
    other_colour = MIN(colour_one, colour_one_enarc);

    if(best_colour < best_straight_colour_1 ||
       (best_colour == best_straight_colour_1 && other_colour < best_straight_colour_2)) {
        // This expansion CAN'T produce a canonical child unless it destroys
        // the parent's shorter/better reductions. Check:
        can_be_canonical = might_destroy_previous_L0_extensions();
    }
}
```

### might_destroy_previous_straight_extension (5291-5313)

Checks if the current expansion's vertices overlap with previously cached
short reductions. If a short reduction is NOT destroyed (no overlap with
expansion patch), then the child still has that short reduction, and a
longer expansion can't be canonical.

This is a key optimization for large graphs: it avoids applying and testing
expansions that can't possibly be canonical.

---

## 11. File Index

### C source
| File | Lines | Description |
|------|-------|-------------|
| `buckygen.c` | 1-14936 | Complete C implementation |

### C code function index (sorted by line number)

| Line | Function | Category |
|------|----------|----------|
| 647 | `code_c20` | Seed data |
| 651 | `code_c28` | Seed data |
| 656 | `code_c30` | Seed data |
| 1184 | `extend_L0` | Expansion operation |
| 1338 | `reduce_L0` | Reduction operation |
| 1395 | `extend_straight` | Expansion operation |
| 1541 | `reduce_straight` | Reduction operation |
| 1602 | `extend_bent_zero` | Expansion operation |
| 1826 | `reduce_bent_zero` | Reduction operation |
| 1886 | `extend_bent` | Expansion operation |
| 2160 | `reduce_bent` | Reduction operation |
| 2246 | `testcanon` | BFS canonical (CCW/->next) |
| 2356 | `testcanon_short` | BFS canonical (short scan) |
| 2459 | `testcanon_mirror` | BFS canonical (CW/->prev) |
| 3267 | `canon` | Full canonicity + automorphism group |
| 4088 | `canon_edge_oriented` | BFS from specific edges |
| 4211 | `canon_edge_oriented_short` | BFS short scan from edges |
| 4868 | `determine_max_pathlength_straight` | Bounding lemmas |
| 4976 | `mark_edges_L0` | Rule 2: mark equivalent L0 edges |
| 5170 | `get_colour_prev_5` | Colour: 5-bit CW walk |
| 5187 | `get_colour_next_5` | Colour: 5-bit CCW walk |
| 5203 | `get_colour_prev_3` | Colour: 3-bit CW walk (IPR) |
| 5221 | `get_colour_next_3` | Colour: 3-bit CCW walk (IPR) |
| 5246 | `get_path_colour_next` | Path colour: 7-bit CCW straight-ahead |
| 5269 | `get_path_colour_prev` | Path colour: 7-bit CW straight-ahead |
| 5291 | `might_destroy_previous_straight_extension` | Lookahead |
| 5325 | `might_destroy_previous_L0_extensions` | Lookahead |
| 5432 | `find_L0_extensions_next` | Rule 2: find L0 expansions (CCW) |
| 5517 | `find_L0_extensions_prev` | Rule 2: find L0 expansions (CW) |
| 5602 | `find_L1_extensions_next` | Rule 2: find L1 expansions |
| 6784 | `find_extensions_fuller` | Rule 2: find all expansions |
| 8990 | `has_L0_reduction` | Validate L0 reduction |
| 9017 | `has_short_straight_reduction` | Find same/shorter straight reductions |
| 9094 | `has_short_straight_reduction_L1` | Optimized for L1 case |
| 10137 | `is_best_third_colour` | Path colour cascade (straight types) |
| 10214 | `is_best_third_colour_bent` | Path colour cascade (bent types) |
| 10292 | `is_best_L0_reduction` | Full canonicity test for L0 |
| 10442 | `is_L0_or_L1_reduction` | Quick check: does L0/L1 exist? |
| 10468 | `has_bent_zero_reduction` | Validate B00 reduction |
| 10506 | `has_bent_zero_reduction_ipr` | B00 validation (IPR) |
| 10561 | `has_B10_reduction` | Check for B10 reduction |
| 10716 | `is_best_bent_zero_reduction` | Full canonicity test for B00 |
| 10848 | `is_best_bent_zero_reduction_ipr` | B00 canonicity (IPR) |
| 11659 | `has_B00_reductions` | Quick check: does B00 exist? |
| 11812 | `is_best_straight_reduction` | Full canonicity test for L_i (i≥1) |
| 11939 | `is_best_straight_reduction_ipr` | L_i canonicity (IPR) |
| 12100 | `scansimple_fuller` | Main generation loop |

### Paper source
| File | Description |
|------|-------------|
| `buckygen_paper_source/fullerene-generation.tex` | Main paper |
| `buckygen_paper_source/figures/` | Figures |

### Paper section index

| Section | Lines | Topic |
|---------|-------|-------|
| 2.1 | 166-273 | Construction algorithm (L, B, F expansions; Theorem 1.1) |
| 2.2 | 275-423 | Isomorphism rejection (Rules 1-2, 6-tuple, BFS) |
| 2.3 | 512-698 | Bounding lemmas (Lemmas 1-6, geometric bound) |
| 3 | 704+ | IPR generation |

### Haskell source
| File | Description | Status |
|------|-------------|--------|
| `Seeds.hs` | DualGraph type, 49 seed graphs | Complete |
| `Expansion.hs` | Types, navigation, all expansion/reduction ops | Complete |
| `Spiral.hs` | Canonical spiral computation (dedup) | Complete |
| `Canonical.hs` | Canonicity test (cascade + BFS) | In development |
| `TestCanonical.hs` | Generation with canonical test + spiral dedup | In development |
| `Generate.hs` | Generation with spiral dedup only (reference) | Complete, verified C160 |
| `Canonical-wrong.hs` | Old broken implementation (kept for reference) | Archived |

---

## 12. Known Issues and Debugging Notes

### Current Canonical.hs test failures

At C40: C26=0/1, C30=1/3, C32=4/6, C34=1/6, C36=13/15, C38=5/17, C40=33/40.
Pattern: systematic under-counting = false negatives in isCanonical.

### Potential bug sources (ordered by likelihood)

1. **allReductions may miss valid reduction sites.**
   The function reuses `expansionsL` and `expansionsB` (via `inverse`) for L_i≥1
   and B_{i,j}. If these enumeration functions don't find all sites (e.g., due
   to self-intersection checks or path validation), the canonical reduction might
   be missing from the list, causing `minimum allOrds` to be wrong.

2. **Inverse identification may be too narrow.**
   `isCanonical` identifies the inverse by: kind matches AND start vertex is one
   of the new vertices. But the C code identifies the inverse via the actual
   edge stored during extension (`last_edge_L0`, `edge_list[nv-3][nv-2]`, etc.).
   The Haskell approach may miss cases where the inverse reduction starts from
   the new vertex but at a different neighbor.

3. **secondColour for bent types returns (0,0).**
   The C code computes actual colours for B00 reductions (using the
   asymmetric `get_colour_next_5`/`get_colour_prev_5` pattern). The Haskell
   `secondColour` returns `(0,0)` for all non-straight types, meaning all B
   reductions tie on x2 and go straight to x3 (which also returns 0 for bent).
   This means the only discriminator for bent reductions is BFS, which is
   expensive and may give different ordering than the proper colour cascade.

4. **Direction mapping errors in colour functions.**
   The current `straightSecondColour` in Canonical.hs uses:
   - DRight → `advanceCW 2` + `colourCW5` + `nextCW`
   - DLeft → `advanceCW (deg-2)` + `colourCCW5` + `prevCW`

   This matches the C code (verified above). But verify the path colour
   (`straightThirdColour`) starting edges:
   - DRight → should use `advanceCW (deg-2)` + `pathColour DRight`
   - DLeft → should use `advanceCW 2` + `pathColour DLeft`

   **Check**: the DRight path colour starts from `test_edge->next->next` in C
   (= 2 CCW = deg-2 CW), with `get_path_colour_prev` (CW straight-ahead).
   The DLeft path colour starts from `test_edge->prev->prev` (= 2 CW),
   with `get_path_colour_next` (CCW straight-ahead).

   **Key insight**: The COLOUR starting point and the PATH COLOUR starting
   point are on OPPOSITE sides of the edge! Colour goes 2 steps in the walk
   direction; path colour goes 2 steps AGAINST the walk direction.

### Verification strategy

1. For small cases (C24, C26), print all reductions found by `allReductions`
   and compare against the C code output.
2. For a specific missing fullerene (e.g., the C26 isomer), trace the expansion
   that should produce it, check why `isCanonical` rejects it.
3. Verify colour functions produce same values as C code on known graphs.
