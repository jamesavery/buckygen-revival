# The Buckygen Algorithm: Complete Specification

**Purpose**: Generate exactly one representative of every fullerene isomer with
a given number of vertices, using McKay's canonical construction path method.

**Scope**: This document is a self-contained algorithmic specification. A reader
should be able to implement the complete generator from this file alone, without
consulting the original paper or C source code.

**Notation**: We work entirely in the **dual** representation. A fullerene with
*n* carbon atoms has *n*/2 + 2 dual vertices. The dual is a triangulation of the
sphere where every vertex has degree 5 or 6, with exactly 12 degree-5 vertices.

---

## 1. Graph Representation

A **dual graph** *G* is stored as:

- An integer *nv* (number of vertices, labelled 0 through *nv* - 1)
- For each vertex *v*, a **cyclic neighbour list** `nbrs(v)` storing its
  neighbours in **clockwise planar order**
- A list `deg5(G)` of the 12 degree-5 vertices (always exactly 12)

Every vertex has degree 5 or 6. The cyclic ordering of neighbours encodes the
planar embedding and is essential for all operations.

### 1.1 Navigation primitives

Given the CW-ordered neighbour list of vertex *u*:

| Primitive | Definition |
|-----------|------------|
| `nextCW(u, v)` | The neighbour of *u* one position **clockwise** after *v* in `nbrs(u)` |
| `prevCW(u, v)` | The neighbour of *u* one position **counter-clockwise** after *v* (= one position before *v* in the CW list) |
| `advanceCW(u, v, k)` | The neighbour of *u* that is *k* positions clockwise from *v* in `nbrs(u)` |
| `straightAhead(dir, u, from)` | The vertex "straight ahead" through *u* when entering from `from`. At a degree-*d* vertex, this is the neighbour *d*/2 positions away in the walk direction. For degree 6 entering from `from`: `advanceCW(u, from, 3)` for either direction. |
| `turnAhead(dir, u, from)` | The vertex reached by a "turn" (2 positions instead of 3): `advanceCW(u, from, 2)` for DRight, `advanceCW(u, from, deg(u)-2)` for DLeft. |

### 1.2 Directions

A **direction** is either **DRight** (clockwise walk) or **DLeft** (counter-clockwise walk).

- **DRight**: Walking around each vertex follows the CW ordering. Navigation steps use `nextCW`. The "side" neighbour (parallel to a path) is obtained by `prevCW`.
- **DLeft**: Walking follows the CCW ordering. Navigation steps use `prevCW`. The "side" neighbour is `nextCW`.

The function `flipDir` maps DRight to DLeft and vice versa.

---

## 2. Seed Graphs

The algorithm starts from three irreducible fullerenes:

| Seed | Dual vertices | Carbon atoms | Symmetry | Role |
|------|---------------|-------------|----------|------|
| C20 (dodecahedron) | 12 | 20 | I_h | Root for most fullerenes |
| C28 | 16 | 28 | T_d | Root for a small family unreachable from C20 |
| C30 | 17 | 30 | D_5h | Root for (5,0) nanotube family (via F expansion) |

**Theorem** (Hasheminezhad et al.): Every fullerene except C28 and the (5,0)
nanotubes can be constructed from C20 by repeatedly applying L and B expansions.
The (5,0) nanotubes are constructed from C30 by repeatedly applying F.

The seed graphs are encoded in **planar code format**: a byte giving *nv*,
followed by the 1-indexed CW neighbour list of each vertex terminated by 0.

---

## 3. Expansion Operations

An **expansion** transforms a dual graph *G* with *nv* vertices into a larger
dual graph *G'* by replacing a small local patch with a larger one. An expansion
is fully described by a triple (*kind*, *edge*, *dir*) where:

- *kind* specifies the expansion type (L_0, L_1, ..., B_{0,0}, B_{1,0}, ..., or F)
- *edge* = (*u*, *v*) is a directed starting edge at a degree-5 vertex *u*
- *dir* is DLeft or DRight

Each expansion has a unique inverse **reduction** that recovers *G* from *G'*.

### 3.1 L_i expansions (straight)

An L_i expansion operates on a "straight path" between two degree-5 vertices.
The path has length *i* + 1 (i.e., *i* + 2 vertices including both endpoints).

**Parameters**:
- *i* >= 0 (L_0 is the shortest; it operates on adjacent degree-5 vertices)
- New vertices added: *i* + 2
- Reduction length: *i* + 1 (= distance between the two 5-vertices)

**Path computation** (`computeStraightPath`):
Starting from directed edge (*u*, *v*) with direction *dir*, compute *numEntries*
= *i* + 3 vertices by following straight-ahead steps:

```
path[0] = u
For k = 1, 2, ..., numEntries-1:
    path[k] = straightAhead(dir, path[k-1], path[k-2])
    // (for k=1, path[k-2] is notional; the step goes from u toward v)

parallelPath[k] = sideNbr(dir, path[k], path[k+1])   for each k
```

where `sideNbr(DRight, a, b) = prevCW(a, b)` and `sideNbr(DLeft, a, b) = nextCW(a, b)`.

**Validity conditions for an L_i expansion site**:
1. `path[0]` (= *u*) has degree 5
2. `parallelPath[last]` has degree 5 (it will gain an edge)
3. The path and parallel path are disjoint (no self-intersection)
4. Each `path[k]` is adjacent to `parallelPath[k-1]` and `parallelPath[k]`
5. Each `parallelPath[k]` is adjacent to `path[k+1]`

**The surgery**: The expansion adds *i* + 2 new vertices (numbered *nv* through
*nv* + *i* + 1) that form a new strip between the main path and the parallel
path. The first and last new vertices become degree-5; the old degree-5 vertices
(`path[0]` and `parallelPath[last]`) become degree-6.

For **L_0** specifically: 2 new vertices are added. The new vertices *nv* and
*nv* + 1 are both degree-5 and adjacent to each other. The old 5-vertices
(`path[0]` and `parallelPath[2]`) become degree-6.

**Exact surgery**: See Section 3.4.1 (L_0) and Section 3.4.2 (L_i, i >= 1).

### 3.2 B_{i,j} expansions (bent)

A B_{i,j} expansion operates on a path between two degree-5 vertices that
includes a **turn** (bend) at an intermediate vertex, instead of going straight.

**Parameters**:
- *i* >= 0, *j* >= 0, with *i* + *j* >= 0 (B_{0,0} is the shortest)
- New vertices added: *i* + *j* + 3
- Reduction length: *i* + *j* + 2

**Path computation** (`computeBentPath` / `computeBentZeroPath`):
Starting from (*u*, *v*) with direction *dir*:

1. **Pre-bend phase**: Follow *i* + 2 straight-ahead steps (same as L_i).
2. **Turn**: At the vertex reached, instead of `straightAhead` (3 positions for
   degree-6), use `turnAhead` (2 positions). This creates the bend.
3. **Post-bend phase**: Follow *j* + 1 more straight-ahead steps.
4. **Final endpoint**: One more straight-ahead step reaches the final vertex.

The main path has *i* + *j* + 5 entries. The parallel path has *i* + *j* + 3
entries (the turn vertex is skipped for the parallel path).

**Validity conditions for a B_{i,j} expansion site**:
1. `path[0]` (= *u*) has degree 5
2. `path[last]` has degree 5
3. The path and parallel path are disjoint
4. Various adjacency conditions hold (each path vertex adjacent to its
   parallel path neighbours, and the bend area adjacencies are correct)

**Exact surgery**: See Section 3.4.3 (B_{0,0}) and Section 3.4.4 (B_{i,j}, i+j > 0).

### 3.3 F expansion (nanotube ring)

The F expansion applies only to (5,0) nanotubes (the family rooted at C30).
It adds 5 new degree-6 vertices around the equatorial ring.

**Detection**: Find a 5-cycle of degree-6 vertices where each consecutive pair
shares a common neighbour not in the ring. The ring is found by following
straight-ahead steps: start from a degree-6 vertex edge, follow straight-ahead
for 5 steps; if the path closes, it's the equatorial ring.

**Surgery**: Insert 5 new degree-6 vertices between the equatorial ring and one
cap, reconnecting edges appropriately. No degree-5 vertices are created or
destroyed, so F-expanded graphs remain in the nanotube family.

F expansions are always canonical (there is only one way to expand a nanotube).
No canonical test is needed.

**Exact surgery**: See Section 3.4.5.

### 3.4 Expansion Surgery Specifications

This section gives the exact neighbour-list modifications for each expansion type.
These are the operations an implementor must perform to produce the child graph
from the parent graph.

**Notation**:
- `p[k]` = `path[k]`, `q[k]` = `parallelPath[k]`
- New vertices are numbered `nv, nv+1, ..., nv+numNew-1` where `nv` is the
  parent's vertex count
- `insertAfter(u, ref, new)`: in `u`'s CW neighbour list, insert `new`
  immediately CW after `ref`, increasing `u`'s degree by 1
- `replaceNbr(u, old, new)`: in `u`'s CW neighbour list, replace `old` with
  `new`, degree unchanged

#### 3.4.1 L_0 Surgery

**Parameters**: path = [p0, p1, p2], par = [q0, q1, q2]. New vertices:
a = nv (deg-5), b = nv+1 (deg-5).

**New vertex CW neighbour lists**:
```
DRight: a = [p0, q0, q1, b, p1]    b = [a, q1, q2, p2, p1]
DLeft:  a = [p0, p1, b, q1, q0]    b = [a, p1, p2, q2, q1]
```

**Existing vertex modifications**:
1. p0 (becomes deg-6): `insertAfter(p0, q0, a)` [DRight] / `insertAfter(p0, p1, a)` [DLeft]
2. q2 (becomes deg-6): `insertAfter(q2, p2, b)` [DRight] / `insertAfter(q2, prevCW(q2, p2), b)` [DLeft]
3. p1: `replaceNbr(p1, q0, a)`, `replaceNbr(p1, q1, b)`
4. p2: `replaceNbr(p2, q1, b)`
5. q0: `replaceNbr(q0, p1, a)`
6. q1: `replaceNbr(q1, p1, a)`, `replaceNbr(q1, p2, b)`
7. Degree-5 list: remove p0, q2; add a, b

#### 3.4.2 L_i (i >= 1) Surgery

**Parameters**: path has i+3 entries, par has i+3 entries. `pathlength` = i+2
new vertices (nv+0 through nv+pathlength-1). First and last new vertices are
deg-5; middle ones are deg-6.

**New vertex CW neighbour lists** (for 0 <= k <= pathlength-1):

k == 0 (deg-5):
```
DRight: [p[0], q[0], q[1], nv+1, p[1]]
DLeft:  [p[0], p[1], nv+1, q[1], q[0]]
```

0 < k < pathlength-1 (deg-6):
```
DRight: [nv+k-1, q[k], q[k+1], nv+k+1, p[k+1], p[k]]
DLeft:  [nv+k-1, p[k], p[k+1], nv+k+1, q[k+1], q[k]]
```

k == pathlength-1 (deg-5):
```
DRight: [nv+k-1, q[k], q[k+1], p[k+1], p[k]]
DLeft:  [nv+k-1, p[k], p[k+1], q[k+1], q[k]]
```

**Existing vertex modifications**:
1. p[0] (becomes deg-6): `insertAfter(p[0], q[0], nv)` [DRight] / `insertAfter(p[0], p[1], nv)` [DLeft]
2. q[pathlength] (becomes deg-6): `insertAfter(q[pl], p[pl], nv+pl-1)` [DRight] / `insertAfter(q[pl], prevCW(q[pl], p[pl]), nv+pl-1)` [DLeft]
3. Path vertices p[k] for 1 <= k <= pathlength: `replaceNbr(p[k], q[k-1], nv+k-1)`; if k < pathlength: `replaceNbr(p[k], q[k], nv+k)`
4. Parallel vertices q[k] for 0 <= k <= pathlength-1: if k > 0: `replaceNbr(q[k], p[k], nv+k-1)`; `replaceNbr(q[k], p[k+1], nv+k)`
5. Degree-5 list: remove p[0], q[pathlength]; add nv+0, nv+pathlength-1

#### 3.4.3 B_{0,0} Surgery

**Parameters**: path = [p0, p1, p2, p3, p4] (5 entries), par = [q0, q1, q2]
(3 entries, turn vertex skipped). New vertices: a = nv (deg-5),
b = nv+1 (deg-6, bend), c = nv+2 (deg-5).

**New vertex CW neighbour lists**:
```
DRight: a = [p0, q0, q1, b, p1]    b = [a, q1, c, p3, p2, p1]    c = [b, q1, q2, p4, p3]
DLeft:  a = [p0, p1, b, q1, q0]    b = [a, p1, p2, p3, c, q1]    c = [b, p3, p4, q2, q1]
```

**Existing vertex modifications**:
1. p0 (becomes deg-6): `insertAfter(p0, q0, a)` [DRight] / `insertAfter(p0, p1, a)` [DLeft]
2. p4 (becomes deg-6): `insertAfter(p4, p3, c)` [DRight] / `insertAfter(p4, prevCW(p4, p3), c)` [DLeft]
3. p1: `replaceNbr(p1, q0, a)`, `replaceNbr(p1, q1, b)`
4. p2: `replaceNbr(p2, q1, b)`
5. p3: `replaceNbr(p3, q1, b)`, `replaceNbr(p3, q2, c)`
6. q0: `replaceNbr(q0, p1, a)`
7. q1: `replaceNbr(q1, p1, a)`, `replaceNbr(q1, p2, b)`, `replaceNbr(q1, p3, c)`
8. q2: `replaceNbr(q2, p3, c)`
9. Degree-5 list: remove p0, p4; add a, c

**Key insight**: The bend vertex b has degree 6 and connects to THREE path
vertices (p1, p2, p3) because the turn causes q1 to be adjacent to all three
in the original graph.

#### 3.4.4 B_{i,j} (i+j > 0) Surgery

**Parameters**: `bentPos` = i, `bentLen` = i+j. Path has bentLen+5 entries, par
has bentLen+3 entries. `numNew` = bentLen+3 vertices. nv+0 = deg-5,
nv+numNew-1 = deg-5, rest deg-6. Bend vertex is nv+bentPos+1.

**New vertex CW neighbour lists** (using k for new vertex index 0..numNew-1):

Before bend (k = 0..bentPos):

k == 0 (deg-5):
```
DRight: [p[0], q[0], q[1], nv+1, p[1]]
DLeft:  [p[0], p[1], nv+1, q[1], q[0]]
```

k > 0 (deg-6):
```
DRight: [nv+k-1, q[k], q[k+1], nv+k+1, p[k+1], p[k]]
DLeft:  [nv+k-1, p[k], p[k+1], nv+k+1, q[k+1], q[k]]
```

Bend vertex (k = bentPos+1, deg-6):
```
DRight: [nv+k-1, q[k], nv+k+1, p[k+2], p[k+1], p[k]]
DLeft:  [nv+k-1, p[k], p[k+1], p[k+2], nv+k+1, q[k]]
```

After bend (k = bentPos+2..numNew-1):

k < numNew-1 (deg-6):
```
DRight: [nv+k-1, q[k-1], q[k], nv+k+1, p[k+2], p[k+1]]
DLeft:  [nv+k-1, p[k+1], p[k+2], nv+k+1, q[k], q[k-1]]
```

k == numNew-1 (deg-5):
```
DRight: [nv+k-1, q[k-1], q[k], p[k+2], p[k+1]]
DLeft:  [nv+k-1, p[k+1], p[k+2], q[k], q[k-1]]
```

**Note**: After the bend, the q-index shifts by -1 relative to the path index.
New vertex nv+k references q[k-1] and q[k] (not q[k] and q[k+1] as in the
pre-bend segment). This is because the turn edge is skipped in the parallel path.

**Existing vertex modifications**:
1. p[0] (becomes deg-6): `insertAfter(p[0], q[0], nv)` [DRight] / `insertAfter(p[0], p[1], nv)` [DLeft]
2. p[bentLen+4] (becomes deg-6): `insertAfter(p[bl+4], p[bl+3], nv+numNew-1)` [DRight] / `insertAfter(p[bl+4], prevCW(p[bl+4], p[bl+3]), nv+numNew-1)` [DLeft]
3. Pre-bend path (k = 1..bentPos+1): `replaceNbr(p[k], q[k-1], nv+k-1)`, `replaceNbr(p[k], q[k], nv+k)`
4. Pre-bend parallel (k = 0..bentPos): if k > 0: `replaceNbr(q[k], p[k], nv+k-1)`; `replaceNbr(q[k], p[k+1], nv+k)`
5. Bend area (bendI = bentPos+2): `replaceNbr(p[bendI], q[bendI-1], nv+bendI-1)`; `replaceNbr(q[bendI-1], p[bendI-1], nv+bendI-2)`; `replaceNbr(q[bendI-1], p[bendI], nv+bendI-1)`; `replaceNbr(q[bendI-1], p[bendI+1], nv+bendI)`
6. Post-bend path (k = bendI+1..bentLen+3): `replaceNbr(p[k], q[k-2], nv+k-2)`, `replaceNbr(p[k], q[k-1], nv+k-1)`
7. Post-bend parallel (k = bentPos+2..bentLen+2): `replaceNbr(q[k], p[k+1], nv+k)`; if k < bentLen+2: `replaceNbr(q[k], p[k+2], nv+k+1)`
8. Degree-5 list: remove p[0], p[bentLen+4]; add nv+0, nv+numNew-1

#### 3.4.5 F (Nanotube Ring) Surgery

**Parameters**: ring = [r0, r1, r2, r3, r4] (5-cycle of deg-6 vertices),
outer = [o0, o1, o2, o3, o4] (cap vertices where o[i] is adjacent to both r[i]
and r[(i+1)%5]). New vertices: nv+0 through nv+4, all degree-6.

**New vertex CW neighbour list** for nv+i:
```
[r[i], nv+(i-1+5)%5, o[i], o[(i+1)%5], nv+(i+1)%5, r[(i+1)%5]]
```

**Existing vertex modifications** (for i = 0..4):
- r[i]: `replaceNbr(r[i], o[(i-1+5)%5], nv+(i-1+5)%5)`, `replaceNbr(r[i], o[i], nv+i)`
- o[i]: `replaceNbr(o[i], r[i], nv+(i-1+5)%5)`, `replaceNbr(o[i], r[(i+1)%5], nv+i)`

Degree-5 list: unchanged (all new vertices are deg-6).

### 3.5 Path Computation Details

#### Straight path computation

```
advance(from, to) = (to, straightAhead(dir, to, from))
edges = iterate advance (u0, v0)
path = [fst(edges[0]), fst(edges[1]), ..., fst(edges[numEntries-1])]
par[k] = sideNbr(dir, path[k], path[k+1])   for k = 0..numEntries-2
```

where:
- `straightAhead(DRight, u, from)` = `advanceCW(u, from, 3)` (for both deg-5 and deg-6)
- `straightAhead(DLeft, u, from)` = `advanceCW(u, from, deg(u)-3)` (= 3 for deg-6, 2 for deg-5)
- `sideNbr(DRight, a, b)` = `prevCW(a, b)`
- `sideNbr(DLeft, a, b)` = `nextCW(a, b)`

#### Bent path computation

For B_{i,j} with `bentPos` = i, `bentLen` = i+j:

```
Phase 1 (pre-bend, bentPos+2 edges):
    preBendEdges = take(bentPos+2, iterate(advance, (u0, v0)))
    prePath = [fst(e) for e in preBendEdges]   -- bentPos+2 vertices

Phase 2 (turn):
    (preLast, turnV) = last(preBendEdges)
    afterTurn = turnAhead(dir, turnV, preLast)

Phase 3 (post-bend, bentLen-bentPos+1 edges):
    postBendEdges = take(bentLen-bentPos+1, iterate(advance, (turnV, afterTurn)))
    postPath = [fst(e) for e in postBendEdges]  -- bentLen-bentPos+1 vertices

Phase 4 (endpoint):
    (lastFrom, lastTo) = last(postBendEdges)
    finalV = straightAhead(dir, lastTo, lastFrom)

Full path = prePath ++ postPath ++ [lastTo, finalV]   -- bentLen+5 entries

Parallel path (turn edge SKIPPED):
    par[k] = sideNbr(dir, path[k], path[k+1])  for k = 0..bentLen+2
```

The key insight: the parallel path has bentLen+3 entries (not bentLen+5),
computed uniformly from the assembled full path. The turn edge contributes to
the full path but is not separately represented in the parallel path — it is
folded into the indexing of the path vertices.

For B_{0,0}: path has 5 entries [p0, p1, p2, p3, p4], par has 3 entries
[q0, q1, q2] where q0 = sideNbr(p0, p1), q1 = sideNbr(p2, p3),
q2 = sideNbr(p3, p4). Note that q1 is NOT sideNbr(p1, p2) — that edge is the
turn edge.

---

## 4. Reductions

A **reduction** is the inverse of an expansion. It is described by the same
triple (*kind*, *edge*, *dir*), where *edge* connects two degree-5 vertices
in the child graph.

**Enumerating all reductions of a graph** (`allReductions`):

For each reduction type, scan the graph for valid reduction sites:

**L_0 reductions**: For each degree-5 vertex *u*, for each neighbour *v* with
degree 5, check both directions. A direction *d* is valid iff the two "flanking"
vertices (2 positions away from the opposite endpoint at each end, in the walk
direction) both have degree 6:

```
isValidL0Direction(G, (u, v), DRight):
    return deg(advanceCW(u, v, 2)) == 6  AND  deg(advanceCW(v, u, 2)) == 6

isValidL0Direction(G, (u, v), DLeft):
    return deg(advanceCW(u, v, deg(u)-2)) == 6  AND  deg(advanceCW(v, u, deg(v)-2)) == 6
```

**L_i reductions (i >= 1)**: For each degree-5 vertex *u* and each neighbour *v*
(with *v* NOT degree-5), follow straight-ahead steps until reaching a degree-5
vertex or exceeding the maximum length. If a degree-5 vertex *w* is reached at
distance *d* (= *i* + 1), check flanking validity at both endpoints (same check
as L_0 but at (*u*, *v*) and (*w*, *prevW*)).

**B_{0,0} reductions**: For each degree-5 vertex *u* and each neighbour *v*,
compute the B_{0,0} path using `computeBentZeroPath`. The site is valid if the
path has no self-intersection, the far endpoint is degree-5, and path/parallel
path are disjoint.

**B_{i,j} reductions (i+j > 0)**: For each degree-5 vertex *u*, neighbour *v*,
direction *d*, and valid (*i*, *j*) values, compute the bent path using
`computeBentPathSafe`. The site is valid if:
1. The path has no self-intersection
2. All intermediate vertices are degree-6
3. The far endpoint is degree-5
4. The flanking vertex at the far endpoint is degree-6 and not on the path

The flanking vertex for bent reductions is on the **opposite side** from
straight reductions: for DLeft, use `advanceCW(w, prevW, 2)` (2 CW from
incoming); for DRight, use `advanceCW(w, prevW, deg(w)-2)` (2 CCW from
incoming).

**F reductions**: Search for a 5-cycle of degree-6 vertices by trying all
degree-6 edges and following `straightAhead` for 5 steps. If the path closes
to form a cycle, verify that each consecutive pair in the ring shares a common
neighbour not in the ring (the "outer" or "cap" vertices), that the 5 outer
vertices are distinct, and that each ring vertex's two outer neighbours are
CW-adjacent in its neighbour list.

---

## 5. The Canonical Construction Path

### 5.1 Core idea

Every reducible fullerene has a unique **canonical reduction** — the reduction
whose "representing triple" has the lexicographically smallest **canonical
ordering tuple**. The fullerene obtained by applying this reduction is the
**canonical parent**.

The generation algorithm:
1. Start from the three seed graphs.
2. For each graph *G*, enumerate expansion sites, apply each expansion to get
   a child *G'*, and accept *G'* if and only if the expansion's inverse is
   the canonical reduction of *G'*.
3. Recurse on accepted children.

This guarantees each isomorphism class is generated exactly once, without any
hash tables or isomorphism testing.

### 5.2 Representing triples

A reduction is represented by a triple (*e*, *x*, *d*) where *e* is a directed
edge at one endpoint, *x* is the type parameters (e.g., L_2 or B_{1,3}), and
*d* is a direction. Since the central path has two endpoints, each reduction has
exactly **two** representing triples — one from each end.

For the two triples of the same reduction:
- **L_i**: Both triples have the **same** direction.
- **B_{i,j}**: The triples have **opposite** directions and **swapped** parameters: if one is (*e*, B_{i,j}, DLeft), the other is (*e'*, B_{j,i}, DRight).

### 5.3 The canonical ordering tuple

Each representing triple is assigned a **5-tuple** (x0, x1, x2, x3, x4),
compared lexicographically. **Smaller is more canonical.** The canonical
reduction is the one with the smallest 5-tuple across all reductions and all
their representing triples.

The tuple components are of increasing discriminating power and cost:

| Component | Name | Cost | Description |
|-----------|------|------|-------------|
| x0 | Reduction length | O(1) | Distance between the two 5-vertices |
| x1 | Longest straight | O(1) | Negative of the longest straight segment |
| x2 | Colour pair | O(1) | 5-bit degree patterns near both endpoints |
| x3 | Path colour | O(1) | 7-bit straight-ahead degree pattern |
| x4 | BFS canonical form | O(n) | Full graph encoding from the edge |

**Lazy evaluation**: x_i is computed only for triples that are tied on
(x0, ..., x_{i-1}). At 300 carbon atoms, x0 through x3 resolve >99.9% of
cases without computing x4.

---

## 6. The 5-Tuple Components in Detail

### 6.0 x0: Reduction length

```
x0(L_i)   = i + 1       (distance between the two 5-vertices)
x0(B_{i,j}) = i + j + 2
```

Shorter reductions are more canonical. An L_0 reduction (length 1) always beats
any L_1 or B_{0,0} (length 2), which always beats any L_2 or B_{1,0} (length 3),
etc.

**Practical consequence**: Before computing x1-x4 for any candidate, first check
if the child has ANY reduction with shorter x0. If so, reject immediately.

### 6.1 x1: Longest straight segment (negated)

```
x1(L_i)     = -(i + 1)     (the entire path is straight)
x1(B_{i,j}) = -(max(i, j) + 1)   (longest segment before or after the bend)
```

Within the same x0, L_i reductions have x1 = -x0 (= -(i+1)), while B_{i,j}
reductions with the same total length have x1 = -(max(i,j)+1) > -x0 (larger =
worse). So **L always beats B of the same length**.

### 6.2 x2: Colour pair (5-bit neighbourhood invariant)

This is a pair of 5-bit integers encoding the degree-5/degree-6 pattern in the
neighbourhood of the reduction edge. The computation differs by reduction type.

#### Colour walk functions

**colourCW5(G, w, start)**: Starting from neighbour `start` at vertex *w*, walk
5 steps clockwise, recording a 5-bit integer where bit *i* (most significant
first) = 1 iff the *i*-th vertex visited has degree 5.

```
colour = 0
nbr = start
for i = 4 downto 0:
    if deg(nbr) == 5: colour |= (1 << i)
    nbr = nextCW(w, nbr)
return colour
```

**colourCCW5(G, w, start)**: Same but walking counter-clockwise (using `prevCW`).

#### Starting edges for STRAIGHT reductions (L_0 and L_i)

Given a representing triple with directed edge (*u* -> *v*) and direction *d*:

**DRight** (CW):
```
w = advanceCW(u, v, 2)                    // 2 CW from v around u
c = colourCW5(G, w, nextCW(w, u))         // 1 CW from u at w, then walk 5 CW
```

**DLeft** (CCW):
```
w = advanceCW(u, v, deg(u) - 2)           // 2 CCW from v around u
c = colourCCW5(G, w, prevCW(w, u))        // 1 CCW from u at w, then walk 5 CCW
```

Compute this for **both** endpoints of the reduction. For L_0, the two endpoints
are (*u*, *v*) and (*v*, *u*). For L_i (i >= 1), the second endpoint is at
the other end of the straight path.

The x2 value is:
```
c1 = colour at first endpoint
c2 = colour at second endpoint
x2 = (-(max(c1, c2)), -(min(c1, c2)))
```

The negation ensures that higher colour = smaller tuple = more canonical.

**Important for L_i (i >= 1)**: The second component of x2 is set to 0 (no
tiebreaker on the min colour). Only the max matters. So `x2 = (-(max(c1, c2)), 0)`.

#### Starting edges for BENT reductions (B_{0,0} and B_{i,j})

The bent colour computation is **asymmetric**: the two endpoints of the bent
path use **different** colour walk functions because the bend reverses the
local orientation.

Given test edges (*u1* -> *v1*) at the first endpoint and (*u2* -> *v2*) at the
second endpoint:

**DLeft** (use_next):
```
c1 = colourCCW5(G, v1, prevCW(v1, u1))    // CCW walk at first endpoint
c2 = colourCW5(G, v2, nextCW(v2, u2))     // CW walk at second endpoint
```

**DRight** (use_prev):
```
c1 = colourCW5(G, v1, nextCW(v1, u1))     // CW walk at first endpoint
c2 = colourCCW5(G, v2, prevCW(v2, u2))    // CCW walk at second endpoint
```

Note the asymmetry: one endpoint uses CCW, the other uses CW (or vice versa).
The starting edge is just one hop from the edge endpoint (`->invers->next` or
`->invers->prev`), NOT two hops like the straight types.

For B_{0,0}: `x2 = (-(max(c1, c2)), -(min(c1, c2)))`.
For B_{i,j} with i+j > 0: `x2 = (-(max(c1, c2)), 0)` (no min tiebreaker).

### 6.3 x3: Path colour (7-bit straight-ahead invariant)

The path colour follows a straight-ahead walk for 7 steps, recording a 7-bit
degree pattern.

**pathColour(G, dir, u, w)**: Starting from vertex *w* in the neighbourhood of
*u*, take 7 straight-ahead steps in direction *dir*, recording degree-5 bits:

```
colour = 0
(atVtx, fromVtx) = (w, u)
for i = 6 downto 0:
    next = straightAhead(dir, atVtx, fromVtx)
    if deg(next) == 5: colour |= (1 << i)
    (atVtx, fromVtx) = (next, atVtx)
return colour
```

#### Starting edges for STRAIGHT path colour

The path colour starting point is on the **OPPOSITE side** of the edge from the
colour starting point. Colour goes 2 steps **in** the walk direction; path
colour goes 2 steps **against** it.

**DLeft** (use_next):
```
w = advanceCW(u, v, 2)                    // 2 CW (opposite from colour's CCW)
pc = pathColour(G, DLeft, u, w)
```

**DRight** (use_prev):
```
w = advanceCW(u, v, deg(u) - 2)           // 2 CCW (opposite from colour's CW)
pc = pathColour(G, DRight, u, w)
```

Compute for both endpoints. `x3 = -(max(pc1, pc2))`.

#### Starting edges for BENT path colour

Also asymmetric, like the bent colour:

**DLeft** (use_next):
```
// First endpoint: CCW path from v1
w1 = advanceCW(v1, u1, deg(v1) - 2)       // 2 CCW from u1 at v1
pc1 = pathColour(G, DLeft, v1, w1)

// Second endpoint: CW path from v2
w2 = advanceCW(v2, u2, 2)                 // 2 CW from u2 at v2
pc2 = pathColour(G, DRight, v2, w2)
```

**DRight** (use_prev):
```
w1 = advanceCW(v1, u1, 2)                 // 2 CW from u1 at v1
pc1 = pathColour(G, DRight, v1, w1)

w2 = advanceCW(v2, u2, deg(v2) - 2)       // 2 CCW from u2 at v2
pc2 = pathColour(G, DLeft, v2, w2)
```

`x3 = -(max(pc1, pc2))`.

### 6.4 x4: BFS canonical form

If x0 through x3 fail to distinguish candidates, the full BFS encoding of the
graph determines the canonical reduction. This is also the mechanism by which the
automorphism group is computed.

#### BFS encoding algorithm

Given a directed edge (*u* -> *v*) and direction *dir*:

```
function bfsCode(G, u, v, dir):
    // Assign BFS numbers starting from 1
    number[u] = 1, number[v] = 2
    refEdge[1] = (u, v), refEdge[2] = (v, u)
    nextNum = 3

    // step(w, x) = next neighbour of w after x in walk direction
    step = if dir == DRight then nextCW else prevCW

    code = []
    for curNum = 1 to nv:
        (w, ref) = refEdge[curNum]
        // Walk deg(w)-1 neighbours starting after ref in the walk direction
        nbr = step(w, ref)
        repeat deg(w)-1 times:
            if nbr has a number:
                append number[nbr] to code
            else:
                append (deg(nbr) + nv + 1) to code    // "colour" = degree + offset
                number[nbr] = nextNum
                refEdge[nextNum] = (nbr, w)
                nextNum += 1
            nbr = step(w, nbr)
        append 0 to code      // separator

    return code
```

The offset `nv + 1` ensures colour values (>= nv + 6) cannot be confused with
vertex numbers (<= nv).

**BFS codes are compared lexicographically.** The smallest BFS code across all
starting edges and both directions defines the canonical form.

---

## 7. The Canonical Test (Rule 1)

### 7.1 isCanonical

After applying expansion *e* = (*kind*, (*u*, *v*), *dir*) to parent *G* with
*nv* vertices, producing child *G'* with *nv'* vertices:

```
function isCanonical(e, parentNV, G'):
    allReds = allReductions(G')
    allOrds = [canonOrd(r, G') for r in allReds]

    // Identify the inverse: a reduction in G' of the same kind,
    // whose starting edge connects TWO vertices that are BOTH newly added,
    // in the SAME direction as the expansion.
    numNew = numNewVertices(kind)
    newRange = [parentNV .. parentNV + numNew - 1]
    inverseReds = [r for r in allReds
                   where r.kind == kind
                     AND r.edge.u IN newRange
                     AND r.edge.v IN newRange
                     AND r.dir == dir]          // ← CRITICAL: direction must match

    if inverseReds is empty: return false

    invOrds = [canonOrd(r, G') for r in inverseReds]
    return min(invOrds) <= min(allOrds)
```

### 7.2 Critical details (from implementation experience)

**Both vertices must be new** (not just one). After an expansion adds vertices
*nv* through *nv* + *k* - 1, the inverse reduction's starting edge must have
BOTH endpoints in that range. A new degree-5 vertex adjacent to an OLD degree-5
vertex creates a spurious L_0 site that is NOT the expansion's inverse.

**Direction must match**. The expansion was applied in a specific direction
(DLeft or DRight). The canonical test must check the inverse reduction in that
SAME direction. If the opposite direction has a better `canonOrd`, the expansion
is rejected — this prevents duplicates when a reversing automorphism maps the
expansion edge to itself (swapping the two directions).

This direction-matching requirement is the single most critical detail omitted
from the original paper. Without it, the generator produces duplicates.

---

## 8. Rule 2: Expansion Orbit Filtering

Before applying Rule 1, we first reduce the number of expansions to try.
**Rule 2**: from each parent *G*, apply only one expansion per equivalence class.

### 8.1 Equivalence of expansions

Two expansions of *G* are **equivalent** if they produce isomorphic children.
The equivalence relation on expansion triples is generated by:

1. **Same reduction, other endpoint**: Each expansion has a "twin" from the other
   end of the central path (the "other triple").
2. **Orientation-preserving automorphism**: An OP automorphism sigma maps
   (*kind*, (*u*, *v*), *d*) to (*kind*, (sigma(*u*), sigma(*v*)), *d*) — same
   direction.
3. **Orientation-reversing automorphism**: An OR automorphism sigma maps
   (*kind*, (*u*, *v*), *d*) to (*kind*, (sigma(*u*), sigma(*v*)), flipDir(*d*))
   — **flipped** direction.

### 8.2 Computing the other triple

For a given expansion (*kind*, (*u*, *v*), *d*):

| Type | Other triple |
|------|-------------|
| L_0 | Follow the straight path from (*u*, *v*) to find `otherU` (= `parallelPath[2]`, the other degree-5 vertex). `otherV` = `sideNbr(flipDir(d), otherU, path[2])`. Direction is **same**: (*L_0*, (*otherU*, *otherV*), *d*). |
| L_i (i >= 1) | Same idea: follow the path to the other endpoint. Direction is **same**. |
| B_{0,0} | Other endpoint of the bent path. Direction is **flipped**: (*B_{0,0}*, (*otherU*, *otherV*), flipDir(*d*)). |
| B_{i,j} | Other endpoint, parameters **swapped** and direction **flipped**: (*B_{j,i}*, (*otherU*, *otherV*), flipDir(*d*)). |
| F | No other triple. |

### 8.3 Computing the automorphism group

The automorphism group is computed as a byproduct of the BFS canonical form. Try
all starting configurations (vertex, neighbour, direction). For each one, compute
the BFS code and numbering. Starting edges that produce the minimum BFS code
define automorphisms.

```
function canonicalBFSAndGroup(G):
    bestCode = infinity
    bestDir = null
    matches = []

    for each vertex u, neighbour v of u, direction d in {DLeft, DRight}:
        (code, numbering) = bfsWithNumbering(G, u, v, d)
        if code < bestCode:
            bestCode = code
            bestDir = d
            matches = [(numbering, d)]
        else if code == bestCode:
            matches.append((numbering, d))

    // Build automorphisms from matches
    refInv = inverse of matches[0].numbering   // BFS number → vertex
    auts = []
    for (numMap, d) in matches:
        perm = {v: refInv[numMap[v]] for each vertex v}
        orientation = Preserving if d == bestDir else Reversing
        auts.append(Automorphism(perm, orientation))

    return (bestCode, auts, count of Preserving auts)
```

### 8.4 Applying automorphisms to expansions

```
function applyAutToExp(aut, expansion):
    sigma = aut.perm
    (kind, (u, v), dir) = expansion
    if aut.orientation == Preserving:
        return (kind, (sigma(u), sigma(v)), dir)
    else:  // Reversing
        return (kind, (sigma(u), sigma(v)), flipDir(dir))
```

### 8.5 The filter algorithm

```
function filterByRule2(G, auts, expansions):
    seen = empty set
    result = []
    for e in expansions:
        key = (e.kind, e.edge.u, e.edge.v, e.dir)
        if key in seen: continue

        // Compute the full orbit
        orbit = {expKey(applyAutToExp(a, e)) for a in auts}
        other = computeOtherTriple(G, e)
        if other is not Nothing:
            orbit = orbit ∪ {expKey(applyAutToExp(a, other)) for a in auts}

        seen = seen ∪ orbit
        result.append(e)

    return result
```

### 8.6 The or_same_edge_found optimization

If an orientation-reversing automorphism maps a directed edge (*u*, *v*) to
itself (i.e., sigma(*u*) = *u* and sigma(*v*) = *v*), then expansion
(*kind*, (*u*, *v*), DLeft) and (*kind*, (*u*, *v*), DRight) are in the same
orbit. The `applyAutToExp` + orbit computation handles this automatically: the
reversing automorphism maps one direction to the other, so both end up in the
same orbit set.

---

## 9. The Generation Loop

```
function generate(maxDualVertices):
    seeds = [(c20, autGroup(c20)), (c28, autGroup(c28)), (c30, autGroup(c30))]
    results = empty map (dualVertices → list of graphs)

    for (seed, auts) in seeds:
        processTree(seed, auts, maxDualVertices, results)

    return results


function processTree(G, auts, maxDV, results):
    if numVertices(G) > maxDV: return
    record G in results[numVertices(G)]
    expandChildren(G, auts, maxDV, results)


function expandChildren(G, auts, maxDV, results):
    nv = numVertices(G)
    maxLen = min(maxDV - nv + 2, 5)    // maximum reduction length to try

    // Step 1: Enumerate all expansions up to maxLen
    exps = expansions(maxLen, G)

    // Step 2: Filter by Rule 2
    expsR2 = filterByRule2(G, auts, exps)

    // Step 3: Apply each expansion, test canonicity, recurse on accepted children
    for e in expsR2:
        G' = applyExpansion(e, G)
        if numVertices(G') > maxDV: continue
        if isCanonical(e, nv, G'):
            childAuts = autGroup(G')
            processTree(G', childAuts, maxDV, results)

    // Step 4: F expansion for nanotubes (always canonical)
    if findNanotubeRing(G) returns (ring, outer) AND nv + 5 <= maxDV:
        G' = applyRing(G, ring, outer)
        childAuts = autGroup(G')
        processTree(G', childAuts, maxDV, results)
```

### 9.1 Correctness argument

**Completeness**: Every reducible fullerene *G'* has a canonical reduction *r*.
The inverse expansion of *r* applied to the canonical parent *G* produces *G'*
(up to isomorphism). By induction from the seeds, every isomorphism class is
reached.

**No duplicates**: Suppose two expansions from the same parent *G* produce
isomorphic children *G'* and *G''*. Since they are isomorphic, they have the
same canonical reduction. If both expansions are canonical, they must be the
inverse of the same canonical reduction — but then they are equivalent under
Rule 2 and only one would have been tried. Contradiction.

---

## 10. Bounding Lemmas (Pruning)

The bounding lemmas restrict the maximum **reduction length** of expansions
tried at each node. Recall that reduction length is the distance between the
two degree-5 vertices:
- L_i has reduction length i + 1
- B_{i,j} has reduction length i + j + 2

The bounding function computes the tightest applicable bound by composing all
lemmas via `min`. Soundness is preserved: a tighter bound only skips expansions
that cannot be canonical, never misses valid ones.

### 10.1 Base bound: size ceiling

The simplest bound: an expansion cannot produce a child larger than the target.
L_i adds i+2 vertices and B_{i,j} adds i+j+3 vertices. So:

```
baseBound = maxDV - numVertices(G) - 1
```

This ensures no expansion overshoots the target vertex count.

### 10.2 Geometric bound

An expansion of reduction length *d* cannot be canonical if the child has fewer
than 12 * f(floor((d-1)/2)) vertices, where f(x) = 1 + (5/2) * x * (x+1).

This gives:

| Reduction length d | min child vertices |
|---|---|
| 1-2 | 12 |
| 3-4 | 72 |
| 5-6 | 192 |
| 7-8 | 372 |

To apply: find the maximum *d* such that nv + d + 1 >= the minimum child
vertex count for that *d*. For small graphs this is very permissive; for large
graphs it is mildly restrictive.

### 10.3 Lemma 1: Non-IPR bound

Every reducible non-IPR dual fullerene (one with adjacent degree-5 vertices) has
a reduction of length at most 2 (i.e., L_0, L_1, or B_{0,0}).

This is not directly used as a per-node bound, but justifies why the algorithm
only needs to consider short reductions for non-IPR graphs.

### 10.4 Lemma 2: Child bound

If *G* has a reduction of length *d* <= 2, every child of *G* has a reduction of
length at most *d* + 2.

**Application**:
- If parent has an L_0 reduction (length 1): max reduction length = 3
- If parent has an L_1 or B_{0,0} reduction (length 2): max reduction length = 4
- If parent has a reduction of length 3 (L_2, B_{1,0}, etc.): max reduction length = 5

### 10.5 Lemma 3: L_0 tightening

If *G* has an L_0 reduction, all **canonical** children of *G* have a reduction
of length at most 2.

**Application**: if parent has an L_0 reduction, max reduction length = 2.
This is strictly tighter than Lemma 2's bound of 3 for L_0 parents.

At C60, this eliminates ~26% of expansion sites (when applied alone). Combined
with other lemmas, the total reduction in expansion sites is ~70%.

### 10.6 Lemma 4: Two independent length-2 reductions

If *G* has at least 2 reductions of length 2 (L_1 or B_{0,0}) with **different
5-vertex pairs** (the sorted pair of degree-5 vertices involved), all canonical
children have a reduction of length at most 3.

**Application**: max reduction length = 3.

**Note**: "Different" means different sorted 5-vertex pairs. Two L_1 reductions
at the same pair of 5-vertices but in different directions count as one pair.

### 10.7 Lemma 5: Three independent length-2 reductions

If *G* has at least 3 reductions of length 2 (L_1 or B_{0,0}) with **pairwise
disjoint 5-vertex sets**, all canonical children have a reduction of length at
most 2.

**Application**: max reduction length = 2.

**Algorithm to check**: Extract the 5-vertex pair (a, b) from each length-2
reduction. Deduplicate by sorted pair. Then search (with backtracking) for 3
pairs that are pairwise vertex-disjoint. The search is:

```
function hasThreeIndependent(pairs):
    return search(pairs, index=0, count=0, usedVertices={})

function search(pairs, index, count, usedVertices):
    if count >= 3: return true
    if index >= len(pairs): return false
    (a, b) = pairs[index]
    if a not in usedVertices and b not in usedVertices:
        // Try taking this pair
        if search(pairs, index+1, count+1, usedVertices ∪ {a, b}):
            return true
    // Skip this pair
    return search(pairs, index+1, count, usedVertices)
```

### 10.8 Lemma 6: Two distant L_0s

If *G* has two L_0 reductions R1 and R2 with BFS distance > 4 between their
degree-5 vertex pairs, all canonical children have an L_0 reduction.

**Application**: max reduction length = 1 (only L_0 expansions).

**Distance computation**: For L_0 reduction at edge (a, b), the vertices are
a and b. The distance between two L_0 reductions (a1, b1) and (a2, b2) is
min(dist(x, y) for x in {a1, b1}, y in {a2, b2}), where dist is BFS distance
in the dual graph.

### 10.9 Extracting 5-vertex pairs from reductions

For the bounding lemmas, you need the pair of degree-5 vertices from each
reduction:

- **L_0**: The edge (u, v) — both u and v are degree-5.
- **L_i** (i >= 1): Starting vertex u (degree-5) and the far endpoint found
  by following straight-ahead i steps from u through v. The far endpoint is
  also degree-5.
- **B_{0,0}**: Starting vertex u (degree-5) and the last vertex of the bent
  zero path (degree-5).
- **B_{i,j}** (i+j > 0): Starting vertex u (degree-5) and the last vertex
  of the bent path (degree-5).

### 10.10 Complete bounding function

```
function maxExpansionLength(G, maxDV, reductions):
    nv = numVertices(G)
    baseBound = maxDV - nv - 1

    // Geometric bound
    geoBound = max d in [1..baseBound] such that
        nv + d + 1 >= 12 * f(floor((d-1)/2))
        where f(x) = 1 + (5/2)*x*(x+1)

    hasL0 = any reduction has kind L_0

    // Length-2 reductions (L_1 and B_{0,0})
    len2Pairs = deduplicate by sorted pair:
        { (min(a,b), max(a,b)) | Red(k, ...) in reductions, reductionLength(k) == 2,
          (a, b) = fiveVertices(Red(k, ...)) }

    // Lemma 3: L_0 → max 2
    lemma3 = if hasL0 then 2 else ∞

    // Lemma 5: ≥3 independent length-2 → max 2
    lemma5 = if hasThreeIndependent(len2Pairs) then 2 else ∞

    // Lemma 4: ≥2 different length-2 → max 3
    lemma4 = if |len2Pairs| >= 2 then 3 else ∞

    // Lemma 2: any length-2 → max 4; any L_0 → max 3
    lemma2 = if |len2Pairs| > 0 then 4
             else if hasL0 then 3
             else ∞

    // Lemma 6: two distant L_0s → max 1
    l0Pairs = deduplicate L_0 five-vertex pairs
    lemma6 = if hasTwoDistantL0s(G, l0Pairs) then 1 else ∞

    return max(1, min(baseBound, geoBound, lemma3, lemma5, lemma4, lemma2, lemma6))
```

### 10.11 Effectiveness

At C60 (32 dual vertices), applying all bounding lemmas reduces:
- Expansion sites from ~1.25M to ~380K (70% reduction)
- After Rule 2: from ~686K to ~181K (74% reduction)
- Total generation time: ~45% faster

At 300 carbon atoms (152 dual vertices), Lemmas 2-6 determine the bound in
93.9% of cases.

---

## 11. Expansion/Reduction Type Summary

| Type | x0 (length) | x1 (neg longest straight) | New vertices | Always non-IPR? |
|------|-------------|---------------------------|-------------|----------------|
| L_0 | 1 | -1 | 2 | Yes |
| L_1 | 2 | -2 | 3 | No |
| B_{0,0} | 2 | -1 | 3 | No |
| L_2 | 3 | -3 | 4 | No |
| B_{1,0} | 3 | -2 | 4 | No |
| B_{0,1} | 3 | -2 | 4 | No |
| L_i | i+1 | -(i+1) | i+2 | No |
| B_{i,j} | i+j+2 | -(max(i,j)+1) | i+j+3 | No |
| F | ∞ (always worst) | 0 | 5 | No |

Within the same x0, L always has a smaller (better) x1 than B.

---

## 12. Worked Example: L_0 Canonical Test

Consider expanding C20 (12 vertices) by L_0 at edge (0, 1) with direction
DRight to produce C24 (14 vertices), with new vertices 12 and 13.

1. Compute `allReductions(C24)`. This finds all L_0, L_i, B_{0,0}, and B_{i,j}
   reduction sites in C24.

2. Identify the inverse: among the L_0 reductions, find those where BOTH edge
   vertices are in {12, 13} and the direction is DRight. This should be exactly
   the reduction `Red(L_0, (12, 13), DRight)` (or `Red(L_0, (13, 12), DRight)`).

3. Compute `canonOrd` for the inverse and for all other reductions.

4. If `min(inverseOrds) <= min(allOrds)`, accept the child.

If C24 happens to have another L_0 reduction (at some other pair of degree-5
vertices) with a smaller `canonOrd`, the expansion is rejected — that other
reduction defines the canonical parent, and C24 should be reached by a different
expansion from a different parent.

---

## 13. Expected Isomer Counts

These are the known counts for validating a correct implementation:

| C_n | Dual vertices | Isomers |
|-----|--------------|---------|
| C20 | 12 | 1 |
| C24 | 14 | 1 |
| C26 | 15 | 1 |
| C28 | 16 | 2 |
| C30 | 17 | 3 |
| C32 | 18 | 6 |
| C34 | 19 | 6 |
| C36 | 20 | 15 |
| C38 | 21 | 17 |
| C40 | 22 | 40 |
| C42 | 23 | 45 |
| C44 | 24 | 89 |
| C46 | 25 | 116 |
| C48 | 26 | 199 |
| C50 | 27 | 271 |
| C52 | 28 | 437 |
| C54 | 29 | 580 |
| C56 | 30 | 924 |
| C58 | 31 | 1205 |
| C60 | 32 | 1812 |

---

## 14. Implementation Notes

### 14.1 Lazy evaluation of the 5-tuple

If your language supports lazy evaluation (e.g., Haskell), represent `canonOrd`
as a lazy tuple. Comparing two tuples lexicographically will automatically
short-circuit: if x0 differs, x1-x4 are never computed.

In strict languages, implement the cascade explicitly: compute x0 for all
candidates, filter to those tied at minimum x0, compute x1 for survivors, etc.

### 14.2 Performance characteristics

At C60 (32 dual vertices):
- ~1.25 million expansion sites before Rule 2
- ~686,000 after Rule 2 (Rule 2 eliminates ~45%)
- ~5,764 pass the canonical test (~1% of post-Rule-2)
- 0 duplicates (the direction-specific test is exact)

### 14.3 Seed graph encoding

The three seeds in planar code (hex bytes, 1-indexed neighbours):

**C20** (12 vertices, all degree-5):
```
0C 02 03 04 05 06 00 01 03 06 07 00 01 02 04 07 08 00 01 03 05 08 09 00
01 04 06 09 0A 00 01 05 02 0A 0B 00 02 01 03 08 0B 00 03 02 04 09 0B 00
03 04 05 0A 0B 00 05 04 06 0B 0C 00 06 05 02 07 0C 00 07 08 09 0A 0C 00
```
(Decode: first byte `0C` = 12 vertices; then each vertex's neighbour list
terminated by 00.)

**C28** (16 vertices):
```
10 02 03 04 05 06 00 01 03 09 06 00 01 02 04 0A 09 00 01 03 05 0B 0A 00
01 04 06 0C 0B 00 01 05 02 08 0C 00 0D 0E 0F 10 00 06 0E 0D 00 02 01 03
0A 0F 0E 00 03 02 04 0B 10 0F 00 04 03 05 0C 10 00 05 04 06 0D 10 00 07
06 08 0E 00 07 06 09 0F 00 09 08 0A 10 00 0B 0A 0C 07 0F 10 00
```

**C30** (17 vertices):
```
11 06 05 04 03 02 00 02 03 07 06 00 01 03 04 09 08 07 00 01 02 04 05 0A
09 00 01 03 05 06 0B 0A 00 01 04 06 0C 0B 00 01 05 07 0D 0C 00 02 01 06
08 0E 0D 00 02 07 09 0F 0E 00 03 02 08 0A 10 0F 00 04 03 09 0B 11 10 00
05 04 0A 0C 11 00 05 06 0B 0D 11 00 07 06 0C 0E 11 00 08 07 0D 0F 11 00
09 08 0E 10 11 00 0A 09 0F 0B 0C 0D 0E 00
```

Note: These are 0-indexed internally (subtract 1 from each neighbour after
decoding). The planar code uses 1-indexed neighbours with 0 as terminator.

### 14.4 Common pitfalls

1. **Not matching direction in the inverse test**. This is the #1 source of
   duplicate generation. The expansion's direction MUST match the inverse
   reduction's direction.

2. **Identifying the inverse by only one new vertex**. Both vertices of the
   inverse reduction's edge must be newly added. A new degree-5 vertex adjacent
   to an existing degree-5 vertex creates a spurious L_0 match.

3. **Using symmetric colours for bent reductions**. The two endpoints of a bent
   path use DIFFERENT colour walk functions (one CW, one CCW) because the bend
   reverses the local orientation.

4. **Forgetting that B_{i,j} other triple flips direction AND swaps i,j**. For
   L-type reductions, the other triple has the same direction. For B-type, the
   direction is flipped and parameters are swapped.

5. **Not checking flanking degree-6 for reduction validity**. An L_0 between two
   degree-5 vertices is only a valid reduction if specific flanking vertices have
   degree 6. Without this check, you'll find spurious reductions.

---

## Appendix A: Straight-Ahead at Degree-5 Vertices

At degree-6 vertices, straight-ahead is unambiguous: the vertex exactly opposite
(3 positions away in either CW or CCW order, since 6/2 = 3).

At degree-5 vertices, straight-ahead depends on the direction:
- **DRight**: `advanceCW(u, from, 3)` — 3 positions CW from the entry edge
- **DLeft**: `advanceCW(u, from, deg(u) - 3)` = `advanceCW(u, from, 2)` — 2
  positions CW from entry (= 3 positions CCW)

This difference matters for paths that pass through degree-5 vertices.

## Appendix B: Turn-Ahead (Bent Path)

The turn at a bent path vertex uses 2 positions instead of 3:
- **DRight**: `advanceCW(u, from, 2)` — 2 positions CW from entry
- **DLeft**: `advanceCW(u, from, deg(u) - 2)` — 2 positions CCW from entry

This creates a "side step" that distinguishes B-type from L-type paths.

## Appendix C: Reduction Length vs. Expansion Length

The **reduction length** of a type is the distance between the two degree-5
vertices *before* the reduction is applied:
- L_i: reduction length = i + 1
- B_{i,j}: reduction length = i + j + 2

The **expansion length** (number of new vertices) differs:
- L_i: i + 2 new vertices
- B_{i,j}: i + j + 3 new vertices

The **pathlength** in the generation loop is the reduction length + 1 (counting
edges, not vertices): L_0 has pathlength 2, L_1 has pathlength 3, B_{0,0} has
pathlength 3, etc.
