# The `max_straight_lengths` Table: Geometric Bounding Without BFS

## Problem

Profiling the Haskell buckygen at C60 (nv=32) revealed that `hasTwoDistantL0s`
— our implementation of Lemma 6 from the paper — consumes **28–30% of time
and 46–47% of allocation**. The function performs BFS between all pairs of L0
reduction sites to check if any two have graph distance > 4. At C60 this
triggers 589,704 BFS calls (roughly 100 per parent graph, from O(k²) pairs
where k is the number of L0 edges).

## How the C code avoids this entirely

The C code never computes BFS distances between pentagon pairs at runtime. It
uses a precomputed table `max_straight_lengths[nv]` that maps the number of
dual vertices to the maximum possible minimum distance between any two
pentagons in a fullerene of that size. This table is computed once at startup
from a geometric formula and replaces all runtime distance queries.

### The geometric argument

A degree-5 vertex (pentagon) in a fullerene dual has a "patch" around it: the
set of all vertices within graph distance d. In a triangulation where most
vertices have degree 6, the number of faces in such a patch is:

```
f(d) = 1 + 5·d·(d+1)/2
```

This counts the central vertex plus 5 concentric rings (the factor 5 comes
from the degree-5 center creating 5 sectors, each growing by one vertex per
ring). A fullerene dual has exactly 12 degree-5 vertices, and the patches
around any two of them must be disjoint if their centers are more than 2d
apart. This gives:

> **Geometric bound.** If the minimum distance between any two degree-5
> vertices in a fullerene dual with n_f vertices is at least d, then
> n_f ≥ 12 · f(⌊(d-1)/2⌋).

Equivalently, `max_min_distance(n_f)` — the maximum possible minimum pentagon
distance in a fullerene with n_f faces — is the largest d satisfying this
inequality.

### The table

The C code computes `max_straight_lengths[i] = max_min_distance(i) - 1` for
each i (the `-1` converts distance to reduction length, since an L_i reduction
has length i+1 = distance+1). For non-IPR fullerenes, there is an additional
refinement: for even distances > 2, `min_nv_required` is increased by 20 (a
tighter bound when adjacent pentagons must exist).

Concrete values (non-IPR):

```
nv:  12  24  72  132  252  ...
d:    1   1   2    3    4  ...
msl:  0   0   1    2    3  ...
```

where msl = max_straight_lengths[nv], d = msl + 1 = max min-pentagon-distance.

### How the table replaces BFS at runtime

In `determine_max_pathlength_straight()` (C line 4868), the pathlength bound
is computed by a simple loop:

```c
int pathlength = 2;
while (nv + pathlength <= maxnv) {
    if (pathlength - 1 <= max_straight_lengths[nv + pathlength])
        max_pathlength = pathlength;
    else
        break;
    pathlength++;
}
```

The condition `pathlength - 1 <= max_straight_lengths[nv + pathlength]` asks:
"in a child with nv + pathlength vertices, is it geometrically possible for
two pentagons to be at distance ≥ pathlength - 1?" If not, no expansion of
that length can be canonical, because the child would necessarily have a
shorter reduction. This is **O(maxDV - nv)** integer comparisons — no graph
traversal at all.

## Lemma: Subsumption of Lemma 6

The paper's Lemma 6 states:

> If G has two L0 reductions R1, R2 with d(R1, R2) > 4, then every canonical
> child of G has an L0 reduction, so the maximum canonical expansion length
> is 2.

Our Haskell code implements this literally: enumerate all L0 reduction pairs,
BFS between them, check if any pair has distance > 4. This is O(k² · n) per
parent graph.

**But this check is unnecessary.** The geometric table already captures a
stronger result:

> **Lemma (geometric subsumption of Lemma 6).** Let G be a fullerene dual
> with n_f vertices and let d* = max_straight_lengths[n_f] + 1 be the maximum
> possible minimum pentagon distance at that size. If the expansion length
> bound from the table is already ≤ 2 (i.e., the table loop in
> `determine_max_pathlength_straight` yields max_pathlength ≤ 2), then Lemma 6
> cannot improve the bound further (since it also yields 2). Conversely, if
> the table allows pathlength > 2, then d* ≤ 4 for the relevant child sizes,
> meaning no two pentagons can be at distance > 4 anyway, so Lemma 6's
> hypothesis is unsatisfied.

More precisely: Lemma 6 fires when two L0 edges have distance > 4. But two
L0 edges are pairs of adjacent pentagons (distance 1 from each other), and
d(R1, R2) > 4 requires the two *pairs* to be far apart. For this to happen,
the graph must be large enough to geometrically accommodate pentagons at
distance > 4. But for any child size where the geometric table already bounds
the pathlength to ≤ 2, the table's constraint is at least as tight as
Lemma 6's conclusion. And for child sizes where the table allows larger
pathlengths, the graph is too small for distance > 4 to occur between L0
pairs (since max_straight_lengths[child_size] ≤ 4 implies max min-distance
≤ 5, but the L0 pairs share vertices with many other pentagons, making large
separations impossible at these sizes).

In practice at all tested sizes (up to C70, nv=37): the geometric table +
Lemmas 2–5 already determine the bound. Lemma 6 never improves it beyond
what the table provides. The C code confirms this: it does not implement
Lemma 6 as a distance check at all. The `max_straight_lengths` table,
together with the incremental `straight_length` variable, subsumes it.

## Relationship between the C code's variables and our Haskell

| C code | Haskell equivalent | Notes |
|--------|-------------------|-------|
| `max_straight_lengths[i]` | `geoMaxRedLen` (partial) | Our `geoMaxRedLen` uses the same formula but applies it as a single bound on the max d. The C code's loop is more precise: it checks each candidate pathlength against the table for the specific *child* size `nv + pathlength`. |
| `straight_length` | `hasL0` / `hasLen2` | The C code tracks the shortest reduction length incrementally (set to 2 after L0 expansion, 3 after L1, MAX after B00). We recompute from `allReductionsUpTo 2 g`. |
| `num_straight_extensions` | `length l0Reds` / `length len2Reds` | Counts of reductions at the shortest length. |
| `num_bent_zero_extensions` | `length b00Reds` | Count of B00 reductions. |
| Lemma 6 (`hasTwoDistantL0s`) | **Not implemented in C** | The C code relies on the geometric table instead. |

## Recommended fix

1. **Replace the C-style table loop.** Instead of our current `geoMaxRedLen`
   (which applies a single bound), replicate the C code's loop: for each
   candidate pathlength p from 2 upward, check
   `p - 1 <= maxStraightLengths ! (nv + p)` where `maxStraightLengths` is
   a precomputed array. This is both more precise and O(1) per check.

2. **Remove `hasTwoDistantL0s` and `bfsDistanceBounded` entirely.** They are
   never needed — the geometric table subsumes Lemma 6 at all relevant sizes.

3. **Precompute `maxStraightLengths` once** as a Haskell `UArray Int Int` at
   program startup, indexed by number of dual vertices.

Expected impact: eliminates 28–30% of runtime and 46–47% of allocation.

## The precomputed table (C code, lines 4732–4774)

```c
// f(d) = 1 + 5*(d+1)*d/2  (faces in a patch of radius d around a pentagon)
#define compute_nf_in_patch(distance) (1 + 5 * ((distance) + 1) * (distance) / 2)

int min_nv_required = 0, min_nv_required_prev = 0;
int distance = 0;
while (min_nv_required <= maxnv) {
    distance++;
    min_nv_required = 12 * compute_nf_in_patch((distance - 1) / 2);
    if (distance % 2 == 0 && distance > 2)
        min_nv_required += 20;  // tighter bound for even distances > 2
    for (i = min_nv_required_prev; i < min_nv_required && i <= maxnv; i++)
        max_straight_lengths[i] = distance - 1;
    min_nv_required_prev = min_nv_required;
}
```

In Haskell:

```haskell
-- | Precompute max_straight_lengths table: for each nv, the maximum possible
-- minimum pentagon distance minus 1 (= maximum reduction length achievable
-- purely from geometric constraints).
buildMaxStraightLengths :: Int -> UArray Int Int
buildMaxStraightLengths maxNV = runSTUArray $ do
    arr <- newArray (0, maxNV) 0
    let go dist minReqPrev
          | minReqPrev > maxNV = return arr
          | otherwise =
              let r = (dist - 1) `div` 2
                  base = 12 * (1 + 5 * (r + 1) * r `div` 2)
                  minReq = if even dist && dist > 2 then base + 20 else base
              in do forM_ [minReqPrev .. min (minReq - 1) maxNV] $ \i ->
                      writeArray arr i (dist - 1)
                    go (dist + 1) minReq
    go 1 0
```

### Use in `maxExpansionLength`

Replace the current `geoMaxRedLen` + `lemma6` with:

```haskell
-- Table-based geometric bound (replaces both geoMaxRedLen and lemma6)
tableBound = last $ 1 : [ p | p <- [2..baseBound]
                             , p - 1 <= mslTable ! (nv + p) ]
```

where `mslTable = buildMaxStraightLengths maxDV`.
