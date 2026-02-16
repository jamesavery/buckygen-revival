# Windup Bug Analysis

## Summary

`windupSpiral` produces invalid graphs for most spirals with 16+ vertices. The bug
is in the cascade logic of `stepK`: `cascFwd` can cascade past the boundary into
nodes that belong to the back end, causing duplicate edges and corrupted adjacency
lists.

## Smallest Failing Case

C28 D2 isomer, RSPI `[1,2,3,4,5,6,8,12,13,14,15,16]`, spiral `[5,5,5,5,5,5,6,5,6,6,6,5,5,5,5,5]`.

C20 (12v), C24 (14v), and C28 Td (16v) all pass. C28 D2 (16v) is the first failure.

## Root Cause

At step k=14 (the last peel step for this 16-vertex graph), the boundary is:

```
[9(0), 10(1), 11(1), 13(0)]
```

The connect-forward decrements 9 to -1 (it was already 0). Then `cascFwd` triggers:
1. Pop 9 (remaining=0 after connect_forward decremented it to -1... actually remaining was 0 before connect_forward)

Wait — let me be more precise. The initial boundary is `[9(1), 10(1), 11(1), 13(1)]`.

After `connect k-front(=9)`: `decAt 0` makes it `[9(0), 10(1), 11(1), 13(1)]`.
After `connect k-back(=13)`: `decLast` makes it `[9(0), 10(1), 11(1), 13(0)]`.

Now `cascFwd` runs:
1. Front is 9 with remaining=0 → pop 9, connect k-10, dec 10: `[10(0), 11(1), 13(0)]`
2. Front is 10 with remaining=0 → pop 10, connect k-11, dec 11: `[11(0), 13(0)]`
3. Front is 11 with remaining=0 → pop 11, connect k-13, dec 13: `[13(-1)]`
4. Front is 13 with remaining=-1 (≠0) → stop.

**The problem**: At step 3, `cascFwd` connects k=14 to vertex 13, but vertex 13 was
already connected to k=14 by the initial `connect k-back`. This creates a **duplicate
edge** (14,13). Vertex 13 ends up with remaining=-1 (negative!), and the final
adjacency list for vertex 13 has 14 listed twice.

## The Invariant Being Violated

The cascade should never pop a node that was the back endpoint (or vice versa).
`cascFwd` and `cascBwd` operate on opposite ends of the boundary deque. When they
converge in the middle, the algorithm must stop — the nodes at the other end have
already been connected.

In the C++ reference implementation (`triangulation.cc`), this is handled by the
loop structure: `connect_forward` and cascading happen first, then `connect_backward`
and cascading. Both track their own end and don't cross. The key difference is that
the C++ code uses a circular buffer with explicit front/back indices, so it naturally
can't cascade past the opposite end.

## How to Fix

The cascade functions must not pop nodes past the opposite end. Options:

1. **Pass the back node to cascFwd and front node to cascBwd** as a stopping
   sentinel. `cascFwd` stops if the next front IS the back node (or if boundary
   length reaches 1).

2. **Check boundary length**: `cascFwd` should stop when `Seq.length ov <= 1`
   (only one node left = it's the back node, don't eat it).

3. **Match the C++ structure more closely**: The C++ code does
   `connect_forward`, `cascade_forward`, `connect_backward`, `cascade_backward`
   in sequence, with each phase only operating on its own end. The issue may also
   be that after `cascFwd` eats a node, `cascBwd` should see the updated boundary
   but currently the Haskell code runs them independently.

Actually, re-reading the trace more carefully: the problem specifically occurs when
both the front AND back nodes reach remaining=0 simultaneously after the initial
two connects. The cascade from one end then eats into the other end's territory.

## Affected Sizes

- C20 (12v): PASS
- C24 (14v): PASS
- C28 Td (16v): PASS (lucky — cascades don't collide for this specific spiral)
- C28 D2 (16v): **FAIL**
- C32+ : All tested RSPIs fail

## Test Data

`WindupTestdata.hs` contains 22 expansion-generated graphs (C20–C40, 2 per size)
with known-correct CCW adjacency lists and their canonical spirals. These can be
used to test windup round-trips:

```haskell
import WindupTestdata
import Spiral (canonicalSpiral)
import Windup (windupSpiral)

-- For each test graph, extract spiral, wind up, extract again, compare
testRoundTrip (name, g, expectedSpiral) =
    let rebuilt = windupSpiral expectedSpiral
        sp2 = canonicalSpiral rebuilt
    in sp2 == Just expectedSpiral
```

`WindupFail.hs` contains a validation harness that checks `fromRSPI` output for
degree violations, asymmetric edges, and duplicate neighbors.
