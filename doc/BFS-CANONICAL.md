# BFS Canonical Form

The BFS canonical form is the ultimate tiebreaker (x5) in the canonicity cascade. It produces a **complete structural fingerprint** of the graph as seen from a specific directed edge. Two reductions produce the same BFS code if and only if there's a graph automorphism mapping one to the other.

## The algorithm

Given a directed edge `(u -> v)` and a walk direction (CW or CCW):

1. **Seed the BFS**: Assign number 1 to `u`, number 2 to `v`. These are the first two vertices in the BFS queue.

2. **Process each vertex in BFS order**: For vertex `w` (discovered via edge from `ref`), walk its neighbors in planar order starting from the neighbor *after* `ref`, for `deg(w) - 1` steps (i.e., all neighbors except `ref`). For each neighbor:
   - If **already numbered**: emit its BFS number (an integer in `1..n`)
   - If **not yet seen**: emit its **colour** (degree + offset, where offset = n+1 so colours can't collide with BFS numbers), then assign it the next available BFS number and record its discovery edge

3. **Terminate** each vertex's neighbor list with `0`.

4. The output is a flat list of integers -- the **BFS code**.

## A concrete example

Imagine processing vertex 1 (= `u`) in a degree-5 graph. Its reference is `v` (= vertex 2). We walk 4 neighbors of `u` starting after `v` in CW order. Say they are vertices `a`, `b`, `c`, `d`:
- `a` is new -> emit `colour(a)`, assign number 3
- `b` is new -> emit `colour(b)`, assign number 4
- `c` is new -> emit `colour(c)`, assign number 5
- `d` is new -> emit `colour(d)`, assign number 6
- Emit `0` (end of vertex 1's list)

Then process vertex 2 (= `v`), walking from the neighbor after `u`. Some of these may already be numbered (e.g., vertex 3 = `a` if `a` is also a neighbor of `v`), so we'd emit `3` instead of a colour.

## Why it works

The BFS code captures the **entire combinatorial structure** of the planar graph relative to a starting directed edge. Two different starting edges on the same graph produce the same BFS code if and only if there exists a planar-embedding-preserving automorphism mapping one starting edge to the other.

The key details that make it work:

- **Planar order matters**: Walking neighbors in CW (or CCW) order means the code encodes the embedding, not just the abstract graph. This is essential because fullerene duality is defined on the embedded (planar) graph.
- **The colour vs. number distinction**: New vertices emit their degree (as a colour), already-visited vertices emit their BFS number. This records both the local structure (what degrees appear) and the global connectivity (which previously-seen vertex is reconnected to).
- **The `0` terminator**: Separates each vertex's neighbor list, making the code unambiguous.

## Why it's expensive (O(n)) but rarely needed

Every vertex and edge is visited exactly once -- linear time. But the cheap filters x0-x4 are just local neighbourhood probes (5-bit and 7-bit bitstrings from a fixed number of hops). At 152 dual vertices, these local probes are enough to distinguish >99.9% of reductions. The BFS code only fires when two reductions happen to have identical local neighbourhoods out to ~7 hops -- which is extremely rare in practice.

## The automorphism connection

If you compute the BFS code from *every* possible starting edge and find that two starting edges produce the *same* code, the vertex-renaming that maps one BFS numbering to the other is a graph automorphism. So the full automorphism group falls out as a byproduct: it's exactly the set of starting edges that achieve the minimum BFS code.
