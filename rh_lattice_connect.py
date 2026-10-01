"""
Connectivity of the interior echo lattice (paper Proposition 9.16).  The closures of the two edges are disjoint as point sets at a
non-resonant support (a common point would put 2a in the lattice), so the component that carries the edge's measure mu_b is the
interior of the right edge's closure.  This script builds that closure by breadth-first search (rh_lattice_band.lattice, the
points keyed by edge), removes the edge vertex, and counts the connected components of the interior among the first N points,
reporting how many of the first-generation points (a - log m) fall in the largest component and the depth at which they join.

Usage: python3 rh_lattice_connect.py [N] [a ...]      (default N = 20000, a = 0.85 0.9 1.0 1.25 1.5 2.0)
"""
import math, sys
from rh_lattice_band import lattice

def components(n, adj, skip):
    comp = [-1]*n; c = 0
    for s in range(n):
        if comp[s] >= 0 or s in skip: continue
        stack = [s]; comp[s] = c
        while stack:
            u = stack.pop()
            for v, _ in adj[u]:
                if v in skip or comp[v] >= 0: continue
                comp[v] = c; stack.append(v)
        c += 1
    return comp, c

if __name__ == "__main__":
    args = sys.argv[1:]; N = int(args[0]) if args and args[0].isdigit() else 20000
    supports = [float(x) for x in (args[1:] if args and args[0].isdigit() else args)] or [0.85, 0.9, 1.0, 1.25, 1.5, 2.0]
    for a in supports:
        xs, gen, adj, steps = lattice(a, N); n = len(xs)
        right = [j for j in range(n) if xs[j] is not None]      # placeholder, refined below
        # xs are positions; the edge vertices are indices 0 (+a) and 1 (-a); points of the right closure are those reached from +a
        from rh_lattice_band import lattice as _l
        comp, c = components(n, adj, skip={0, 1})
        # which component holds the right edge's first generation (neighbours of vertex 0)?
        first = [v for v, _ in adj[0]]
        comps_first = {comp[v] for v in first}
        sizes = {}
        for j in range(n):
            if j in (0, 1): continue
            sizes[comp[j]] = sizes.get(comp[j], 0) + 1
        big = max(sizes, key=sizes.get)
        print(f"a = {a}: {n} points (chains <= {max(gen)}), interior components {c}, largest {sizes[big]} points; first-generation points of +a: {len(first)}, "
              f"in {len(comps_first)} component(s) {sorted(comps_first)}; -a's first generation in components {sorted({comp[v] for v, _ in adj[1]})}", flush=True)
