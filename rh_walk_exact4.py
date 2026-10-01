"""
Exact evaluation of the fourth-order (n = 4) polytope integrals of the confined walk: J_4 (distinct-prime groups of m_6) and
L_6 (groups of m_6 with one cancelling pair), by the same method as rh_walk_exact.py in four dimensions: exact rational vertices,
a recursive fan triangulation (cone from the centroid over the facets, each facet fanned recursively), exact monomial integration
over the simplices.  kappa_6 = J_4/12 + 2 L_6.  Slow in pure Python (hours); run in the background.

Usage: python3 rh_walk_exact4.py [J4|L6|both]
"""
import sys, itertools, math
from fractions import Fraction as Fr
from math import factorial
from rh_walk_exact import constraints, solve, simplex_monomial_integral

def vertices(rows, n):
    rows = list({(tuple(a), b) for a, b in rows})
    verts = set()
    for combo in itertools.combinations(rows, n):
        v = solve([list(a) for a, _ in combo], [b for _, b in combo])
        if v is None: continue
        if all(sum(ai*vi for ai, vi in zip(a, v)) <= b for a, b in rows): verts.add(tuple(v))
    return [list(v) for v in verts], rows

def orthobasis(normals, n):
    """orthonormal basis (floats) of the orthogonal complement of the given normal vectors in R^n"""
    basis = []
    cand = [[1.0 if i == j else 0.0 for i in range(n)] for j in range(n)]
    N = [[float(c) for c in v] for v in normals]
    def proj_out(v, B):
        for b in B:
            d = sum(x*y for x, y in zip(v, b)); v = [x - d*y for x, y in zip(v, b)]
        return v
    Nn = []
    for v in N:
        v = proj_out(v, Nn); nr = math.sqrt(sum(x*x for x in v))
        if nr > 1e-9: Nn.append([x/nr for x in v])
    for v in cand:
        v = proj_out(v, Nn + basis); nr = math.sqrt(sum(x*x for x in v))
        if nr > 1e-9: basis.append([x/nr for x in v])
    return basis

def fan(V, active_rows, all_rows, n, dim):
    """simplices (lists of exact points) triangulating the convex hull of V, which lies in the affine subspace of dimension dim
    cut out by the active constraints (equalities); recursive cone from the centroid over the facets within that subspace"""
    if len(V) < dim + 1: return []
    cen = [sum(v[j] for v in V)/len(V) for j in range(n)]
    if dim == 1:
        return [[V[0], V[-1]]] if len(V) >= 2 else []
    if dim == 2:
        normals = [a for a, _ in active_rows]; B = orthobasis(normals, n)
        fc = [float(c) for c in cen]
        def ang(v):
            d = [float(v[j]) - fc[j] for j in range(n)]
            return math.atan2(sum(d[j]*B[1][j] for j in range(n)), sum(d[j]*B[0][j] for j in range(n)))
        V = sorted(V, key=ang)
        return [[cen, V[i], V[(i+1) % len(V)]] for i in range(len(V))]
    out = []
    seen = set()
    for a, b in all_rows:
        if any(a == aa and b == bb for aa, bb in active_rows): continue
        face = [v for v in V if sum(ai*vi for ai, vi in zip(a, v)) == b]
        if len(face) < dim: continue
        key = frozenset(tuple(v) for v in face)
        if key in seen: continue                 # two hyperplanes can cut the same face of this face: count it once
        seen.add(key)
        # a lower-dimensional face passing the vertex-count test gives degenerate simplices, which integrate to zero
        for simp in fan(face, active_rows + [(a, b)], all_rows, n, dim - 1):
            out.append([cen] + simp)
    return out

def polytope_integral(rows, n):
    V, rows = vertices(rows, n)
    if len(V) < n + 1: return Fr(0)
    total = Fr(0)
    for simp in fan(V, [], rows, n, n):
        if len(simp) == n + 1: total += simplex_monomial_integral(simp, n)
    return total

def J4():
    n = 4; total = Fr(0); perms = list(itertools.permutations(range(n)))
    for eps in itertools.product((1, -1), repeat=n):
        valid = [p for p in perms if eps[p[0]] == 1]
        for p in valid:
            for q in valid:
                total += polytope_integral(constraints(n, p, eps) + constraints(n, q, eps), n)
        print(f"  pattern {eps} done, running total {float(total):.6f}", flush=True)
    return total

def L6():
    """L_6 = (1/2) sum_{e1,e2} int x y [int z N(x,y,z) dz]^2: variables (x, y, z, z'), steps {z+, z-, x^e1, y^e2} ordered by p on
    (x,y,z) and by q on (x,y,z'); the monomial x y z z' over the product polytope."""
    n = 4; perms = list(itertools.permutations(range(4))); total = Fr(0)
    for e1 in (1, -1):
        for e2 in (1, -1):
            eps = (e1, e2, 1, -1)            # variables: 0 = x, 1 = y, 2 = z(+), 3 = z(-)  for the walk with z, z'
            # first step in; an ordering that begins z+, z- returns to the edge itself, which is not a point of the echo set
            # (the partial sum is then identically 0 and the strict inequality fails): such orderings are excluded
            valid = [p for p in perms if eps[p[0]] == 1 and set(p[:2]) != {2, 3}]
            for p in valid:
                for q in valid:
                    # p acts on (x, y, z, zbar) with the cancelling pair z (coordinates 2, 3 both = z); q likewise with z'
                    rows = []
                    for order, zi in ((p, 2), (q, 3)):
                        # build constraints in 4 variables (x, y, z, z'): the walk uses step lengths x, y, z_var, z_var
                        part = [Fr(0)]*4
                        for k in order:
                            part = part[:]
                            if k in (2, 3): part[zi] += Fr(eps[k])
                            else: part[k] += Fr(eps[k])
                            rows.append(([-c for c in part], Fr(0))); rows.append((part[:], Fr(1)))
                    for i in range(4):
                        a = [Fr(0)]*4; a[i] = Fr(-1); rows.append((a, Fr(0)))
                        a = [Fr(0)]*4; a[i] = Fr(1); rows.append((a, Fr(1)))
                    total += polytope_integral(rows, 4)
            print(f"  (e1,e2) = ({e1},{e2}) done, running total {float(total):.6f}", flush=True)
    return total/2

if __name__ == "__main__":
    what = sys.argv[1] if len(sys.argv) > 1 else "both"
    j4 = Fr(3329, 1680) if what == "L6" else None
    if what in ("J4", "both"):
        j4 = J4(); print(f"J_4 = {j4} = {float(j4):.9f};  J_4/12 = {j4/12} = {float(j4/12):.9f}", flush=True)
    if what in ("L6", "both"):
        l6 = L6(); print(f"L_6 = {l6} = {float(l6):.9f};  2 L_6 = {2*l6} = {float(2*l6):.9f}", flush=True)
    if what in ("both", "L6"):
        k6 = j4/12 + 2*l6; print(f"kappa_6 = J_4/12 + 2 L_6 = {k6} = {float(k6):.9f};  kappa_6/kappa_2^3 = {k6/Fr(5,12)**3} = {float(k6/Fr(5,12)**3):.6f}", flush=True)
