"""
Exact evaluation of the polytope integrals of the confined walk (paper, Proposition 9.39).

J_n = int_{[0,1]^n} x_1...x_n S_n(x) dx with S_n = sum over the 2^n direction patterns of N(x,pattern)^2, N the number of orderings
of the n signed steps whose first step goes in and whose partial sums stay in (0,1).  Expanding N^2 = sum_{ord,ord'} 1_P 1_P',
J_n is a sum of integrals of the monomial x_1...x_n over the polytopes P_ord & P_ord', each cut out by inequalities with
coefficients in {-1,0,1}; every vertex is rational, so the integral is computed exactly: vertices by exact 3x3 solves (fractions),
a Delaunay triangulation of the vertex set (floating point, combinatorics only), and the monomial integrated exactly over each
simplex by the standard-simplex formula.  J_2 = 5/12 is the check.  Also the exact one-dimensional integral 11/60 of the
cancelling-pair groups and the resulting kappa_4 = J_3/3 + 11/60 and kurtosis kappa_4/kappa_2^2.

Usage: python3 rh_walk_exact.py [nmax]        (default nmax = 3; n = 4 is slow in pure Python)
"""
import sys, itertools
from fractions import Fraction as Fr
from math import factorial
import math

def constraints(n, order, eps):
    """inequalities a.x <= b for the ordering 'order' of the signed steps: partial sums in (0,1), cube [0,1]^n"""
    rows = []
    for i in range(n):
        a = [Fr(0)]*n; a[i] = Fr(-1); rows.append((a, Fr(0)))      # x_i >= 0
        a = [Fr(0)]*n; a[i] = Fr(1); rows.append((a, Fr(1)))       # x_i <= 1
    part = [Fr(0)]*n
    for k in order:
        part = part[:]; part[k] += Fr(eps[k])
        rows.append(([-c for c in part], Fr(0)))                    # partial sum >= 0  (strict a.e.)
        rows.append((part[:], Fr(1)))                               # partial sum <= 1
    return rows

def solve(A, b):
    """exact Gaussian elimination for a square system; None if singular"""
    n = len(A); M = [row[:] + [bb] for row, bb in zip(A, b)]
    for c in range(n):
        piv = next((r for r in range(c, n) if M[r][c] != 0), None)
        if piv is None: return None
        M[c], M[piv] = M[piv], M[c]
        for r in range(n):
            if r != c and M[r][c] != 0:
                f = M[r][c]/M[c][c]; M[r] = [x - f*y for x, y in zip(M[r], M[c])]
    return [M[i][n]/M[i][i] for i in range(n)]

def vertices(rows, n):
    rows = list({(tuple(a), b) for a, b in rows})                   # dedupe
    verts = set()
    for combo in itertools.combinations(rows, n):
        v = solve([list(a) for a, _ in combo], [b for _, b in combo])
        if v is None: continue
        if all(sum(ai*vi for ai, vi in zip(a, v)) <= b for a, b in rows): verts.add(tuple(v))
    return [list(v) for v in verts]

def simplex_monomial_integral(V, n):
    """exact int over the simplex with vertices V (n+1 points) of x_1 ... x_n"""
    v0 = V[0]; D = [[V[i+1][j] - v0[j] for j in range(n)] for i in range(n)]   # rows: edge vectors
    # determinant (exact)
    def det(M):
        M = [r[:] for r in M]; d = Fr(1); m = len(M)
        for c in range(m):
            piv = next((r for r in range(c, m) if M[r][c] != 0), None)
            if piv is None: return Fr(0)
            if piv != c: M[c], M[piv] = M[piv], M[c]; d = -d
            d *= M[c][c]
            for r in range(c+1, m):
                f = M[r][c]/M[c][c]; M[r] = [x - f*y for x, y in zip(M[r], M[c])]
        return d
    jac = abs(det(D))
    if jac == 0: return Fr(0)
    # x_j(lambda) = v0_j + sum_i lambda_i D[i][j]; the monomial prod_j x_j is a polynomial in lambda: represent as dict exps->coef
    polys = []
    for j in range(n):
        p = {tuple([0]*n): v0[j]}
        for i in range(n):
            e = [0]*n; e[i] = 1; p[tuple(e)] = p.get(tuple(e), Fr(0)) + D[i][j]
        polys.append(p)
    prod = {tuple([0]*n): Fr(1)}
    for p in polys:
        new = {}
        for e1, c1 in prod.items():
            for e2, c2 in p.items():
                e = tuple(a + b for a, b in zip(e1, e2)); new[e] = new.get(e, Fr(0)) + c1*c2
        prod = new
    total = Fr(0)
    for e, c in prod.items():
        if c == 0: continue
        num = 1
        for a in e: num *= factorial(a)
        total += c*Fr(num, factorial(n + sum(e)))
    return jac*total

def polytope_integral(rows, n):
    """exact integral of x_1...x_n over the polytope {a.x <= b}: fan triangulation from the centroid through the facets
    (n = 3: each facet is a convex polygon of the vertices on an active plane, fanned from one of its vertices)"""
    V = vertices(rows, n)
    if len(V) < n + 1: return Fr(0)
    rows = list({(tuple(a), b) for a, b in rows})
    cen = [sum(v[j] for v in V)/len(V) for j in range(n)]
    total = Fr(0)
    if n == 2:
        fc = [float(c) for c in cen]
        V.sort(key=lambda v: math.atan2(float(v[1]) - fc[1], float(v[0]) - fc[0]))
        for i in range(len(V)):
            total += simplex_monomial_integral([cen, V[i], V[(i+1) % len(V)]], n)
        return total
    assert n == 3
    for a, b in rows:
        face = [v for v in V if sum(ai*vi for ai, vi in zip(a, v)) == b]
        if len(face) < 3: continue
        # 2D coordinates in the plane for the cyclic order (floats, ordering only)
        af = [float(c) for c in a]; nrm = math.sqrt(sum(c*c for c in af)); af = [c/nrm for c in af]
        u = [1.0, 0.0, 0.0] if abs(af[0]) < 0.9 else [0.0, 1.0, 0.0]
        u = [u[i] - af[i]*sum(u[j]*af[j] for j in range(3)) for i in range(3)]; un = math.sqrt(sum(c*c for c in u)); u = [c/un for c in u]
        w = [af[1]*u[2] - af[2]*u[1], af[2]*u[0] - af[0]*u[2], af[0]*u[1] - af[1]*u[0]]
        fc = [sum(float(v[j]) for v in face)/len(face) for j in range(3)]
        def ang(v):
            d = [float(v[j]) - fc[j] for j in range(3)]
            return math.atan2(sum(d[j]*w[j] for j in range(3)), sum(d[j]*u[j] for j in range(3)))
        face.sort(key=ang)
        # drop collinear duplicates is unnecessary: a zero-volume tetrahedron contributes zero
        for i in range(1, len(face) - 1):
            total += simplex_monomial_integral([cen, face[0], face[i], face[i+1]], n)
    return total

def J(n):
    total = Fr(0)
    perms = list(itertools.permutations(range(n)))
    for eps in itertools.product((1, -1), repeat=n):
        valid = [p for p in perms if eps[p[0]] == 1]
        for p in valid:
            for q in valid:
                total += polytope_integral(constraints(n, p, eps) + constraints(n, q, eps), n)
    return total

if __name__ == "__main__":
    nmax = int(sys.argv[1]) if len(sys.argv) > 1 else 3
    for n in range(2, nmax + 1):
        Jn = J(n); print(f"J_{n} = {Jn} = {float(Jn):.9f}   kappa_{2*(n-1)} (distinct-prime part) = 2 J_{n}/{n}! = {2*Jn/factorial(n)} = {float(2*Jn/factorial(n)):.9f}", flush=True)
        if n == 3:
            K4 = Fr(11, 60); k4 = 2*Jn/6 + K4; k2 = Fr(5, 12)
            print(f"kappa_4 = J_3/3 + 11/60 = {k4} = {float(k4):.9f};   kurtosis kappa_4/kappa_2^2 = {k4/k2**2} = {float(k4/k2**2):.6f}", flush=True)
