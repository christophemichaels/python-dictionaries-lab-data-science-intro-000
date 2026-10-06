"""
The echo tree of a support (paper Lemma 8.13, Proposition 8.14, Corollary 8.15, Computations 8.16 and 8.16).

The entries d = log m < 2a copy the minimizer onto itself shifted by +-d. Starting from the edges +-a, the points reached by
walks with steps +-d that stay in [-a, a] form the echo set S; its interior points S° are the singular points of the
minimizer inside the window (the echoes of the edges), and the partition of the window by S is invariant under the entries.
The right edge sees the interior points through the weighted adjacency C (C_pq = c_d if |p - q| = d, c_d = Lambda(m)/sqrt m)
and the couplings b_p = c_d if a - p = d.  The effective symbol of the right edge is the Schur complement
    sigma_eff = sigma~ - b^T (sigma~ I - C)^{-1} b,
a rational function of the archimedean symbol sigma~ = sigma~_inf(t) (a finite continued fraction along a chain), and the
echo amplitudes are rho = (sigma~ I - C)^{-1} b.  The tree is finite as long as the walks are confined: this script computes
S for a range of supports, the size of the tree, the longest chain, the confinement of the walks from arbitrary interior
points, the spectra of C and of C with the edge adjoined (whose eigenvalues above min sigma~ are the real poles and zeros of
sigma_eff), and prints the tree at a = 0.6.

Usage: python3 rh_echo_tree.py [a_min a_max step]
"""
import math, sys, numpy as np

def vm(n):
    q = min(p for p in range(2, n + 1) if n % p == 0); m = n
    while m % q == 0: m //= q
    return math.log(q) if m == 1 else 0.0

def entries(a):
    """(d, c_d) for the prime powers m with log m < 2a."""
    out = []; m = 2
    while math.log(m) < 2*a:
        if vm(m): out.append((math.log(m), vm(m)/math.sqrt(m), m))
        m += 1
    return out

def walk_closure(starts, a, ds, cap=20000, tol=1e-9):
    """Points of [-a, a] reachable from `starts` by steps +-d (d in ds) staying in [-a, a]; None if more than `cap` points."""
    pts = []
    def find(x):
        for i, p in enumerate(pts):
            if abs(p - x) < tol: return i
        return -1
    todo = list(starts)
    for s in todo:
        if find(s) < 0: pts.append(s)
    k = 0
    while k < len(pts):
        p = pts[k]; k += 1
        for d in ds:
            for q in (p - d, p + d):
                if -a - tol <= q <= a + tol and find(q) < 0:
                    pts.append(q)
                    if len(pts) > cap: return None
    return sorted(pts)

def tree(a):
    """The echo set S, its interior S°, the matrices C (on S°), b, bbar, and the chain lengths (graph distance from a)."""
    ent = entries(a); ds = [d for d, _, _ in ent]; cd = {round(d, 9): c for d, c, _ in ent}
    S = walk_closure([a, -a], a, ds)
    if S is None: return None
    tol = 1e-9
    So = [p for p in S if abs(p - a) > tol and abs(p + a) > tol]
    n = len(So); C = np.zeros((n, n)); b = np.zeros(n); bb = np.zeros(n)
    for i, p in enumerate(So):
        for j, q in enumerate(So):
            key = round(abs(p - q), 9)
            for d in ds:
                if abs(abs(p - q) - d) < tol: C[i, j] = cd[round(d, 9)]
        for d in ds:
            if abs(a - p - d) < tol: b[i] = cd[round(d, 9)]
            if abs(p + a - d) < tol: bb[i] = cd[round(d, 9)]
    # graph distance from the right edge (chain length of the shortest chain), BFS on S
    dist = {a: 0}; todo = [a]
    while todo:
        p = todo.pop(0)
        for q in S:
            if q not in dist and any(abs(abs(p - q) - d) < tol for d in ds):
                dist[q] = dist[p] + 1; todo.append(q)
    # are the trees of the two edges connected inside the window? (then the cross-return bbar^T (sigma~ - C)^{-1} b is not zero)
    comp = {a}; todo = [a]
    while todo:
        p = todo.pop()
        for q in S:
            if q not in comp and any(abs(abs(p - q) - d) < tol for d in ds): comp.add(q); todo.append(q)
    return dict(S=S, So=So, C=C, b=b, bb=bb, dist=dist, ds=ds, cd=cd, ent=ent, connected=(-a in comp))

def confinement(a, ds, ngrid=400, cap=20000):
    """Largest orbit of a walk started anywhere in the window (None if some orbit exceeds the cap)."""
    worst = 0
    for p in np.linspace(-a, a, ngrid + 2)[1:-1]:
        O = walk_closure([p], a, ds, cap=cap)
        if O is None: return None
        worst = max(worst, len(O))
    return worst

def sigma_eff_rational(a):
    """Eigenvalues of C and of C with the edge adjoined (b as the coupling): sigma_eff = det(s - C_full)/det(s - C)."""
    T = tree(a)
    if T is None: return None
    C, b = T["C"], T["b"]; n = len(b)
    Cf = np.zeros((n + 1, n + 1)); Cf[:n, :n] = C; Cf[:n, n] = b; Cf[n, :n] = b
    return np.linalg.eigvalsh(C) if n else np.array([]), np.linalg.eigvalsh(Cf), T

if __name__ == "__main__":
    a0, a1, st = (float(sys.argv[1]), float(sys.argv[2]), float(sys.argv[3])) if len(sys.argv) > 3 else (0.35, 1.5, 0.025)
    print("the echo set S of the window: |S°| interior points, the longest shortest chain from the edge, the largest orbit of a walk from")
    print("any interior point (confinement), the top eigenvalue of C (poles of sigma_eff where sigma~ = eigenvalue) and of C with the edge")
    print("adjoined (zeros of sigma_eff), against min sigma~ ~ 1 and sigma~(T*) = 2a - lambda, sigma~(e T*) = 2a + 1 - lambda")
    print("    a    entries                    |S°|   chain   orbit    max eig C   max eig C+edge   |b|^2 = sum c_d^2   2a     trees of the two edges")
    for a in np.arange(a0, a1 + 1e-9, st):
        a = round(float(a), 6)
        ent = entries(a)
        if any(abs(math.log(m) - 2*a) < 1e-6 for _, _, m in ent) or any(abs(math.log(m) - 2*a) < 1e-6 for m in range(2, 40)):
            continue
        T = tree(a)
        if T is None:
            print(f"  {a:5.3f}  {[m for _,_,m in ent]!s:26s}  > cap: the walks are not confined"); continue
        n = len(T["So"]); chain = max(T["dist"].values())
        orb = confinement(a, T["ds"], ngrid=200, cap=5000)
        ev = np.linalg.eigvalsh(T["C"]) if n else np.array([0.0])
        Cf = np.zeros((n + 1, n + 1)); Cf[:n, :n] = T["C"]; Cf[:n, n] = T["b"]; Cf[n, :n] = T["b"]; evf = np.linalg.eigvalsh(Cf)
        b2 = float(np.sum(T["b"]**2))
        print(f"  {a:5.3f}  {[m for _,_,m in ent]!s:26s}  {n:4d}   {chain:4d}   {('%d' % orb) if orb is not None else '>5000':>6s}   {ev.max():9.4f}   {evf.max():13.4f}     {b2:8.4f}        {2*a:5.3f}   {'connected' if T['connected'] else 'disjoint'}")
    # the tree at a = 0.6
    a = 0.6; T = tree(a)
    print(f"\nthe tree at a = {a}: entries {[(m, round(d, 4), round(c, 4)) for d, c, m in T['ent']]}")
    print("  S° =", [round(p, 4) for p in T["So"]], " chain lengths", {round(p, 4): k for p, k in T["dist"].items()})
    print("  C =\n", np.array2string(T["C"], precision=4), "\n  b =", np.array2string(T["b"], precision=4), " bbar =", np.array2string(T["bb"], precision=4))
    ev = np.linalg.eigvalsh(T["C"]); n = len(T["b"]); Cf = np.zeros((n + 1, n + 1)); Cf[:n, :n] = T["C"]; Cf[:n, n] = T["b"]; Cf[n, :n] = T["b"]
    print("  eigenvalues of C:", np.array2string(ev, precision=4), " of C with the edge:", np.array2string(np.linalg.eigvalsh(Cf), precision=4))
    for s in (1.0, 1.5, 2.0, 2.8, 3.5, 4.2, 5.0):
        rho = np.linalg.solve(s*np.eye(n) - T["C"], T["b"]); R = float(T["b"] @ rho)
        print(f"  sigma~ = {s:4.1f}: R = {R:.4f} (first order |b|^2/sigma~ = {np.sum(T['b']**2)/s:.4f}), sigma_eff = {s - R:.4f}, rho = {np.array2string(rho, precision=4)}, "
              f"mean |E|^2 = {1 + float(rho @ rho):.4f}, M_chi - 1 = {(R + s*float(rho @ rho))/(s - R):.4f} (first order 2|b|^2/sigma~^2 = {2*np.sum(T['b']**2)/s**2:.4f})")
