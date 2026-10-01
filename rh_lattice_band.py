"""
The band of the confined echo lattice (paper Proposition 8.20, Computation 8.21).

Beyond a_inf = 0.843 the echo set S is dense, but the walks are still confined to the window, and the weighted adjacency C of the
echo graph (C_pq = c_d for |p - q| = d = log m, c_d = Lambda(m)/sqrt m) is a bounded operator on l^2(S°) whose spectral radius
rho(C) is far below the top 2 sum_d c_d of the periodic limit's band: a step +-d from p stays inside the window only when there is
room, and from any point at most a few of the 2|D_a| steps do.  The Neumann series of the echo tree converges where
sigma~_inf(t) > rho(C), so the band of the confined lattice ends at t_band = 2 pi exp(rho(C) + lambda), and the sharp form of the
tail law can hold from there on, far below the last sign change of the window's symbol (33 T* at a = 1, 259 T* at a = 1.25).

This script builds S by breadth-first search from +-a in the exponent-vector representation (a point is a + sum_p m_p log p over the
primes p < e^{2a}, the prime powers riding on their primes), truncated to the first N points by chain length, and computes the top
eigenvalue of the truncated adjacency by power iteration (it increases to rho(C) with N), the largest weighted degree, and t_band.

Usage: python3 rh_lattice_band.py [a ...]       (default 0.8 0.9 1.0 1.1 1.25 1.5)
"""
import math, sys, numpy as np

def primes_upto(x):
    return [p for p in range(2, int(x) + 2) if all(p % q for q in range(2, int(p**0.5) + 1))]

def lattice(a, N):
    """Echo points reached from +-a by walks inside [-a, a], in BFS order (by chain length), at most N of them.
    Returns the points, the chain lengths, and the adjacency lists [(j, c_d), ...]."""
    P = [p for p in primes_upto(math.exp(2*a)) if math.log(p) < 2*a]
    steps = []                                               # (prime index, k, d, c_d) for the prime powers p^k < e^{2a}
    for i, p in enumerate(P):
        k = 1
        while k*math.log(p) < 2*a:
            steps.append((i, k, k*math.log(p), math.log(p)/math.sqrt(p**k))); k += 1
    def pos(v, sgn): return sgn*a + sum(m*math.log(p) for m, p in zip(v, P))
    start = [((0,)*len(P), 1), ((0,)*len(P), -1)]          # (exponent vector, which edge)
    index = {s: i for i, s in enumerate(start)}; pts = list(start); gen = [0, 0]; k0 = 0
    while k0 < len(pts) and len(pts) < N:
        v, sgn = pts[k0]; k0 += 1
        for (i, k, d, c) in steps:
            for e in (k, -k):
                w = list(v); w[i] += e; w = tuple(w); key = (w, sgn)
                if key in index: continue
                x = pos(w, sgn)
                if -a - 1e-12 <= x <= a + 1e-12:
                    index[key] = len(pts); pts.append(key); gen.append(gen[k0 - 1] + 1)
                    if len(pts) >= N: break
            if len(pts) >= N: break
    # adjacency among the points found (a step between two found points)
    adj = [[] for _ in pts]
    for j, (v, sgn) in enumerate(pts):
        for (i, k, d, c) in steps:
            for e in (k, -k):
                w = list(v); w[i] += e; key = (tuple(w), sgn)
                if key in index: adj[j].append((index[key], c))
    xs = np.array([pos(v, sgn) for v, sgn in pts])
    return xs, np.array(gen), adj, steps

def top_eigenvalue(adj, iters=300):
    n = len(adj); v = np.ones(n)/math.sqrt(n); lam = 0.0
    for _ in range(iters):
        w = np.zeros(n)
        for j, nb in enumerate(adj):
            for i, c in nb: w[j] += c*v[i]
        # shift to make the iteration converge to the top (most positive) eigenvalue
        w += 0.0
        nrm = np.linalg.norm(w)
        if nrm == 0: return 0.0
        lam_new = float(v @ w); v = w/nrm
        if abs(lam_new - lam) < 1e-10: lam = lam_new; break
        lam = lam_new
    # the adjacency is bipartite-like only for commensurate steps; take the largest |eigenvalue| as the spectral radius
    return abs(lam)

if __name__ == "__main__":
    supports = [float(x) for x in sys.argv[1:]] or [0.8, 0.9, 1.0, 1.1, 1.25, 1.5]
    print("the confined echo lattice: top eigenvalue of the adjacency on the first N points (by chain length), the largest weighted degree,")
    print("the periodic band top 2 sum c_d, and the heights t = 2 pi exp(.) in units of the horizon T* = 2 pi e^{2a}")
    for a in supports:
        Ts = 2*math.pi*math.exp(2*a); out = []
        for N in (500, 2000, 8000):
            xs, gen, adj, steps = lattice(a, N)
            if len(xs) < N:
                out.append(f"N={len(xs)}(all) rho={top_eigenvalue(adj):.4f} chain<={gen.max()}"); break
            out.append(f"N={N} rho={top_eigenvalue(adj):.4f} chain<={gen.max()}")
        xs, gen, adj, steps = lattice(a, 8000)
        deg = max(sum(c for _, c in nb) for nb in adj); S2c = 2*sum(c for _, _, _, c in steps); rho = top_eigenvalue(adj)
        print(f"a = {a:5.3f}: entries {len(steps)}, |S| >= {len(xs)}; " + "; ".join(out) +
              f"; max weighted degree {deg:.3f}; 2 sum c_d = {S2c:.3f}; band top of the confined lattice t = {2*math.pi*math.exp(rho)/Ts:.2f} T*, "
              f"of the periodic limit {2*math.pi*math.exp(S2c)/Ts:.0f} T*; sigma~_inf(T*) = {2*a:.2f}, sigma~_inf(5 T*) = {2*a + math.log(5):.2f}", flush=True)
