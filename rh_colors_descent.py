"""Colors inside the subpower bound (COLORS_AND_DESCENT.md). The B-energy of any charge vector c is int_0^N (sum_{n>x} c_n)^2 dx
(min(a,b) = int 1_{x<a} 1_{x<b}), so every Gram entry below is an O(N) tail product. Sections: (A) the orthogonal boundary compression
of the directive (two-moment projection on the mesh x_{r+1} = x_r + floor(N^{-1/6} x_r^{2/3})): R = E_P + E_res, orthogonality, the
baseline row (N = 16384: m = 382, E_P = 1.588343108, E_res = 0.015097113), E_res as the within-interval variance of h, E_P to 1e7;
(B) the colored Gram G_kl = p_k^T B p_l by factor count omega, the parity direction, the small-prime refinement; (C) the prime-packet
identity c = prod (I - D_p) f_C and the subset Gram; (D) the six-term reconstruction (14) and its Gram; (E) the signed work E_P vs D_N;
(F) the three controls.  Usage: python3 rh_colors_descent.py"""
import math, numpy as np
from itertools import combinations
rng = np.random.default_rng(5)
X = 10**7
spf = np.zeros(X + 1, dtype=np.int64)
for p in range(2, int(X**0.5) + 1):
    if spf[p] == 0: spf[p*p::p][spf[p*p::p] == 0] = p
s = np.ones(X + 1, dtype=bool); s[:2] = False
for p in range(2, int(X**0.5) + 1):
    if s[p]: s[p*p::p] = False
primes = np.nonzero(s)[0]
MU = np.ones(X + 1, dtype=np.int8)
for p in primes: MU[p::p] *= -1; MU[p*p::p*p] = 0
MU[0] = 0
OM = np.zeros(X + 1, dtype=np.int8)
for p in primes: OM[p::p] += 1
def energy(c):
    """c^T B c for coefficient vector c indexed n = 1..N: int_0^N tail(x)^2 dx, tail constant on [x, x+1)."""
    tail = c[::-1].cumsum()[::-1]; return float(np.sum(tail**2))     # tail[x] = sum_{n > x} c_n for x = 0..N-1
def gram(c1, c2):
    t1 = c1[::-1].cumsum()[::-1]; t2 = c2[::-1].cumsum()[::-1]; return float(np.dot(t1, t2))
def mesh(N):
    xs = [0]; x = 0; d = N**(-1/6)
    while x < N: x = min(N, x + max(1, int(d * x**(2/3)))); xs.append(x)
    return xs
def project(c, xs, N):
    """the two-moment boundary projection p = T_N c on the mesh; coefficient at n = 1 kept (initial singleton)."""
    nn = np.arange(1, N + 1); p = np.zeros(N); p[0] = c[0]
    for a, b in zip(xs[1:-1], xs[2:]):
        seg = (nn > a) & (nn <= b); Q = float(np.sum(c[seg]*nn[seg])); H = float(np.sum(c[seg]))
        p[a-1] += (b*H - Q)/(b - a); p[b-1] += (Q - a*H)/(b - a)
    return p
print("== A. orthogonal boundary compression ==")
for N in (256, 1024, 4096, 16384):
    nn = np.arange(1, N + 1); c = MU[1:N+1]/nn; xs = mesh(N); p = project(c, xs, N); z = c - p
    R = energy(c); EP = energy(p); Eres = energy(z); orth = gram(z, p)
    M = np.cumsum(MU[1:N+1].astype(np.int64)); h = np.cumsum(c); hN = h[-1]
    u = lambda x: (M[x-1] if x >= 1 else 0) + x*(hN - (h[x-1] if x >= 1 else 0))
    EP_u = sum((u(b) - u(a))**2/(b - a) for a, b in zip(xs[:-1], xs[1:]))
    hx = np.concatenate([[0.0], h[:-1]])                        # h(x) for x = 0..N-1 (constant on [x, x+1)): h(x) = sum_{n<=x}
    var = sum(float(np.sum((hx[a:b] - np.mean(hx[a:b]))**2)) for a, b in zip(xs[:-1], xs[1:]))
    print(f"   N = {N:5d}: m = {len(xs)-1:4d} (<= 13 sqrt N + 2 = {int(13*math.sqrt(N))+2});  R = {R:.9f} = E_P {EP:.9f} + E_res {Eres:.9f};  z^T B p = {orth:.1e};  E_P from u = {EP_u:.9f};  E_res = within-interval variance of h: {var:.9f};  E_res <= 2^(-4/3) = {2**(-4/3):.4f}")
print("   E_P and E_res by the h-form at large N:")
for N in (10**5, 10**6, 10**7):
    c = MU[1:N+1]/np.arange(1, N + 1); xs = mesh(N); M = np.cumsum(MU[1:N+1].astype(np.int64)); h = np.cumsum(c); hN = h[-1]
    u = lambda x: (M[x-1] if x >= 1 else 0) + x*(hN - (h[x-1] if x >= 1 else 0))
    EP = sum((u(b) - u(a))**2/(b - a) for a, b in zip(xs[:-1], xs[1:])); R = energy(c); D = float(np.sum(MU[1:N+1].astype(float)**2/np.arange(1, N + 1)))
    print(f"   N = {N:8d}: m = {len(xs)-1:6d};  R = {R:.6f}  E_P = {EP:.6f}  E_res = {R-EP:.6f};  D_N = sum mu^2/n = {D:.4f};  E_P/D_N = {EP/D:.3f};  M(N)^2/N = {M[-1]**2/N:.4f}")
print("== B. the colored Gram by factor count omega, projected on the mesh ==")
for N in (1024, 4096, 16384):
    nn = np.arange(1, N + 1); xs = mesh(N); sq = MU[1:N+1] != 0; om = OM[1:N+1]; K = int(om[sq].max())
    pk = [project(np.where(sq & (om == k), 1.0/nn, 0.0), xs, N) for k in range(K + 1)]
    G = np.array([[gram(pk[k], pk[l]) for l in range(K + 1)] for k in range(K + 1)]); sgn = np.array([(-1)**k for k in range(K + 1)])
    EP = float(sgn @ G @ sgn); diag = float(np.trace(G)); mixed = EP - float(np.sum(np.diag(G)))
    Gfull = np.array([[gram(np.where(sq & (om == k), 1.0/nn, 0.0), np.where(sq & (om == l), 1.0/nn, 0.0)) for l in range(K + 1)] for k in range(K + 1)])
    print(f"   N = {N:5d}, colors k = 0..{K}: E_P = (-1)^(k+l) G_kl summed = {EP:.6f};  sum of same-color G_kk = {diag:.2f};  mixed-color signed part = {mixed:.2f};  angular average (8) = {diag:.2f};  unprojected: parity energy {float(sgn @ Gfull @ sgn):.6f}, diag {float(np.trace(Gfull)):.2f}")
    print("      G_kl (projected):\n" + "\n".join("        " + " ".join(f"{G[k,l]:11.3f}" for l in range(K + 1)) for k in range(K + 1)))
    if N == 16384:
        # refinement: classes by (2|n, 3|n) and omega outside {2,3}
        o23 = om - (nn % 2 == 0) - (nn % 3 == 0); labels = []; vecs = []
        for e2 in (0, 1):
            for e3 in (0, 1):
                for k in range(int(o23[sq].max()) + 1):
                    sel = sq & ((nn % 2 == 0) == e2) & ((nn % 3 == 0) == e3) & (o23 == k)
                    if sel.any(): labels.append((e2, e3, k)); vecs.append(project(np.where(sel, 1.0/nn, 0.0), xs, N))
        Gr = np.array([[gram(a, b) for b in vecs] for a in vecs]); sg = np.array([(-1)**(e2+e3+k) for e2, e3, k in labels])
        EPr = float(sg @ Gr @ sg); within = {}
        for (e2, e3, k), v in zip(labels, vecs): within.setdefault((e2, e3), []).append(((-1)**(e2+e3+k), v))
        parts = {key: energy(sum(s_*v for s_, v in lst)) for key, lst in within.items()}
        print(f"      refinement by (2|n, 3|n) x omega outside {{2,3}}: {len(labels)} classes, parity energy {EPr:.6f};  energy of each (2,3)-signature's own signed combination: " + ", ".join(f"({e2},{e3}): {parts[(e2,e3)]:.3f}" for (e2, e3) in sorted(parts)) + f";  sum of the four = {sum(parts.values()):.3f} (the rest is cross-signature)")
print("== C. prime packets: c = prod_{p in C} (I - D_p) f_C, subset Gram ==")
for N in (4096, 16384):
    nn = np.arange(1, N + 1); xs = mesh(N); c = MU[1:N+1]/nn
    for C in ((2, 3), (2, 3, 5)):
        qC = int(np.prod(C)); f = np.where(np.gcd(nn, qC) == 1, c, 0.0)
        def dil(v, q):
            out = np.zeros(N); K = N // q; out[q*np.arange(1, K + 1) - 1] = v[:K]/q; return out
        rec = np.zeros(N); terms = {}
        for k in range(len(C) + 1):
            for S in combinations(C, k):
                v = f.copy()
                for p in S: v = dil(v, p)
                terms[S] = v; rec += (-1)**k * v
        err = np.max(np.abs(rec - c)); keys = list(terms); P = {S: project(terms[S], xs, N) for S in keys}
        G = np.array([[gram(P[S], P[T]) for T in keys] for S in keys]); sg = np.array([(-1)**len(S) for S in keys])
        print(f"   N = {N:5d}, C = {C}: identity (9) error {err:.1e};  E(f_C) = {energy(f):.4f} (R(N) = {energy(c):.4f});  subset Gram diag: " + ", ".join(f"{S}:{G[i,i]:.3f}" for i, S in enumerate(keys)) + f";  signed sum = E_P = {float(sg @ G @ sg):.6f};  sum of diagonal {float(np.trace(G)):.3f}")
print("== D. the six-term reconstruction (14) at scale K = N^(1/6) ==")
def dconv(a, b, N):
    out = np.zeros(N + 1)
    for i in np.nonzero(a)[0]:
        if i == 0: continue
        out[i::i][:N//i] += a[i]*b[1:N//i+1]
    return out
for N in (4096, 15625, 16384):
    K = int(round(N**(1/6))) if N == 15625 else int(N**(1/6)); nn = np.arange(1, N + 1); xs = mesh(N)
    A = np.zeros(N + 1); A[1:K+1] = MU[1:K+1]; one = np.zeros(N + 1); one[1:] = 1.0
    Aj = [None, A.copy()]; 
    for j in range(2, 7): Aj.append(dconv(Aj[-1], A, N))
    onej = [np.zeros(N + 1), one.copy()]; onej[0][1] = 1.0
    for j in range(2, 6): onej.append(dconv(onej[-1], one, N))
    terms = []; rec = np.zeros(N + 1)
    for j in range(1, 7):
        t = dconv(Aj[j], onej[j-1], N) if j > 1 else Aj[1].copy(); coef = (-1)**(j-1)*math.comb(6, j); terms.append(coef*t); rec += coef*t
    err = np.max(np.abs(rec[1:N+1] - MU[1:N+1])); P = [project(t[1:N+1]/nn, xs, N) for t in terms]
    G = np.array([[gram(a, b) for b in P] for a in P])
    print(f"   N = {N:5d}, K = {K}: identity (14) error {err:.1e};  projected energies of the six signed terms: " + " ".join(f"{G[j,j]:.3e}" for j in range(6)) + f";  total E_P = {float(np.sum(G)):.6f};  sum of diagonal {float(np.trace(G)):.3e};  S(K) = max R(n), n <= K: {max(energy(MU[1:n+1]/np.arange(1, n+1)) for n in range(1, K+1)):.4f}")
print("== E. signed work: E_P against D_N, and the work W_N = sum mu(n) M(n-1)/n ==")
for N in (10**3, 10**4, 10**5, 10**6, 10**7):
    nn = np.arange(1, N + 1); c = MU[1:N+1]/nn; M = np.cumsum(MU[1:N+1].astype(np.int64)); D = float(np.sum(MU[1:N+1].astype(float)**2/nn))
    W = float(np.sum(MU[2:N+1]*M[:-1]/nn[1:])); R = energy(c)
    print(f"   N = {N:8d}: D_N = {D:.4f};  R = {R:.4f};  R - D_N = {R - D:.4f} = 2 W_N = {2*W:.4f};  trial E_P <= D_N holds with margin D_N - R = {D - R:.3f}")
print("== F. controls on the same mesh at N = 16384: E_P, same-color and mixed-color parts ==")
N = 16384; nn = np.arange(1, N + 1); xs = mesh(N); sq = MU[1:N+1] != 0; om = OM[1:N+1]; K = int(om[sq].max())
def shuffled(block):
    sgn = MU[1:N+1].astype(float).copy()
    for a in range(0, N, block):
        idx = np.nonzero(sq[a:a+block])[0] + a; sgn[idx] = rng.permutation(sgn[idx])
    return sgn
controls = {"Moebius": MU[1:N+1].astype(float), "positive": sq.astype(float), "random signs": sq*rng.choice([-1.0, 1.0], size=N), "shuffled in blocks of 64": shuffled(64), "shuffled in blocks of 1024": shuffled(1024)}
for name, sgn in controls.items():
    pk = [project(np.where(sq & (om == k), sgn/nn, 0.0), xs, N) for k in range(K + 1)]
    G = np.array([[gram(pk[k], pk[l]) for l in range(K + 1)] for k in range(K + 1)])
    print(f"   {name:27s}: E_P = {float(np.sum(G)):12.4f};  same-color {float(np.trace(G)):12.4f};  mixed-color {float(np.sum(G)) - float(np.trace(G)):+12.4f};  R = {energy(sgn/nn):.4f};  M(N) = {int(sgn.sum())}")
