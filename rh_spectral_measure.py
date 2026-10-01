"""
The band the edge sees: the spectral measure of the confined echo lattice at the edge's coupling (paper Proposition 8.23, Computation 8.24).

The edge of the window couples to the echo lattice only through b (b_p = c_d at the first-generation points a - d), so the return
R(s) = b^T (s - C)^{-1} b = int dmu_b(nu)/(s - nu) is the Stieltjes transform of the spectral measure mu_b of the adjacency C at b,
with total mass |b|^2.  The effective symbol sigma_eff = s - R(s) at the frozen height s = sigma~_inf(t) is real where s lies above the
support of mu_b, and acquires the imaginary part pi mu_b'(s) (the edge radiating into the lattice) inside it.  The band the edge sees
is therefore the support of mu_b, and the height from which the sharp form can hold is where the weight of mu_b above s is small,
not where s exceeds the spectral radius rho(C), which mu_b may not reach.

The lattice is truncated to the first N points by chain length (rh_lattice_band.py); the moments <b, C^k b> of mu_b for k up to twice the
largest chain length are then exact, so the truncation resolves mu_b like a Lanczos iteration of that many steps.  For each support the
script prints rho_N, the top of the support of mu_b, the weight of mu_b above s = sigma~_inf(kappa T*) for several kappa, the smoothed
density of mu_b there (the radiative width), and R(s), sigma_eff(s)/s and the first-order 1 - |b|^2/s^2 at those heights.

Usage: python3 rh_spectral_measure.py [N] [a ...]      (default N = 6000, a = 1.0 1.25 1.5 1.75 2.0)
"""
import math, sys, numpy as np
from rh_lattice_band import lattice

def spectral_measure(a, N):
    xs, gen, adj, steps = lattice(a, N)
    n = len(xs); C = np.zeros((n, n))
    for j, nb in enumerate(adj):
        for i, c in nb: C[j, i] = c
    interior = np.array([k for k in range(n) if k not in (0, 1)])          # the points other than the two edges
    b = C[0, interior]                                                        # the coupling of the right edge to the interior
    Ci = C[np.ix_(interior, interior)]
    nu, V = np.linalg.eigh(Ci)
    w = (V.T @ b)**2                                                           # weights of mu_b at the eigenvalues
    return dict(n=n, gmax=int(gen.max()), nu=nu, w=w, b2=float(b @ b), rho=float(np.abs(nu).max()), steps=steps)

def lanczos_measure(a, N, m):
    """Gauss quadrature of mu_b with m nodes from m Lanczos steps started at b on the lattice truncated to N points (full
    reorthogonalization); exact moments of mu_b to order 2m-1 as long as the walks of that length stay inside the truncation."""
    xs, gen, adj, steps = lattice(a, N)
    n = len(xs); interior = np.array([k for k in range(n) if k not in (0, 1)]); pos = -np.ones(n, dtype=int); pos[interior] = np.arange(len(interior))
    rows = [[] for _ in interior]; cols = []; vals = []
    for j, nb in enumerate(adj):
        if j in (0, 1): continue
        for i, c in nb:
            if i not in (0, 1): rows[pos[j]].append((pos[i], c))
    b = np.zeros(len(interior))
    for i, c in adj[0]:
        if i not in (0, 1): b[pos[i]] = c
    def Cv(v):
        w = np.zeros_like(v)
        for j, nb in enumerate(rows):
            for i, c in nb: w[j] += c*v[i]
        return w
    # sparse matvec via index arrays (faster than the loop)
    I = np.concatenate([[j]*len(nb) for j, nb in enumerate(rows)]).astype(int); J = np.array([i for nb in rows for i, _ in nb], dtype=int); W = np.array([c for nb in rows for _, c in nb])
    def Cv(v):
        return np.bincount(I, weights=W*v[J], minlength=len(v))
    b2 = float(b @ b); q = b/math.sqrt(b2); Q = [q]; alpha = []; beta = []
    for k in range(m):
        w = Cv(Q[-1]); al = float(Q[-1] @ w); w = w - al*Q[-1] - (beta[-1]*Q[-2] if beta else 0)
        for qq in Q: w -= float(qq @ w)*qq                                      # full reorthogonalization
        be = float(np.linalg.norm(w)); alpha.append(al)
        if be < 1e-12 or k == m - 1: break
        beta.append(be); Q.append(w/be)
    T = np.diag(alpha) + np.diag(beta, 1) + np.diag(beta, -1)
    nu, V = np.linalg.eigh(T); w = b2*V[0, :]**2
    return dict(n=n, gmax=int(gen.max()), nu=nu, w=w, b2=b2, steps=steps, m=len(alpha))

def moments_direct(a):
    """The first two moments of mu_b from the entries alone: m_1 |b|^2 = <b, C b> = sum over ordered pairs of entries (m, m') with
    m/m' a prime power of c_m c_m' c_{m/m'} (closed walks of length 3 through the edge), and |b|^2 = sum c_m^2."""
    ent = {}
    m = 2
    while math.log(m) < 2*a:
        q = min(p for p in range(2, m + 1) if m % p == 0); r = m
        while r % q == 0: r //= q
        if r == 1: ent[m] = math.log(q)/math.sqrt(m)
        m += 1
    b2 = sum(c*c for c in ent.values())
    m1 = 2*sum(ent[m]*ent[mp]*ent[m//mp] for m in ent for mp in ent if m != mp and m % mp == 0 and (m//mp) in ent)   # both orders of the pair
    return b2, m1/b2, sorted(ent)

if __name__ == "__main__" and len(sys.argv) > 1 and sys.argv[1] == "moments":
    N = int(sys.argv[2]); m = int(sys.argv[3]); supports = [float(x) for x in sys.argv[4:]] or [1.0, 1.25, 1.5, 1.75, 2.0, 2.25, 2.5]
    print("the moments m_k = <b, C^k b>/|b|^2 of the edge's spectral measure: closed walks of length k+2 through the edge, i.e. multiplicative")
    print("relations m_1^{e_1} ... m_{k+2}^{e_{k+2}} = 1 among prime powers below e^{2a} whose partial products stay in the window; m_k^{1/k} -> top of supp mu_b")
    for a in supports:
        M = lanczos_measure(a, N, m); nu, w, b2 = M["nu"], M["w"], M["b2"]
        b2d, m1d, ent = moments_direct(a)
        mk = [float(np.sum(w*nu**k))/b2 for k in range(1, 15)]
        roots = [abs(v)**(1.0/k) if v != 0 else 0.0 for k, v in enumerate(mk, 1)]
        print(f"a = {a}: entries {len(ent)}, |b|^2 = {b2:.4f} (direct {b2d:.4f}), m_1 = {mk[0]:.4f} (direct, pairs m/m' a prime power: {m1d:.4f}), depth {M['gmax']} (moments exact to order {2*M['gmax']-2});"
              f" 2a = {2*a:.2f}, 2a + log 5 = {2*a + math.log(5):.2f}", flush=True)
        print("   k:      " + "  ".join(f"{k:6d}" for k in range(1, 15)))
        print("   m_k:    " + "  ".join(f"{v:6.3f}" if abs(v) < 1e3 else f"{v:6.0f}" for v in mk))
        print("   m_k^1/k:" + "  ".join(f"{v:6.3f}" for v in roots), flush=True)
    sys.exit(0)

if __name__ == "__main__" and len(sys.argv) > 1 and sys.argv[1] == "lanczos":
    N = int(sys.argv[2]); m = int(sys.argv[3]); supports = [float(x) for x in sys.argv[4:]] or [1.5, 1.75, 2.0, 2.25, 2.5]
    for a in supports:
        Ts = 2*math.pi*math.exp(2*a); M = lanczos_measure(a, N, m); nu, w, b2 = M["nu"], M["w"], M["b2"]
        cum = np.cumsum(w[::-1])[::-1]
        print(f"a = {a}: N = {M['n']}, chains <= {M['gmax']}, Lanczos steps {M['m']}, |b|^2 = {b2:.4f}, top node of the Gauss quadrature of mu_b = {nu.max():.4f} "
              f"(weight {w[-1]/b2:.2e}); sigma~_inf(T*) = {2*a:.2f}", flush=True)
        print("   kappa   s=sigma~(kT*)   weight above s     R(s)     sigma_eff/s   1-|b|^2/s^2")
        for k in (1, 2, 3, 5, 7, 10, 20):
            s = 2*a + math.log(k); above = float(w[nu >= s].sum())/b2
            Rr = float(np.sum(w/(s - nu))) if above == 0 else float('nan')
            print(f"   {k:4d}   {s:9.3f}        {above:9.4f}      {Rr:8.4f}   {1 - Rr/s:10.4f}   {1 - b2/s**2:10.4f}", flush=True)
        out = []
        for qt in (0.5, 0.2, 0.1, 0.05, 0.01, 0.001):
            idx = np.nonzero(cum/b2 <= qt)[0]; sq = nu[idx[0]] if len(idx) else nu[-1]
            out.append(f"{100*qt:g}%: s = {sq:.3f} (t = {math.exp(sq - 2*a):.2f} T*)")
        print("   weight of mu_b above s drops below  " + ";  ".join(out), flush=True)
    sys.exit(0)

if __name__ == "__main__":
    args = sys.argv[1:]; N = int(args[0]) if args and args[0].isdigit() else 6000
    supports = [float(x) for x in (args[1:] if args and args[0].isdigit() else args)] or [1.0, 1.25, 1.5, 1.75, 2.0]
    kappas = (1, 2, 3, 5, 7, 10, 20)
    for a in supports:
        Ts = 2*math.pi*math.exp(2*a); M = spectral_measure(a, N); nu, w, b2 = M["nu"], M["w"], M["b2"]
        order = np.argsort(nu); nu, w = nu[order], w[order]; cum = np.cumsum(w[::-1])[::-1]                  # weight at or above nu
        top = nu[w > 1e-8*b2].max()
        # mean degree lower bound and the first moments
        m1 = float(nu @ w)/b2; m2 = float((nu*nu) @ w)/b2
        print(f"a = {a}: N = {M['n']}, chains <= {M['gmax']}, |b|^2 = {b2:.4f}, rho_N = {M['rho']:.4f}, top of supp mu_b = {top:.4f}, "
              f"mean of mu_b/|b|^2 = {m1:.4f}, second moment {m2:.4f} (= (b^T C^2 b)/|b|^2); sigma~_inf(T*) = {2*a:.2f}", flush=True)
        print("   kappa   s=sigma~(kT*)   weight above s   density at s (eta=0.1)   R(s)      sigma_eff/s   1-|b|^2/s^2   Im R / pi (radiative width)")
        for k in kappas:
            s = 2*a + math.log(k)                                                 # sigma~_inf(kappa T*) ~ log(kappa T*/2pi) (lambda ~ 0)
            above = float(w[nu >= s].sum())/b2
            eta = 0.1; dens = float(np.sum(w*eta/np.pi/((s - nu)**2 + eta**2)))/b2
            Rc = complex(np.sum(w/(s - nu + 1j*eta)))                          # limiting absorption with eta
            Rr = float(np.sum(w/(s - nu))) if above == 0 else Rc.real
            print(f"   {k:4d}   {s:9.3f}        {above:9.4f}        {dens:12.4f}          {Rr:8.4f}   {1 - Rr/s:10.4f}   {1 - b2/s**2:10.4f}      {-Rc.imag/math.pi:10.4f}", flush=True)
        # the spectral weight quantiles: the height at which the weight above s is 10%, 5%, 1% of |b|^2
        out = []
        for q in (0.5, 0.2, 0.1, 0.05, 0.01):
            idx = np.nonzero(cum/b2 <= q)[0]
            sq = nu[idx[0]] if len(idx) else nu[-1]
            out.append(f"{int(100*q)}%: s = {sq:.3f} (t = {math.exp(sq - 2*a):.2f} T*)")
        print("   weight of mu_b above s drops below  " + ";  ".join(out), flush=True)
