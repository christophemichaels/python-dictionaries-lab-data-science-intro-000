"""
The continuum constants of the edge's spectral measure (paper, the local law after Proposition 9.38).

In the limit a -> inf the squared weights of the entries have the density dB(y) = y dy on (0, 2a) (Mertens), the closed walks through
the edge of length 2k+2 are grouped by the rational they land on, and for k+1 steps on powers of k+1 DISTINCT primes all the walks of a
group have the same weight, so that  m_{2k} |b|^2 = sum_groups (N_group c_1 ... c_{k+1})^2  up to the single-prime corrections, where
N_group is the number of orderings of the k+1 signed steps that stay inside the window.  Scaling the step lengths to units of 2a:

    kappa_{2k} := lim m_{2k}/(2a)^{2k} = 2 J_{k+1}/(k+1)!,   J_n = int_{[0,1]^n} x_1...x_n S_n(x) dx,   S_n = sum_{sign patterns} N(x, pattern)^2,

N(x, pattern) = number of orderings of the n signed steps (in = towards the far edge, out = back) whose first step is 'in' and whose
partial sums stay strictly in (0, 1), the landing point included.  J_2 = 5/12 exactly (Proposition 9.33(v)); the script checks it and
computes J_3, J_4, J_5 by Monte Carlo with standard errors, hence kappa_4, kappa_6, kappa_8 and the limiting kurtosis kappa_4/kappa_2^2.

Usage: python3 rh_walk_constants.py [samples] [nmax]     (default 2e6 samples, nmax = 5)
"""
import sys, itertools, numpy as np
S = int(float(sys.argv[1])) if len(sys.argv) > 1 else 2_000_000
nmax = int(sys.argv[2]) if len(sys.argv) > 2 else 5
rng = np.random.default_rng(1)

def Jn(n, S, chunk=200_000):
    perms = list(itertools.permutations(range(n)))
    pats = list(itertools.product((1, -1), repeat=n))          # +1 = in, -1 = out, attached to the elements
    tot = 0.0; tot2 = 0.0; cnt = 0
    while cnt < S:
        m = min(chunk, S - cnt); x = rng.random((m, n))
        Ssum = np.zeros(m)
        for pat in pats:
            eps = np.array(pat, dtype=float)
            N = np.zeros(m)
            for perm in perms:
                if pat[perm[0]] != 1: continue                     # first step must go in
                pos = np.zeros(m); ok = np.ones(m, dtype=bool)
                for i in perm:
                    pos = pos + eps[i]*x[:, i]
                    ok &= (pos > 0) & (pos < 1)
                N += ok
            Ssum += N*N
        w = np.prod(x, axis=1)*Ssum
        tot += w.sum(); tot2 += (w*w).sum(); cnt += m
    mean = tot/cnt; var = tot2/cnt - mean*mean
    return mean, np.sqrt(var/cnt)

import math

def valid_orderings(x, eps, perms):
    """number of orderings of the signed steps (eps_i x_i) with first step in and all partial sums in (0,1); x: (m,n) array"""
    m = x.shape[0]; N = np.zeros(m)
    for perm in perms:
        if eps[perm[0]] != 1: continue
        pos = np.zeros(m); ok = np.ones(m, dtype=bool)
        for i in perm:
            pos = pos + eps[i]*x[:, i]; ok &= (pos > 0) & (pos < 1)
        N += ok
    return N

def K4_exact():
    """the backtracking groups of m_4: walks {m+, m-, m'+} landing on the first-generation point a - log m'; group weight
    c_{m'} sum_m c_m^2 [2 1(d+d'<2a) + 1(d<d')] -> c_{m'} (2a)^2 [ (1-x')^2 + x'^2/2 ]; contribution 2 int_0^1 x [(1-x)^2 + x^2/2]^2 dx."""
    xs = np.linspace(0, 1, 200001); f = xs*((1 - xs)**2 + xs**2/2)**2
    return 2*np.trapezoid(f, xs)                       # = 2 * 11/120 = 11/60

def L6(S, zgrid=400, chunk=50_000):
    """the two-cancellation groups of m_6: {m+, m-, m'^e', m''^e''} landing on the rational m'^e' m''^e'';
    group weight c_m' c_m'' sum_m c_m^2 N(m) -> (2a)^2 int_0^1 z N(x,y,z) dz; contribution to kappa_6: 2 * (1/2) int int x y sum_eps [int z N dz]^2."""
    perms = list(itertools.permutations(range(4))); zs = (np.arange(zgrid) + 0.5)/zgrid
    tot = 0.0; tot2 = 0.0; cnt = 0
    while cnt < S:
        m = min(chunk, S - cnt); xy = rng.random((m, 2)); acc = np.zeros(m)
        for e1 in (1, -1):
            for e2 in (1, -1):
                inner = np.zeros(m)
                for z in zs:
                    X = np.column_stack([np.full(m, z), np.full(m, z), xy[:, 0], xy[:, 1]])
                    N = valid_orderings(X, np.array([1.0, -1.0, e1, e2]), perms)
                    inner += z*N/zgrid
                acc += inner**2
        w = xy[:, 0]*xy[:, 1]*acc
        tot += w.sum(); tot2 += (w*w).sum(); cnt += m
    mean = tot/cnt; var = tot2/cnt - mean*mean
    return 0.5*mean, 0.5*np.sqrt(var/cnt)

print(f"Monte Carlo with {S} samples per order")
print("  n   J_n (± s.e.)        kappa_{2(n-1)} = 2 J_n/n!      exact J_2 = 5/12 = 0.416667")
res = {}
for n in range(2, nmax + 1):
    J, se = Jn(n, S)
    k = 2*J/math.factorial(n); res[n] = k
    print(f"  {n}   {J:.6f} ± {se:.6f}     {k:.6f}", flush=True)
k2 = res[2]
if 3 in res: print(f"limiting kurtosis kappa_4/kappa_2^2 = {res[3]/k2**2:.4f}   (semicircle 2, Gaussian 3)")
if 4 in res: print(f"kappa_6/kappa_2^3 = {res[4]/k2**3:.4f}   (semicircle 5, Gaussian 15)")
if 5 in res: print(f"kappa_8/kappa_2^4 = {res[5]/k2**4:.4f}   (semicircle 14, Gaussian 105)")
K4 = K4_exact(); print(f"backtracking groups of m_4 (exact integral): 11/60 = {K4:.6f};  kappa_4 = J_3/3 + 11/60 = {res[3]/1 + K4 if 3 in res else float('nan'):.6f}  (17/72 + 11/60 = 151/360 = {151/360:.6f} if J_3 = 17/24)")
if 3 in res: print(f"limiting kurtosis with the backtracking groups: {(res[3] + K4)/k2**2:.4f}")
if 4 in res:
    l6, se6 = L6(min(S, 20000)); print(f"two-cancellation groups of m_6: 2 L_6 = {2*l6:.5f} ± {2*se6:.5f};  kappa_6 = J_4/12 + 2 L_6 = {res[4] + 2*l6:.5f};  kappa_6/kappa_2^3 = {(res[4] + 2*l6)/k2**3:.3f}")
print("exact moments for comparison (edge_moments.log): m_4/(2a)^4 = 0.243, 0.287, 0.308, 0.324 at a = 2, 2.5, 3, 3.5; m_2/(2a)^2 = 0.322, 0.348, 0.359, 0.368")
