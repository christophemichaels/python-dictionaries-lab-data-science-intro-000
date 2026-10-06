"""Decay inside the growth (PARITY_DECAY.md). Exact mesh (integer sixth root), the projection (3); the parity defect delta_N = E_P/A_N
with the two terms of (10) from the slope representation; the Gram split (11)-(12); the completed-packet remainder for C = {2,3}:
components of (17)-(18), the kernel (19), and the scan of I_C(N) to 1e7; the energy records, first passages and the retained set (21);
the profile curvature (30); the conditional-shuffle expectation (31); the controls.  Usage: python3 rh_parity_decay.py"""
import math, numpy as np
rng = np.random.default_rng(9)
X = 10**7
s = np.ones(X + 1, dtype=bool); s[:2] = False
for p in range(2, int(X**0.5) + 1):
    if s[p]: s[p*p::p] = False
primes = np.nonzero(s)[0]
MU = np.ones(X + 1, dtype=np.int8)
for p in primes: MU[p::p] *= -1; MU[p*p::p*p] = 0
MU[0] = 0
OM = np.zeros(X + 1, dtype=np.int8)
for p in primes: OM[p::p] += 1
def iroot6(v):
    r = int(round(v**(1/6)))
    while r**6 > v: r -= 1
    while (r+1)**6 <= v: r += 1
    return r
def mesh(N):
    xs = [0]; x = 0
    while x < N: x = min(N, x + max(1, iroot6(x**4 // N))); xs.append(x)
    return xs
def tails(c): return c[::-1].cumsum()[::-1]
def energy(c): t = tails(c); return float(np.dot(t, t))
def gram(a, b): return float(np.dot(tails(a), tails(b)))
def project(c, xs, N):
    nn = np.arange(1, N + 1); p = np.zeros(N); p[0] = c[0]
    for a, b in zip(xs[1:-1], xs[2:]):
        seg = slice(a, b); w = c[seg]; n = nn[seg]
        p[a-1] += float(np.sum(w*(b - n)))/(b - a); p[b-1] += float(np.sum(w*(n - a)))/(b - a)
    return p
def slopes(c, xs, N):
    """Euclidean slope representation v_r = Delta u / sqrt(Delta x) of the projected source: ||v||^2 = E_P."""
    nn = np.arange(1, N + 1); M = np.cumsum(c*nn); h = np.cumsum(c); hN = h[-1]
    u = lambda x: (M[x-1] if x >= 1 else 0.0) + x*(hN - (h[x-1] if x >= 1 else 0.0))
    return np.array([(u(b) - u(a))/math.sqrt(b - a) for a, b in zip(xs[:-1], xs[1:])])
print("== 1. reconciliation at N = 16384 (expected m = 382, R = 1.6034402209313143, E_P = 1.5883431079039165, E_res = 0.015097113027358108) ==")
N = 16384; nn = np.arange(1, N + 1); c = MU[1:N+1]/nn; xs = mesh(N); p = project(c, xs, N)
print(f"   m = {len(xs)-1};  R = {energy(c):.16f};  E_P = {energy(p):.16f} (slopes: {float(np.sum(slopes(c, xs, N)**2)):.16f});  E_res = {energy(c - p):.16f};  z^T B p = {gram(c - p, p):.1e};  first mesh points {xs[:12]}")
print("== 2. the parity defect across horizons: E_P, A_N, delta_N, N delta_N, the two terms of (10), the terminal ==")
Ns = sorted(set([1024, 16384, 65536, 262144, 10**6] + [int(2**(k/2)) for k in range(20, 47)]))
rows = []
for N in Ns:
    nn = np.arange(1, N + 1); sq = MU[1:N+1] != 0; om = OM[1:N+1]; xs = mesh(N)
    ce = np.where(sq & (om % 2 == 0), 1.0/nn, 0.0); co = np.where(sq & (om % 2 == 1), 1.0/nn, 0.0)
    ve = slopes(ce, xs, N); vo = slopes(co, xs, N); a = math.sqrt(float(np.dot(ve, ve))); b = math.sqrt(float(np.dot(vo, vo)))
    EP = float(np.sum((ve - vo)**2)); A = a*a + b*b; cos = float(np.dot(ve, vo))/(a*b); amp = (a - b)**2; dirn = a*b*float(np.sum((ve/a - vo/b)**2))
    M = int(np.sum(MU[1:N+1])); rows.append((N, EP, A, EP/A, N*EP/A, amp, dirn, M*M/N))
    print(f"   N = {N:8d}: E_P = {EP:.9f}  A_N = {A:.6f}  delta = {EP/A:.6e}  N delta = {N*EP/A:.6f}  (a-b)^2 = {amp:.6f}  2ab(1-cos) = {dirn:.6f}  M^2/N = {M*M/N:.6f}  A_N/N = {A/N:.4f}")
print("   local exponents of N delta_N between consecutive horizons: " + " ".join(f"{math.log(rows[i+1][4]/rows[i][4])/math.log(rows[i+1][0]/rows[i][0]):+.3f}" for i in range(len(rows)-1)))
print("== 3. the Gram split (11)-(12) at N = 16384 ==")
N = 16384; nn = np.arange(1, N + 1); sq = MU[1:N+1] != 0; om = OM[1:N+1]; xs = mesh(N); K = int(om[sq].max())
V = np.array([slopes(np.where(sq & (om == k), 1.0/nn, 0.0), xs, N) for k in range(K + 1)]).T       # V[r, k]
q = np.array([int(np.sum(sq & (om == k))) for k in range(K + 1)], dtype=float); w = np.array([math.sqrt(b - a) for a, b in zip(xs[:-1], xs[1:])])
G = V.T @ V; C = (V - np.outer(w, q)/N).T @ (V - np.outer(w, q)/N); sgn = np.array([(-1)**k for k in range(K + 1)])
print(f"   G - qq^T/N - C: max error {np.max(np.abs(G - np.outer(q, q)/N - C)):.1e};  min eigenvalue of C = {np.linalg.eigvalsh(C).min():.2e};  E_P = {float(sgn @ G @ sgn):.9f} = M^2/N {float((sgn @ q)**2/N):.9f} + s^T C s {float(sgn @ C @ sgn):.9f};  same-colour {np.trace(G):.9f}, signed mixed {float(sgn @ G @ sgn) - np.trace(G):.9f}")
print("== 4. the completed-packet remainder, C = {2,3} ==")
def packet_parts(N, C=(2, 3)):
    q = int(np.prod(C)); divs = [d for d in range(1, q + 1) if q % d == 0]; nn = np.arange(1, N + 1); c = MU[1:N+1]/nn
    par = np.arange(1, N//q + 1); par = par[np.gcd(par, q) == 1]; g = np.zeros(N)
    for d in divs: g[d*par - 1] += MU[par]*MU[d]/(d*par)
    t = c - g; kappa = sum(MU[d]*MU[e]/max(d, e) for d in divs for e in divs)
    diag = kappa*float(np.sum(MU[par].astype(float)**2/par)); Eg = energy(g); overlap = Eg - diag; cross = 2*gram(g, t); unf = energy(t)
    return kappa, diag, overlap, cross, unf, Eg, energy(c)
kappa, diag, overlap, cross, unf, Eg, R = packet_parts(16384)
print(f"   N = 16384: kappa = {kappa:.6f};  diagonal {diag:.6f}  local overlap {overlap:.6f}  completed-unfinished cross {cross:.6f}  unfinished square {unf:.6f};  sum {diag+overlap+cross+unf:.6f} = R {R:.6f};  I_C = R - diagonal = {R - diag:.6f}")
# kernel (19): K_C(m,n) = (mn)^(-1/2) k_C(log(n/m)), k_C(xi) = sum mu(d)mu(e)(de)^(-1/2) exp(-|xi + log(e/d)|/2); compact support
def KC(m, n, C=(2, 3)):
    q = int(np.prod(C)); divs = [d for d in range(1, q + 1) if q % d == 0]; return sum(MU[d]*MU[e]/max(d*m, e*n) for d in divs for e in divs)
def kC(xi, C=(2, 3)):
    q = int(np.prod(C)); divs = [d for d in range(1, q + 1) if q % d == 0]; return sum(MU[d]*MU[e]/math.sqrt(d*e)*math.exp(-abs(xi + math.log(e/d))/2) for d in divs for e in divs)
print("   kernel (19): " + "; ".join(f"K({m},{n}) = {KC(m, n):.6f} vs (mn)^-1/2 k(log n/m) = {kC(math.log(n/m))/math.sqrt(m*n):.6f}" for m, n in ((5, 7), (7, 25), (5, 29), (5, 31), (1, 6), (1, 7))) + f";  k_C(xi) at xi = log 6: {kC(math.log(6)):.1e}, at 1.9: {kC(1.9):.1e}, at 0: {kC(0):.6f} (= kappa)")
# scan of I_C(N) = R(N) - kappa * sum_{m <= N/6, (m,6)=1} mu(m)^2/m, to 1e7
n_ = np.arange(X + 1, dtype=float); M = np.cumsum(MU.astype(np.int64))
inc = np.zeros(X + 1); inc[1:] = (MU[1:].astype(float)**2 + 2*MU[1:]*np.concatenate([[0], M[1:-1]]))/n_[1:]; Rn = np.cumsum(inc)
cop = (MU != 0) & (np.gcd(np.arange(X + 1), 6) == 1); Dcop = np.cumsum(np.where(cop, 1.0, 0.0)/np.maximum(n_, 1))
IC = Rn[1:] - (2/3)*Dcop[(np.arange(1, X + 1))//6]
pos = np.nonzero(IC > 0)[0] + 1
print(f"   I_C(N) for N <= 1e7: positive at {len(pos)} values of N, the largest being N = {pos.max() if len(pos) else None};  max I_C = {IC.max():.6f} at N = {int(np.argmax(IC))+1};  I_C at 1e3..1e7: " + ", ".join(f"{IC[10**k-1]:.4f}" for k in range(3, 8)))
print("== 5. energy records, first passages, the retained set (21), through 1e7 ==")
rec = [1]; best = Rn[1]
for n in range(2, X + 1):
    if Rn[n] > best + 1e-12: rec.append(n); best = Rn[n]
absM = np.abs(M); tau = {}
for n in range(1, X + 1):
    h = int(absM[n])
    if h not in tau: tau[h] = n
fp = all(tau.get(int(absM[n])) == n for n in rec)
S = np.maximum.accumulate(Rn)
retained = []
for n in rec[1:]:
    h = int(absM[n]); a = tau[h-1]; thr = math.sqrt(n)/(1 + S[int(n**(1/6))])**2.5
    if n - a < thr: retained.append((n, h, a, thr))
print(f"   strict energy records: {len(rec)} (first ten {rec[:10]}; last {rec[-1]}, height {int(absM[rec[-1]])});  every record is a first passage of |M|: {fp};  retained noninitial records by (21): {len(retained)}" + (f" -> {retained[:5]}" if retained else " (empty)"))
print(f"   gaps N - tau_(h-1) at the last five records: {[ (n, n - tau[int(absM[n])-1], round(math.sqrt(n)/(1 + S[int(n**(1/6))])**2.5, 2)) for n in rec[-5:] ]}")
print("== 6. profile curvature (30) and width ==")
for N in (16384, 10**5, 10**6):
    nn = np.arange(1, N + 1); sq = MU[1:N+1] != 0; om = OM[1:N+1]; xs = mesh(N); K = int(om[sq].max())
    Vk = [slopes(np.where(sq & (om == k), 1.0/nn, 0.0), xs, N) for k in range(K + 1)]
    pN = sum((-1)**k*Vk[k] for k in range(K + 1)); H1 = sum(k*(-1)**k*Vk[k] for k in range(K + 1)); H2 = sum(k*k*(-1)**k*Vk[k] for k in range(K + 1))
    E = lambda t: float(np.sum(sum(math.cos(k*t)*Vk[k] for k in range(K + 1))**2) + np.sum(sum(math.sin(k*t)*Vk[k] for k in range(K + 1))**2))
    curv = 2*float(np.dot(H1, H1)) - 2*float(np.dot(pN, H2)); Epi = E(math.pi); dE = (E(math.pi) - E(math.pi - 1e-4))/1e-4
    w = next((x for x in np.linspace(0, 0.5, 5001) if E(math.pi - x) >= 2*Epi), None)
    print(f"   N = {N:7d}: E(pi) = {Epi:.6f};  E'(pi) numerically {dE:.2e};  E''(pi) = 2||H1||^2 - 2<p,H2> = {curv:.4f} (numerical {(E(math.pi) - 2*E(math.pi-1e-3) + E(math.pi-2e-3))/1e-6:.4f});  width to 2E(pi): pi - t = {w:.4f};  curvature/E(pi) = {curv/Epi:.2f}")
print("== 7. the conditional-shuffle expectation (31) at N = 16384 (expected 2.035225552 = 1.473291335 + 0.561934217) ==")
N = 16384; nn = np.arange(1, N + 1); sq = MU[1:N+1] != 0; xs = mesh(N); c = MU[1:N+1]/nn
abar = np.zeros(N); V = 0.0
for j in range(0, 15):
    lo, hi = 2**j, min(2**(j+1), N + 1); idx = np.arange(lo, hi)[sq[lo-1:hi-1]]; k = len(idx)
    if k == 0: continue
    sj = float(np.sum(MU[idx])); abar[idx-1] = sj/k
    if k >= 2:
        tr = 0.0
        for n in idx:
            # <T e_n, T e_n>_B: T e_n = alpha e_a + beta e_b on its interval
            r = np.searchsorted(xs, n); a, b = xs[r-1], xs[r]
            if a == 0: tr += 1.0/n**2; continue
            al, be = (b - n)/(b - a), (n - a)/(b - a); tr += (al*al*a + 2*al*be*a + be*be*b)/n**2
        one = energy(project(np.where(np.isin(nn, idx), 1.0/nn, 0.0), xs, N))
        V += (1 - (sj/k)**2)/(k - 1)*(k*tr - one)
EPbar = energy(project(abar/nn, xs, N)); mc = np.mean([energy(project(np.concatenate([rng.permutation(MU[1:N+1][sq]) if False else MU[1:N+1]])[0:0] if False else 0, xs, N)) for _ in range(0)]) if False else None
def shuffle_once():
    sgn = MU[1:N+1].astype(float).copy()
    for j in range(0, 15):
        lo, hi = 2**j, min(2**(j+1), N + 1); idx = np.arange(lo, hi)[sq[lo-1:hi-1]] - 1
        if len(idx) > 1: sgn[idx] = rng.permutation(sgn[idx])
    return energy(project(sgn/nn, xs, N))
mcs = [shuffle_once() for _ in range(200)]
print(f"   E_P(abar) = {EPbar:.9f};  V_N = {V:.9f};  expectation {EPbar + V:.9f};  Monte Carlo over 200 shuffles: {np.mean(mcs):.4f} +/- {np.std(mcs)/math.sqrt(200):.4f};  actual E_P = {energy(project(c, xs, N)):.9f};  D_N = {float(np.sum(sq/nn)):.4f}")
print("== 8. controls at N = 16384 ==")
aC = np.where(sq, (-1.0)**((nn % 2 == 0).astype(int) + (nn % 3 == 0).astype(int)), 0.0); alpha = (6/math.pi**2)*(1/3)*(2/4)
for name, src in (("all positive", np.ones(N)), ("squarefree positive", sq.astype(float)), ("random signs on squarefree", sq*rng.choice([-1.0, 1.0], size=N)), ("conditional shuffle (one draw)", None), ("a_{2,3}(n) = mu^2 (-1)^{#{p in {2,3}: p | n}}", aC), ("Moebius", MU[1:N+1].astype(float))):
    e = shuffle_once() if src is None else energy(project(src/nn, xs, N))
    print(f"   {name:42s}: E_P = {e:12.4f}" + (f"   (2 alpha_C^2 N = {2*alpha**2*N:.1f})" if name.startswith("a_") else ""))
