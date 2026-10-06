"""The Rubik's-cube continuation (CUBE_MOVES.md).  Fixed-ratio increment Delta_6 R(N) = R(N) - R(floor(N/6)).
Section 1: exact rational checks of the moves (5)-(8), (9)-(13), (19) at small horizons.
Section 2: Delta_6 R(N) for every N <= 1e7 from the stable prefix sums; peaks, declared/held-out ranges, chain budgets.
Section 3: the complete fixed-packet difference (13) with every signed component; the newly-completed-parent band.
Section 4: the separated family of Section 5 of the directive: G_N, the proved allowance, (16)-(17), the calibration, increments.
Section 5: controls on one realization across all horizons.  Section 6: the returned finite budget.
Usage: python3 rh_cube_moves.py   (writes data/cube_moves.json)"""
import math, json, numpy as np
from fractions import Fraction as Fr
X = 10**7
s = np.ones(X + 1, dtype=bool); s[:2] = False
for p in range(2, int(X**0.5) + 1):
    if s[p]: s[p*p::p] = False
PR = np.nonzero(s)[0]
MU = np.ones(X + 1, dtype=np.int8)
for p in PR: MU[p::p] *= -1; MU[p*p::p*p] = 0
MU[0] = 0
def tails(c): return c[::-1].cumsum()[::-1]
def energy(c): t = tails(c); return float(np.dot(t, t))
def gram(a, b): return float(np.dot(tails(a), tails(b)))
def energies_from_counts(Mv):
    """R(n) for every n <= len(Mv)-1 from the prefix charges Mv[k] = sum_{j<=k} v_j j (Mv[0] = 0): the h-form (1)."""
    n = np.arange(len(Mv), dtype=float); n[0] = 1.0
    inner = np.cumsum(Mv.astype(float)**2/(n*(n + 1)))                 # sum_{k<=n} M(k)^2/(k(k+1))
    R = np.zeros(len(Mv)); R[1:] = inner[:-1] + Mv[1:].astype(float)**2/n[1:]
    return R
out = {}
# ---------------------------------------------------------------- 1. exact rational checks
print("== 1. exact rational checks of the moves, N <= 120 ==")
def Rfrac(c):  # c: dict n -> Fraction (coefficients), B = min
    T = {}; tot = Fr(0)
    for n in sorted(c, reverse=True): tot += c[n]; T[n] = tot
    ks = sorted(c); e = Fr(0); prev = 0
    for k in ks: e += (k - prev)*T[k]**2; prev = k       # sum_{j=1}^{N} tail(j)^2 with tail constant between support points
    return e
def gramfrac(a, b):
    sup = sorted(set(a) | set(b)); Ta = {}; Tb = {}; ta = Fr(0); tb = Fr(0)
    for n in sorted(sup, reverse=True): ta += a.get(n, Fr(0)); tb += b.get(n, Fr(0)); Ta[n] = ta; Tb[n] = tb
    e = Fr(0); prev = 0
    for k in sup: e += (k - prev)*Ta[k]*Tb[k]; prev = k
    return e
mu = [int(MU[n]) for n in range(121)]; Mx = np.cumsum(mu)
def cvec(N): return {n: Fr(mu[n], n) for n in range(1, N + 1) if mu[n]}
Rex = {0: Fr(0)}
for N in range(1, 121):
    Rex[N] = sum(Fr(int(Mx[k])**2, k*(k + 1)) for k in range(1, N)) + Fr(int(Mx[N])**2, N)
ok = True; worst12 = Fr(0); d23 = {}
for N in range(1, 121):
    cN = cvec(N); assert Rfrac(cN) == Rex[N]
    if N <= 60: assert sum(Fr(mu[a]*mu[b], max(a, b)) for a in range(1, N + 1) for b in range(1, N + 1)) == Rex[N]
    # (9): completed {2,3} family and D23
    par = [m for m in range(1, N//6 + 1) if mu[m] and math.gcd(m, 6) == 1]
    g = {}
    for m in par:
        for d in (1, 2, 3, 6): g[d*m] = g.get(d*m, Fr(0)) + Fr(mu[m]*mu[d], d*m)
    g = {n: v for n, v in g.items() if v}; t = {n: cN.get(n, Fr(0)) - g.get(n, Fr(0)) for n in set(cN) | set(g)}
    D = Fr(2, 3)*sum(Fr(1, m) for m in par); Eg = Rfrac(g) if g else Fr(0); I = Rex[N] - D
    assert (Eg - D) + 2*gramfrac(g, t) + Rfrac(t) == I; d23[N] = (D, I)
for N in range(6, 121):
    L = N//6; cN = cvec(N); cL = cvec(L); ML = int(Mx[L])
    bL = {L: Fr(ML, L)} if ML else {}
    rL = {n: cL.get(n, Fr(0)) - bL.get(n, Fr(0)) for n in set(cL) | set(bL)}; rL = {n: v for n, v in rL.items() if v}
    u = {n: cN[n] for n in cN if n > L}; w = dict(u); 
    if ML: w[L] = w.get(L, Fr(0)) + Fr(ML, L)
    assert sum(n*v for n, v in rL.items()) == 0                                  # zero physical total
    assert gramfrac(rL, u) == 0 and (not bL or gramfrac(rL, bL) == 0)           # B-orthogonality
    assert Rfrac(rL) + Rfrac(w) == Rex[N] and Rfrac(rL) + Fr(ML*ML, L) == Rex[L]   # (5)
    d6 = Rfrac(w) - Fr(ML*ML, L)                                                 # (6)
    d7 = 2*ML*sum(Fr(mu[n], n) for n in range(L + 1, N + 1)) + sum(Fr(mu[a]*mu[b], max(a, b)) for a in range(L + 1, N + 1) for b in range(L + 1, N + 1))
    d8 = sum(Fr(int(Mx[k])**2, k*(k + 1)) for k in range(L, N)) + Fr(int(Mx[N])**2, N) - Fr(ML*ML, L)
    assert d6 == d7 == d8 == Rex[N] - Rex[L]
    DN, IN = d23[N]; DL, IL = d23[L]; assert Rex[N] - Rex[L] == (DN - DL) + (IN - IL)   # (13)
    assert 0 <= DN - DL; worst12 = max(worst12, DN - DL)                                # (12)
    chain = []; n = N
    while n >= 6: chain.append(n); n //= 6
    assert Rex[N] == Rex[n] + sum(Rex[m] - Rex[m//6] for m in chain)                       # (19)
print(f"   (1) two forms of R agree for N <= 120 (double sum to 60); (5) orthogonality and both energy splits exact for 6 <= N <= 120;")
print(f"   (6) = (7) = (8) = R(N) - R(L) exactly; (9)-(10) components sum to I_23 exactly; (13) exact; (19) telescopes exactly;")
print(f"   (12): max over N <= 120 of D23(N) - D23(L) = {float(worst12):.6f} < 2/3 (1 + log 6) = {2/3*(1 + math.log(6)):.6f}")
# ---------------------------------------------------------------- 2. the increment to 1e7
print("== 2. Delta_6 R(N) for every 6 <= N <= 1e7 ==")
n_ = np.arange(X + 1, dtype=float); M = np.cumsum(MU.astype(np.int64)); R = energies_from_counts(M)
assert all(abs(R[N] - float(Rex[N])) < 1e-13 for N in range(1, 121))
S = np.maximum.accumulate(R)
idx = np.arange(X + 1); D6 = np.full(X + 1, np.nan); D6[6:] = R[6:] - R[idx[6:]//6]
mN = np.cumsum(MU.astype(float)/np.maximum(n_, 1))                                 # m(x) = sum mu(n)/n
def parts7(N):
    L = N//6; nn = np.arange(1, N + 1); c = MU[1:N+1]/nn; u = c.copy(); u[:L] = 0
    w = u.copy(); w[L-1] += M[L]/L
    cross = 2*M[L]*(mN[N] - mN[L]); tail = energy(u)
    return energy(w) - M[L]**2/L, cross + tail, cross, tail
for N in (16384, 10**5, 10**6, 10**7):
    d6, d7, cr, tl = parts7(N)
    print(f"   N = {N:8d}: (6) {d6:.9f}  (7) {d7:.9f}  (8) {D6[N]:.9f};  cross 2M(L)(m(N)-m(L)) = {cr:+.6f}  tail square ||u||^2 = {tl:.6f};  M(L) = {M[N//6]}, M(N) = {M[N]}")
valid = D6[6:]; Ns = idx[6:]
def stats(lo, hi):
    sel = (Ns >= lo) & (Ns <= hi); v = valid[sel]; k = int(np.argmax(v)); Nk = int(Ns[sel][k])
    return dict(lo=int(lo), hi=int(hi), max=float(v[k]), argmax=Nk, M_N=int(M[Nk]), M_L=int(M[Nk//6]), mean=float(v.mean()), frac_neg=float(np.mean(v < 0)), q99=float(np.quantile(v, 0.99)), min=float(v.min()))
rows2 = [stats(max(6, 10**k), min(X, 10**(k+1))) for k in range(0, 7)]
for r in rows2:
    print(f"   [{r['lo']:8d},{r['hi']:8d}]: max Delta_6 R = {r['max']:.6f} at N = {r['argmax']} (M(N) = {r['M_N']}, M(L) = {r['M_L']});  mean {r['mean']:+.5f}  q99 {r['q99']:.4f}  min {r['min']:+.4f}  negative fraction {r['frac_neg']:.3f}")
dec = stats(6, 10**6); held = stats(10**6 + 1, X)
print(f"   declared range [6, 1e6]: max {dec['max']:.6f} at {dec['argmax']};  held-out (1e6, 1e7]: max {held['max']:.6f} at {held['argmax']};  overall mean {valid.mean():+.6f} = slope {valid.mean()/math.log(6):.5f} per e-fold")
# peaks: top local maxima separated by a factor 1.5
order = np.argsort(-valid); peaks = []
for k in order:
    N = int(Ns[k])
    if all(not (N/1.5 < P < N*1.5) for P, _ in peaks): peaks.append((N, float(valid[k])))
    if len(peaks) == 8: break
print("   eight largest separated peaks (N, Delta_6 R, M(N), M(L), M(N)^2/N - M(L)^2/L, cross, tail):")
for N, v in peaks:
    _, _, cr, tl = parts7(N); print(f"      N = {N:8d}  Delta = {v:.6f}  M(N) = {M[N]:5d}  M(L) = {M[N//6]:5d}  terminal diff {M[N]**2/N - M[N//6]**2/(N//6):+.4f}  cross {cr:+.4f}  tail {tl:.4f}")
# completion horizons N = 6m, (m,6)=1
comp = (Ns % 6 == 0) & (np.gcd(Ns//6, 6) == 1)
print(f"   packet-completion horizons N = 6m, (m,6) = 1: count {int(comp.sum())}, max Delta {valid[comp].max():.6f}, mean {valid[comp].mean():+.5f};  other horizons: max {valid[~comp].max():.6f}, mean {valid[~comp].mean():+.5f}")
# chain budgets B+(N) = sum_j [Delta_6 R(N_j)]_+ along N -> N//6 -> ..., for every N
pos = np.where(np.isnan(D6), 0.0, np.maximum(D6, 0.0)); neg = np.where(np.isnan(D6), 0.0, np.maximum(-D6, 0.0))
Bp = np.zeros(X + 1); Bm = np.zeros(X + 1); lo = 6
while lo <= X:
    hi = min(X + 1, lo*6); par = idx[lo:hi]//6
    Bp[lo:hi] = Bp[par] + pos[lo:hi]; Bm[lo:hi] = Bm[par] + neg[lo:hi]; lo *= 6
# exact check of (19) for all N: R(N) = R(N_J) + Bp - Bm, with N_J the chain end (< 6)
end = idx.copy(); lo = 6
while lo <= X:
    hi = min(X + 1, lo*6); end[lo:hi] = end[idx[lo:hi]//6]; lo *= 6
assert np.max(np.abs(R - (R[end] + Bp - Bm))) < 1e-9
steps = np.zeros(X + 1, dtype=int); lo = 6
while lo <= X:
    hi = min(X + 1, lo*6); steps[lo:hi] = steps[idx[lo:hi]//6] + 1; lo *= 6
kB = int(np.argmax(Bp)); kBm = int(np.argmax(Bm)); ratio = Bp[6:]/np.maximum(R[6:], 1e-9)
print(f"   (19) holds for every N <= 1e7 (max error < 1e-9).  Chain budgets: max B+(N) = {Bp[kB]:.6f} at N = {kB} (R(N) = {R[kB]:.6f}, J = {steps[kB]} steps, B- = {Bm[kB]:.6f});")
print(f"      max B-(N) = {Bm[kBm]:.6f} at N = {kBm};  B+(1e7) = {Bp[X]:.6f} = R(1e7) {R[X]:.6f} - R(N_J) {R[end[X]]:.6f} + B- {Bm[X]:.6f} over J = {steps[X]} steps;  max B+/R = {ratio.max():.3f} at N = {int(np.argmax(ratio)) + 6}")
print(f"      max_N B+(N)/J(N) (mean positive cost per step along the worst chain) = {np.max(Bp[6:]/steps[6:]):.6f}")
out['increment'] = dict(ranges=rows2, declared=dec, held_out=held, peaks=peaks, Bplus_max=[kB, float(Bp[kB]), float(R[kB]), int(steps[kB])], Bplus_1e7=[float(Bp[X]), float(Bm[X]), float(R[X]), int(steps[X])])
# log-binned envelope for the figure
edges = np.unique(np.round(np.exp(np.linspace(math.log(6), math.log(X), 241))).astype(int)); env = []
for a, b in zip(edges[:-1], edges[1:]):
    v = D6[a:b+1]; env.append([int(a), int(b), float(np.nanmin(v)), float(np.nanmean(v)), float(np.nanmax(v))])
out['envelope'] = env
# ---------------------------------------------------------------- 3. the complete fixed-packet difference (13)
print("== 3. the fixed-packet difference (13), every signed component, N against L = floor(N/6) ==")
G6 = np.gcd(idx, 6); copsf = (MU != 0) & (G6 == 1)
Dcum = np.cumsum(np.where(copsf, 1.0, 0.0)/np.maximum(n_, 1))                      # sum_{m<=x,(m,6)=1} mu^2/m
def packet23(N):
    nn = idx[1:N+1]; c = MU[1:N+1]/nn; g = np.where(nn//G6[1:N+1] <= N//6, c, 0.0); t = c - g
    D = (2/3)*Dcum[N//6]; Eg = energy(g); gc = gram(g, c); RN = R[N]
    return dict(N=N, D=D, kern=Eg - D, cross=2*(gc - Eg), unf=RN - 2*gc + Eg, Eg=Eg, I=RN - D, R=RN)
grid = sorted(set([6**k for k in range(3, 9)] + [10**k for k in range(3, 8)] + [16384, 4*10**6]))
rows3 = []
for N in grid:
    a = packet23(N); b = packet23(N//6)
    d = {k: a[k] - b[k] for k in ('D', 'kern', 'cross', 'unf', 'I', 'R', 'Eg')}
    assert abs(a['kern'] + a['cross'] + a['unf'] - a['I']) < 1e-9 and abs(d['R'] - D6[N]) < 1e-9 and abs(d['D'] + d['I'] - d['R']) < 1e-9
    rows3.append(dict(N=N, L=N//6, at_N=a, at_L=b, diff=d))
    print(f"   N = {N:8d} L = {N//6:7d}: Delta D = {d['D']:+.6f}  Delta kern = {d['kern']:+.6f}  Delta cross = {d['cross']:+.6f}  Delta unf = {d['unf']:+.6f}  | Delta I = {d['I']:+.6f}  Delta_6 R = {d['R']:+.6f}  (completed-family energy Delta ||g||^2 = {d['Eg']:+.6f})")
print(f"   at N = 16384: components {rows3[[r['N'] for r in rows3].index(16384)]['at_N']['kern']:.6f} {rows3[[r['N'] for r in rows3].index(16384)]['at_N']['cross']:.6f} {rows3[[r['N'] for r in rows3].index(16384)]['at_N']['unf']:.6f} (previous note: -0.537728, -0.024996, 0.112661)")
def kern_band(N):
    """2 sum_{L/6 < n <= N/6} sum_{n/6 < m < n} mu(m)mu(n) K23(m,n) over (mn,6)=1: the newly completed parents' overlap."""
    L = N//6; tot = 0.0; per = []
    ms = idx[1:N//6 + 1]; ms = ms[copsf[ms]]
    for n in ms[ms > L//6]:
        mm = ms[(ms > n/6) & (ms < n)]
        if len(mm) == 0: continue
        K = np.zeros(len(mm))
        for d in (1, 2, 3, 6):
            for e in (1, 2, 3, 6): K += MU[d]*MU[e]/np.maximum(d*mm, e*n)
        v = 2*float(MU[n]*np.dot(MU[mm], K)); tot += v; per.append(v)
    return tot, np.array(per)
for N in (7776, 16384, 46656):
    tot, per = kern_band(N); d = packet23(N)['kern'] - packet23(N//6)['kern']
    print(f"   newly-completed-parent band at N = {N}: 2 sum over (L/6, N/6] x (n/6, n) = {tot:+.6f} = Delta kern {d:+.6f};  {len(per)} parents, per-parent contributions in [{per.min():+.4f}, {per.max():+.4f}], sum |.| = {np.abs(per).sum():.4f}")
out['packet23'] = rows3
# ---------------------------------------------------------------- 4. the separated family (directive Section 5)
print("== 4. the separated family: G_N, its allowance, (16)-(17), calibration at 4e6, increments ==")
Astar = 18/17*(17 + 12*math.sqrt(2))*(7 + 4*math.sqrt(3)); print(f"   A* = {Astar:.9f}")
def sep_family(N):
    l0 = int(N**(5/6))//6; l = int(PR[np.searchsorted(PR, l0, side='right')]); K = N//(6*l)
    nn = idx[1:N+1]; c = MU[1:N+1]/nn
    h = np.where(nn//G6[1:N+1] <= K, c, 0.0)                                  # (I-D2)(I-D3) f_N: packets of parents <= K
    Dl = np.zeros(N); Dl[l*idx[1:N//l+1] - 1] = h[:N//l]/l                    # D_l h
    g = h - Dl; t = c - g
    Eh = energy(h); Eg = energy(g); gc = gram(g, c)
    return dict(N=N, l=l, K=K, hDh=gram(h, Dl), Eh=Eh, G=Eg, G_formula=(1 + 1/l)*Eh, cross=2*(gc - Eg), unf=R[N] - 2*gc + Eg, W=R[N] - Eg, R=R[N], S16=float(S[int(N**(1/6))]), gc=gc, h=h, K6=6*K)
def V_of(x_int, K):  # V(x) = sum_{m<=K,(m,6)=1} mu(m) P23(x/m) on integers 0..x_int
    V = np.zeros(x_int + 1); xs = np.arange(x_int + 1)
    for m in range(1, K + 1):
        if MU[m] and math.gcd(m, 6) == 1:
            y = xs/m; V += MU[m]*((y >= 1).astype(float) - (y >= 2) - (y >= 3) + (y >= 6))
    return V
rows4 = []
for N in (46656, 10**6, 4*10**6, 10**7):
    f = sep_family(N); l, K = f['l'], f['K']
    # (16): <g,c> = sum over fine partition k/l, k = l .. 6Kl-1
    k = np.arange(l, 6*K*l); V = V_of(6*K, K); Vk = V[k//l]
    g16 = float(np.sum(Vk*(M[k//l] - M[k]/l)*l/(k*(k + 1.0))))
    assert np.max(np.abs(V - np.cumsum(idx[:6*K+1]*np.concatenate([[0.0], f['h'][:6*K]])))) < 1e-9
    # (17)
    kk = np.arange(1, N); Vfull = np.zeros(N + 1); Vfull[:6*K+1] = V
    first = float(np.sum((M[kk] - Vfull[kk] + Vfull[kk//l])**2/(kk*(kk + 1.0)))) + M[N]**2/N
    last = 2*float(np.sum(Vk*(M[k//l] - M[k]/l - (1 + 1/l)*Vk)*l/(k*(k + 1.0))))
    W17 = first + last
    print(f"   N = {N:8d}: l_N = {l}, K_N = {K} (parents {[m for m in range(1, K+1) if MU[m] and math.gcd(m,6)==1]});  <h, D_l h> = {f['hDh']:.1e};  G_N = {f['G']:.9f} = (1+1/l)||h||^2 {f['G_formula']:.9f};  A* S(N^(1/6)) = {Astar*f['S16']:.2f} (S({int(N**(1/6))}) = {f['S16']:.4f}, ratio {Astar*f['S16']/f['G']:.0f})")
    print(f"      R = {f['G']:.9f} {f['cross']:+.9f} + {f['unf']:.9f} = {f['G'] + f['cross'] + f['unf']:.9f} (R(N) = {f['R']:.9f});  (16): {g16:.9f} vs <g,c> {f['gc']:.9f};  (17): W = {W17:.9f} vs R - G = {f['W']:.9f}")
    rows4.append(dict(N=N, l=l, K=K, G=f['G'], cross=f['cross'], unf=f['unf'], W=f['W'], R=f['R'], allowance=Astar*f['S16']))
print("   increments (18) for the pairs (N, L): G_N - G_L, W_N - W_L, Delta_6 R; the families at N and L differ in prime and parent cutoff")
for N in (46656, 10**6, 4*10**6, 10**7):
    a = sep_family(N); b = sep_family(N//6)
    print(f"   N = {N:8d} (l = {a['l']}, K = {a['K']}), L = {N//6:7d} (l = {b['l']}, K = {b['K']}): G_N - G_L = {a['G'] - b['G']:+.6f}  W_N - W_L = {a['W'] - b['W']:+.6f}  sum {a['G'] - b['G'] + a['W'] - b['W']:+.6f} = Delta_6 R {D6[N]:+.6f};  Delta cross {a['cross'] - b['cross']:+.6f}  Delta unf {a['unf'] - b['unf']:+.6f}")
# how often does l_N change? count distinct l over N <= 1e7 at sample
ls = sorted(set(int(PR[np.searchsorted(PR, int(N**(5/6))//6, side='right')]) for N in range(216, X + 1, 997)))
print(f"   the dilation prime l_N takes about {len(ls)} distinct values on a sample of horizons to 1e7 (one per prime gap above N^(5/6)/6); K_N <= N^(1/6) <= {int(X**(1/6))}")
out['separated'] = rows4
# ---------------------------------------------------------------- 5. controls on a single realization across all horizons
print("== 5. controls (one realization, both horizons of every pair) ==")
rng = np.random.default_rng(2026)
def control_stats(name, Mv):
    Rv = energies_from_counts(Mv); Dv = Rv[6:] - Rv[idx[6:]//6]
    Bpv = np.zeros(X + 1); lo = 6
    posv = np.concatenate([np.zeros(6), np.maximum(Dv, 0)])
    while lo <= X:
        hi = min(X + 1, lo*6); Bpv[lo:hi] = Bpv[idx[lo:hi]//6] + posv[lo:hi]; lo *= 6
    envv = []
    for a, b in zip(edges[:-1], edges[1:]):
        v = Dv[a-6:b-5]; envv.append([int(a), int(b), float(v.min()), float(v.mean()), float(v.max())])
    k = int(np.argmax(Dv)); row = dict(name=name, envelope=envv, max=float(Dv[k]), argmax=int(k + 6), mean=float(Dv.mean()), frac_neg=float(np.mean(Dv < 0)), R_X=float(Rv[X]), Bp_X=float(Bpv[X]), Bp_max=float(Bpv.max()),
        by_decade=[float(np.max(Dv[max(6, 10**j) - 6: min(X, 10**(j+1)) - 5])) for j in range(1, 7)])
    print(f"   {name:34s}: max Delta_6 = {row['max']:.5g} at N = {row['argmax']}, mean {row['mean']:+.5g}, negative fraction {row['frac_neg']:.3f}, R(1e7) = {row['R_X']:.5g}, B+(1e7) = {row['Bp_X']:.5g};  decade maxima {['%.4g' % v for v in row['by_decade']]}")
    return row
rows5 = [control_stats("Mobius", M)]
rows5.append(control_stats("all positive (M = k)", idx.astype(np.int64)))
rows5.append(control_stats("squarefree positive (M = Q)", np.cumsum((MU != 0).astype(np.int64))))
eps = np.where(MU != 0, rng.choice([-1, 1], size=X + 1), 0).astype(np.int64); eps[0] = 0
rows5.append(control_stats("independent signs, squarefree support", np.cumsum(eps)))
sh = MU.astype(np.int64).copy(); lo = 1
while lo <= X:                                                                 # shuffle nonzero values inside each block (6^j, 6^(j+1)]
    hi = min(X, lo*6); blk = idx[lo+1:hi+1]; nz = blk[MU[blk] != 0]; sh[nz] = rng.permutation(MU[nz]); lo *= 6
rows5.append(control_stats("block-conditioned shuffle (hexadic)", np.cumsum(sh)))
print(f"   all-positive check: Delta_6 R = 2(N - L) - (H_N - H_L) at N = 1e7: {2*(X - X//6) - (np.sum(1/n_[X//6+1:X+1])):.6f} vs computed {rows5[1]['max'] if rows5[1]['argmax']==X else float('nan'):.6f}")
out['controls'] = rows5
# ---------------------------------------------------------------- 6. the returned finite budget
print("== 6. return (19)-(21) with the measured budget: a finite statement on [6, 1e7] ==")
Cmax = float(np.nanmax(D6)); J = steps[X]
print(f"   measured C = max_{{6<=N<=1e7}} Delta_6 R(N) = {Cmax:.6f};  (19) with this C: R(1e7) <= max_{{n<6}} R(n) + J C = {R[:6].max():.4f} + {J} x {Cmax:.4f} = {R[:6].max() + J*Cmax:.4f}, actual R(1e7) = {R[X]:.6f}, actual B+(1e7) = {Bp[X]:.6f}")
print(f"   per-step mean along the chain to 1e7: {(R[X] - R[end[X]])/J:+.5f};  S(1e7) = {S[X]:.6f} at N = {int(np.argmax(R)) if R.max()==S[X] else -1}")
print(f"   the sixth-root factor on this range: S(N^(1/6)) <= S({int(X**(1/6))}) = {S[int(X**(1/6))]:.6f}, constant to within the range, so (A) on [6, 1e7] cannot separate gamma from C_eps")
out['budget'] = dict(Cmax=Cmax, J=int(J), R_X=float(R[X]), S_X=float(S[X]), S16=float(S[int(X**(1/6))]))
print("== 7. sharpness of the local bound M(N)^2/N - M(L)^2/L <= Delta_6 R(N) <= (1 + log(N/L)) max_{L<=k<=N} M(k)^2/k ==")
M2k = M[1:].astype(float)**2/n_[1:]                                               # M(k)^2/k at index k-1
BS = 1024; nb = (X + BS - 1)//BS; bmax = np.array([M2k[i*BS:(i+1)*BS].max() for i in range(nb)])
def winmax(L, N):                                                                  # max_{L<=k<=N} M(k)^2/k
    a, b = L - 1, N - 1; ia, ib = a//BS, b//BS
    if ia == ib: return float(M2k[a:b+1].max())
    v = max(float(M2k[a:(ia+1)*BS].max()), float(M2k[ib*BS:b+1].max()))
    if ib > ia + 1: v = max(v, float(bmax[ia+1:ib].max()))
    return v
sample = sorted(set(list(range(100, X + 1, 491)) + [N for N, _ in peaks if N >= 100] + [42968, 300551, 1065673, 6481601, 10**7]))
rat_up = []; rat_lo = []
for N in sample:
    L = N//6; wm = winmax(L, N); d = D6[N]
    rat_up.append((d/((1 + math.log(N/L))*wm), N)); rat_lo.append((d - (M[N]**2/N - M[L]**2/L), N))
ru = max(rat_up); print(f"   over {len(sample)} horizons N >= 100: max Delta_6 R / [(1 + log(N/L)) max M^2/k] = {ru[0]:.4f} at N = {ru[1]};  mean ratio {np.mean([r for r, _ in rat_up]):.4f};  min of Delta_6 R - (M(N)^2/N - M(L)^2/L) = {min(rat_lo)[0]:+.2e} (>= 0 always: the difference is the block sum of M(k)^2/(k(k+1)) over L <= k < N)")
for N in (42968, 300551, 1065673, 6481601, 10**7):
    L = N//6; print(f"      N = {N:8d}: Delta_6 R = {D6[N]:.6f}, lower M(N)^2/N - M(L)^2/L = {M[N]**2/N - M[L]**2/L:+.6f}, window max M^2/k = {winmax(L, N):.6f}, upper = {(1 + math.log(N/L))*winmax(L, N):.6f}")
S100 = float(S[100]); C100 = float(np.nanmax(D6[101:])); J100 = int(np.sum([(X // 6**j) > 100 for j in range(0, 12)]))
print(f"   two-piece finite return: S(100) = {S100:.6f}, max_{{100 < N <= 1e7}} Delta_6 R = {C100:.6f} at N = {int(np.nanargmax(D6[101:])) + 101}; steps of the chain from 1e7 above 100: {J100}; bound R(1e7) <= {S100:.4f} + {J100} x {C100:.4f} = {S100 + J100*C100:.4f} (actual {R[X]:.6f}); for every N <= 1e7: R(N) <= S(100) + ceil(log(N/100)/log 6) x {C100:.4f}")
out['sharpness'] = dict(max_ratio=ru, S100=S100, C100=C100, J100=J100)
json.dump(out, open('data/cube_moves.json', 'w'), indent=1, default=float)
print("wrote data/cube_moves.json")
