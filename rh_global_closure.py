"""Global closure and Connes geometry (GLOBAL_CLOSURE.md).  Checks, each beside its derivation in the note:
1 checkpoints; 2 the local bridge (19)-(22) with the delta jumps and the finite-interval identity (21); 3 the scalar circle
model (23)-(25), the energy scaling under translation, the two-period density; 4 the packet action (46) in potential
coordinates, the two-prime assembly with every cross term, the all-prime product; 5 the orbit-sum identities and the fold
onto C_p (divisor on the unfinished packets, constants from complete packets, convergence iff M(N) = 0, slope denominators);
6 the spectral identity (6) and the compactified metric (26); 7 the chain-refined ledger (37) and the return (42);
8 the sharp diagonal 4/5; 9 controls.   Usage: python3 rh_global_closure.py  (writes data/global_closure.json)"""
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
idx = np.arange(X + 1); n_ = idx.astype(float); n_[0] = 1.0
M = np.cumsum(MU.astype(np.int64))
inner = np.cumsum(M.astype(float)**2/(n_*(n_ + 1))); R = np.zeros(X + 1); R[1:] = inner[:-1] + M[1:].astype(float)**2/n_[1:]
hN = np.cumsum(MU.astype(float)/n_)                                  # h(x) = sum_{n<=x} mu(n)/n
def tails(c): return c[::-1].cumsum()[::-1]
def energy(c): t = tails(c); return float(np.dot(t, t))
def gram(a, b): return float(np.dot(tails(a), tails(b)))
def iroot6(v):
    r = int(round(v**(1/6)))
    while r**6 > v: r -= 1
    while (r+1)**6 <= v: r += 1
    return r
def mesh(N):
    xs = [0]; x = 0
    while x < N: x = min(N, x + max(1, iroot6(x**4 // N))); xs.append(x)
    return xs
def project(c, xs, N):
    nn = np.arange(1, N + 1); p = np.zeros(N); p[0] = c[0]
    for a, b in zip(xs[1:-1], xs[2:]):
        seg = slice(a, b); w = c[seg]; n = nn[seg]
        p[a-1] += float(np.sum(w*(b - n)))/(b - a); p[b-1] += float(np.sum(w*(n - a)))/(b - a)
    return p
out = {}
print("== 1. checkpoints ==")
for N in (1000, 10**4, 10**5, 10**6, 10**7, 16384): print(f"   R({N}) = {R[N]:.10f}")
N = 16384; c = MU[1:N+1]/np.arange(1, N + 1); xs = mesh(N); p = project(c, xs, N)
print(f"   original mesh at 16384: m = {len(xs)-1}, E_P = {energy(p):.16f}, E_res = {energy(c - p):.16f}, 2^(-4/3) = {2**(-4/3):.6f}")
# ---------------------------------------------------------------- 2. the local bridge
print("== 2. the local bridge (19)-(22) ==")
def V(N, x):      # Mobius potential V_N(x) = sum_{n<=N} mu(n)/n min(x,n), vectorised in x
    x = np.asarray(x, dtype=float); k = np.minimum(np.floor(x), N).astype(int)          # breakpoints at integers
    return (M[k] + x*(hN[N] - hN[k]))*(x > 0)                                              # V = M(k) + x (h(N) - h(k)) on [k, k+1)
def phi_green(N, u):   # phi_N(u) = sum mu(n)/sqrt(n) exp(-|u - log n|/2)   (Green sum, (4))
    u = np.asarray(u, dtype=float); ns = np.arange(1, N + 1); w = MU[1:N+1]/np.sqrt(ns)
    return np.array([float(np.sum(w*np.exp(-np.abs(uu - np.log(ns))/2))) for uu in u])
def dphi_green(N, u):
    u = np.asarray(u, dtype=float); ns = np.arange(1, N + 1); w = MU[1:N+1]/np.sqrt(ns)
    return np.array([float(np.sum(w*(-0.5)*np.sign(uu - np.log(ns))*np.exp(-np.abs(uu - np.log(ns))/2))) for uu in u])
N = 50; u = np.linspace(-3, 6, 901); e19 = np.max(np.abs(V(N, np.exp(u)) - np.exp(u/2)*phi_green(N, u)))
e7 = np.max(np.abs(phi_green(N, u) - (np.exp(-u/2)*M[np.minimum(np.floor(np.exp(u)), N).astype(int)] + np.exp(u/2)*(hN[N] - hN[np.minimum(np.floor(np.exp(u)), N).astype(int)]))))
jumps = [(n, n*((hN[N] - hN[n]) - (hN[N] - hN[n-1]))) for n in range(1, N + 1)]            # f' = x V'(x): jump at log n is n (V'(n+) - V'(n-))
print(f"   N = 50: max |V_N(e^u) - e^(u/2) phi_N(u)| = {e19:.1e} (19);  max |Green sum - two-tail form (7)| = {e7:.1e};  f'-jumps at log n equal -mu(n): max error {max(abs(j + MU[n]) for n, j in jumps):.1e} (20)")
# (22) exact in rationals: R(N) = int |V'|^2 dx = h(N)^2 + sum_{k<N} (h(N)-h(k))^2 for N = 50
hq = [Fr(0)]
for n in range(1, 51): hq.append(hq[-1] + Fr(int(MU[n]), n))
Rq = sum(Fr(int(M[k])**2, k*(k + 1)) for k in range(1, 50)) + Fr(int(M[50])**2, 50)
print(f"   (22) exact at N = 50: int |V'|^2 dx = h(N)^2 + sum (h(N)-h(k))^2 = {float(hq[50]**2 + sum((hq[50]-hq[k])**2 for k in range(1, 50))):.12f} = R(50) = {float(Rq):.12f}: {hq[50]**2 + sum((hq[50]-hq[k])**2 for k in range(1, 50)) == Rq}")
# (21) on [a,b]: LHS by Gauss-Legendre on the Green-sum field per cell, RHS from the potential
gx, gw = np.polynomial.legendre.leggauss(24)
def lhs21(N, a, b):
    bps = sorted(set([a, b] + [math.log(n) for n in range(1, N + 1) if a < math.log(n) < b])); tot = 0.0
    for lo, hi in zip(bps[:-1], bps[1:]):
        uu = (hi - lo)/2*gx + (hi + lo)/2; tot += (hi - lo)/2*float(np.sum(gw*(dphi_green(N, uu)**2 + phi_green(N, uu)**2/4)))
    return tot
def rhs21(N, a, b):
    xa, xb = math.exp(a), math.exp(b); bps = sorted(set([xa, xb] + [n for n in range(1, N + 1) if xa < n < xb])); tot = 0.0
    for lo, hi in zip(bps[:-1], bps[1:]):
        k = min(int(math.floor((lo + hi)/2)), N); tot += (hN[N] - hN[k])**2*(hi - lo)            # int |V'|^2 dx on the cell
    return tot - 0.5*(math.exp(-b)*float(V(N, xb))**2 - math.exp(-a)*float(V(N, xa))**2)
for a, b in ((-1.0, 2.3), (0.7, 3.5), (-4.0, 8.0)):
    print(f"   (21) on [{a}, {b}] at N = 50: LHS {lhs21(50, a, b):.10f}  RHS {rhs21(50, a, b):.10f}")
print(f"   (21) on [-30, 30] (full line to 1e-13): LHS {lhs21(50, -30, 30):.10f} vs R(50) = {float(Rq):.10f}")
print(f"   exterior tails at N = 50: h(N)^2/2 = {hN[50]**2/2:.6f}, M(N)^2/(2N) = {M[50]**2/100:.6f}")
# ---------------------------------------------------------------- 3. the scalar circle model
print("== 3. the scalar circle model (23)-(25); translation scaling; two periods ==")
rng = np.random.default_rng(5); p = 5; l = math.log(p); q = p**-0.5
co = rng.normal(size=(2, 7))
def fper(u, c):  # period-l trigonometric polynomial and derivative
    ks = np.arange(1, 4); kap = 2*np.pi*ks/l
    f = c[0] + np.sum(c[1:4]*np.cos(kap*u) + c[4:7]*np.sin(kap*u)); df = np.sum(-c[1:4]*kap*np.sin(kap*u) + c[4:7]*kap*np.cos(kap*u)); return f, df
def psi(u, c): f, df = fper(u, c); return math.exp(-u/2)*f, math.exp(-u/2)*(df - f/2)
ps0, dps0 = psi(0.0, co[0]); ch0, dch0 = psi(0.0, co[1]); psl, dpsl = psi(l, co[0]); chl, dchl = psi(l, co[1])
print(f"   (23): psi(l)/psi(0) = {psl/ps0:.9f}, psi'(l)/psi'(0) = {dpsl/dps0:.9f}, p^(-1/2) = {q:.9f}")
print(f"   (24): [psi chi' - psi' chi]_0^l = {psl*dchl - dpsl*chl - (ps0*dch0 - dps0*ch0):.9f};  (q^2 - 1)(psi(0)chi'(0) - psi'(0)chi(0)) = {(q*q - 1)*(ps0*dch0 - dps0*ch0):.9f}")
kap = 2*np.pi*2/l; z = complex(-0.5, kap); print(f"   (25): for psi_k = e^(-u/2) e^(i kappa u), -psi'' + psi/4 = (-(z^2) + 1/4) psi with z = -1/2 + i kappa: coefficient {-(z*z) + 0.25:.9f} = kappa^2 + i kappa = {complex(kap*kap, kap):.9f}")
uu = np.linspace(0, l, 20001); fv = np.array([fper(x, co[0])[1] for x in uu]); e0 = np.trapezoid(np.exp(-uu)*fv**2, uu); e1 = np.trapezoid(np.exp(-(uu + l))*fv**2, uu)
print(f"   raw psi-energy int e^(-u)|f'|^2 over [0,l] = {e0:.9f}, over [l,2l] = {e1:.9f}, ratio {e1/e0:.9f} = 1/p;  half-line [0,inf) total = cell x p/(p-1) = {e0*p/(p-1):.9f}")
# two periods: continued-fraction convergents of log3/log2
a_, b_ = math.log(3), math.log(2); r = a_/b_; cf = []; x = r
for _ in range(14): ai = math.floor(x); cf.append(ai); x = 1/(x - ai)
h0, h1, k0, k1 = 1, cf[0], 0, 1; conv = [(h1, k1)]
for ai in cf[1:]: h0, h1 = h1, ai*h1 + h0; k0, k1 = k1, ai*k1 + k0; conv.append((h1, k1))
print("   two periods: |a log 2 - b log 3| for convergents (a, b) = " + ", ".join(f"({a},{b}): {abs(a*b_ - b*a_):.2e}" for a, b in conv[3:9]) + ";  the set Z log 2 + Z log 3 is dense (log 3/log 2 irrational since 2^a != 3^b)")
# ---------------------------------------------------------------- 4. the packet action in potential coordinates; two primes with every cross term; the all-prime product
print("== 4. packet action (46) in potential coordinates; two-prime assembly; all-prime product ==")
def VQ(Q, N, x):   # potential of the squarefree integers n <= N with primes in Q
    qs = [n for n in range(1, N + 1) if MU[n] and all(r in Q for r in prime_factors(n))]
    x = np.asarray(x, dtype=float); return sum(MU[n]/n*np.minimum(x, n) for n in qs)
def prime_factors(n):
    f = []; d = 2
    while d*d <= n:
        if n % d == 0: f.append(d)
        while n % d == 0: n //= d
        d += 1
    if n > 1: f.append(n)
    return f
xg = np.linspace(0, 1200, 4801); Q = (2, 3); pnew = 5; N = 1000
lhs = VQ((2, 3, 5), N, xg); rhs = VQ(Q, N, xg) - VQ(Q, N//pnew, xg/pnew)
print(f"   (46) in potential coordinates, Q = {{2,3}}, p = 5, N = 1000: max |V_(Q+p,N)(x) - [V_(Q,N)(x) - V_(Q,N/p)(x/p)]| = {np.max(np.abs(lhs - rhs)):.1e} (coefficient -1, pure dilation: the p^(-1/2) of (46) is the half-density e^(-u/2))")
def coefQ(Q, N):
    c = np.zeros(N); 
    for n in range(1, N + 1):
        if MU[n] and all(r in Q for r in prime_factors(n)): c[n-1] = MU[n]/n
    return c
def coefQ_fast(Qset, N):
    nn = np.arange(1, N + 1); keep = (MU[1:N+1] != 0); rest = nn.copy()
    for pp in Qset:
        rest = np.where(rest % pp == 0, rest//pp, rest)
    keep &= (rest == 1); return np.where(keep, MU[1:N+1]/nn, 0.0)
rows4 = []
for N in (10**4, 10**7):
    Qs = []; prev = None
    for pp in [int(v) for v in PR[:8]]:
        Qs.append(pp); cQ = coefQ_fast(Qs, N); EQ = energy(cQ)
        if prev is None: row = dict(N=N, p=pp, E=EQ, E_prev=None, child=None, cross=None, complete=bool(int(np.prod(Qs)) <= N))
        else:
            cprev, Eprev = prev; child = pp**-1*energy(coefQ_fast(Qs[:-1], N//pp)); cross = EQ - Eprev - child
            row = dict(N=N, p=pp, E=EQ, E_prev=Eprev, child=child, cross=cross, complete=bool(int(np.prod(Qs)) <= N))
        rows4.append(row); prev = (cQ, EQ)
        print(f"   N = {N:8d}, Q = {Qs}: R_Q = {EQ:.9f}" + ("" if row['child'] is None else f" = R_Q' {row['E_prev']:.9f} + p^-1 R_Q'(N/p) {row['child']:.9f} + cross {row['cross']:+.9f}") + f"  [{'complete' if row['complete'] else 'unfinished'}: q_Q = {int(np.prod(Qs))}]")
out['two_prime'] = rows4
# ---------------------------------------------------------------- 5. orbit sums and the fold onto C_p
print("== 5. orbit-sum identities; the fold onto C_p ==")
def Vq(N, x):   # exact rational V_N at rational x
    k = min(int(math.floor(x)), N); return sum(Fr(int(MU[n]), n)*min(x, n) for n in range(1, N + 1)) if N <= 80 else None
N = 60; pts = [Fr(j, 2) for j in range(0, 2*N + 3)]
err_p = {}
for pp in (2, 3):
    err_p[pp] = max(abs(sum(Vq(N//pp**k, x/pp**k) for k in range(0, 7) if N//pp**k >= 1) - sum(Fr(int(MU[n]), n)*min(x, n) for n in range(1, N + 1) if n % pp)) for x in pts)
err_all = max(abs(sum(Vq(N//d, x/d) for d in range(1, N + 1)) - min(x, 1)) for x in pts)
print(f"   exact at N = 60 on half-integers: sum_k V_(N/p^k)(x/p^k) - V^(p)_N(x) = {err_p[2]}, {err_p[3]} (p = 2, 3);  sum_d V_(N/d)(x/d) - min(x,1) = {err_all}")
N = 10**5; xg = np.linspace(0, 1.2*N, 20001)
orbit = sum(V(N//d, xg/d) for d in range(1, N + 1) if N//d >= 1)
print(f"   float at N = 1e5: max |sum_d V_(N/d)(x/d) - min(x,1)| = {np.max(np.abs(orbit - np.minimum(xg, 1))):.1e};  energy of the unit potential: 1")
# the fold onto C_p for N with M(N) = 0 and p = 5, N = 39
def fold(N, pp, x, kmin=-60):
    x = np.asarray(x, dtype=float); kmax = int(math.ceil(math.log(N*1.0/max(x.min(), 1e-300))/math.log(pp))) + 2
    return sum(V(N, pp**k*x) for k in range(kmin, kmax + 1))
N, pp = 39, 5; assert M[N] == 0
xf = np.exp(np.linspace(math.log(1.0), math.log(1.0) + math.log(pp), 4001)); F = fold(N, pp, xf)
per = np.max(np.abs(fold(N, pp, xf*pp) - F)); unf = [m for m in range(1, N + 1) if MU[m] and m % pp and pp*m > N]; comp = [m for m in range(1, N + 1) if MU[m] and m % pp and pp*m <= N]
# slope jumps of the fold at the classes of all breakpoints
def fold_slope(N, pp, x, eps=1e-7): return (fold(N, pp, np.array([x + eps]))[0] - fold(N, pp, np.array([x - eps]))[0])/(2*eps)
cls = {}                                                                                   # classes [m] = m p^Z, representative in [1, p)
for m in range(1, N + 1):
    if MU[m]: cls.setdefault(m*pp**(-math.floor(math.log(m)/math.log(pp))), []).append(m)
ord_by_class = {lam: sum(-int(MU[m])*(1 if pp*m > N else 0) for m in ms if m % pp) for lam, ms in cls.items()}      # parents m' coprime to p: order -mu(m')[p m' > N]
meas = {lam: lam*(fold(N, pp, np.array([lam + 1e-6]))[0] - 2*fold(N, pp, np.array([lam]))[0] + fold(N, pp, np.array([lam - 1e-6]))[0])/1e-6 for lam in cls}   # order = lambda x (slope jump), second difference
errs = max(abs(meas[lam] - ord_by_class[lam]) for lam in cls)
const = sum(MU[m] for m in comp)
print(f"   fold of V_39 onto C_5 (M(39) = 0): periodic to {per:.1e}; unfinished parents {unf} (charge {sum(MU[m] for m in unf)}), complete parents {comp} (constant {const} = M^(5)(7));  measured orders at every class vs -mu(m)[5m > N]: max error {errs:.1e};  degree = {sum(ord_by_class.values())}")
# slopes of the fold: denominators
c5 = 1
for n in range(1, N + 1):
    if n % pp: c5 = c5*n//math.gcd(c5, n)
hq = [Fr(0)]
for n in range(1, N + 1): hq.append(hq[-1] + Fr(int(MU[n]), n))
def fold_slope_exact(x):   # x in (1, p): F'(x) = h(N)/(p-1) + sum_{k>=0, p^k x < N} p^k (h(N) - h(floor(p^k x)))
    v = hq[N]/(pp - 1); k = 0
    while pp**k*x < N: v += pp**k*(hq[N] - hq[int(math.floor(pp**k*x))]); k += 1
    return v
sl_ex = [fold_slope_exact(Fr(x0)) for x0 in (Fr(13, 10), Fr(12, 5), Fr(37, 10), Fr(47, 10))]
den = [((pp - 1)*c5*v).denominator for v in sl_ex]
print(f"   exact slopes of the fold at x = 1.3, 2.4, 3.7, 4.7: {['%.9f' % float(v) for v in sl_ex]} (finite differences: {['%.9f' % float(fold_slope(N, pp, x0)) for x0 in (1.3, 2.4, 3.7, 4.7)]});  denominators of (p-1) c_5 F' with c_5 = 5'-part of lcm(1..39): {den} (powers of 5: slopes of (p-1) c F lie in Z[1/p]);  denominators of c_5 F' alone: {[(c5*v).denominator for v in sl_ex]}")
h2 = Fr(1) + Fr(-1, 2); print(f"   the factor p - 1 is needed in general: N = 2 (M(2) = 0), p = 5, x in (1, 2): F' = h(2)/4 + (h(2) - h(1)) = {h2/4 + (h2 - 1)}, c_5(2) = 2, c_5 F' = {2*(h2/4 + (h2 - 1))} not in Z[1/5], (p-1) c_5 F' = {8*(h2/4 + (h2 - 1))}")
# divergence when M(N) != 0
N2 = 38; parts = [sum(float(V(N2, pp**k*2.0)) for k in range(-40, K)) for K in (5, 10, 20)]
print(f"   N = 38 (M = {M[38]}): partial folds with k < 5, 10, 20 at x = 2: {['%.3f' % v for v in parts]} (grow by M(N) per unit of k: no global section)")
# energy not circle invariant
ug = np.linspace(0.0, l, 4001); Fu = fold(N, pp, np.exp(ug)); dF = np.gradient(Fu, ug); e_cell = np.trapezoid(np.exp(-ug)*dF**2, ug); Fu2 = fold(N, pp, np.exp(ug + l)); dF2 = np.gradient(Fu2, ug); e_cell2 = np.trapezoid(np.exp(-(ug + l))*dF2**2, ug)
print(f"   psi-energy of the fold on one period [0, log 5]: {e_cell:.6f}; on the next period: {e_cell2:.6f}; ratio {e_cell2/e_cell:.6f} = 1/5")
# ---------------------------------------------------------------- 6. the spectral identity (6) and the compactified metric (26)
print("== 6. (6) and (26) ==")
N = 20; T = 4000.0; t = np.arange(-T, T, 0.002); ns = np.arange(1, N + 1); DN = np.zeros(len(t), dtype=complex)
for n in ns: DN += MU[n]*n**(-0.5 - 1j*t)
I6 = np.trapezoid(np.abs(DN)**2/(t*t + 0.25), t)/(2*np.pi); tail = float(np.sum(MU[1:N+1]**2/ns))/(np.pi*T)
print(f"   (6) at N = 20: (1/2pi) int_{{|t|<{T:.0f}}} |D_N(1/2+it)|^2/(t^2+1/4) dt = {I6:.6f}, tail beyond T about {tail:.2e}, sum {I6 + tail:.6f} vs R(20) = {R[20]:.6f}")
N = 10; vb = [math.tanh(math.log(k)/2) for k in range(1, N + 1)]; vb = [-1.0] + vb + [1.0]; tot = 0.0
gx64, gw64 = np.polynomial.legendre.leggauss(200)
for lo, hi in zip(vb[:-1], vb[1:]):
    vv = (hi - lo)/2*gx64 + (hi + lo)/2; uu = 2*np.arctanh(vv); ph = phi_green(N, uu); dph = dphi_green(N, uu)*2/(1 - vv**2)
    tot += (hi - lo)/2*float(np.sum(gw64*((1 - vv**2)/2*dph**2 + ph**2/(2*(1 - vv**2)))))
print(f"   (26) at N = 10: v-integral = {tot:.9f} vs R(10) = {R[10]:.9f}")
# ---------------------------------------------------------------- 7. the chain-refined ledger (37) and the return (42)
print("== 7. the chain-refined projection: (37) and (42) ==")
G6 = np.gcd(idx, 6); Dcum = np.cumsum(np.where((MU != 0) & (G6 == 1), 1.0, 0.0)/n_); Dcum[0] = 0.0
def D23(N): return (2/3)*Dcum[N//6]
rows7 = []
for N in (10**6, 10**7):
    chain = []; n = N
    while n >= 6: chain.append(n); n //= 6
    chain.append(n); refined = sorted(set(mesh(N)) | set(chain)); EP = {}; ER = {}
    for Nj in chain:
        if Nj < 1: EP[Nj] = 0.0; ER[Nj] = 0.0; continue
        xs = [x for x in refined if x <= Nj]; c = MU[1:Nj+1]/np.arange(1, Nj + 1); pj = project(c, xs, Nj); EP[Nj] = energy(pj); ER[Nj] = energy(c - pj)
        assert abs(EP[Nj] + ER[Nj] - R[Nj]) < 1e-8
    eps = [ER[a] - ER[b] for a, b in zip(chain[:-1], chain[1:])]; dhat = [EP[a] - EP[b] for a, b in zip(chain[:-1], chain[1:])]
    d23 = [D23(a) - D23(b) for a, b in zip(chain[:-1], chain[1:])]; B23 = sum(max(x - y, 0.0) for x, y in zip(dhat, d23))
    bound = 43/30 + D23(N) + B23 + 2**(-4/3)
    rows7.append(dict(N=N, chain=chain, EP=[EP[a] for a in chain], ER=[ER[a] for a in chain], eps=eps, dhat=dhat, d23=d23, B23=B23, bound=bound, R=float(R[N]), m=len(refined) - 1))
    print(f"   N = {N}: refined mesh {len(refined)-1} intervals; chain {chain}")
    print(f"      eps_j = E_res(N_j) - E_res(N_j+1) = {['%.5f' % e for e in eps]} (all >= 0: {all(e >= -1e-12 for e in eps)}), sum {sum(eps):.6f} = E_res(N) - E_res(N_J) <= C0 = {2**(-4/3):.6f}")
    print(f"      Delta-hat_j = {['%+.5f' % d for d in dhat]};  d23(N_j) = {['%.5f' % d for d in d23]};  B-hat_23(N) = {B23:.6f}")
    print(f"      (42): R(N) = {R[N]:.6f} <= 43/30 + D23(N) {D23(N):.6f} + B-hat {B23:.6f} + C0 = {bound:.6f};  R(N_J) = {R[chain[-1]]:.6f} <= 43/30 = {43/30:.6f}")
out['ledger'] = rows7
# the all-horizon positive budget B23(N) = sum_j [Delta_6 R(N_j) - d23(N_j)]_+ >= B-hat_23(N), computable for every N
NN = idx[6:]; d23all = (2/3)*(Dcum[NN//6] - Dcum[(NN//6)//6])
D6 = np.zeros(X + 1); D6[6:] = R[6:] - R[idx[6:]//6]; d23v = np.zeros(X + 1); d23v[6:] = d23all
posb = np.maximum(D6 - d23v, 0.0); posb[:6] = 0.0; B23 = np.zeros(X + 1); lo = 6
while lo <= X:
    hi = min(X + 1, lo*6); B23[lo:hi] = B23[idx[lo:hi]//6] + posb[lo:hi]; lo *= 6
kb = int(np.argmax(B23)); D23v = (2/3)*Dcum[idx//6]; slack = 43/30 + D23v + B23 + 2**(-4/3) - R; ks = int(np.argmin(slack[6:])) + 6
print(f"   all-horizon budget: max_N B23(N) = {B23[kb]:.6f} at N = {kb} (chain {[int(v) for v in [kb//6**j for j in range(0, 12) if kb//6**j >= 1]]}); B23(N) = 0 for {np.mean(B23[6:] == 0)*100:.2f} % of N <= 1e7; positive parts come only from chain elements N_j <= {int(np.max(idx[posb > 0])) if np.any(posb > 0) else 0}")
print(f"      (42) with B23 for every N <= 1e7: min slack of 43/30 + D23(N) + B23(N) + C0 - R(N) is {slack[ks]:.6f} at N = {ks};  max_N [R(N) - D23(N)] = {np.max(R[6:] - D23v[6:]):.6f} at N = {int(np.argmax(R[6:] - D23v[6:])) + 6}")
out['B23'] = dict(max=float(B23[kb]), argmax=int(kb), frac_zero=float(np.mean(B23[6:] == 0)), min_slack=float(slack[ks]), argmin_slack=int(ks))
# exact complete-packet energies R_{Q_P} = sum_{d,e | q} mu(d)mu(e)/max(d,e)
def Rpacket(primes):
    divs = [1]
    for p_ in primes: divs = divs + [d*p_ for d in divs]
    return sum(Fr(int(MU[d])*int(MU[e]), max(d, e)) for d in divs for e in divs)
print("   exact complete-packet energies: " + "; ".join(f"Q={[int(v) for v in PR[:k]]}: {Rpacket([int(v) for v in PR[:k]])} = {float(Rpacket([int(v) for v in PR[:k]])):.6f}" for k in range(1, 7)))
# ---------------------------------------------------------------- 8. the sharp diagonal 4/5
print("== 8. the diagonal increment d23(N) <= 4/5 ==")
NN = idx[6:]; d23all = (2/3)*(Dcum[NN//6] - Dcum[(NN//6)//6]); k = int(np.argmax(d23all))
print(f"   max over 6 <= N <= 1e7 of d23(N) = {d23all[k]:.12f} at N = {NN[k]} (= 4/5: {abs(d23all[k] - 0.8) < 1e-12}); attained exactly on N in {[int(N) for N in NN[np.abs(d23all - 0.8) < 1e-12]]}")
print(f"   max for N >= 396 (a = floor(L/6) >= 11): {d23all[NN >= 396].max():.6f} <= analytic (2/3)[(1/3) log(6(a+1)/(a-5)) + 1/(a+1)] at a = 11: {(2/3)*((1/3)*math.log(6*12/6) + 1/12):.6f}")
out['d23'] = dict(max=float(d23all[k]), argmax=int(NN[k]), max_above_396=float(d23all[NN >= 396].max()))
# ---------------------------------------------------------------- 9. controls
print("== 9. controls ==")
Rp = np.zeros(X + 1); Q = np.cumsum((MU != 0).astype(np.int64)); innerp = np.cumsum(Q.astype(float)**2/(n_*(n_ + 1))); Rp[1:] = innerp[:-1] + Q[1:].astype(float)**2/n_[1:]
print(f"   squarefree positive: R+(1e6)/(2 (6/pi^2)^2 N) = {Rp[10**6]/(2*(6/math.pi**2)**2*10**6):.6f}, at 1e7: {Rp[X]/(2*(6/math.pi**2)**2*X):.6f}")
aC = np.where(MU != 0, 1, 0)*np.where((idx % 2 == 0), -1, 1)*np.where((idx % 3 == 0), -1, 1); aC[0] = 0; MC = np.cumsum(aC.astype(np.int64))
innerC = np.cumsum(MC.astype(float)**2/(n_*(n_ + 1))); RC = np.zeros(X + 1); RC[1:] = innerC[:-1] + MC[1:].astype(float)**2/n_[1:]
alpha = 6/math.pi**2*(1/3)*(1/2)
for N in (10**6, 10**7): print(f"   fixed signs C = {{2,3}}: Delta_6 R_C({N}) = {RC[N] - RC[N//6]:.2f} vs (5/3) alpha_C^2 N = {(5/3)*alpha**2*N:.2f} (ratio {(RC[N] - RC[N//6])/((5/3)*alpha**2*N):.4f});  mean of a_C to N: {MC[N]/N:.6f} vs alpha_C = {alpha:.6f}")
print(f"   Mobius: Delta_6 R(1e7) = {R[X] - R[X//6]:.6f}; R(1e7) = {R[X]:.6f}")
json.dump(out, open('data/global_closure.json', 'w'), indent=1, default=float); print("wrote data/global_closure.json")
