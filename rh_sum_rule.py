"""
The tail identity from the Euler-Lagrange equation (PROGRAM.md item 2), tested on an edge-FEM minimizer with primes.

With K the operator of symbol Psi (archimedean kernel -J_Gamma(y-z) off the diagonal, plus the prime shifts) and
f the minimizer, (K - 2|s><s|) f = lambda f on the window with s = sinh(y/2) 1_window, i.e. K f = lambda f + 2P sinh(y/2)
there, and K f = g_out outside it, so (Psi - lambda) F = G + 2P Sigma and, f and g_out having disjoint supports, for every T
   (1/2pi) int_{|t|>T} |F|^2 (Psi - lambda) dt  =  (2P/2pi) int_{|t|>T} conj(F) Sigma dt  -  <f_T, g_out>,
f_T the low-pass of f at T.  All three pieces are computed here for the FEM minimizer; the identity holds up to the
Galerkin residual.  The last piece is the boundary overlap: the leak of f_T just outside the edge, ~ C L^{-1/2},
against the singular part of g_out ~ -C L^{1/2}, which has no logarithm and is the tail law; the overlap with the
smooth part of g_out (regular part of the kernel, left half-window, reflected primes) and the polar piece cancel.

Usage: python3 rh_sum_rule.py data/fem_a0.5_full.json [T/T* list] [lambda'] [C beta]      or: python3 rh_sum_rule.py jl   or: python3 rh_sum_rule.py envelope FEM.json lambda' primes(0/1)
"""
import sys, json, math, numpy as np
from numpy.polynomial.legendre import leggauss
import mpmath as mp

def JL(L):
    """J(L) of paper Lemma 7.11: the edge overlap is 2C^2/(pi T) J(L(1/T)); no 1/L term, J = 1 + O(1/L^2)."""
    def gL(w, n=4000):
        sg = (np.arange(n) + 0.5)/n
        num = L + np.log(1/(w*(1-sg))); den = L + np.log(1/(w*sg))
        return np.mean(np.sqrt(np.clip(num, 1e-9, None)/np.clip(den, 1e-9, None)))
    W = 0.5*math.exp(L); eps = 5.0/W
    w = np.linspace(1e-6, W, int(min(4e6, 200*W)) + 1); wc = w[::max(1, len(w)//2000)]
    g = np.interp(w, wc, np.array([gL(x) for x in wc]))
    return np.trapezoid(np.sin(w)*np.exp(-eps*w)*g, w)/np.trapezoid(np.sin(w)*np.exp(-eps*w), w)

def envelope(fn, lam_prime, primes):
    """Computation 7.18: <t^2 (Psi - lambda) |F|^2> over a window of width 4 pi/a at T = k T*, divided by -lambda'."""
    d = json.load(open(fn)); a = d["a"]; lam = d["lambda"]
    dl = np.array(d["deltas"]); f = np.array(d["f"]); x = a - dl; o = np.argsort(x); x, f = x[o], f[o]
    if x[0] > 0: x = np.concatenate([[0.0], x]); f = np.concatenate([[0.0], f])
    B = np.diff(f)/np.diff(x); fa = f[-1]
    def F2(t):
        t = np.atleast_1d(t)[:, None]
        S = (-fa*np.cos(t[:, 0]*a)/t[:, 0]) + (B[None, :]*(np.sin(t*x[1:][None, :]) - np.sin(t*x[:-1][None, :]))).sum(1)/t[:, 0]**2
        return 4*S**2
    tg = np.concatenate([np.linspace(0.02, 60, 3000), np.linspace(60, 40000, 8000)[1:]])
    psi_g = np.array([float(mp.re(mp.psi(0, mp.mpf(1)/4 + 1j*mp.mpf(float(t))/2))) for t in tg])
    def vm(n):
        q = min(p for p in range(2, n+1) if n % p == 0); m = n
        while m % q == 0: m //= q
        return math.log(q) if m == 1 else 0.0
    pr = [n for n in range(2, 60) if math.log(n) < 2*a and vm(n)] if primes else []
    Psi = lambda t: np.interp(t, tg, psi_g) - math.log(math.pi) - sum(2*vm(n)/math.sqrt(n)*np.cos(t*math.log(n)) for n in pr)
    Ts = 2*math.pi*math.exp(2*a); W = 4*math.pi/a
    for k in (1, 2, 3, 5, 10, 20, 50, 100, 300):
        T = k*Ts; tt = np.linspace(T - W/2, T + W/2, 20001)
        print(f"  T/T*={k:4d}: <t^2 (Psi - lambda) |F|^2> / (-lambda') = {np.trapezoid(tt**2*(Psi(tt) - lam)*F2(tt), tt)/W/(-lam_prime):.4f}")

if len(sys.argv) > 1 and sys.argv[1] == "envelope":     # python3 rh_sum_rule.py envelope data/fem_a0.5_primefree.json -2.25815 0|1
    envelope(sys.argv[2], float(sys.argv[3]), bool(int(sys.argv[4]))); sys.exit(0)

if len(sys.argv) > 1 and sys.argv[1] == "jl":
    for L in (3, 4, 5, 6, 7, 8): print(f"L={L}: J(L) = {JL(L):.5f}   1 + pi^2/(24 L^2) = {1 + math.pi**2/(24*L*L):.5f}")
    sys.exit(0)

fn = sys.argv[1]; ks = [float(k) for k in sys.argv[2].split(",")] if len(sys.argv) > 2 else [5, 10, 20, 50, 100]
lam_prime = float(sys.argv[3]) if len(sys.argv) > 3 else 0.0     # exact derivative at the support (data/dilation_*.log), for the law
Cfit, bfit = (float(sys.argv[4]), float(sys.argv[5])) if len(sys.argv) > 5 else (None, None)   # edge-law fit (Computation 3.4), for the model edge overlap
d = json.load(open(fn)); a = d["a"]; lam = d["lambda"]
dl = np.array(d["deltas"]); f = np.array(d["f"]); x = a - dl; o = np.argsort(x); x, f = x[o], f[o]
if x[0] > 0: x = np.concatenate([[0.0], x]); f = np.concatenate([[0.0], f])
B = np.diff(f)/np.diff(x); A = f[:-1] - B*x[:-1]                      # f = A + B z on [x_j, x_{j+1}]
fa = f[-1]                                                            # jump at the edge (last node)
primes = [n for n in range(2, 60) if math.log(n) < 2*a and any(n == p**k for p in range(2, n+1) for k in range(1, 7) if all(p % q for q in range(2, p)))]
def vm(n):
    q = min(p for p in range(2, n+1) if n % p == 0); m = n
    while m % q == 0: m //= q
    return math.log(q) if m == 1 else 0.0
Ts = 2*math.pi*math.exp(2*a)
# --- f on the line (odd, piecewise linear, zero outside) ---
def fval(y):
    y = np.atleast_1d(y); s = np.sign(y); u = np.abs(y); out = np.zeros_like(u)
    m = u < a; idx = np.searchsorted(x, u[m], side="right") - 1; idx = np.clip(idx, 0, len(B)-1)
    out[m] = A[idx] + B[idx]*u[m]; return s*out
# --- transform: F(t) = 2i S(t), S(t) = int_0^a f sin(tx) ---
def S(t):
    t = np.atleast_1d(t)[:, None]
    return (-fa*np.cos(t[:, 0]*a)/t[:, 0]) + (B[None, :]*(np.sin(t*x[1:][None, :]) - np.sin(t*x[:-1][None, :]))).sum(1)/t[:, 0]**2
# --- symbol Psi(t) ---
def repsi(t):   # Re psi(1/4 + i t/2), vectorized through mpmath at modest precision
    return np.array([float(mp.re(mp.psi(0, mp.mpf(1)/4 + 1j*mp.mpf(float(tt))/2))) for tt in t])
def Psi(t):
    return repsi(t) - math.log(math.pi) - sum(2*vm(n)/math.sqrt(n)*np.cos(t*math.log(n)) for n in primes)
# --- polar: P and Sigma(t) = int_{-a}^{a} sinh(y/2) e^{ity} dy = 2i int_0^a sinh(y/2) sin(ty) dy ---
gx, gw = leggauss(200); ys = a*(gx+1)/2; ws = a*gw/2
P = 2*np.sum(ws*fval(ys)*np.sinh(ys/2))
def Sig_im(t):  # Sigma = 2i * Sig_im
    t = np.atleast_1d(t)[:, None]; return (np.sinh(ys/2)[None, :]*np.sin(t*ys[None, :])*ws[None, :]).sum(1)
# --- g_out(y) for y > a: -int_{-a}^{a} J(y-z) f(z) dz - sum_n Lambda(n)/sqrt(n) (f(y - log n) + f(y + log n)) ---
def J(s): return np.exp(-s/2)/(-np.expm1(-2*s))
def Jsmooth(s): return J(s) - 1/(2*s)                                  # smooth remainder; J = 1/(2s) + Jsmooth
gq, gwq = leggauss(24)
def gout(y, split=False):
    y = np.atleast_1d(y); out = np.zeros_like(y); sing = np.zeros_like(y)
    for i, yy in enumerate(y):
        tot = 0.0; ts = 0.0
        # right half [0, a]: segments; singular part exactly, smooth part by GL
        for j in range(len(B)):
            z0, z1 = x[j], x[j+1]; Aj, Bj = A[j], B[j]
            # int (A + B z)/(2(y - z)) dz, y > z: with s = y - z: (A + B y - B s)/(2 s)
            tot += 0.5*((Aj + Bj*yy)*math.log((yy - z0)/(yy - z1)) - Bj*(z1 - z0)); ts += 0.5*((Aj + Bj*yy)*math.log((yy - z0)/(yy - z1)) - Bj*(z1 - z0))
            zz = z0 + (z1 - z0)*(gq + 1)/2; wz = (z1 - z0)*gwq/2
            tot += np.sum(wz*(Aj + Bj*zz)*Jsmooth(yy - zz))
        # left half [-a, 0]: f(z) = -f(-z), J(y - z) smooth (y - z > y > a)
        for j in range(len(B)):
            z0, z1 = -x[j+1], -x[j]; zz = z0 + (z1 - z0)*(gq + 1)/2; wz = (z1 - z0)*gwq/2
            tot += np.sum(wz*fval(zz)*J(yy - zz))
        out[i] = -tot - sum(vm(n)/math.sqrt(n)*float(fval(yy - math.log(n))[0] + fval(yy + math.log(n))[0]) for n in primes); sing[i] = -ts
    return (out, sing) if split else out
# --- low-pass leak: f_T(y) = (2/pi) int_0^T S(t) sin(ty) dt ---
def fT(y, T, npts=None):
    n = npts or int(max(20000, 40*T*(y.max() + a)))
    t = np.linspace(1e-9, T, n); St = S(t)
    return np.array([(2/math.pi)*np.trapezoid(St*np.sin(t*yy), t) for yy in np.atleast_1d(y)])
print(f"{fn}: a={a} lambda_FEM={lam:.8e} primes in window {primes} P={P:.6f} f(a-)={fa:.5f} T*={Ts:.3f}")
Tmax = 400*Ts
for k in ks:
    T = k*Ts
    # (A) lhs = (1/2pi) int_{|t|>T} |F|^2 (Psi - lambda) = (1/pi) int_T^inf 4 S^2 (Psi - lambda)
    tt = np.concatenate([np.linspace(T, 4*T, 20001), np.linspace(4*T, Tmax, 40001)[1:]])
    F2 = 4*S(tt)**2; Ps = Psi(tt)
    lhs = np.trapezoid(F2*(Ps - lam), tt)/math.pi + (-lam_prime/(math.pi*Tmax) if lam_prime else 4*fa**2*(np.log(Tmax/(2*math.pi)) + 1)/(2*math.pi*Tmax))  # tail beyond Tmax: the tail law (Computation 7.6), else the jump estimate
    rhoT = np.trapezoid(F2*np.log(tt/(2*math.pi))/(2*math.pi), tt)/math.pi                                                   # (1/pi) int_T |F|^2 rho: the tail law's object
    # (B) polar = (4P/2pi) int_{|t|>T} conj(F) Sigma = (4P/pi) int_T^inf (2 S)(2 Sig_im) -> conj(2iS)(2i Sig_im) = 4 S Sig_im
    pol = (2*P/math.pi)*np.trapezoid(4*S(tt)*Sig_im(tt), tt)
    # (C) overlap = <f_T, g_out> = 2 int_a^inf f_T(y) g_out(y) dy   (f, g_out odd)
    yy = np.concatenate([a + np.logspace(-7, -2, 60), a + np.linspace(0.011, 6, 1200)])
    ft = fT(yy, T); go = gout(yy); ov = 2*np.trapezoid(ft*go, yy)
    gmodel = -Cfit*np.sqrt(np.log(a/(yy - a)) + bfit) if Cfit else np.zeros_like(yy)        # the edge-law force -C L(delta)^{1/2}, valid for small delta only
    m = (yy - a) < 0.1*a; ovs = 2*np.trapezoid((ft*gmodel)[m], yy[m])
    law = -lam_prime/(math.pi*T) if lam_prime else float('nan')
    print(f"  T/T*={k:6.1f}: lhs = (1/2pi)int_{{|t|>T}}|F|^2(Psi-lam) = {lhs:.5e}  [2 int_T |F|^2 rho = {2*math.pi*rhoT:.5e}; law -lambda'/(pi T) = {law:.5e}]   rhs = polar {pol:.3e} + overlap {-ov:.5e} = {pol - ov:.5e}   rhs/lhs = {(pol-ov)/lhs:.4f}", flush=True)
    print(f"           model edge overlap -<f_T, -C L^(1/2)> on delta < a/10 = {-ovs:.5e} = {-ovs/law:.4f} of the law (2C^2/(pi T) = {2*Cfit**2/(math.pi*T) if Cfit else float('nan'):.5e});  rest of the overlap + polar = {-(ov-ovs) + pol:.3e} = {(-(ov-ovs)+pol)/law:+.4f} of the law", flush=True)
    i = np.argmin(abs(yy - (a + 1/T))); print(f"           leak f_T(a + 1/T) = {ft[i]:.4e}   g_out(a + 1/T) = {go[i]:.4e}   g_out(a+1e-6) = {go[0]:.4e}   f_T(a+1e-6) = {ft[0]:.4e}")
