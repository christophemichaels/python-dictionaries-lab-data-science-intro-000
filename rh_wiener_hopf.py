"""
The Wiener-Hopf phase of the edge (paper Conjecture 7.19, Computation 7.20).

For a critical point f of the window [-a, a] with symbol sigma = Psi - lambda, the half-line problem at the right edge
gives F(t) e^{-ita} = A/(t sigma_-(t)) + O(t^{-2}), where sigma = sigma_+ sigma_- is the canonical factorization
(sigma_+ analytic and non-vanishing in the upper half-plane, sigma_- = conj(sigma_+) on the line, theta = arg sigma_+ odd)
and A is purely imaginary.  Hence, with S(t) = int_0^a f sin(tx) dx (F = 2iS),
    S(t) = -|A| cos(ta + theta(t)) / (t sigma(t)^{1/2}) + O(t^{-2}),      |A|^2 = -lambda'/2  (envelope law),
    theta(t) = (1/2) H[log sigma](t) = (t/pi) int_0^inf (log sigma(tau) - log sigma(t)) / (t^2 - tau^2) dtau.
The real zeros +-t_j of sigma are divided out first (sigma~ = sigma prod_j (t^2 + t_j^2)/(t^2 - t_j^2) > 0); they, the
polar term and the far edge contribute one real constant kappa in a phase kappa/t.  Two tests on the FEM critical point:
  (i)  the nulls t_k of S:  t_k a = pi/2 + k pi - theta(t_k) + kappa/t_k   (no free phase; kappa fitted, or 0);
  (ii) S and its Wiener-Hopf form agree in L^2 over each interference period, with an error decaying like 1/t.

Usage: python3 rh_wiener_hopf.py data/fem_a0.5_primefree.json -2.25815 [kmin_T* kmax_T*]
"""
import sys, json, math, numpy as np, mpmath as mp

fn = sys.argv[1]; lam_prime = float(sys.argv[2])
kmin, kmax = (float(sys.argv[3]), float(sys.argv[4])) if len(sys.argv) > 4 else (3.0, 120.0)
d = json.load(open(fn)); a = d["a"]; lam = d["lambda"]
dl = np.array(d["deltas"]); f = np.array(d["f"]); x = a - dl; o = np.argsort(x); x, f = x[o], f[o]
if x[0] > 0: x = np.concatenate([[0.0], x]); f = np.concatenate([[0.0], f])
B = np.diff(f)/np.diff(x); fa = f[-1]; Ts = 2*math.pi*math.exp(2*a)
def vm(n):
    q = min(p for p in range(2, n+1) if n % p == 0); m = n
    while m % q == 0: m //= q
    return math.log(q) if m == 1 else 0.0
primes = [n for n in range(2, 60) if math.log(n) < 2*a and vm(n)] if d.get("primes") else []

def S(t):                                        # S(t) = int_0^a f sin(tx) dx, f piecewise linear with a jump fa at the edge
    t = np.atleast_1d(np.asarray(t, float)); out = np.empty_like(t)
    for i0 in range(0, len(t), 4000):
        tt = t[i0:i0+4000, None]
        out[i0:i0+4000] = -fa*np.cos(tt[:, 0]*a)/tt[:, 0] + (B[None, :]*(np.sin(tt*x[1:][None, :]) - np.sin(tt*x[:-1][None, :]))).sum(1)/tt[:, 0]**2
    return out

def repsi(t):                                    # Re psi(1/4 + i t/2): mpmath below 60, the Stirling series above (error < 1e-14)
    t = np.atleast_1d(np.asarray(t, float)); out = np.empty_like(t); sm = t < 60
    out[sm] = [float(mp.re(mp.psi(0, mp.mpf(1)/4 + 1j*mp.mpf(float(v))/2))) for v in t[sm]]
    z = 0.25 + 0.5j*t[~sm]; s = np.log(z) - 1/(2*z)
    for k, b in enumerate([1/6, -1/30, 1/42, -1/30, 5/66, -691/2730], 1): s -= b/(2*k*z**(2*k))
    out[~sm] = s.real; return out
def sigma(t):
    t = np.atleast_1d(np.asarray(t, float))
    return repsi(t) - math.log(math.pi) - sum(2*vm(n)/math.sqrt(n)*np.cos(t*math.log(n)) for n in primes) - lam

# --- real zeros of the symbol ---
tz = np.arange(0.005, 3000, 0.01); sz = sigma(tz); zeros = []
for i in np.nonzero(np.sign(sz[:-1]) != np.sign(sz[1:]))[0]:
    lo, hi = tz[i], tz[i+1]; slo = sz[i]
    for _ in range(60):
        mid = 0.5*(lo + hi); sm_ = sigma(mid)[0]
        if np.sign(sm_) == np.sign(slo): lo = mid
        else: hi = mid
    zeros.append(0.5*(lo + hi))
zeros = np.array(zeros)
def u_of(t):                                     # log sigma~, the positive symbol with the real zeros divided out
    t = np.atleast_1d(np.asarray(t, float)); s = sigma(t)
    for tj in zeros: s = s*(t*t + tj*tj)/(t*t - tj*tj)
    return np.log(s)
# --- master grid for the Hilbert transform ---
tau = np.concatenate([np.arange(0.0, 8000, 0.01), np.geomspace(8000, 1e11, 40000)[1:]])
U = u_of(tau); dtau = np.diff(tau); w = np.zeros_like(tau); w[:-1] += dtau/2; w[1:] += dtau/2      # trapezoid weights
def theta(t):                                    # (1/2) H[u](t) = (t/pi) int_0^inf (u(tau) - u(t)) / (t^2 - tau^2) dtau, odd in t
    t = np.atleast_1d(np.asarray(t, float)); ut = u_of(t); out = np.empty_like(t)
    for i, (ti, ui) in enumerate(zip(t, ut)):
        den = ti*ti - tau*tau; h = (U - ui)/np.where(np.abs(den) < 1e-12, 1.0, den)
        near = np.abs(tau - ti) < 1e-6
        if near.any():
            up = (u_of(ti + 1e-3)[0] - u_of(ti - 1e-3)[0])/2e-3; h[near] = -up/(2*ti)
        out[i] = ti/math.pi*np.dot(w, h)
    return out
def theta_local(t):                              # the adiabatic value -(pi/4) t u'(t)
    t = np.atleast_1d(np.asarray(t, float)); return -(math.pi/4)*t*(u_of(t*(1 + 1e-4)) - u_of(t*(1 - 1e-4)))/(2e-4*t)

# --- self-test of the Hilbert transform on log((tau^2+1)/(tau^2+4)) -> arctan(1/t) - arctan(2/t) ---
Usave = U.copy(); U = np.log((tau**2 + 1)/(tau**2 + 4))
_u = u_of; u_of = lambda t: np.log((np.atleast_1d(np.asarray(t, float))**2 + 1)/(np.atleast_1d(np.asarray(t, float))**2 + 4))
tt_ = np.array([0.7, 3.0, 25.0]); err = np.abs(theta(tt_) - (np.arctan(1/tt_) - np.arctan(2/tt_))).max()
U = Usave; u_of = _u
print(f"{fn}: a = {a}, lambda = {lam:.6f}, lambda' = {lam_prime}, primes {primes}, T* = {Ts:.3f}; nodes {len(f)}")
print(f"real zeros of the symbol: {np.array2string(zeros, precision=3)}  (t/T*: {np.array2string(zeros/Ts, precision=3)})")
print(f"Hilbert-transform self-test: max error {err:.2e}")

# --- (i) the nulls of S ---
kk = np.arange(int((kmin*Ts*a - math.pi/2)/math.pi), int((kmax*Ts*a - math.pi/2)/math.pi) + 1)
tk = []; ks = []
for k in kk:
    lo, hi = (k*math.pi + math.pi/2 - 0.6)/a, (k*math.pi + math.pi/2 + 1.4)/a; slo, shi = S(lo)[0], S(hi)[0]
    if np.sign(slo) == np.sign(shi): print(f"  no sign change for k = {k}"); continue
    for _ in range(64):
        mid = 0.5*(lo + hi); sm_ = S(mid)[0]
        if np.sign(sm_) == np.sign(slo): lo = mid
        else: hi = mid
    tk.append(0.5*(lo + hi)); ks.append(k)
tk = np.array(tk); ks = np.array(ks); r = tk*a - math.pi/2 - ks*math.pi                     # raw null shift
th = theta(tk); thl = theta_local(tk); thn = -math.pi/(4*np.log(tk/(2*math.pi)))
sel = tk >= 10*Ts
kappa = np.sum((r + th)[sel]/tk[sel])/np.sum(1/tk[sel]**2)                                    # r + theta = kappa/t on t >= 10 T*
kappa0 = np.sum(r[sel]/tk[sel])/np.sum(1/tk[sel]**2)                                            # alternative: r = kappa0/t alone
c_log = np.sum(r[sel]/np.log(tk[sel]/(2*math.pi)))/np.sum(1/np.log(tk[sel]/(2*math.pi))**2)    # alternative: r = c/log(t/2pi)
print(f"\n(i) nulls of S(t) at t_k a = pi/2 + k pi + r_k, {len(tk)} nulls on {tk[0]/Ts:.1f}-{tk[-1]/Ts:.1f} T*; fits on t >= 10 T*:")
print(f"    kappa = {kappa:.4f} (r + theta = kappa/t);  kappa0 = {kappa0:.3f} (r = kappa0/t);  c = {c_log:.4f} (r = c/log(t/2pi); Wiener-Hopf: -pi/4 = {-math.pi/4:.4f})")
print("    t_k/T*      r_k    -theta   -theta_loc  -theta_naive   r+theta   r+theta-kappa/t   r-kappa0/t   r-c/log")
for k in [3, 4, 5, 7, 10, 15, 20, 30, 50, 70, 100, 120]:
    i = np.argmin(np.abs(tk/Ts - k))
    if abs(tk[i]/Ts - k) > 0.6: continue
    print(f"   {tk[i]/Ts:7.2f}  {r[i]:8.4f}  {-th[i]:8.4f}  {-thl[i]:8.4f}  {-thn[i]:9.4f}     {r[i]+th[i]:8.4f}  {r[i]+th[i]-kappa/tk[i]:11.5f}   {r[i]-kappa0/tk[i]:9.4f}  {r[i]-c_log*np.log(tk[i]/(2*math.pi))**-1:8.4f}")
rms = lambda v: math.sqrt(np.mean(v*v))
if primes:                                       # the smooth alternative: the phase of the prime-free symbol at the same lambda, plus kappa/t
    primes_save = primes; primes = []; zeros_save = zeros
    tz = np.arange(0.005, 3000, 0.01); sz = sigma(tz); zeros = []
    for i in np.nonzero(np.sign(sz[:-1]) != np.sign(sz[1:]))[0]:
        lo, hi = tz[i], tz[i+1]; slo = sz[i]
        for _ in range(60):
            mid = 0.5*(lo + hi); sm_ = sigma(mid)[0]
            if np.sign(sm_) == np.sign(slo): lo = mid
            else: hi = mid
        zeros.append(0.5*(lo + hi))
    zeros = np.array(zeros); zeros_inf = zeros; U = u_of(tau); th_inf = theta(tk); primes = primes_save; zeros = zeros_save; U = u_of(tau)
    kap_inf = np.sum((r + th_inf)[sel]/tk[sel])/np.sum(1/tk[sel]**2)
    print(f"    smooth alternative (phase of the prime-free symbol at the same lambda, zeros {np.array2string(zeros_inf, precision=3)}, plus kappa/t, kappa = {kap_inf:.3f}):"
          f" rms r+theta_inf-kappa/t = {rms((r+th_inf-kap_inf/tk)[sel]):.4f} (t >= 10 T*), {rms(r+th_inf-kap_inf/tk):.4f} (t >= 3 T*)")
print(f"    rms over t >= 10 T*:  r: {rms(r[sel]):.4f}   r+theta: {rms((r+th)[sel]):.4f}   r+theta-kappa/t: {rms((r+th-kappa/tk)[sel]):.5f}"
      f"   r-kappa0/t: {rms((r-kappa0/tk)[sel]):.4f}   r-c/log: {rms((r-c_log/np.log(tk/(2*math.pi)))[sel]):.4f}")
print(f"    rms over t >= 3 T*:   r: {rms(r):.4f}   r+theta: {rms(r+th):.4f}   r+theta-kappa/t: {rms(r+th-kappa/tk):.5f}   r-kappa0/t: {rms(r-kappa0/tk):.4f}   r-c/log: {rms(r-c_log/np.log(tk/(2*math.pi))):.4f}")

# --- (ii) S against its Wiener-Hopf form over one interference period, L^2 ---
Aabs = math.sqrt(-lam_prime/2); Wd = 4*math.pi/a
print(f"\n(ii) relative L^2 distance of S from -|A| cos(ta + theta - kappa/t)/(t sigma^{{1/2}}), |A| = (-lambda'/2)^(1/2) = {Aabs:.5f}, over one period 4 pi/a:")
print("    T/T*   dist(kappa)   dist(kappa=0)   dist(theta=0)   dist*T/T*   sign")
for k in (3, 5, 10, 20, 50, 100):
    T = k*Ts; tt = np.linspace(T - Wd/2, T + Wd/2, 801); s = S(tt); sg = sigma(tt); thv = theta(tt); amp = Aabs/(tt*np.sqrt(sg))
    best = None
    for sign in (1, -1):
        swh = -sign*amp*np.cos(tt*a + thv - kappa/tt); dist = math.sqrt(np.trapezoid((s - swh)**2, tt)/np.trapezoid(s*s, tt))
        if best is None or dist < best[0]: best = (dist, sign)
    sign = best[1]
    d0 = math.sqrt(np.trapezoid((s + sign*amp*np.cos(tt*a + thv))**2, tt)/np.trapezoid(s*s, tt))
    dn = math.sqrt(np.trapezoid((s + sign*amp*np.cos(tt*a))**2, tt)/np.trapezoid(s*s, tt))
    print(f"   {k:5d}   {best[0]:10.5f}   {d0:12.5f}   {dn:12.5f}   {best[0]*k:9.4f}   {sign:+d}")
