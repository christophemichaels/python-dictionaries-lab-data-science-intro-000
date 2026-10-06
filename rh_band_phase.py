"""
The K-mode minimizers above the horizon: the edge's transform, and what the polynomial does not resolve (paper Computation 9.18).

Proposition 8.18 gives, at every support, F~ = A_a(t)/(t sigma_{a,-}) + e^{-2ita} Phi_-/sigma~_a with the window's symbol sigma_a, and
Theorem 9.12 says A_a(t) -> A_0 above the periodic band, with the mirror chains (walks from the far edge across the window) as the leading
oscillating correction.  The K-mode minimizer is a polynomial of degree 2K-1 on [-a, a]: it resolves the edge to delta ~ a/(2K-1)^2, i.e.
heights of order (2K-1)^2/a (several hundred horizons), but not the interior echo cusps above heights of order (2K-1)/a (a few horizons), so
above a few horizons its transform is the edge's transform alone: the nulls follow the archimedean Hilbert phase, the amplitude is
(-lambda'/2)^(1/2), and neither the echo phase of Corollary 8.15(iii) nor the window's symbol phase is present (block (iv) tests the
echo phase directly: its regression coefficients vanish, where the FEM at a = 0.6 gives c_2, c_3; rh_mirror_chains.py).  Blocks (i)-(iii)
therefore measure the edge of the minimizer, at heights up to several hundred horizons, against the sharp form with the archimedean symbol
and against the window's symbol with the asymptotic amplitude.

Usage: python3 rh_band_phase.py [a] [kmin_T* kmax_T*] [coefs.json] [band of the echo regression klo khi]      (default a = 1.0, 3-400 T*, data/tail_law_kmode/coefs_a.json, 10-30 T*)
"""
import json, math, sys, numpy as np, mpmath as mp
from rh_band import S_kmode, repsi, vm

a_s = sys.argv[1] if len(sys.argv) > 1 else "1.0"
kmin, kmax = (float(sys.argv[2]), float(sys.argv[3])) if len(sys.argv) > 3 else (3.0, 400.0)
cf = sys.argv[4] if len(sys.argv) > 4 else f"data/tail_law_kmode/coefs_{a_s}.json"                  # optional coefficient file
d = json.load(open(cf)); a = float(d["a"]); K = d["K"]; coef = d["c"]
klo, khi = (float(sys.argv[5]), float(sys.argv[6])) if len(sys.argv) > 6 else (10.0, 30.0)         # the band of the echo regression
lam = float(d["lambda"]); lam_prime = float(d["a_lambda_prime"])/a; Ts = 2*math.pi*math.exp(2*a)
ent = [(math.log(m), vm(m)/math.sqrt(m), m) for m in range(2, 400) if vm(m) and math.log(m) < 2*a]
S2c = 2*sum(c for _, c, _ in ent); t_band = 2*math.pi*math.exp(S2c + lam)
def S(t): return S_kmode(coef, a, np.atleast_1d(np.asarray(t, float)))
def sigma(t):
    t = np.atleast_1d(np.asarray(t, float))
    return repsi(t) - math.log(math.pi) - sum(2*c*np.cos(t*dd) for dd, c, _ in ent) - lam
def sigma_inf(t):
    t = np.atleast_1d(np.asarray(t, float)); return repsi(t) - math.log(math.pi) - lam

# --- the Hilbert phase of a real even symbol (as in rh_wiener_hopf.py) ---
tau = np.concatenate([np.arange(0.0, 8000, 0.01), np.geomspace(8000, 1e11, 40000)[1:]])
dtau = np.diff(tau); w = np.zeros_like(tau); w[:-1] += dtau/2; w[1:] += dtau/2
def real_zeros(sig, tmax=3000.0):
    tz = np.arange(0.005, tmax, 0.01); sz = sig(tz); zs = []
    for i in np.nonzero(np.sign(sz[:-1]) != np.sign(sz[1:]))[0]:
        lo, hi = tz[i], tz[i+1]; slo = sz[i]
        for _ in range(60):
            mid = 0.5*(lo + hi); sm_ = sig(mid)[0]
            if np.sign(sm_) == np.sign(slo): lo = mid
            else: hi = mid
        zs.append(0.5*(lo + hi))
    return np.array(zs)
def hilbert_phase(sig):
    zs = real_zeros(sig)
    def u_of(t):
        t = np.atleast_1d(np.asarray(t, float)); s = sig(t)
        for tj in zs: s = s*(t*t + tj*tj)/(t*t - tj*tj)
        return np.log(s)
    U = u_of(tau)
    def theta(t):
        t = np.atleast_1d(np.asarray(t, float)); ut = u_of(t); out = np.empty_like(t)
        for i, (ti, ui) in enumerate(zip(t, ut)):
            den = ti*ti - tau*tau; h = (U - ui)/np.where(np.abs(den) < 1e-12, 1.0, den)
            near = np.abs(tau - ti) < 1e-6
            if near.any():
                up = (u_of(ti + 1e-3)[0] - u_of(ti - 1e-3)[0])/2e-3; h[near] = -up/(2*ti)
            out[i] = ti/math.pi*np.dot(w, h)
        return out
    return zs, u_of, theta

print(f"a = {a}, K = {K}, lambda = {lam:.4e}, lambda' = {lam_prime:.4e}, T* = {Ts:.2f}, entries {[m for _,_,m in ent]}, 2 sum c_d = {S2c:.3f}, periodic band top {t_band/Ts:.1f} T*", flush=True)
zeros_a, u_a, theta_a = hilbert_phase(sigma); zeros_i, u_i, theta_i = hilbert_phase(sigma_inf)
print(f"real zeros of the window's symbol: {len(zeros_a)}, the last at {zeros_a.max()/Ts:.1f} T*; of the archimedean symbol: {np.array2string(zeros_i, precision=3)}", flush=True)

# --- the nulls of S ---
kk = np.arange(int((kmin*Ts*a - math.pi/2)/math.pi), int((kmax*Ts*a - math.pi/2)/math.pi) + 1)
tk = []; ks = []; miss = 0
for k in kk:
    t0 = (k*math.pi + math.pi/2)/a; lo, hi = t0 - 1.2/a, t0 + 1.2/a
    slo, shi = S(lo)[0], S(hi)[0]
    if np.sign(slo) == np.sign(shi): miss += 1; continue
    for _ in range(50):
        mid = 0.5*(lo + hi); sm_ = S(mid)[0]
        if np.sign(sm_) == np.sign(slo): lo = mid
        else: hi = mid
    tk.append(0.5*(lo + hi)); ks.append(k)
tk = np.array(tk); ks = np.array(ks); r = tk*a - math.pi/2 - ks*math.pi
th_a = theta_a(tk); th_i = theta_i(tk)
rms = lambda v: math.sqrt(np.mean(v*v)) if len(v) else float('nan')
print(f"\n(i) {len(tk)} nulls on {tk[0]/Ts:.1f}-{tk[-1]/Ts:.1f} T* ({miss} brackets without a sign change); r_k = t_k a - pi/2 - k pi", flush=True)
for lo_k, hi_k in ((3, 10), (10, 30), (30, 47), (47, 100), (100, 200), (200, 400)):
    m = (tk >= lo_k*Ts) & (tk < hi_k*Ts)
    if m.sum() < 5: continue
    ka = np.sum((r + th_a)[m]/tk[m])/np.sum(1/tk[m]**2); ki = np.sum((r + th_i)[m]/tk[m])/np.sum(1/tk[m]**2)
    print(f"   {lo_k:4d}-{hi_k:3d} T* ({m.sum():4d} nulls): rms r = {rms(r[m]):.4f};  window's symbol: rms r + theta_a - kappa/t = {rms((r + th_a - ka/tk)[m]):.4f} (kappa = {ka:7.2f});"
          f"  archimedean symbol: rms r + theta_inf - kappa/t = {rms((r + th_i - ki/tk)[m]):.4f} (kappa = {ki:7.2f})", flush=True)

# --- the pointwise form over windows of width 4 T* with |A| = (-lambda'/2)^(1/2) ---
Aabs = math.sqrt(-lam_prime/2)
mb = tk >= t_band
if mb.sum() < 5: mb = tk >= tk[len(tk)//2]
ka = np.sum((r + th_a)[mb]/tk[mb])/np.sum(1/tk[mb]**2); ki = np.sum((r + th_i)[mb]/tk[mb])/np.sum(1/tk[mb]**2)
print(f"\n(ii) relative L^2 distance of S from -|A| cos(ta + theta - kappa/t)/(t sigma~^(1/2)) over [T - 2T*, T + 2T*] (step 0.25), |A| = (-lambda'/2)^(1/2) = {Aabs:.4e}, kappa from the fit on the nulls above the band top ({mb.sum()} nulls): window's symbol kappa = {ka:.2f}, archimedean kappa = {ki:.2f}")
print("    T/T*   dist(window's symbol)   dist(archimedean symbol)   <t^2 sigma_a |F|^2>/(-lambda') over the window")
for k in (5, 10, 20, 30, 47, 70, 100, 200, 300):
    if k*Ts > kmax*Ts: continue
    T = k*Ts; tt = np.arange(T - 2*Ts, T + 2*Ts, 0.25); s_ = S(tt)
    sa = np.exp(u_a(tt)); si = np.exp(u_i(tt))
    swa = -Aabs*np.cos(tt*a + theta_a(tt) - ka/tt)/(tt*np.sqrt(sa)); swi = -Aabs*np.cos(tt*a + theta_i(tt) - ki/tt)/(tt*np.sqrt(si))
    da = min(math.sqrt(np.trapezoid((s_ - sg*swa)**2, tt)/np.trapezoid(s_*s_, tt)) for sg in (1, -1))
    di = min(math.sqrt(np.trapezoid((s_ - sg*swi)**2, tt)/np.trapezoid(s_*s_, tt)) for sg in (1, -1))
    env = np.mean(tt*tt*sigma(tt)*4*s_*s_)/(-lam_prime)
    print(f"   {k:5d}   {da:14.4f}   {di:18.4f}   {env:14.4f}", flush=True)

# --- the local amplitude: |A(t)|^2/|A_0|^2 = <4 t^2 sigma~_a S^2>_period/(-lambda') on one period 4 pi/a per horizon, its mean and scatter by height band ---
print(f"\n(iii) the local amplitude |A(t)|^2/(-lambda'/2) = <4 t^2 sigma~_a S^2> over one period 4 pi/a (101 points), one period per horizon; mean and standard deviation by height band:")
Wd = 4*math.pi/a; loc = []
for k in range(int(kmin), int(kmax) + 1):
    T = k*Ts; tt = np.linspace(T - Wd/2, T + Wd/2, 101); s_ = S(tt); sa = np.exp(u_a(tt))
    loc.append((T, float(np.trapezoid(4*tt*tt*sa*s_*s_, tt)/Wd/(-lam_prime))))   # the period mean of 4 t^2 sigma~ S^2 is 2|A|^2 = -lambda' under the sharp form
loc = np.array(loc)
for lo_k, hi_k in ((3, 10), (10, 30), (30, 47), (47, 100), (100, 200), (200, 400)):
    m = (loc[:, 0] >= lo_k*Ts) & (loc[:, 0] < hi_k*Ts)
    if m.sum() < 3: continue
    v = loc[m, 1]; print(f"   {lo_k:4d}-{hi_k:3d} T* ({m.sum():4d} periods): mean {v.mean():.4f}, std {v.std():.4f}, min {v.min():.4f}, max {v.max():.4f}", flush=True)

# --- (iv) the first-generation echo phase: regression of r + theta_inf - kappa/t on (alpha_d sin(td) + beta_d cos(td))/sigma~_inf on 10-30 T*;
#     the echo form (Corollary 8.15(iii)) has alpha_d = c_d, beta_d = 0, which the FEM at a = 0.6 reproduces (rh_mirror_chains.py) ---
m = (tk >= klo*Ts) & (tk < khi*Ts); s_ = np.exp(u_i(tk[m]))
cols = [1/tk[m]] + [np.sin(tk[m]*dd)/s_ for dd, _, _ in ent] + [np.cos(tk[m]*dd)/s_ for dd, _, _ in ent]
Xm = np.array(cols).T; cf, *_ = np.linalg.lstsq(Xm, (r + th_i)[m], rcond=None); k0 = np.sum((r + th_i)[m]/tk[m])/np.sum(1/tk[m]**2)
print(f"\n(iv) the echo phase on {klo:g}-{khi:g} T* ({m.sum()} nulls): rms r + theta_inf - kappa/t = {rms((r + th_i)[m] - k0/tk[m]):.4f}; regression on (alpha_d sin(td) + beta_d cos(td))/sigma~_inf:")
for j, (dd, c, mm) in enumerate(ent): print(f"      m = {mm}: alpha = {cf[1+j]:+.4f} (echo form: c_d = {c:.4f}), beta = {cf[1+len(ent)+j]:+.4f}")
print(f"   rms after the regression {rms((r + th_i)[m] - Xm @ cf):.4f}; the polynomial resolves the interior to about (2K-1)/a = {(2*K-1)/a/Ts:.1f} T* and the edge to about (2K-1)^2/a = {(2*K-1)**2/a/Ts:.0f} T*", flush=True)
