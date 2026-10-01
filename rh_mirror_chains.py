"""
The mirror chains in the window's symbol amplitude, and what the K-mode minimizers resolve (paper Theorem 9.8, Computation 9.10).

Proposition 8.17 represents the transform with the symbol of the window at every support, F~ = A_a(t)/(t sigma_{a,-}) + ...; Theorem 9.8
says that above the periodic band A_a(t) = A_0 + (mirror chains) + o(1), the mirror chains being the walks from the far edge that cross
the window, oscillating at the frequencies v - 2a (v a lattice point above 2a, the shortest words of length two: v = d + d') with the
amplitudes rho_v(s) = c_d c_d'/s^2 (times 2 for d != d'), s = sigma~_inf(t).  The echo form of Corollary 8.14(iii), exact at a support
with a finite tree, carries instead the echo phase arg E(t) = arg(1 + sum_p rho_p e^{-it(a-p)}), of first order c_d/s.

This script uses the edge FEM at a = 0.6 (entries 2, 3; data/fem_a0.6_full.json, the interior resolved to about 34 T*), and the K-mode
minimizer at the same support (data/tail_law_kmode/coefs_0.6.json, a polynomial of degree 2K-1 = 95): (i) the two transforms agree to
0.1% below 5 T* and depart above (2K-1)/a = 7.6 T*, where the polynomial stops resolving the interior cusps; (ii) the first-generation
echo phase is in the FEM nulls with the amplitudes c_2, c_3 and absent from the K-mode nulls; (iii) on the FEM nulls, the residual of
the window's symbol phase is removed by the two-step mirror chains.

Usage: python3 rh_mirror_chains.py
"""
import json, math, sys, numpy as np, mpmath as mp
from rh_band import S_kmode, repsi, vm
from rh_echo_tree import tree

a_s = "0.6"; dF = json.load(open("data/fem_a0.6_full.json")); a = float(dF["a"]); lam = float(dF["lambda"]); Ts = 2*math.pi*math.exp(2*a)
dK = json.load(open(f"data/tail_law_kmode/coefs_{a_s}.json")); K = dK["K"]; coef = dK["c"]; lamK = float(dK["lambda"]); lpK = float(dK["a_lambda_prime"])/a
lam_prime = -4.1545e-5                                                     # lambda'_64 at a = 0.6 (dilation engine, Computation 8.16)
ent = [(math.log(m), vm(m)/math.sqrt(m), m) for m in range(2, 400) if vm(m) and math.log(m) < 2*a]
X = 2*sum(c for _, c, _ in ent); t_band = 2*math.pi*math.exp(X + lam)
dl = np.array(dF["deltas"]); fF = np.array(dF["f"]); x = a - dl; o = np.argsort(x); x, fF = x[o], fF[o]
if x[0] > 0: x = np.concatenate([[0.0], x]); fF = np.concatenate([[0.0], fF])
B = np.diff(fF)/np.diff(x); fa = fF[-1]
def S_fem(t):
    t = np.atleast_1d(np.asarray(t, float)); out = np.empty_like(t)
    for i0 in range(0, len(t), 2000):
        tt = t[i0:i0+2000, None]
        out[i0:i0+2000] = -fa*np.cos(tt[:, 0]*a)/tt[:, 0] + (B[None, :]*(np.sin(tt*x[1:][None, :]) - np.sin(tt*x[:-1][None, :]))).sum(1)/tt[:, 0]**2
    return out
def S_K(t): return S_kmode(coef, a, np.atleast_1d(np.asarray(t, float)))
def sigma_a(t):
    t = np.atleast_1d(np.asarray(t, float)); return repsi(t) - math.log(math.pi) - sum(2*c*np.cos(t*dd) for dd, c, _ in ent) - lam
def sigma_inf(t):
    t = np.atleast_1d(np.asarray(t, float)); return repsi(t) - math.log(math.pi) - lam

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
            den = ti*ti - tau*tau; h = (U - ui)/np.where(np.abs(den) < 1e-12, 1.0, den); near = np.abs(tau - ti) < 1e-6
            if near.any():
                up = (u_of(ti + 1e-3)[0] - u_of(ti - 1e-3)[0])/2e-3; h[near] = -up/(2*ti)
            out[i] = ti/math.pi*np.dot(w, h)
        return out
    return zs, u_of, theta
def nulls(S, kmin, kmax):
    kk = np.arange(int((kmin*Ts*a - math.pi/2)/math.pi), int((kmax*Ts*a - math.pi/2)/math.pi) + 1); tk = []; ks = []
    for k in kk:
        t0 = (k*math.pi + math.pi/2)/a; lo, hi = t0 - 1.2/a, t0 + 1.2/a; slo, shi = S(lo)[0], S(hi)[0]
        if np.sign(slo) == np.sign(shi): continue
        for _ in range(50):
            mid = 0.5*(lo + hi); sm_ = S(mid)[0]
            if np.sign(sm_) == np.sign(slo): lo = mid
            else: hi = mid
        tk.append(0.5*(lo + hi)); ks.append(k)
    tk = np.array(tk); ks = np.array(ks); return tk, tk*a - math.pi/2 - ks*math.pi
rms = lambda v: math.sqrt(np.mean(v*v))
def fit_kappa(res0, tk): kap = np.sum(res0/tk)/np.sum(1/tk**2); return kap, res0 - kap/tk

print(f"a = {a}, entries {[m for _,_,m in ent]}, FEM nodes {len(fF)}, lambda = {lam:.4e}; K-mode K = {K}, degree {2*K-1}, lambda_K = {lamK:.4e}; T* = {Ts:.2f}; "
      f"(2K-1)/a = {(2*K-1)/a/Ts:.1f} T*, (2K-1)^2/a = {(2*K-1)**2/a/Ts:.0f} T*; X = 2 sum c_d = {X:.3f}, periodic band top {t_band/Ts:.1f} T*", flush=True)

# (i) FEM against K-mode
print("\n(i) relative L^2 distance of the K-mode transform from the FEM transform over [T - T*, T + T*], and the ratio of their mean squares:")
for k in (2, 3, 5, 7, 10, 15, 20, 30):
    T = k*Ts; tt = np.arange(T - Ts, T + Ts, 0.2); sk = S_K(tt); sf = S_fem(tt)
    print(f"   T = {k:3d} T*: distance {math.sqrt(np.trapezoid((sk - sf)**2, tt)/np.trapezoid(sf*sf, tt)):.4f}, <S_K^2>/<S_fem^2> = {np.trapezoid(sk*sk, tt)/np.trapezoid(sf*sf, tt):.4f}", flush=True)

# (ii) the first-generation echo phase in the nulls
zeros_i, u_i, theta_i = hilbert_phase(sigma_inf); zeros_a, u_a, theta_a = hilbert_phase(sigma_a)
print(f"\nreal zeros of the window's symbol: {len(zeros_a)}, the last at {zeros_a.max()/Ts:.2f} T*; of the archimedean symbol: {np.array2string(zeros_i, precision=3)}")
print("\n(ii) the nulls on 10-30 T*: residual r + theta_inf - kappa/t and its regression on (alpha_d sin(td) + beta_d cos(td))/sigma~_inf; the echo form predicts alpha_d = c_d, beta_d = 0:")
for name, S in (("FEM", S_fem), ("K-mode", S_K)):
    tk, r = nulls(S, 10, 30); th = theta_i(tk); s_ = np.exp(u_i(tk)); kap, res = fit_kappa(r + th, tk)
    cols = [1/tk] + [np.sin(tk*dd)/s_ for dd, _, _ in ent] + [np.cos(tk*dd)/s_ for dd, _, _ in ent]
    Xm = np.array(cols).T; cf, *_ = np.linalg.lstsq(Xm, r + th, rcond=None)
    print(f"   {name:7s}: {len(tk)} nulls, rms {rms(res):.4f} (kappa = {kap:.1f}); " + ", ".join(f"m = {m}: alpha = {cf[1+j]:+.4f} (c_d = {c:.4f}), beta = {cf[1+len(ent)+j]:+.4f}" for j, (dd, c, m) in enumerate(ent))
          + f"; rms after {rms(r + th - Xm @ cf):.4f}", flush=True)

# (iii) the window's symbol phase on the FEM nulls and the two-step mirror chains
Tr = tree(a); So = np.array(Tr["So"]); Cm = Tr["C"]; bv = Tr["b"]; nS = len(So)
st_inf = lambda t: np.exp(u_i(t))
def rho_of(t):
    t = np.atleast_1d(np.asarray(t, float)); sv = st_inf(t); M = sv[:, None, None]*np.eye(nS)[None] - Cm[None]
    return np.linalg.solve(M, np.broadcast_to(bv, (len(t), nS))[..., None])[..., 0]
def sigma_tree(t): t = np.atleast_1d(np.asarray(t, float)); return st_inf(t) - rho_of(t) @ bv
def E_tree(t): t = np.atleast_1d(np.asarray(t, float)); return 1 + (rho_of(t)*np.exp(-1j*np.outer(t, a - So))).sum(1)
_, _, theta_tree = hilbert_phase(sigma_tree)
words = []                                                                 # ordered two-step words with v = d + d' > 2a
for d1, c1, m1 in ent:
    for d2, c2, m2 in ent:
        if d1 + d2 > 2*a: words.append((d1 + d2, c1*c2, m1, m2))
freqs = sorted(set(round(v - 2*a, 10) for v, _, _, _ in words)); wt = {om: sum(cc for v, cc, _, _ in words if round(v - 2*a, 10) == om) for om in freqs}
print(f"\n(iii) FEM nulls on 5-34 T*: the window's symbol phase theta_a, and the regression of its residual on (p sin + q cos)(t (v - 2a))/sigma~_inf^2 over the two-step mirror words v = d + d' > 2a:")
print("   mirror frequencies v - 2a and chain weights sum c_d c_d': " + "; ".join(f"{om:.4f} ({wt[om]:.4f}: " + ",".join(f"{m1}x{m2}" for v, _, m1, m2 in words if round(v - 2*a, 10) == om) + ")" for om in freqs))
tk, r = nulls(S_fem, 5, 34); s_ = st_inf(tk)
tha = theta_a(tk); kap_a, res_a = fit_kappa(r + tha, tk)
tht = theta_tree(tk) + np.angle(E_tree(tk)); kap_t, res_t = fit_kappa(r + tht, tk)
thi = theta_i(tk); kap_i, res_i = fit_kappa(r + thi, tk)
cols = [1/tk] + [np.sin(tk*om)/s_**2 for om in freqs] + [np.cos(tk*om)/s_**2 for om in freqs]
Xm = np.array(cols).T; cf, *_ = np.linalg.lstsq(Xm, r + tha, rcond=None); fit = Xm @ cf
print(f"   {len(tk)} nulls; rms r + phase - kappa/t:  window's symbol {rms(res_a):.4f} (kappa = {kap_a:.1f}),  tree {rms(res_t):.4f} (kappa = {kap_t:.1f}),  archimedean {rms(res_i):.4f} (kappa = {kap_i:.1f})")
print("   window's symbol residual regressed on the mirror frequencies: " + "; ".join(f"omega = {om:.4f}: amplitude {math.hypot(cf[1+j], cf[1+len(freqs)+j]):.4f} (weight {wt[om]:.4f})" for j, om in enumerate(freqs)) + f"; rms after {rms(r + tha - fit):.4f}", flush=True)
# the same regression on the tree residual (control: the tree form has no mirror terms when the trees are disjoint)
cf_t, *_ = np.linalg.lstsq(Xm, r + tht, rcond=None)
print("   control, the tree's residual regressed on the same frequencies: " + "; ".join(f"omega = {om:.4f}: amplitude {math.hypot(cf_t[1+j], cf_t[1+len(freqs)+j]):.4f}" for j, om in enumerate(freqs)) + f"; rms after {rms(r + tht - Xm @ cf_t):.4f}", flush=True)
# by height band
for lo_k, hi_k in ((5, 10), (10, 20), (20, 34)):
    m = (tk >= lo_k*Ts) & (tk < hi_k*Ts)
    print(f"   {lo_k:3d}-{hi_k:3d} T* ({m.sum():3d} nulls): window's symbol rms {rms(res_a[m]):.4f}, after the mirror fit {rms((r + tha - fit)[m]):.4f}, tree {rms(res_t[m]):.4f}, archimedean {rms(res_i[m]):.4f};  (X/s)^2 at the band's centre = {(X/st_inf(0.5*(lo_k+hi_k)*Ts)[0])**2:.4f}")
