"""
The sampling defect of the zeros above a height is the prime sum at the entries (paper Proposition 7.29, Corollary 7.30, Computation 7.31; the echo-tree form of the constant, Corollary 8.13).

Under RH, Weil's explicit formula for the even test function h = |F|^2 chi, chi an entire even cutoff, reads
    sum_gamma h(gamma) = (1/2pi) int h Psi_inf - sum_n Lambda(n) n^{-1/2} [h^(log n) + h^(-log n)] + h(i/2) + h(-i/2),
h^ = (1/2pi) int h(t) e^{itx} dt = (g * chi^)(x) the smoothed autocorrelation of f.  With chi = chi_T a high-pass at height T
the left side is the (smoothed) tail of the zero sum, the first term on the right the density integral of (S_kappa), and the
prime sum is the only thing between them: it lives on the entries, where the autocorrelation of the minimizer is singular.
This script evaluates every term on the edge-FEM minimizers with primes, with the first 6000 zeros, for the sharp tail
(the definition of M(kappa) in the paper) and for the entire cutoff chi_T(t) = 1 - exp(-(t/T)^4), and compares the defect with
the leading form of the prime term, 2 c^2 / sigma(T)^2 for one entry inside the window (c = Lambda(n) n^{-1/2}).

Usage: python3 rh_sampling_defect.py [2m]   (the cutoff 1 - exp(-(t/T)^{2m}), m even, default 2m = 4; its polar term is 2P^2(1 - exp(-(2T)^{-2m})))
"""
import json, math, re, sys, numpy as np, mpmath as mp
EXP = int(sys.argv[1]) if len(sys.argv) > 1 else 4
assert EXP % 4 == 0, 'the exponent 2m must have m even, so that chi_T(i/2) = 1 - exp(-(2T)^{-2m}) is positive'

def vm(n):
    q = min(p for p in range(2, n+1) if n % p == 0); m = n
    while m % q == 0: m //= q
    return math.log(q) if m == 1 else 0.0
def repsi(t):
    t = np.atleast_1d(np.asarray(t, float)); out = np.empty_like(t); sm = t < 60
    out[sm] = [float(mp.re(mp.psi(0, mp.mpf(1)/4 + 1j*mp.mpf(float(v))/2))) for v in t[sm]]
    z = 0.25 + 0.5j*t[~sm]; s = np.log(z) - 1/(2*z)
    for k, b in enumerate([1/6, -1/30, 1/42, -1/30, 5/66, -691/2730], 1): s -= b/(2*k*z**(2*k))
    out[~sm] = s.real; return out

zeros = np.array([float(l.split()[1]) for l in open("data/zeros_6000.txt")])
tg = np.concatenate([np.arange(0.005, 3000, 0.01), np.arange(3000, 1.2e5, 0.05)])         # integration grid (oscillation period pi/a ~ 6-8)
Psi_g = repsi(tg) - math.log(math.pi)
from rh_echo_tree import tree
for a, fn in [(0.4, "data/fem_a0.4_full.json"), (0.45, "data/fem_a0.45_full.json"), (0.5, "data/fem_a0.5_full.json"), (0.6, "data/fem_a0.6_full.json")]:
    d = json.load(open(fn)); lam = d["lambda"]; dl = np.array(d["deltas"]); f = np.array(d["f"]); x = a - dl; o = np.argsort(x); x, f = x[o], f[o]
    if x[0] > 0: x = np.concatenate([[0.0], x]); f = np.concatenate([[0.0], f])
    B = np.diff(f)/np.diff(x); fa = f[-1]
    txt = open(f"data/dilation_a{a}_full.log").read(); lp = float(re.search(r"a lambda' = \[([-0-9.e+]+)", txt).group(1))/a
    def F2(t):
        t = np.atleast_1d(np.asarray(t, float)); out = np.empty_like(t)
        for i0 in range(0, len(t), 4000):
            tt = t[i0:i0+4000, None]
            S = -fa*np.cos(tt[:, 0]*a)/tt[:, 0] + (B[None, :]*(np.sin(tt*x[1:][None, :]) - np.sin(tt*x[:-1][None, :]))).sum(1)/tt[:, 0]**2
            out[i0:i0+4000] = 4*S*S
        return out
    from numpy.polynomial.legendre import leggauss
    gx, gw = leggauss(400); ys = a*(gx+1)/2; ws = a*gw/2
    fy = np.interp(ys, x, f); P = 2*np.sum(ws*fy*np.sinh(ys/2))
    Ts = 2*math.pi*math.exp(2*a); F2g = F2(tg); F2z = F2(zeros); rho_g = np.log(tg/(2*math.pi))/(2*math.pi)
    entries = [n for n in range(2, 400) if vm(n) > 0]; c2 = {n: vm(n)/math.sqrt(n) for n in entries}
    inside = [n for n in entries if math.log(n) < 2*a]
    print(f"{fn.split('/')[-1]}: a = {a}, lambda = {lam:.6e}, lambda' = {lp:.6e}, P = {P:.5f}, T* = {Ts:.3f}, entries inside the window {inside}, "
          f"first 6000 zeros to {zeros[-1]:.1f} = {zeros[-1]/Ts:.0f} T*; 2 sum_6000 |F|^2 / lambda = {2*F2z.sum()/lam:.5f}")
    Tr = tree(a); Cm = Tr["C"]; bv = Tr["b"]; nS = len(bv)
    def M_tree(sg):                                   # paper Corollary 8.14(vi): M_chi - 1 = [R + sigma~(|rho|^2 - 1)]/sigma_eff, rho = (sigma~ - C)^{-1} b, R = b.rho
        rho = np.linalg.solve(sg*np.eye(nS) - Cm, bv); R = float(bv @ rho); return 1 + (R + sg*float(rho @ rho))/(sg - R)
    # the tree's prediction averaged over the cutoff (paper Corollary 8.14(vi)): M_chi - 1 = <e^T C_a e / sigma~> / <sigma_eff / sigma~>,
    # <g> = int chi_T g t^-2 dt / int chi_T t^-2 dt, with sigma~ = sigma~_inf(t) on the grid (the real zero t_lambda divided out)
    f_sig = lambda t: float(repsi(np.array([t]))[0]) - math.log(math.pi) - lam
    lo, hi = 3.0, 12.0
    for _ in range(60):
        mid = 0.5*(lo + hi)
        if f_sig(mid) < 0: lo = mid
        else: hi = mid
    tl = 0.5*(lo + hi)
    mask_t = tg > 1.5*tl
    sig_g = (Psi_g[mask_t] - lam)*(tg[mask_t]**2 + tl*tl)/(tg[mask_t]**2 - tl*tl); tgm = tg[mask_t]
    sv = np.linspace(sig_g.min(), sig_g.max(), 400); eCe_s = np.empty_like(sv); seff_s = np.empty_like(sv)
    for i, sg_ in enumerate(sv):
        rho = np.linalg.solve(sg_*np.eye(nS) - Cm, bv); R = float(bv @ rho); eCe_s[i] = R + sg_*float(rho @ rho); seff_s[i] = sg_ - R
    eCe_g = np.interp(sig_g, sv, eCe_s); seff_g = np.interp(sig_g, sv, seff_s); b2 = float(bv @ bv)
    def M_avg(T):
        chi = 1 - np.exp(-(tgm/T)**EXP); w = chi/tgm**2
        return 1 + np.trapezoid(w*eCe_g/sig_g, tgm)/np.trapezoid(w*seff_g/sig_g, tgm), 1 + 2*b2*np.trapezoid(w/sig_g**2, tgm)/np.trapezoid(w, tgm)
    print(f"   echo tree: S° = {np.array2string(np.array(Tr['So']), precision=4)}, |b|^2 = sum_(log m < 2a) Lambda(m)^2/m = {float(bv @ bv):.4f}, t_lambda = {tl:.4f}")
    print(f"   kappa    M sharp      | smooth cutoff 1-exp(-(t/T)^{EXP}):   Z_chi       D_chi      prime sum (inside | all)     identity residual   M_chi = D/Z    1 + 2|b|^2/sigma(T)^2   tree at T   tree averaged   2|b|^2<sigma~^-2> averaged")
    for k in (1, 1.5, 2, math.e, 3, 4, 6, 10, 20):
        T = k*Ts
        # sharp
        m = zeros > T; Z = 2*F2z[m].sum() + 2*np.trapezoid((F2g*rho_g)[tg > zeros[-1]], tg[tg > zeros[-1]])          # tail beyond the last zero at the density
        D = 2*np.trapezoid((F2g*rho_g)[tg > T], tg[tg > T])
        # smooth
        chi_g = 1 - np.exp(-(tg/T)**EXP); chi_z = 1 - np.exp(-(zeros/T)**EXP)
        Zc = 2*(F2z*chi_z).sum() + 2*np.trapezoid((F2g*rho_g)[tg > zeros[-1]], tg[tg > zeros[-1]])
        Dc = 2*np.trapezoid(F2g*chi_g*Psi_g, tg)/(2*math.pi)
        g_chi = lambda xx: 2*np.trapezoid(F2g*chi_g*np.cos(tg*xx), tg)/(2*math.pi)
        pr_in = sum(2*c2[n]*g_chi(math.log(n)) for n in inside); pr_all = pr_in + sum(2*c2[n]*g_chi(math.log(n)) for n in entries if n not in inside)
        polar = -2*P*P*(1 - math.exp(-(2*T)**(-EXP)))
        resid = Zc - (Dc - pr_all + polar)
        sig = math.log(T/(2*math.pi)) - lam; pred = 1 + 2*sum(c2[n]**2 for n in inside)/sig**2
        print(f"   {k:5.3f}   {D/Z:8.4f}      |   {Zc:.5e}  {Dc:.5e}   {pr_in:+.3e} | {pr_all:+.3e}    {resid:+.2e} ({resid/Zc:+.4f} of Z)    {Dc/Zc:8.4f}        {pred:.4f}            {M_tree(sig):.4f}      {M_avg(T)[0]:.4f}          {M_avg(T)[1]:.4f}")
