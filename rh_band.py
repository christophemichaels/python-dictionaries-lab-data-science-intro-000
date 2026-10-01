"""
The sharp form inside the band of the echo lattice (paper Proposition 8.18, Computation 8.20).

Beyond a_inf = 0.843 the echo set is dense and the window's symbol sigma_a = Psi_inf - lambda - sum_{log m < 2a} 2 Lambda(m) m^{-1/2} cos(t log m)
is negative on a set of t reaching far above the horizon (its last sign change is at 33 T* for a = 1 and 259 T* for a = 1.25;
the periodic limit's band top, sigma_inf = 2 sum c_d, is at 47 and 412 T*).  The tail law of the zeros was nevertheless measured
there from ~5 T* (Computation 7.6).  This script measures, on the K-mode minimizers of data/tail_law_kmode (odd Legendre modes,
coefficients at 600-800 bits; F in closed form through spherical Bessel functions, int P_n(x/a) e^{itx} dx = 2a i^n j_n(ta),
evaluated with mpmath at 60 digits because |F| is 10^-13 of the coefficients), the local envelope of the sharp form
    <t^2 (Psi - lambda) |F|^2> / (-lambda')      over windows of width 4 T* centred at T = kappa T*,
together with the fraction of the window where sigma_a < 0, the parts of the envelope carried by the positive and the negative
set, and the test of the pointwise form |F|^2 ~ |A|^2/(t^2 |sigma_a|) (the window-symbol representation's edge piece, whose
modulus is set by |sigma_a|): the correlation of t^2 |F|^2 with 1/|sigma_a| and the envelope <t^2 |sigma_a| |F|^2>.

Usage: python3 rh_band.py [a ...]      (default 1.0 1.25 1.5)
"""
import json, math, sys, numpy as np, mpmath as mp

def vm(n):
    q = min(p for p in range(2, n + 1) if n % p == 0); m = n
    while m % q == 0: m //= q
    return math.log(q) if m == 1 else 0.0
def repsi(t):
    t = np.atleast_1d(np.asarray(t, float)); out = np.empty_like(t); sm = t < 60
    out[sm] = [float(mp.re(mp.psi(0, mp.mpf(1)/4 + 1j*mp.mpf(float(v))/2))) for v in t[sm]]
    z = 0.25 + 0.5j*t[~sm]; s = np.log(z) - 1/(2*z)
    for k, b in enumerate([1/6, -1/30, 1/42, -1/30, 5/66, -691/2730], 1): s -= b/(2*k*z**(2*k))
    out[~sm] = s.real; return out

def S_kmode(coef, a, ts, dps=80):
    """S(t) = int_0^a f sin(tx) dx = sqrt(a) sum_i c_i (-1)^i j_{2i+1}(ta) for f = sum c_i P_{2i+1}(x/a)/sqrt(a) (the K-mode normalization
    of rh_zero_side.py, the stored coefficients being the unit eigenvector in the orthonormalized basis: F = 2i S).  j_n by the upward recurrence at `dps` digits where z = ta > 2K (the second solution grows like
    exp(2n log(2n/(e z))) relative to j_n, within the working precision there), and by mpmath's Bessel function below."""
    K = len(coef); out = np.empty(len(ts)); nmax = 2*K - 1
    with mp.workdps(dps):
        c = [mp.mpf(x)*mp.sqrt(mp.mpf(4*i + 3)/2) for i, x in enumerate(coef)]; sa = mp.sqrt(mp.mpf(a))     # the stored c are the unit eigenvector; P_{2i+1} orthonormalized
        for k, t in enumerate(ts):
            z = mp.mpf(float(t))*mp.mpf(a)
            if z < 1.05*nmax:
                acc = sum(c[i]*(-1)**i*mp.sqrt(mp.pi/(2*z))*mp.besselj(2*i + mp.mpf(3)/2, z) for i in range(K))
            else:
                s_, co = mp.sin(z), mp.cos(z); j0 = s_/z; j1 = s_/z**2 - co/z
                acc = c[0]*j1; jm, jn = j0, j1
                for n in range(1, nmax):                      # j_{n+1} from j_n, j_{n-1}
                    jp = (2*n + 1)/z*jn - jm; jm, jn = jn, jp
                    if (n + 1) % 2 == 1:
                        i = n//2
                        if i < K: acc += c[i]*(-1)**i*jp
            out[k] = float(sa*acc)
    return out

if __name__ == "__main__":
    supports = sys.argv[1:] or ["1.0", "1.25", "1.5"]
    for a_s in supports:
        d = json.load(open(f"data/tail_law_kmode/coefs_{a_s}.json")); a = float(d["a"]); K = d["K"]; coef = d["c"]
        lam = float(d["lambda"]); lp = float(d["a_lambda_prime"])/a; Ts = 2*math.pi*math.exp(2*a)
        ent = [(math.log(m), vm(m)/math.sqrt(m), m) for m in range(2, 400) if vm(m) and math.log(m) < 2*a]
        S2c = 2*sum(c for _, c, _ in ent)
        print(f"a = {a}, K = {K}, lambda = {lam:.4e}, lambda' = {lp:.4e}, T* = {Ts:.2f}, entries {[m for _,_,m in ent]}, 2 sum c_d = {S2c:.3f} "
              f"(band top sigma_inf = 2 sum c at {2*math.pi*math.exp(S2c + lam)/Ts:.0f} T*)", flush=True)
        print("   T/T*   neg.frac   <t^2 (Psi-lam)|F|^2>/(-lam')   on sigma_a>0   on sigma_a<0   <t^2|sigma_a||F|^2>/(-lam')   corr(t^2|F|^2, 1/|sigma_a|)   <t^2 sigma_inf |F|^2>/(-lam')   1-2 neg.frac")
        for kap in (2, 3, 5, 7, 10, 15, 20, 30, 50, 100):
            T = kap*Ts; step = 0.05 if T*a > 1.05*(2*K - 1) else 0.25; tt = np.arange(max(T - 2*Ts, 0.6*T), T + 2*Ts, step)      # a window of width 4 T* (narrower on the left at kappa = 2)
            S = S_kmode(coef, a, tt); F2 = 4*S*S
            Psi_inf = repsi(tt) - math.log(math.pi); osc = sum(2*c*np.cos(tt*dd) for dd, c, _ in ent)
            sig = Psi_inf - osc - lam
            w = tt*tt*F2
            env = np.mean(w*sig)/(-lp); pos = np.mean(w*sig*(sig > 0))/(-lp); neg = np.mean(w*sig*(sig < 0))/(-lp)
            envabs = np.mean(w*np.abs(sig))/(-lp); envinf = np.mean(w*(Psi_inf - lam))/(-lp); nf = (sig < 0).mean()
            corr = np.corrcoef(w, 1/np.abs(sig))[0, 1] if np.abs(sig).min() > 0 else float('nan')
            print(f"   {kap:4d}   {nf:7.4f}   {env:12.4f}            {pos:9.4f}      {neg:9.4f}        {envabs:12.4f}                {corr:8.3f}                 {envinf:10.4f}            {1-2*nf:8.4f}", flush=True)
