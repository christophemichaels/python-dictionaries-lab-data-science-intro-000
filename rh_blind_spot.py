"""
The blind spot of the window (paper Theorem 9.31, Computation 9.32), and the variation of the K-mode minimizers (Proposition 9.28).

The quartet of an off-line zero rho_0 = 1/2 + delta + i gamma contributes 16 Re S(gamma - i delta)^2 to the floor, against 16 S(gamma)^2
for a double zero on the line at the same height: the difference is -16 delta^2 (S'(gamma)^2 + S(gamma) S''(gamma)) + O(delta^4), even in
delta.  So the floor at support a is a function of the configuration of the zeros below its horizon with the single-zero sensitivity
  delta_c(gamma) = (lambda |f|^2 / (16 S'(gamma)^2))^{1/2}
(the displacement of a double zero at a null of the minimizer's transform that would annul the floor), tiny where the minimizer hides
from the zeros (the mass region, below a few horizons) and above 1/2 where the transform is the edge's cosine (above five horizons):
the window resolves the zeros below its horizon to e^{-cT*/2} and is blind above it.  This script prints, on the K-mode minimizers
(unit norm, transform through spherical Bessel functions at 60 digits, rh_plane.py), (i) the second-order formula against the exact
quartet difference at the first zero and at a zero near five horizons, (ii) delta_c at the first zeros and at zeros near 1/2, 1, 2, 3,
5, 10, 20 horizons, and (iii) sup|f| and the total variation of the minimizer (the bound S(f,f) <= M_a^2 |f'|_1^2 / pi on the Suzuki norm).

Usage: python3 rh_blind_spot.py [a ...]      (default 0.6 0.8 1.0)
"""
import json, math, sys, numpy as np
from numpy.polynomial import legendre as L
from rh_plane import S_complex

if __name__ == "__main__":
    supports = sys.argv[1:] or ["0.6", "0.8", "1.0"]
    zeros = np.array([float(l.split()[1]) for l in open("data/zeros_6000.txt")])
    for a_s in supports:
        d = json.load(open(f"data/tail_law_kmode/coefs_{a_s}.json")); a = float(d["a"]); K = d["K"]; coef = d["c"]; coefF = [float(c) for c in coef]                                     # strings: full precision for the transform
        lam = float(d["lambda"]); lp = float(d["a_lambda_prime"])/a; Ts = 2*math.pi*math.exp(2*a)
        # (iii) the minimizer in position space: f(u) = sum_i c_i sqrt((4i+3)/(2a)) P_{2i+1}(u/a), unit norm
        cl = np.zeros(2*K); cl[1::2] = [c*math.sqrt((4*i + 3)/(2*a)) for i, c in enumerate(coefF)]
        u = np.linspace(0, a, 200001); fu = L.legval(u/a, cl); norm2 = 2*np.trapezoid(fu*fu, u)
        tv = 2*float(np.sum(np.abs(np.diff(fu)))); sup = float(np.abs(fu).max()); argmax = u[np.argmax(np.abs(fu))]
        print(f"a = {a}, K = {K}, lambda = {lam:.4e}, lambda' = {lp:.4e}, T* = {Ts:.2f}; |f|^2 = {norm2:.6f}, sup|f| = {sup:.4f} at u = {argmax:.4f} = a - {a - argmax:.4f}, "
              f"total variation {tv:.4f} = {tv/sup:.3f} sup|f|", flush=True)
        # (i) the second-order formula at the first zero and near 5 T*
        def Sv(z): return S_complex(coef, a, np.array([z]))[0]
        for g in (zeros[0], zeros[np.argmin(np.abs(zeros - 5*Ts))]):
            h = 1e-4; S0 = Sv(g + 0j).real; Sp = (Sv(g + h + 0j).real - Sv(g - h + 0j).real)/(2*h); Spp = (Sv(g + h + 0j).real - 2*S0 + Sv(g - h + 0j).real)/h**2
            row = []
            for dl in (1e-3, 1e-2, 1e-1):
                exact = 16*(Sv(g - 1j*dl)**2).real - 16*S0*S0; second = -16*dl*dl*(Sp*Sp + S0*Spp)
                row.append(f"delta={dl:g}: exact {exact:+.4e}, second order {second:+.4e}, ratio {exact/second if second else float('nan'):.4f}")
            print(f"   (i) gamma = {g:.3f} ({g/Ts:.2f} T*): S = {S0:+.3e}, S' = {Sp:+.3e}, S'' = {Spp:+.3e};  " + ";  ".join(row), flush=True)
        # (ii) the single-zero sensitivity delta_c(gamma) = (lambda/(16 S'^2))^{1/2}
        picks = list(range(5)) + [int(np.argmin(np.abs(zeros - k*Ts))) for k in (0.5, 1, 2, 3, 5, 10, 20) if k*Ts < zeros[-1]]
        print("   (ii)  gamma      gamma/T*     S(gamma)        S'(gamma)       16 S'^2 / lambda     delta_c = (lambda/16S'^2)^(1/2)")
        for j in picks:
            g = zeros[j]; h = 1e-4; S0 = Sv(g + 0j).real; Sp = (Sv(g + h + 0j).real - Sv(g - h + 0j).real)/(2*h)
            dc = math.sqrt(lam/(16*Sp*Sp)) if Sp != 0 else float('inf')
            print(f"        {g:9.3f}   {g/Ts:8.3f}   {S0:+.4e}   {Sp:+.4e}   {16*Sp*Sp/lam:14.4e}   {dc:12.4e}", flush=True)
