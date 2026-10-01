"""
The floor on the plane (paper Theorem 9.3, Corollary 9.4, Computation 9.5).

On the zero side the floor is lambda = sum_rho q(gamma_rho) with q(z) = F(z)F(-z) = -F(z)^2 = 4 S(z)^2, S(z) = int_0^a f sin(zu) du, the
quadratic functional of the minimizer's transform evaluated at every zero rho = 1/2 + i gamma_rho of zeta, gamma_rho complex when the zero
is off the line (a quartet rho, 1-rho, conj: gamma = +-gamma_0 +- i delta contributes -4 Re F(gamma_0 - i delta)^2 = 8 Re S(gamma_0 - i delta)^2).
On the line q >= 0; off it, to second order in the displacement y, q(x + iy) = 4(S(x) + iy S'(x))^2 has real part 4 S^2 - 4 y^2 S'^2,
negative near the nulls of S, which is where the zeros are.  This script computes, on the K-mode minimizers:
  (i)  the stiffness of the floor against the zeros leaving the line, kappa_2 = sum_(gamma>0) 16 S'(gamma)^2 (per unit y^2), and the critical
       displacement y_c = (lambda/kappa_2)^{1/2} at which a uniform displacement of the zeros would annul the floor to second order;
  (ii) the growth of |q| off the line, r(y) = <|q(x+iy)|>/<q(x)> over windows of heights, against the pure-edge prediction cosh(2ya):
       above the horizon the transform is the edge's cosine of frequency a and r follows cosh(2ya); in the mass region the transform
       hides from the zeros (dense nulls, S' >> a S) and r is larger;
  (iii) the fraction of the heights in [T*, 5 T*] at which Re q(x + iy) < 0 for y = 0.05 ... 0.45.

Usage: python3 rh_plane.py [a ...]      (default 0.6 0.8 1.0)
"""
import json, math, sys, numpy as np, mpmath as mp

def S_complex(coef, a, zs, dps=60):
    """S(z) for complex z (the K-mode transform, as in rh_band.py), by the Bessel series below 1.05(2K-1) and the recurrence above."""
    K = len(coef); nmax = 2*K - 1; out = []
    with mp.workdps(dps):
        c = [mp.mpf(x)*mp.sqrt(mp.mpf(4*i + 3)/2) for i, x in enumerate(coef)]; sa = mp.sqrt(mp.mpf(a))
        for z in zs:
            zz = mp.mpc(z.real, z.imag)*mp.mpf(a)
            if abs(zz) < 1.05*nmax:
                acc = sum(c[i]*(-1)**i*mp.sqrt(mp.pi/(2*zz))*mp.besselj(2*i + mp.mpf(3)/2, zz) for i in range(K))
            else:
                s_, co = mp.sin(zz), mp.cos(zz); j0 = s_/zz; j1 = s_/zz**2 - co/zz
                acc = c[0]*j1; jm, jn = j0, j1
                for n in range(1, nmax):
                    jp = (2*n + 1)/zz*jn - jm; jm, jn = jn, jp
                    if (n + 1) % 2 == 1:
                        i = n//2
                        if i < K: acc += c[i]*(-1)**i*jp
            out.append(complex(sa*acc))
    return np.array(out)

if __name__ == "__main__":
    supports = sys.argv[1:] or ["0.6", "0.8", "1.0"]
    zeros = np.array([float(l.split()[1]) for l in open("data/zeros_6000.txt")])
    for a_s in supports:
        d = json.load(open(f"data/tail_law_kmode/coefs_{a_s}.json")); a = float(d["a"]); K = d["K"]; coef = d["c"]
        lam = float(d["lambda"]); lp = float(d["a_lambda_prime"])/a; Ts = 2*math.pi*math.exp(2*a)
        # (i) stiffness: S'(gamma) by a central difference, zeros up to 6000 (the tail beyond at the density is small for the derivative sum)
        h = 1e-4
        Sp = S_complex(coef, a, zeros + h); Sm = S_complex(coef, a, zeros - h); S0 = S_complex(coef, a, zeros + 0j)
        dS = (Sp - Sm).real/(2*h); k2 = float(np.sum(16*dS*dS)); lam_z = float(np.sum(4*S0.real**2))*2      # a quartet at gamma -+ iy replaces the pair's 8 S^2 by -16 y^2 S'^2
        yc = math.sqrt(lam/k2) if k2 > 0 else float('inf')
        # where the stiffness sits: shares by height
        shares = [(T, float(np.sum(16*dS[zeros <= T]**2))/k2) for T in (Ts, 2*Ts, 3*Ts, 5*Ts, 10*Ts)]
        print(f"a = {a}, K = {K}, lambda = {lam:.4e}, lambda' = {lp:.4e}, T* = {Ts:.2f}; 2 sum_6000 4S^2 = {lam_z:.4e} ({lam_z/lam:.4f} of lambda)", flush=True)
        print(f"   stiffness kappa_2 = sum_(gamma>0) 16 S'(gamma)^2 = {k2:.4e};  kappa_2/lambda = {k2/lam:.3e};  critical displacement y_c = (lambda/kappa_2)^(1/2) = {yc:.4e}"
              f"  (a y_c = {a*yc:.4f});  shares of kappa_2 below " + ", ".join(f"{T/Ts:.0f}T*: {s:.3f}" for T, s in shares), flush=True)
        # (ii) growth off the line: r(y) = <|q(x+iy)|>/<q(x)> over x in a window, against the pure-edge prediction cosh(2ya)
        #      (a cosine of frequency a in x grows like cosh(2ya) in modulus squared); above the horizon the transform is the edge's and
        #      r follows cosh(2ya); in the mass region the transform hides from the zeros (dense nulls, S' >> a S) and r is larger
        print("   (ii) growth off the line: r(y) = <|q(x+iy)|>/<q(x)> over x in [kT*, (k+1)T*] (step 0.5) against cosh(2ya) [in brackets], and the fraction of x with Re q(x+iy) < 0")
        for k in (1, 2, 3, 5, 10):
            xs = np.arange(k*Ts, (k + 1)*Ts, 0.5); row = []
            S_line = S_complex(coef, a, xs + 0j); q_line = 4*S_line.real**2
            for y in (0.05, 0.1, 0.2, 0.3, 0.45):
                S_off = S_complex(coef, a, xs + 1j*y); q_off = 4*S_off**2
                ratio = np.mean(np.abs(q_off))/np.mean(q_line); neg = float(np.mean(q_off.real < 0))
                row.append(f"y={y}: r={ratio:7.4f} [{math.cosh(2*y*a):6.4f}] neg {neg:.2f}")
            print(f"      x in [{k}T*, {k+1}T*]:  " + "  ".join(row), flush=True)
