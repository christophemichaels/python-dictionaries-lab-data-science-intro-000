"""
The tail law on the zero side (paper Section 7.1, Computation 7.6): the fraction of the floor carried by the zeros above
height T is theta(T) = Phi'(a)/(pi T) beyond a few horizons, i.e. Phi' = pi T_eff with T_eff = lim T theta(T).

  python3 rh_tail_law.py fem data/fem_a0.5_full.json data/zeros_6000.txt LAMBDA_PRIME
      edge-FEM minimizer on the zeros: F in closed form for the odd piecewise-linear f (jump at the edge included),
      2 sum |F(gamma)|^2 against lambda_FEM, theta(T) at fixed T/T* against Phi'/(pi T), median height of the mass.
      LAMBDA_PRIME is the exact derivative at the support (from data/dilation_*.log: a lambda' / a).
  python3 rh_tail_law.py coefs a K bits out.json
      eigenvector of the K-mode form (two inverse-iteration steps), lambda_K and a lambda_K', saved as decimal strings.
  python3 rh_tail_law.py shares coefs.json zeros.txt n0 n1 dps out.csv
      per-zero shares 2|F_K(gamma)|^2/lambda_K for zeros n0..n1 (F_K through spherical Bessel functions at dps digits;
      the working precision must exceed the digits lost in the cancellation to |F_K(gamma)| ~ sqrt(lambda_K), and the
      zeros must be accurate to |F_K(gamma)|/|F_K'(gamma)| in the plunge region: 36 digits at a = 1, 50 at a = 1.5).
  python3 rh_tail_law.py shares_arb coefs.json zeros.txt n0 n1 bits out.csv
      the same shares with arb: j_n(a gamma) for all n < 2K by Miller's backward recurrence (started above
      max(2K, a gamma (1+eps)) and normalized by sum (2n+1) j_n^2 = 1; midpoints, since the ball radii of an
      oscillatory recurrence grow geometrically), the sum over the modes at `bits` bits; bits must exceed the
      cancellation to |F_K(gamma)| ~ sqrt(lambda_K) (600 bits at a = 2).  Not rigorous (nor is the eigenvector).
  python3 rh_tail_law.py overlay DIR
      theta(T/T*) for every coefs_*.json / shares_*.csv pair in DIR: the universality overlay and the medians.
"""
import sys, os, json, glob, math
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

def fem(argv):
    import numpy as np
    fn, zfile, lam_prime = argv[0], argv[1], float(argv[2])
    d = json.load(open(fn)); a = d["a"]; lam = d["lambda"]
    dl = np.array(d["deltas"]); f = np.array(d["f"]); x = a - dl; o = np.argsort(x); x, f = x[o], f[o]      # centre -> edge
    if x[0] > 0: x = np.concatenate([[0.0], x]); f = np.concatenate([[0.0], f])
    B = np.diff(f)/np.diff(x)                                     # slopes on the segments
    zeros = np.array([float(l.split()[1]) for l in open(zfile)])
    def S(t):  # int_0^a f(x) sin(tx) dx for piecewise-linear f, f(0)=0, jump f(a) at the edge
        t = t[:, None]
        return (-f[-1]*np.cos(t[:, 0]*a)/t[:, 0]) + (B[None, :]*(np.sin(t*x[1:][None, :]) - np.sin(t*x[:-1][None, :]))).sum(1)/t[:, 0]**2
    F2 = 4*S(zeros)**2                                            # |F|^2 = |2i S|^2
    Ts = 2*np.pi*np.exp(2*a); Phip = -lam_prime/lam
    share = 2*F2/lam; cum = np.cumsum(share)
    C2 = -lam_prime/2
    print(f"{fn.split('/')[-1]}: a={a:.5f} lambda_FEM={lam:.8e} lambda'={lam_prime:.6e} Phi'/T*={Phip/Ts:.4f} f(a)={f[-1]:.5f} C^2=-lambda'/2={C2:.4e}")
    print(f"  2 sum_{{6000 zeros}} |F|^2 / lambda = {cum[-1]:.6f}   (tail beyond gamma_6000 by the law Phi'/(pi T): {Phip/(np.pi*zeros[-1]):.6f})  ->  total {cum[-1] + Phip/(np.pi*zeros[-1]):.6f}")
    nb = int((zeros < Ts).sum()); print(f"  zeros below T*={Ts:.3f}: {nb}, share {cum[nb-1] if nb else 0:.3e}")
    print("   T/T*   theta(T)   T theta/T*   Phi'/(pi T*)=%.4f   ratio" % (Phip/(np.pi*Ts)))
    for k in (1, 1.5, 2, 2.5, 3, 4, 5, 7, 10, 15, 20, 30, 50, 100):
        i = np.searchsorted(zeros, k*Ts); th = 1 - cum[i-1] if i > 0 else 1.0
        print(f"  {k:5.1f}   {th:.5f}   {th*k:.4f}                          {th*k/(Phip/(np.pi*Ts)):.4f}")
    i50 = np.searchsorted(cum, 0.5); print(f"  median height of the floor's mass: gamma_{i50+1} = {zeros[i50]:.3f} = {zeros[i50]/Ts:.3f} T*   (e = 2.718)")

def coefs(argv):
    import time
    from flint import arb, arb_mat, acb_mat, ctx
    from rh_weil_arb import OddWeilArb
    a_s, K, prec, out = argv[0], int(argv[1]), int(argv[2]), argv[3]
    ctx.prec = prec; a_arb = arb(a_s); t0 = time.time()
    W = OddWeilArb(K=K, prec=prec); Q = W.matrix(a_arb)
    E = acb_mat(Q).eig(nonstop=True); E = sorted(E, key=lambda z: z.real.mid()); lam0 = E[0].real
    Qm = arb_mat(K, K, [arb(Q[i, j].mid()) for i in range(K) for j in range(K)])
    sigma = lam0.mid()*(1 - arb(2)**(-20))
    for i in range(K): Qm[i, i] = Qm[i, i] - sigma
    v = Qm.solve(arb_mat(K, 1, [arb(1) for _ in range(K)]), nonstop=True); v = Qm.solve(v, nonstop=True)
    nrm = sum(v[i, 0]*v[i, 0] for i in range(K)).sqrt(); c = arb_mat(K, 1, [v[i, 0]/nrm for i in range(K)])
    def quad(X): return (c.transpose()*X*c)[0, 0]
    h = arb(10)**(-30); Ap, Pp, sp = W.parts(a_arb + h); Am, Pm, sm = W.parts(a_arb - h)
    lam_prime = quad((Ap - Am)/(2*h)) + sum(quad((Pp[n] - Pm[n])/(2*h)) for n in Pp) - 2*quad((sp*sp.transpose() - sm*sm.transpose())/(2*h))
    json.dump({"a": a_s, "K": K, "prec": prec, "lambda": lam0.mid().str(60, radius=False), "a_lambda_prime": (a_arb*lam_prime).mid().str(60, radius=False),
               "c": [c[i, 0].mid().str(120, radius=False) for i in range(K)]}, open(out, "w"))
    print(f"a={a_s} K={K}: lambda_K={lam0.str(10)} a lambda'={(a_arb*lam_prime).str(10)} ({time.time()-t0:.0f}s)")

def shares(argv):
    import mpmath as mp
    cf, zfile, n0, n1, dps, out = argv[0], argv[1], int(argv[2]), int(argv[3]), int(argv[4]), argv[5]
    d = json.load(open(cf)); mp.mp.dps = dps
    a = mp.mpf(d["a"]); K = d["K"]; lamK = mp.mpf(d["lambda"])
    coef = [mp.mpf(d["c"][i])*mp.sqrt(mp.mpf(4*i+3)/2) for i in range(K)]
    def G(t):
        z = t*a; return 2*mp.sqrt(a)*sum(coef[i]*(-1)**i*mp.sqrt(mp.pi/(2*z))*mp.besselj(2*i + mp.mpf(3)/2, z) for i in range(K))
    zeros = {}
    for l in open(zfile):
        p = l.split()
        if len(p) == 2 and n0 <= int(p[0]) <= n1: zeros[int(p[0])] = mp.mpf(p[1])
    with open(out, "w") as fh:
        for n in sorted(zeros):
            g = zeros[n]; fh.write(f"{n},{mp.nstr(g, 20)},{mp.nstr(2*G(g)**2/lamK, 12)}\n"); fh.flush()

def shares_arb(argv):
    from flint import arb, ctx
    cf, zfile, n0, n1, bits, out = argv[0], argv[1], int(argv[2]), int(argv[3]), int(argv[4]), argv[5]
    d = json.load(open(cf)); ctx.prec = bits
    a = arb(d["a"]); K = d["K"]; lamK = arb(d["lambda"])
    coef = [arb(d["c"][i])*(arb(4*i + 3)/2).sqrt()*(1 if i % 2 == 0 else -1) for i in range(K)]
    ra2 = 2*a.sqrt(); digits = bits*0.30103
    def spherical(z, nmax):
        """j_0..j_nmax at z > 0 by Miller's backward recurrence."""
        zf = float(z.mid())
        eps = 0.5*(3.45*digits/zf)**(2.0/3) if zf > 0 else 10.0
        N = max(nmax + 1, int(zf*(1 + eps))) + 20
        jp, jc = arb(0), arb(2)**(-bits)                      # j_{N+1}, j_N (tiny seed)
        vals = [None]*(N + 1); vals[N] = jc
        zm = arb(z.mid())
        for n in range(N, 0, -1):                             # floating point (midpoints): the ball radii of an
            jm = arb((jc*(2*n + 1)/zm - jp).mid())            # oscillatory three-term recurrence grow like 3^N
            jp, jc = jc, jm
            vals[n - 1] = jc
        norm = sum(((2*n + 1)*vals[n]*vals[n] for n in range(N + 1)), arb(0)).sqrt()
        sgn = 1 if (vals[0]*(z.sin()/z)) > 0 else -1
        return [v*sgn/norm for v in vals[:nmax + 1]]
    zeros = {}
    for l in open(zfile):
        p = l.split()
        if len(p) == 2 and n0 <= int(p[0]) <= n1: zeros[int(p[0])] = arb(p[1])
    with open(out, "w") as fh:
        for n in sorted(zeros):
            g = zeros[n]; z = a*g
            j = spherical(z, 2*K - 1)
            F = ra2*sum((coef[i]*j[2*i + 1] for i in range(K)), arb(0))
            share = 2*F*F/lamK
            fh.write(f"{n},{g.mid().str(20, radius=False)},{share.mid().str(12, radius=False)}\n"); fh.flush()

def overlay(argv):
    rows = []
    for cf in sorted(glob.glob(argv[0] + "/coefs_*.json")):
        d = json.load(open(cf)); a = float(d["a"]); lam = float(d["lambda"]); alp = float(d["a_lambda_prime"]); Phip = -alp/a/lam
        Ts = 2*math.pi*math.exp(2*a); pts = []
        for f in sorted(glob.glob(cf.replace("coefs_", "shares_").replace(".json", "_*.csv"))):
            for l in open(f):
                n, g, s = l.split(","); pts.append((int(n), float(g), float(s)))
        if not pts: continue
        pts.sort(); cum = 0; th = {}; med = None; grid = [1, 1.5, 2, 2.5, 3, 4, 5, 7, 10, 15]
        gi = 0
        for n, g, s in pts:
            while gi < len(grid) and g >= grid[gi]*Ts: th[grid[gi]] = 1 - cum; gi += 1
            cum += s
            if med is None and cum >= 0.5: med = g/Ts
        tail = Phip/(math.pi*Ts)
        print(f"a={a:<5} K={d['K']:<4} zeros to {pts[-1][1]/Ts:5.1f} T*  Phi'/T*={Phip/Ts:.4f}  Phi'/(pi T*)={tail:.4f}  cum(all)={cum:.4f}  median={med if med else float('nan'):.3f} T*")
        print("   T/T*:   " + "  ".join(f"{k:>6}" for k in grid if k in th))
        print("   theta:  " + "  ".join(f"{th[k]:6.4f}" for k in grid if k in th))
        print("   Tth/T*: " + "  ".join(f"{th[k]*k:6.3f}" for k in grid if k in th))
        print("   ratio:  " + "  ".join(f"{th[k]*k/tail:6.3f}" for k in grid if k in th))

if __name__ == "__main__":
    {"fem": fem, "coefs": coefs, "shares": shares, "shares_arb": shares_arb, "overlay": overlay}[sys.argv[1]](sys.argv[2:])
