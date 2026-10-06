"""
Bounded variation of the minimizers (paper Computation 5.8, memo section 2.8): the hypothesis (H_BV) behind the
Lipschitz continuity of the floor.

  python3 rh_bv_check.py fem [data/fem_*.json ...]
      sup|f|, total variation over the window, Var/(4 sup) and the Hardy integral H(d0) = int_{d0}^{a/10} f^2/delta
      between d0 = 1e-3, 1e-6, 1e-9, 1e-12 for edge-FEM minimizers (rh_edge_fem.py --out files; default: data/).
  python3 rh_bv_check.py kmode a K1,K2,... [bits]
      the same for the K-mode minimizers of rh_weil_arb.py (two steps of inverse iteration), evaluated on 44,000
      points of the half window graded to 1e-12 a; the variation includes the jump f_K(a) at the edge.
"""
import sys, os, glob, json, time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

def fem(files):
    print("file                              a       primes ||f||    sup|f|   Var(f)    Var/(4sup)  H(1e-3)  increments per three decades")
    for fn in files:
        d = json.load(open(fn)); dl = np.array(d["deltas"]); f = np.array(d["f"]); a = d["a"]
        x = a - dl; o = np.argsort(x); xs, fs = x[o], f[o]                     # centre to edge
        nrm = np.sqrt(2*np.trapezoid(fs*fs, xs))
        fv = np.concatenate([[0.0], fs, [0.0]]) if xs[0] > 0 else np.concatenate([fs, [0.0]])
        var = 2*np.abs(np.diff(fv)).sum(); sup = np.abs(fs).max()
        od = np.argsort(dl); ds, fd = dl[od], f[od]                           # edge outwards, delta increasing
        H = []
        for dm in (1e-3, 1e-6, 1e-9, 1e-12):
            m = (ds >= dm) & (ds <= 0.1*a)
            H.append(np.trapezoid(fd[m]**2/ds[m], ds[m]) if (m.sum() > 2 and ds[ds > 0].min() <= 2*dm) else float("nan"))
        inc = ", ".join(f"{H[i+1]-H[i]:.4f}" for i in range(3))
        print(f"{os.path.basename(fn):32s} {a:.5f} {d['primes']:6d} {nrm:.5f}  {sup:.6f} {var:.6f}  {var/(4*sup):.4f}    {H[0]:.4f}   {inc}")

def kmode(a, Ks, prec=500):
    from flint import arb, arb_mat, acb_mat, ctx
    from rh_weil_arb import OddWeilArb
    from numpy.polynomial.legendre import legval
    ctx.prec = prec; a_arb = arb(str(a))
    u = np.unique(np.concatenate([np.linspace(0, 1 - 1e-3, 40001), 1 - np.logspace(-3, -12, 4000)]))
    for K in Ks:
        t0 = time.time(); W = OddWeilArb(K=K, prec=prec); Q = W.matrix(a_arb)
        E = acb_mat(Q).eig(nonstop=True); E = sorted(E, key=lambda z: z.real.mid()); lam0 = E[0].real
        Qm = arb_mat(K, K, [arb(Q[i, j].mid()) for i in range(K) for j in range(K)])
        sigma = lam0.mid()*(1 - arb(2)**(-20))
        for i in range(K): Qm[i, i] = Qm[i, i] - sigma
        v = Qm.solve(arb_mat(K, 1, [arb(1) for _ in range(K)]), nonstop=True)
        v = Qm.solve(v, nonstop=True)
        nrm = sum(v[i, 0]*v[i, 0] for i in range(K)).sqrt(); c = [float((v[i, 0]/nrm).mid()) for i in range(K)]
        coef = np.zeros(2*K)
        for i in range(K): coef[2*i+1] = c[i]*np.sqrt((4*i+3)/2)
        f = legval(u, coef)/np.sqrt(a)
        fv = np.concatenate([f, [0.0]]); var = 2*np.abs(np.diff(fv)).sum(); sup = np.abs(f).max()
        print(f"a={a} K={K:3d} lambda_K={float(lam0.mid()):.4e} sup|f_K|={sup:.6f} f_K(a)={f[-1]:+.3e} Var f_K={var:.6f} Var/(4sup)={var/(4*sup):.4f} sup*Var={sup*var:.4f}  ({time.time()-t0:.0f}s)", flush=True)

if __name__ == "__main__":
    if len(sys.argv) < 2 or sys.argv[1] == "fem":
        files = sys.argv[2:] or sorted(glob.glob(os.path.join(os.path.dirname(os.path.abspath(__file__)), "data", "fem_*.json")))
        fem(files)
    elif sys.argv[1] == "kmode":
        kmode(float(sys.argv[2]), [int(k) for k in sys.argv[3].split(",")], int(sys.argv[4]) if len(sys.argv) > 4 else 500)
