"""
Soft kink at the entry of 3.

Finite differences of the edge-FEM floor across a_3 = (1/2) log 3 at scales eps = 1e-3 ... 1e-7:
    D(eps) = [lambda(a_3+eps) - 2 lambda(a_3) + lambda(a_3-eps)] / eps,
against the first-order entering energy on the a_3 minimizer,
    P(eps)/eps,  P(eps) = 2 Lambda(3) 3^{-1/2} int_0^{2 eps} f(a-d) f(a-(2eps-d)) dd,
and the asymptotic law K/(log(1/eps) + beta), K = 4 Lambda(3) 3^{-1/2} C^2, from the edge fit
    f(a-d) = C (log(a/d) + beta)^{-1/2}.
A hard kink would give D(eps) -> const; the exact form gives D(eps) ~ K/log(1/eps).

Usage: python3 rh_soft_kink.py N_INT [two] [eps1,eps2,...]     ("two" selects the two-stage graded mesh)
"""
import sys, os, json, math, time, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from rh_edge_fem import EdgeFEM, mesh
a3 = math.log(3)/2
n_int, rho, r, dmin = int(sys.argv[1]), 0.1, 0.5, 1e-11
r1, dsw = (0.9, 1e-4) if len(sys.argv) > 2 and sys.argv[2] == 'two' else (None, None)
if r1: rho = 0.05
def solve(a):
    fem = EdgeFEM(a, mesh(a, n_int, rho, r, dmin, r1, dsw), primes=True); fem.assemble(); lam, f = fem.solve(); return lam, fem
t0 = time.time()
lam0, fem0 = solve(a3)
print(f"lambda(a3) = {lam0:.12e}  nodes={fem0.n} ({time.time()-t0:.0f}s)", flush=True)
# first-order entering energy from the a3 minimizer: P(eps) = 4 Lam(3) 3^{-1/2} int_0^{2eps} f(a-d) f(a-(2eps-d)) dd
d = fem0.d[1:]; f = fem0.f                # d decreasing to 0 (edge), nodal values
fe = lambda dd: np.interp(-dd, -d, f)     # f at distance dd from the edge (linear interpolation on the mesh)
lam3 = math.log(3); K3 = 4*lam3/math.sqrt(3)
def P(eps):
    x = np.linspace(0, 2*eps, 20001)
    return (K3/2)*np.trapezoid(fe(x)*fe(2*eps-x), x)   # P = 2 Lam(3) 3^{-1/2} int_0^{2eps} f f
# fit C, beta on 8 < t < 20
t = np.log(a3/d[(d>0)&(d<rho*a3)]); ff = f[(d>0)&(d<rho*a3)]
m = (t>8)&(t<20); sl, ic = np.polyfit(t[m], 1/ff[m]**2, 1); C = 1/math.sqrt(sl); beta = ic/sl
print(f"edge fit: 1/f^2 = {sl:.4f} (t + {beta:.3f})  ->  C = {C:.5f},  K = 4 Lam(3) 3^-1/2 C^2 = {K3*C*C:.4e}", flush=True)
out = {"lambda0": lam0, "C": C, "beta": beta, "rows": []}
for eps in ([1e-3, 1e-4, 1e-5, 1e-6, 1e-7] if len(sys.argv) < 4 else [float(e) for e in sys.argv[3].split(",")]):
    lp, _ = solve(a3+eps); lm, _ = solve(a3-eps)
    D = (lp - 2*lam0 + lm)/eps
    Pe = P(eps)/eps; L = math.log(1/eps)
    row = {"eps": eps, "lam+": lp, "lam-": lm, "D": D, "P/eps": Pe, "K/L": K3*C*C/L, "K/(L+beta)": K3*C*C/(L+beta)}
    out["rows"].append(row)
    print(f"eps={eps:.0e}  lam+={lp:.10e} lam-={lm:.10e}  D(eps)={D:.4e}  P/eps={Pe:.4e}  K/L={K3*C*C/L:.4e}  K/(L+beta)={K3*C*C/(L+beta):.4e}  ({time.time()-t0:.0f}s)", flush=True)
json.dump(out, open(f"softkink_{n_int}{sys.argv[2] if len(sys.argv)>2 else ''}.json", "w"))
