# Analysis of the outputs of rh_mobius_green.py and rh_mobius_zeros.py (Sections 13-14 of mobius_modifier_v2.tex).
import json, numpy as np
g = json.load(open("green.json")); z = json.load(open("ngconst.json"))
returns = np.load("returns.npy"); xg, Mg = np.load("Mgrid.npy"); zd = np.load("zeros_dz.npy")
N = g["N"]

print("=== Green energy ===")
print(f"N={N}: M(N)={g['M_N']}, R(N)={g['R_N']:.5f}, R/logN={g['R_over_logN']:.5f}")
print("checkpoint   M(N)      R(N)        R/log N")
for c, (Mc, Rc, r) in sorted(g["checkpoints"].items(), key=lambda kv: int(kv[0])):
    print(f"{int(c):>10d} {Mc:>7d} {Rc:>10.5f} {r:>10.5f}")
print(f"C_M partial sum from {len(zd)} zeros (T={z['T']:.1f}): {z['S2']:.6f}; Gonek-tail estimate 2*(3/pi^3)/T = {2*3/np.pi**3/z['T']:.2e}")
print(f"Cramer sum partial: {z['S1']:.6f} vs 2+gamma-log 4pi = {2+0.5772156649-np.log(4*np.pi):.6f}")

print("\n=== Blocks (a0=0) ===")
print(f"completed blocks through N: {g['n_returns']}, sign changes of M: {g['sign_changes']}, max|M|: {g['max_abs_M']}")
print("return counts:", g["return_counts"])
cs = sorted(int(c) for c in g["return_counts"]); cnt = [g["return_counts"][str(c)] for c in cs]
lc, ln = np.log(cs[2:]), np.log(cnt[2:])
slope = np.polyfit(lc, ln, 1)[0]
print(f"log-log slope of #returns vs N over N>=1e4: {slope:.3f}  (random walk 0.5, heuristic 0.25)")
print("block length quantiles:", g["block_length_quantiles"])
print("longest blocks:")
for b in g["longest_blocks"]:
    print(f"  ({b['a']}, {b['b']}]  L={b['L']}  L/b={b['L']/b['b']:.3f}  H={b['H']}  H/sqrt(b)={b['H']/np.sqrt(b['b']):.3f}  sign={b['sign']}  E={b['E']:.4f}")

print("\n=== Explicit-formula reconstruction on the log grid ===")
u = np.log(xg.astype(np.float64)); y = Mg/np.sqrt(xg.astype(np.float64))
gam = zd[:,0]; dz = zd[:,1] + 1j*zd[:,2]; rho = 0.5 + 1j*gam
coef = 1/(rho*dz)
def recon(K, uu):
    r = np.zeros_like(uu)
    for kk in range(0, K, 500):
        blk = slice(kk, min(kk+500, K))
        r += 2*np.real(np.exp(1j*np.outer(uu, gam[blk])) @ coef[blk])
    return r - 2*np.exp(-uu/2)
for K in (100, 1000, len(gam)):
    sel = u >= np.log(1e4)
    r = recon(K, u[sel]); yy = y[sel]
    corr = np.corrcoef(r, yy)[0,1]; rms = np.sqrt(np.mean((r-yy)**2)); rmsy = np.sqrt(np.mean(yy**2))
    print(f"K={K:5d} zeros: corr={corr:.4f}  rms error={rms:.4f}  rms signal={rmsy:.4f}  mean sq of M/sqrt(x) on grid={np.mean(yy**2):.5f}")
# sign changes of M on the grid vs reconstruction (all zeros)
sel = u >= np.log(1e5); uu = u[sel]; yy = y[sel]; rr = recon(len(gam), uu)
def crossings(v, uu):
    s = np.sign(v); idx = np.nonzero(s[1:]*s[:-1] < 0)[0]; return uu[idx]
ca, cr = crossings(yy, uu), crossings(rr, uu)
print(f"sign changes of M for x>=1e5 on grid: {len(ca)}; of reconstruction: {len(cr)}")
for c in ca:
    j = np.argmin(np.abs(cr - c)); print(f"  M changes sign at x={np.exp(c):.4e} (u={c:.3f}); nearest reconstructed crossing at u={cr[j]:.3f}, |du|={abs(cr[j]-c):.3f}")
