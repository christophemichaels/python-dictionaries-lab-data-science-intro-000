# Relay transition at the entry of 3 versus the mode count K (24, 32, 40), one-sided second-order differences and Richardson extrapolation.
# Relay transition at the entry of 3, retaining the changing minimizer, as a function of the mode count K.
import mpmath as mp, json, time
from rh_weil_odd import OddWeil
a3 = mp.log(3)/2; h = mp.mpf('0.00025')
out = {}
for K in (24, 32, 40):
    W = OddWeil(K=K); t0 = time.time()
    lam = {}; fa = {}
    for m in (-4, -2, -1, 0, 1, 2, 4):
        l, f, _ = W.floor(a3 + m*h); lam[m] = l; fa[m] = f
    def dm(s):   # one-sided second-order derivatives with step s*h (s=1 or 2)
        L = lambda k: lam[k*s]
        return ((3*L(0) - 4*L(-1) + L(-2))/(2*s*h), (-3*L(0) + 4*L(1) - L(2))/(2*s*h))
    dm1, dp1 = dm(1); dm2, dp2 = dm(2)
    # Richardson (one-sided second-order -> fourth-order)
    dmR = (4*dm1 - dm2)/3; dpR = (4*dp1 - dp2)/3
    pred = 4*mp.log(3)/mp.sqrt(3)*fa[0]**2
    rec = {"lambda_a3": float(lam[0]), "f_a3": float(fa[0]), "pred_jump": float(pred),
           "jump_h": float(dp1-dm1), "jump_2h": float(dp2-dm2), "jump_richardson": float(dpR-dmR),
           "ratio_h": float((dp1-dm1)/pred), "ratio_2h": float((dp2-dm2)/pred), "ratio_richardson": float((dpR-dmR)/pred),
           "dlam_minus": float(dmR), "dlam_plus": float(dpR), "seconds": round(time.time()-t0)}
    out[K] = rec
    print(K, json.dumps(rec), flush=True)
json.dump(out, open("kinkK.json","w"), indent=1)
print("done")
