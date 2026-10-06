# Relay transition at the entry of 3 versus the step size h at K=24 (edit `base` for finer steps).
# Step-size study of the relay jump at the entry of 3, K=24: h = 1.25e-4 and 6.25e-5 (plus reuse of coarser points).
import mpmath as mp, json
from rh_weil_odd import OddWeil
W = OddWeil(K=24); a3 = mp.log(3)/2
base = mp.mpf('0.0000625')
lam = {}; fa0 = None
for m in (-8, -4, -2, -1, 0, 1, 2, 4, 8):
    l, f, _ = W.floor(a3 + m*base); lam[m] = l
    if m == 0: fa0 = f
pred = 4*mp.log(3)/mp.sqrt(3)*fa0**2
print("f(a3) =", mp.nstr(fa0, 10), " predicted jump =", mp.nstr(pred, 8))
for s in (1, 2, 4, 8):
    h = s*base
    dm = (3*lam[0] - 4*lam[-s] + lam[-2*s])/(2*h) if -2*s in lam else None
    dp = (-3*lam[0] + 4*lam[s] - lam[2*s])/(2*h) if 2*s in lam else None
    if dm is not None and dp is not None:
        print(f"h={float(h):.3e}: jump={mp.nstr(dp-dm,8)}  ratio={mp.nstr((dp-dm)/pred,6)}", flush=True)
# first-order one-sided differences as a cross-check
for s in (1, 2, 4, 8):
    h = s*base
    dm1 = (lam[0]-lam[-s])/h; dp1 = (lam[s]-lam[0])/h
    print(f"first-order h={float(h):.3e}: ratio={mp.nstr((dp1-dm1)/pred,6)}")
json.dump({str(k): float(v) for k, v in lam.items()}, open("kinkh.json", "w"))
print("done")
