# Endpoint value f(a3) and floor lambda(a3) of the K-mode minimizer at the entry of 3, for larger K.
import mpmath as mp, json, time
from rh_weil_odd import OddWeil
a3 = mp.log(3)/2; out = {}
for K in (48, 56, 64):
    t0 = time.time(); W = OddWeil(K=K)
    lam, fa, E = W.floor(a3); Es = sorted(E)
    out[K] = {"lambda": float(lam), "f_a3": float(fa), "lambda1": float(Es[1]), "seconds": round(time.time()-t0)}
    print(K, json.dumps(out[K]), flush=True)
json.dump(out, open("endpointK.json", "w")); print("done")
