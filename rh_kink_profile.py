# Profile of the K-mode minimizer near the endpoint at the entry of 3: is there a shrinking boundary layer?
import mpmath as mp, json, time
from rh_weil_odd import OddWeil, legendre_all
a3 = mp.log(3)/2
us = ['0.5', '0.8', '0.9', '0.95', '0.98', '0.99', '0.995', '0.999', '1.0']
out = {}
for K in (24, 40, 64):
    t0 = time.time(); W = OddWeil(K=K)
    Q = W.matrix(a3); E, V = mp.eigsy(Q)
    k0 = min(range(K), key=lambda i: E[i]); c = V[:, k0]
    sgn = 1 if sum(c[i]*W.N[i] for i in range(K)) > 0 else -1     # fix the sign so that f(a) > 0
    prof = {}
    for u in us:
        P = legendre_all(mp.mpf(u), 2*K-1)
        prof[u] = float(sgn*sum(c[i]*W.N[i]*P[2*i+1] for i in range(K))/mp.sqrt(a3))
    out[K] = prof
    print(K, json.dumps(prof), f"({time.time()-t0:.0f}s)", flush=True)
json.dump(out, open("profileK.json", "w"))
print("done")
