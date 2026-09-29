"""
Floor of the odd Weil form on a grid of supports, for the relay conjecture.

For each a: the three lowest eigenvalues of the K-mode form, the endpoint value f_K(a) of the normalized
minimizer, and the number of negative eigenvalues.  Grid step 0.0125 on [0.30, 1.20] plus the prime-power
entries a_n = (1/2) log n with a_n +- 0.004, so that the drop of Phi' = -lambda'/lambda at each entry is resolved.

Usage: python3 rh_floor_grid.py CHUNK NCHUNKS [K] [dps]   -> writes floor_grid_CHUNK.json
"""
import sys, json, math, time
import mpmath as mp
import rh_weil_odd as eng

chunk, nch = int(sys.argv[1]), int(sys.argv[2])
K = int(sys.argv[3]) if len(sys.argv) > 3 else 40
dps = int(sys.argv[4]) if len(sys.argv) > 4 else 60
mp.mp.dps = dps

grid = [round(0.30 + 0.0125*i, 6) for i in range(73)]
for n in (2, 3, 4, 5, 7, 8, 9, 11):
    an = math.log(n)/2
    grid += [round(an - 0.004, 8), round(an, 8), round(an + 0.004, 8)]
grid = sorted(set(grid))
mine = grid[chunk::nch]

W = eng.OddWeil(K=K)
out = []
t0 = time.time()
for a in mine:
    Q = W.matrix(mp.mpf(a))
    E, V = mp.eigsy(Q)
    idx = sorted(range(K), key=lambda i: E[i])
    k0 = idx[0]; c = V[:, k0]
    fa = sum(c[i]*W.N[i] for i in range(K))/mp.sqrt(a)
    row = {"a": a, "lam0": mp.nstr(E[idx[0]], 20), "lam1": mp.nstr(E[idx[1]], 12), "lam2": mp.nstr(E[idx[2]], 12),
           "fa": mp.nstr(fa, 12), "neg": sum(1 for e in E if e < 0)}
    out.append(row)
    print(f"a={a:.4f} lam0={row['lam0']} lam1={row['lam1']} f(a)={row['fa']} neg={row['neg']} ({time.time()-t0:.0f}s)", flush=True)
    json.dump({"K": K, "dps": dps, "rows": out}, open(f"floor_grid_{chunk}.json", "w"))
print("done", flush=True)
