"""
Floors of the odd Weil form beyond a = 1 with the arb engine, for the drift test of the relay conjecture.

Jobs are (a, K, prec) triples; Phi' at a centre c is taken from the pair c +- 0.005.  Centres avoid the
prime-power entries a_n = (1/2) log n.  Convergence checks repeat two centres at a larger K.

Usage: python3 rh_arb_grid.py CHUNK NCHUNKS   -> appends JSON lines to arb_grid_CHUNK.jsonl
"""
import sys, json, time
from rh_weil_arb import OddWeilArb

chunk, nch = int(sys.argv[1]), int(sys.argv[2])
plan = []
for c, K in ((1.06, 110), (1.10, 110), (1.15, 110), (1.25, 140), (1.33, 140), (1.40, 180), (1.45, 180), (1.50, 180)):
    plan += [(round(c - 0.005, 6), K, 500), (round(c + 0.005, 6), K, 500)]
plan += [(1.245, 180, 500), (1.255, 180, 500), (1.495, 230, 560), (1.505, 230, 560)]   # convergence checks
plan.sort(key=lambda t: t[1]*t[1]*t[1])                                                 # cheap first
mine = plan[chunk::nch]
engines = {}
import os
done = set()
if os.path.exists(f"arb_grid_{chunk}.jsonl"):
    for l in open(f"arb_grid_{chunk}.jsonl"):
        if l.strip(): r = json.loads(l); done.add((round(r["a"], 6), r["K"]))
for a, K, prec in mine:
    if (round(a, 6), K) in done: continue
    t0 = time.time()
    if (K, prec) not in engines: engines[(K, prec)] = OddWeilArb(K=K, prec=prec)
    W = engines[(K, prec)]
    lam0, lam1, fa = W.floor(a)
    row = {"a": a, "K": K, "prec": prec, "lam0": lam0.mid().str(30, radius=False), "lam0_rad": lam0.rad().str(3),
           "lam1": lam1.mid().str(15, radius=False), "fa": fa.mid().str(15, radius=False), "seconds": round(time.time() - t0)}
    open(f"arb_grid_{chunk}.jsonl", "a").write(json.dumps(row) + "\n")
    print(json.dumps(row), flush=True)
print("done", flush=True)
