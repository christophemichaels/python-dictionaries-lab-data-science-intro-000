# Green energy R(N) of the Moebius vector, block-closure statistics, and M(x) on a log grid (Section 13 of mobius_modifier_v2.tex).
# Green energy R(N) = sum_{k<N} M(k)^2/(k(k+1)) + M(N)^2/N of the Moebius vector, block-closure statistics.
import numpy as np, time, json
N = 50_000_000
t0 = time.time()
sieve = np.ones(N+1, dtype=bool); sieve[:2] = False
for p in range(2, int(N**0.5)+1):
    if sieve[p]: sieve[p*p::p] = False
primes = np.nonzero(sieve)[0]; del sieve
mu = np.ones(N+1, dtype=np.int8); mu[0] = 0
for p in primes:
    mu[p::p] *= -1
    if p*p <= N: mu[p*p::p*p] = 0
print(f"mu sieved to {N} in {time.time()-t0:.0f}s, primes {len(primes)}", flush=True)
# chunked cumulative sums
chunk = 5_000_000; M_prev = 0; R = 0.0; checkpoints = {}
returns = []          # k with M(k) = 0
sign_changes = 0; last_sign = 0
maxabs = 0
cps = [10**k for k in range(2, 8)] + [2*10**7, 5*10**7]
for start in range(1, N+1, chunk):
    end = min(start+chunk-1, N)
    k = np.arange(start, end+1, dtype=np.int64)
    M = M_prev + np.cumsum(mu[start:end+1].astype(np.int64))
    Mf = M.astype(np.float64)
    # sum over k<N of M(k)^2/(k(k+1)); handle k=N separately
    inner = Mf**2/(k*(k+1.0))
    if end == N: inner[-1] = 0.0
    R += inner.sum()
    zeros = k[M == 0]; returns.append(zeros)
    s = np.sign(M); nz = s[s != 0]
    if nz.size:
        if last_sign and nz[0] != last_sign: sign_changes += 1
        sign_changes += int(np.count_nonzero(np.diff(nz)))
        last_sign = nz[-1]
    maxabs = max(maxabs, int(np.abs(M).max()))
    for c in cps:
        if start <= c <= end:
            Mc = int(M[c-start]); Rc = R - inner[c-start+1:].sum() + Mc*Mc/c
            checkpoints[c] = (Mc, Rc, Rc/np.log(c))
    M_prev = int(M[-1])
returns = np.concatenate(returns)
np.save("returns.npy", returns)
gaps = np.diff(returns)
# M sampled on a log grid for the explicit-formula comparison
ug = np.arange(np.log(100.0), np.log(N), 0.0005); xg = np.floor(np.exp(ug)).astype(np.int64)
Mfull = np.cumsum(mu.astype(np.int64)); Mg = Mfull[xg]
np.save("Mgrid.npy", np.stack([xg, Mg]))
retcounts = {str(c): int(np.count_nonzero(returns <= c)) for c in cps}
top = np.argsort(gaps)[-10:][::-1]
longest = []
for i in top:
    a, b = int(returns[i]), int(returns[i+1])
    seg = Mfull[a:b+1]
    H = int(np.abs(seg - Mfull[a]).max()); sgn = int(np.sign(seg[1:-1].sum())) if b > a+1 else 0
    E = float(np.sum((seg[1:-1].astype(np.float64))**2/(np.arange(a+1, b, dtype=np.float64)*np.arange(a+2, b+1, dtype=np.float64))))
    longest.append({"a": a, "b": b, "L": b-a, "H": H, "sign": sgn, "E": E})
del Mfull
out = {"N": N, "M_N": M_prev, "R_N": R + M_prev**2/N, "R_over_logN": (R + M_prev**2/N)/np.log(N),
       "checkpoints": {str(c): v for c, v in checkpoints.items()},
       "n_returns": int(returns.size), "sign_changes": sign_changes, "max_abs_M": maxabs,
       "longest_blocks": longest,
       "block_length_quantiles": {q: float(np.quantile(gaps, q)) for q in [0.5, 0.9, 0.99, 0.999]},
       "return_counts": retcounts}
json.dump(out, open("green.json", "w"), indent=1, default=float)
print(json.dumps(out, indent=1, default=float))
print(f"done in {time.time()-t0:.0f}s")
