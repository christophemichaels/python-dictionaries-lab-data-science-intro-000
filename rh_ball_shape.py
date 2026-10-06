"""The shape of the final estimate, in ball arithmetic (BALL_SHAPE.md).  python-flint arb/acb balls, 200-bit working
precision; every reported number is a ball [mid +/- rad] and the radius is printed.  Objects: R(N) = sum_{k<N} M(k)^2/(k(k+1))
+ M(N)^2/N, delta_N = R(N) - R(floor(N/6)), B(X) = 1 + S(X), c_N = [delta_N]_+ / B(N^(1/6))^5 (Hypothesis 1.1 of the
manuscript), the two-mode signed readout (10.10), the Hardy coefficients gamma_{k,N} = sum mu(n)/n P_k(log n) with the
generating function exp(-x z/(1-z)) (9.6), the three-step modified increment (2.31), the two-phase circle readout (8.2)-(8.4).
Usage: python3 rh_ball_shape.py   (writes data/ball_shape.json; about 8 minutes to 1e7)"""
import math, json, time, numpy as np
from fractions import Fraction as Fr
from flint import arb, acb, fmpq, ctx, arb_series
ctx.prec = 200
X = 10**7
t0 = time.time()
s = np.ones(X + 1, dtype=bool); s[:2] = False
for p in range(2, int(X**0.5) + 1):
    if s[p]: s[p*p::p] = False
PR = np.nonzero(s)[0]
MU = np.ones(X + 1, dtype=np.int8)
for p in PR: MU[p::p] *= -1; MU[p*p::p*p] = 0
MU[0] = 0
idx = np.arange(X + 1); n_ = idx.astype(float); n_[0] = 1.0
M = np.cumsum(MU.astype(np.int64))
inner_f = np.cumsum(M.astype(float)**2/(n_*(n_ + 1))); Rf = np.zeros(X + 1); Rf[1:] = inner_f[:-1] + M[1:].astype(float)**2/n_[1:]
D6f = np.full(X + 1, np.nan); D6f[6:] = Rf[6:] - Rf[idx[6:]//6]
out = {}
def ball(x, d=12): return x.str(d, radius=True)
# exact small values: R(n), S(n), B(n) for n <= 14, as rationals
Rq = {0: Fr(0)}
for n in range(1, 15): Rq[n] = sum(Fr(int(M[k])**2, k*(k + 1)) for k in range(1, n)) + Fr(int(M[n])**2, n)
Sq = {n: max([Fr(1)] + [Rq[k] for k in range(1, n + 1)]) for n in range(1, 15)}
def B5(N):   # B(N^(1/6))^5 exactly: X = floor(N^(1/6)) (<= 14 on the range), B = 1 + S(X)
    x = int(round(N**(1/6)))
    while x**6 > N: x -= 1
    while (x + 1)**6 <= N: x += 1
    return (1 + Sq[x])**5, x
print(f"== 0. exact base values: R(n) for n <= 14 and the envelope B(X)^5 = (1 + S(X))^5 ==")
print("   R(n): " + ", ".join(f"R({n}) = {Rq[n]} = {float(Rq[n]):.6f}" for n in range(1, 15)))
print("   S(X), B(X)^5 for X = 1..14: " + ", ".join(f"X={x}: S={float(Sq[x]):.5f}, B^5={float((1+Sq[x])**5):.3f}" for x in range(1, 15)))
print(f"   43/30 = S(5) = {Sq[5]};  B >= 2 so B^5 >= 32")
# ---------------------------------------------------------------- 1. certified energies along chains
def chain(N):
    c = [N]
    while c[-1] >= 6: c.append(c[-1]//6)
    return c
heads = [10**7, 8388608, 4*10**6, 1679616, 10**6, 10**5, 16384, 10**4]
need = set()
for N in heads: need |= set(chain(N))
dec = []
for j in range(1, 7):
    lo, hi = 10**j, min(X, 10**(j+1)); seg = D6f[lo:hi+1]; k = int(np.nanargmax(seg)) + lo; dec.append((lo, hi, k)); need |= {k, k//6}
    cvals = seg/np.array([B5(n)[0] for n in range(lo, hi + 1)], dtype=float) if hi - lo < 2000 else None
peaks = [13, 20, 110, 2838, 42968, 300551, 1065673, 6481601, 7108798]
for N in peaks: need |= {N, N//6}
need = sorted(n for n in need if n >= 1)
print(f"== 1. certified pass to 1e7 ({len(need)} recorded horizons) ==")
t1 = time.time(); inner = arb(0); h = arb(0); rec = {}; needset = set(need); kk = 1
for n in range(1, X + 1):
    mu = int(MU[n]); Mn = int(M[n])
    if mu: h += arb(mu)/n
    if n in needset: rec[n] = (inner + arb(Mn*Mn)/n, h, Mn)
    inner += arb(Mn*Mn)/(n*(n + 1))
print(f"   pass done in {time.time()-t1:.0f} s; R(1e7) = {ball(rec[X][0], 20)}, h(1e7) = {ball(rec[X][1], 15)}; float R(1e7) = {Rf[X]:.12f} (float error bound {X*2.3e-16*2:.1e})")
def two_mode(N):
    """exact two-mode signed readout (10.10): p0 = R(L)/B^5, q0 = ||t||^2/B^5, z0 = <v,t>/B^5 = M(L)(h(N)-h(L))/B^5."""
    L = N//6; RN, hN_, MN = rec[N]; RL, hL_, ML = rec[L]; b5 = arb(fmpq(B5(N)[0].numerator, B5(N)[0].denominator))
    z = arb(ML)*(hN_ - hL_); dlt = RN - RL; tt = dlt - 2*z
    p0 = RL/b5; q0 = tt/b5; z0 = z/b5; Dl = q0 + 2*z0; om = p0*q0 - z0*z0; rt = (Dl*Dl + 4*om).sqrt()
    return dict(N=N, L=L, R=RN, RL=RL, delta=dlt, cross=2*z, tail2=tt, b5=b5, c=(dlt/b5), p0=p0, q0=q0, z0=z0, nu_plus=(Dl + rt)/2, nu_minus=(Dl - rt)/2, omega=om, cos=z0/(p0*q0).sqrt())
rows1 = []
for N0 in heads:
    print(f"   chain from {N0}:")
    for N in chain(N0)[:-1]:
        r = two_mode(N); rows1.append({k: (ball(v, 15) if isinstance(v, arb) else v) for k, v in r.items()})
        print(f"      N = {N:8d} L = {r['L']:7d}: R(N) = {ball(r['R'], 12)}  delta = {ball(r['delta'], 10)}  cross 2M(L)(h(N)-h(L)) = {ball(r['cross'], 8)}  ||t||^2 = {ball(r['tail2'], 8)}  B^5 = {float(r['b5']):.3f}  c_N = {ball(r['c'], 8)}  nu+ = {ball(r['nu_plus'], 7)}  nu- = {ball(r['nu_minus'], 7)}  cos(a,d) = {ball(r['cos'], 5)}")
out['chains'] = rows1
# ---------------------------------------------------------------- 2. the coefficient c_N: certified decade maxima and the shape
print("== 2. the normalized coefficient c_N = [delta_N]_+ / B(N^(1/6))^5: certified maxima (argmax from the float pass, values certified) ==")
rows2 = []
for lo, hi, k in dec:
    r = two_mode(k); rows2.append(dict(lo=lo, hi=hi, argmax=k, delta=ball(r['delta'], 12), c=ball(r['c'], 12), b5=float(r['b5'])))
    print(f"   [{lo:8d},{hi:8d}]: max delta_N = {ball(r['delta'], 12)} at N = {k} (B^5 = {float(r['b5']):.3f}, X = {B5(k)[1]}), c_N there = {ball(r['c'], 10)}")
# the global maximum of c_N: delta_13/32 exactly
c13 = (Rq[13] - Rq[2])/32; print(f"   global maximum of c_N (N >= 6): N = 13, delta_13 = R(13) - R(2) = {Rq[13] - Rq[2]} = {float(Rq[13]-Rq[2]):.9f}, B^5 = 32, c* = {c13} = {float(c13):.9f}")
# maxima of c over N > 100 and over N >= 1e6 (float ranking, certified value at the argmax)
cf = np.full(X + 1, np.nan)
b5arr = np.zeros(X + 1)
xs6 = np.floor(n_**(1/6)).astype(int); xs6 = np.clip(xs6, 1, 14)
# fix rounding at exact sixth powers
for x in range(1, 15):
    if x**6 <= X: xs6[x**6] = x
    if x**6 - 1 >= 1 and x**6 - 1 <= X: xs6[x**6 - 1] = min(xs6[x**6 - 1], x - 1) if x > 1 else 1
b5tab = np.array([0.0] + [float((1 + Sq[x])**5) for x in range(1, 15)]); b5arr = b5tab[xs6]
cf[6:] = np.maximum(D6f[6:], 0)/b5arr[6:]
for lo in (6, 101, 10**4 + 1, 10**6 + 1):
    k = int(np.nanargmax(cf[lo:])) + lo; r = two_mode(k) if k in rec else None
    print(f"   max c_N over N >= {lo}: at N = {k}, float {cf[k]:.9f}" + (f", certified {ball(r['c'], 10)}" if r else f" (R not recorded in the ball pass; float error < 1e-8)"))
# the envelope B(X) itself along decades (float S from the pass; certified at 1e7)
Sf = np.maximum.accumulate(Rf)
print("   B(X) = 1 + S(X) at X = 10^j: " + ", ".join(f"{1 + Sf[10**j]:.5f}" for j in range(1, 8)) + f";  log B(1e7) = {math.log(1 + Sf[X]):.5f};  F(t) = log B(e^t) is flat: F(log 1e7)/F(log 1e7 / 6) = {math.log(1+Sf[X])/math.log(1+Sf[int(round(X**(1/6)))]):.4f}")
# explicit return of Theorem 11.2 in balls
alpha = arb(5).log()/arb(6).log(); cstar = arb(fmpq(c13.numerator, c13.denominator)); Cstar = arb(fmpq(73, 30)) + cstar/arb(6).log()
Kc = 5*(arb(fmpq(5, 4))*arb(2).log() + Cstar.log()/4 + arb(fmpq(5, 16))*arb(6).log())
print(f"   Theorem 11.2 with c* = {float(c13):.6f}: alpha = {ball(alpha, 10)}, C* = {ball(Cstar, 10)}, K(c*) = {ball(Kc, 10)};  bound log B(X) <= K (log X)^alpha: at X = 1e7: {ball(Kc*(arb(10**7).log())**alpha, 8)} (actual log B = {math.log(1+Sf[X]):.4f}); at 1e10: {ball(Kc*(arb(10**10).log())**alpha, 8)}; at 1e100: {ball(Kc*(arb(10**100).log())**alpha, 8)}")
print(f"   with the coefficient above 100 only (c = {cf[110]:.6f}): K = {ball(5*(arb(fmpq(5,4))*arb(2).log() + (arb(fmpq(73,30)) + arb(cf[110])/arb(6).log()).log()/4 + arb(fmpq(5,16))*arb(6).log()), 8)}")
out['decades'] = rows2; out['cstar'] = [c13.numerator, c13.denominator]
# ---------------------------------------------------------------- 3. Hardy coefficients across the horizon (series exponential, balls)
print("== 3. Hardy coefficients gamma_{k,N} (k <= 12) at the recorded horizons; head, tail, increments ==")
Kmax = 12; ctx.cap = Kmax + 1; z = arb_series([0, 1]); gam = arb_series([0]); grec = {}; t1 = time.time()
sq = np.nonzero(MU)[0]
for n in sq:
    n = int(n)
    if n == 1: gam = gam + arb_series([1])
    else: gam = gam + ((-arb(n).log())*z/(1 - z)).exp()*(arb(int(MU[n]))/n)
    if n in needset or n + 1 in needset:
        pass
    # record at the largest needed horizon <= next squarefree: handled below
    if n in needset: grec[n] = [gam[k] for k in range(Kmax + 1)]
# horizons that are not squarefree: gamma at N equals gamma at the last squarefree n <= N
for N in need:
    if N not in grec:
        m = N
        while MU[m] == 0: m -= 1
        grec[N] = grec[m] if m in grec else None
print(f"   series pass done in {time.time()-t1:.0f} s")
# fill any missing (horizon whose predecessor squarefree was not recorded): recompute directly for those
missing = [N for N in need if grec.get(N) is None]
if missing:
    for N in missing:
        g = arb_series([0])
        for n in range(1, N + 1):
            if MU[n]: g = g + (arb_series([1]) if n == 1 else ((-arb(n).log())*z/(1 - z)).exp()*(arb(int(MU[n]))/n))
        grec[N] = [g[k] for k in range(Kmax + 1)]
def Hhead(N, K): return sum(grec[N][k]*grec[N][k] for k in range(K + 1))
def allowance_99(N, K):  # (9.9): C_K (1 + log N)^(2K)
    CK = 1 + 4*sum(sum(arb(math.comb(k-1, j-1))/math.factorial(j) for j in range(1, k + 1))**2 for k in range(1, K + 1)); return CK*(1 + arb(N).log())**(2*K)
def allowance_910(N, K): x = arb(N).log(); return 1 + 4*K*x*x*(4*(K*x).sqrt()).exp()
rows3 = []
for N0 in (10**7, 10**6):
    print(f"   chain from {N0} (K_N = floor((log N)^(3/4)) at the head of each pair; same K at both ends):")
    ch = chain(N0)
    for N, L in zip(ch[:-1], ch[1:]):
        if L < 2: continue
        K = int(math.floor(math.log(N)**0.75)); HN = Hhead(N, K); HL = Hhead(L, K); TN = rec[N][0] - HN; TL = rec[L][0] - HL
        row = dict(N=N, L=L, K=K, HK_N=ball(HN, 10), HK_L=ball(HL, 10), TK_N=ball(TN, 10), TK_L=ball(TL, 10), tail_inc=ball(TN - TL, 10), delta=ball(rec[N][0] - rec[L][0], 10), allow99=ball(allowance_99(N, K), 6), allow910=ball(allowance_910(N, K), 6), gam=[ball(grec[N][k], 10) for k in range(K + 1)])
        rows3.append(row)
        print(f"      N = {N:8d}, L = {L:7d}, K = {K:2d}: H_K(N) = {ball(HN, 9)}  H_K(L) = {ball(HL, 9)}  T_K(N) - T_K(L) = {ball(TN - TL, 9)}  (delta_N = {ball(rec[N][0] - rec[L][0], 9)});  head allowance (9.9) = {ball(allowance_99(N, K), 4)}, (9.10) = {ball(allowance_910(N, K), 4)}")
    print(f"      gamma_k(N0), k = 0..12: " + ", ".join(ball(grec[N0][k], 8) for k in range(Kmax + 1)))
out['hardy'] = rows3
# convergence of the Hardy energy at N = 1000: gamma_k to k = 800
print("   Theorem 9.1 at N = 1000: partial sums H_K(1000) against R(1000)")
ctx.prec = 400; ctx.cap = 801; z = arb_series([0, 1]); g = arb_series([0])
for n in range(1, 1001):
    if MU[n]: g = g + (arb_series([1]) if n == 1 else ((-arb(n).log())*z/(1 - z)).exp()*(arb(int(MU[n]))/n))
R1000 = sum(arb(int(M[k])**2)/(k*(k + 1)) for k in range(1, 1000)) + arb(int(M[1000])**2)/1000
Hk = []; acc = arb(0); conv = []
for k in range(801):
    acc += g[k]*g[k]
    if k in (0, 1, 2, 5, 10, 20, 50, 100, 200, 400, 800): conv.append((k, ball(acc, 12), ball(R1000 - acc, 6)))
print("      " + "; ".join(f"K={k}: H_K = {hk}, tail = {tl}" for k, hk, tl in conv) + f";  R(1000) = {ball(R1000, 14)}")
print(f"      |gamma_k(1000)| at k = 1, 10, 100, 400, 800: " + ", ".join(ball(abs(g[k]), 5) for k in (1, 10, 100, 400, 800)))
out['hardy_conv'] = conv
ctx.prec = 200
# ---------------------------------------------------------------- 4. the three-step modified increment (2.31) along the chain from 1e7
print("== 4. three-step grouping (2.23)-(2.31): Q3_N = E(a~) - R(L) against delta_N; allowance 2 log(N/L) ==")
rows4 = []
for N, L in zip(chain(10**7)[:-1], chain(10**7)[1:]):
    if N - L < 2: continue
    m = N - L; q3, r = divmod(m, 3)
    widths = [3]*q3 if r == 0 else ([3]*q3 + [2] if r == 2 else [3]*(q3 - 1) + [2, 2])
    ends = [L]; 
    for w in widths: ends.append(ends[-1] + w)
    assert ends[-1] == N
    qj = [int(M[b] - M[a]) for a, b in zip(ends[:-1], ends[1:])]
    # energy of the compressed vector: cumulative charge = M(k) for k <= L, then M(L) + sum_{t_j <= k} q_j
    inc = np.zeros(N + 1, dtype=np.int64); inc[np.array(ends[1:])] = np.array(qj, dtype=np.int64)   # block masses at the block ends
    cumc = np.array(M[:N+1], dtype=np.int64).copy(); cumc[L+1:] = M[L]; cumc = cumc + np.cumsum(inc)      # cumulative charge of a~
    assert cumc[N] == M[N]
    # E(a~) in balls: sum_{k<N} cumc(k)^2/(k(k+1)) + cumc(N)^2/N  (exact integers in the numerators)
    E = arb(0)
    for k in range(1, N): E += arb(int(cumc[k])**2)/(k*(k + 1))
    E += arb(int(cumc[N])**2)/N
    Q3 = E - rec[L][0]; dlt = rec[N][0] - rec[L][0]
    rows4.append(dict(N=N, L=L, Q3=ball(Q3, 10), delta=ball(dlt, 10), diff=ball(dlt - Q3, 10), allow=2*math.log(N/L)))
    print(f"   N = {N:8d}, L = {L:7d}: Q3_N = {ball(Q3, 9)}  delta_N = {ball(dlt, 9)}  delta - Q3 = {ball(dlt - Q3, 6)}  (allowance 2 log(N/L) = {2*math.log(N/L):.4f}; c^(3) = [Q3]_+/B^5 = {ball(Q3/arb(fmpq(B5(N)[0].numerator, B5(N)[0].denominator)), 6)})")
out['three_step'] = rows4
# ---------------------------------------------------------------- 5. the two-phase circle readout with rigorous tails (small N)
print("== 5. circle readouts C_{p,theta} (8.2) with certified tails; (8.3) and the two-phase isometry (8.4) ==")
def circle(N, Kc=200000):
    p = int(PR[np.searchsorted(PR, (2*N + 1)**2)]); ell = arb(p).log(); sp = arb(p).sqrt()
    ns = [n for n in range(1, N + 1) if MU[n]]; A = sum(arb(1)/arb(n).sqrt() for n in ns)   # |D| <= A
    res = {}
    for name, th in (("0", arb(0)), ("pi", arb.pi())):
        tot = arb(0)
        for k in range(-Kc, Kc + 1):
            kap = (2*arb.pi()*k + th)/ell; sv = acb(arb(fmpq(1, 2)), kap); D = acb(0)
            for n in ns: D += int(MU[n])*acb(n)**(-sv)
            tot += abs(D)**2/(kap*kap + arb(fmpq(1, 4)))
        tot = tot/ell
        tail = A*A*ell/(arb.pi()**2*(2*Kc - 1))                      # sum_{|k|>K} 1/kappa_k^2 <= ell^2/(pi^2 (2K-1))
        res[name] = tot.union(tot + tail)                                 # the omitted tail lies in [0, tail]
    RN = sum(arb(int(M[k])**2)/(k*(k + 1)) for k in range(1, N)) + arb(int(M[N])**2)/N
    Mh = arb(int(M[N]))*sum(arb(int(MU[n]))/n for n in ns)
    d0 = 1/(sp - 1); dpi = -1/(sp + 1); w0 = (sp - 1)/(2*sp); wpi = (sp + 1)/(2*sp)
    print(f"   N = {N}, p = {p}, ell = {ball(ell, 8)}, K = {Kc}: C_p,0 = {ball(res['0'], 8)} vs R + 2 d_0 M h = {ball(RN + 2*d0*Mh, 10)};  C_p,pi = {ball(res['pi'], 8)} vs R + 2 d_pi M h = {ball(RN + 2*dpi*Mh, 10)};  w0 C_0 + wpi C_pi = {ball(w0*res['0'] + wpi*res['pi'], 8)} vs R(N) = {ball(RN, 10)};  |C_pi - R| = {ball(abs(res['pi'] - RN), 4)} <= 1")
    return dict(N=N, p=p, C0=ball(res['0'], 10), Cpi=ball(res['pi'], 10), R=ball(RN, 12))
rows5 = [circle(5), circle(12, 150000)]
out['circle'] = rows5
json.dump(out, open('data/ball_shape.json', 'w'), indent=1, default=str)
print(f"wrote data/ball_shape.json; total {time.time()-t0:.0f} s")
