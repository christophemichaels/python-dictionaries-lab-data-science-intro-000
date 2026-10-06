"""The subpower closure attempt through the horizon modes (SUBPOWER_CLOSURE.md). Objects of the closure prompt
(received/Michaels_Arithmophysics_Subpower_Closure_Prompt.pdf): B_N = [min(a,b)], its sine modes v_j with gains lambda_j,
c_N(n) = mu(n)/n, b = V^T c, S(J) = sum_{j<=J} |b_j|^2, K_N = max_J lambda_J S(J), the prime induction (B)-(D) with the
transfer T_p = V_N^T E_p V_K.  Computes: (A) two-sided; (C) for several p; the induction over all primes <= N with the
squared and cross terms of (D) at the maximising J; the transfer matrix's diagonal and aliases; mode populations against
the random-sign benchmark (3 pi / 2) J^2 / N^2 and the trivial bound; dyadic mode bands against the dyadic historical
integral.  Usage: python3 rh_subpower_closure.py SUBPOWER_CLOSURE.png [N=2000]"""
import sys, math, json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
out = sys.argv[1]; N = int(sys.argv[2]) if len(sys.argv) > 2 else 2000
rng = np.random.default_rng(7)
X = N
s = np.ones(X + 1, dtype=bool); s[:2] = False
for p in range(2, int(X**0.5) + 1):
    if s[p]: s[p*p::p] = False
primes = np.nonzero(s)[0]
mu = np.ones(X + 1, dtype=np.int8)
for p in primes: mu[p::p] *= -1; mu[p*p::p*p] = 0
mu[0] = 0
M = np.cumsum(mu.astype(np.int64))
nn = np.arange(1, N + 1); jj = np.arange(1, N + 1)
theta = (2*jj - 1) * math.pi / (2*N + 1); lam = 1 / (4 * np.sin(theta/2)**2)
V = (2 / math.sqrt(2*N + 1)) * np.sin(np.outer(nn, theta))            # V[n-1, j-1] = v_j(n)
def R_direct(Nn): return float(np.sum(M[1:Nn].astype(float)**2 / (nn[:Nn-1] * (nn[:Nn-1] + 1))) + M[Nn]**2 / Nn)
c = mu[1:N+1] / nn; b = V.T @ c; S = np.cumsum(b**2)
R = float(np.sum(lam * b**2)); KN = float(np.max(lam * S)); Jstar = int(np.argmax(lam * S)) + 1
mult = 1 + math.log(lam[0] / lam[-1])
print(f"N = {N}: R(N) direct {R_direct(N):.9f}, mode resolution {R:.9f};  K_N = {KN:.9f} at J* = {Jstar};  (A): K_N <= R <= (1 + log(l1/lN)) K_N  ->  {KN:.6f} <= {R:.6f} <= {mult * KN:.6f}  multiplier {mult:.3f}")
print(f"   orthonormality max|V^T V - I| = {np.max(np.abs(V.T @ V - np.eye(N))):.1e};  max|B v_1 - l_1 v_1| = {np.max(np.abs(np.minimum.outer(nn, nn) @ V[:, 0] - lam[0]*V[:, 0])):.1e}")
# (C) for several primes
for p in (2, 3, 7, 101):
    K = N // p; thK = (2*np.arange(1, K+1) - 1) * math.pi / (2*K + 1); VK = (2 / math.sqrt(2*K + 1)) * np.sin(np.outer(np.arange(1, K+1), thK))
    lamK = 1 / (4 * np.sin(thK/2)**2)
    T = V[p*np.arange(1, K+1) - 1, :].T @ VK                                   # T = V_N^T E V_K  (N x K)
    e1 = np.max(np.abs(T.T @ T - np.eye(K))); e2 = np.max(np.abs(T.T @ (lam[:, None] * T) - p * np.diag(lamK)))
    diag = np.array([T[j, j] for j in range(min(K, 30))]); r = N - p*K
    off = np.abs(T.copy()); np.fill_diagonal(off, 0)
    big = np.unravel_index(np.argsort(off.ravel())[-3:][::-1], off.shape)
    print(f"(C) p = {p:3d}, K = {K:4d}, N = pK + {r}: max|T^T T - I| = {e1:.1e}, max|T^T Lam T - p Lam_K| = {e2:.1e};  T_jj for j = 1..5: " + " ".join(f"{d:.5f}" for d in diag[:5])
          + f"  (p^(-1/2) = {p**-0.5:.5f});  largest off-diagonal |T_jl|: " + ", ".join(f"({big[0][i]+1},{big[1][i]+1}) {off[big[0][i], big[1][i]]:.4f}" for i in range(3)))
    col = T[:, 0]; top = np.argsort(np.abs(col))[::-1][:p]
    pred = sorted(set(int(round(((2*N+1)*(2*m - sgn/(2*K+1))/p + 1)/2)) for m in range(1, (p+1)//2 + 1) for sgn in (1, -1)))
    print(f"      column l = 1 of T: squared mass {float(np.sum(col**2)):.6f} = 1; its {p} largest entries (j, T_j1, T_j1^2 * p): " + ", ".join(f"({j+1}, {col[j]:+.4f}, {p*col[j]**2:.3f})" for j in top)
          + f";  predicted alias rows j = ((2N+1)(2m -/+ 1/(2K+1))/p + 1)/2: {pred[:p]}")
# the prime induction (B), (D) with the terms at J* and the running K
cQ = np.zeros(N); cQ[0] = 1.0; bQ = V.T @ cQ
Jlist = [1, 10, 50, Jstar]
rows = []; stage = None
for p in primes:
    K = N // p
    EcK = np.zeros(N); EcK[p*np.arange(1, K+1) - 1] = cQ[:K]                 # E_p c_{Q,K}
    tb = V.T @ EcK                                                           # T_p b_{Q,K} = V_N^T E c_{Q,K}
    rec = {"p": int(p)}
    for Jt in Jlist:
        sq = float(np.sum(tb[:Jt]**2)) / p**2; cr = 2.0 / p * float(np.dot(bQ[:Jt], tb[:Jt]))
        cs = 2.0 / p * math.sqrt(float(np.sum(bQ[:Jt]**2)) * float(np.sum(tb[:Jt]**2)))
        rec[Jt] = (float(np.sum(bQ[:Jt]**2)), sq, cr, cs)
    cQ = cQ - EcK / p; bQ_new = V.T @ cQ
    rec["errB"] = float(np.max(np.abs(bQ_new - (bQ - tb / p))))
    for Jt in Jlist: rec[Jt] = rec[Jt] + (float(np.sum(bQ_new[:Jt]**2)),)
    Sq = np.cumsum(bQ_new**2); rec["KQ"] = float(np.max(lam * Sq)); rec["JQ"] = int(np.argmax(lam * Sq)) + 1; rec["S1"] = float(bQ_new[0]**2)
    if p <= N / 2: stage = dict(rec)
    rows.append(rec); bQ = bQ_new
errD = max(abs(r[Jt][0] + r[Jt][1] - r[Jt][2] - r[Jt][4]) for r in rows for Jt in Jlist); errB = max(r["errB"] for r in rows)
print(f"induction over the {len(primes)} primes <= N: max error in (B) {errB:.1e}, in (D) {errD:.1e}")
for Jt in Jlist:
    sq_sum = sum(r[Jt][1] for r in rows); cr_sum = sum(r[Jt][2] for r in rows); base = float(np.sum(V[0, :Jt]**2))
    neg = sum(1 for r in rows if r[Jt][2] < 0); med = np.median([abs(r[Jt][2])/r[Jt][3] for r in rows if r[Jt][3] > 0])
    print(f"   J = {Jt:5d}: base {base:.3e} + squared {sq_sum:.3e} - cross {cr_sum:+.3e} = S(J) {base + sq_sum - cr_sum:.3e} (direct {S[Jt-1]:.3e});  target scale J^2/N^2 = {Jt**2/N**2:.3e};  cross < 0 at {neg}/{len(rows)} primes; median |cross|/CS {med:.3f}")
    print("        dyadic prime packets [2^k,2^(k+1)): squared, cross, net:  " + "  ".join(f"k{k}: {sum(r[Jt][1] for r in rows if 2**k <= r['p'] < 2**(k+1)):.1e},{sum(r[Jt][2] for r in rows if 2**k <= r['p'] < 2**(k+1)):+.1e},{sum(r[Jt][1]-r[Jt][2] for r in rows if 2**k <= r['p'] < 2**(k+1)):+.1e}" for k in range(1, int(math.log2(N)) + 1) if any(2**k <= r['p'] < 2**(k+1) for r in rows)))
# the squared terms are the slices of the source by largest prime factor: p^{-1} E_p c_{Q_p,K} = -c 1_{P+(n) = p}
Pplus = np.zeros(N + 1, dtype=np.int64)
for q in primes: Pplus[q::q] = q
slice_check = 0.0; Mp_sum = 0.0; sq1 = 0.0
cQ2 = np.zeros(N); cQ2[0] = 1.0
for p in primes:
    K = N // p; EcK = np.zeros(N); EcK[p*np.arange(1, K+1) - 1] = cQ2[:K]
    cslice = np.where(Pplus[1:N+1] == p, c, 0.0); slice_check = max(slice_check, float(np.max(np.abs(EcK / p + cslice))))
    m = np.arange(1, K + 1); smooth = Pplus[1:K+1] < p; Mp = float(np.sum(mu[1:K+1][smooth]))
    Mp_sum += Mp**2; sq1 += float(np.dot(V[:, 0], EcK))**2 / p**2; cQ2 = cQ2 - EcK / p
print(f"   slices: max|p^(-1) E_p c_(Q_p,K) + c 1_(P+(n)=p)| = {slice_check:.1e};  sum_p sq_p(1) = {sq1:.4e};  (4 pi^2/(2N+1)^3) sum_p M_p(N/p)^2 = {4*math.pi**2/(2*N+1)**3*Mp_sum:.4e}  (M_p(x) = Moebius sum over p-smooth m <= x);  J^2/N^2 = {1/N**2:.4e};  R(N)/N^2 = {R/N**2:.4e}")
imax = int(np.argmax([r["KQ"] for r in rows]))
print(f"   running K_Q(N) along the induction: final {rows[-1]['KQ']:.4f} at J = {rows[-1]['JQ']};  maximum {rows[imax]['KQ']:.4f} at J = {rows[imax]['JQ']} after the prime {rows[imax]['p']};  ratio max/final {rows[imax]['KQ']/rows[-1]['KQ']:.1f}")
print(f"   stage Q = primes <= N/2 (after p = {stage['p']}): K_Q = {stage['KQ']:.4f} at J = {stage['JQ']};  lambda_1 S_Q(1) = {lam[0]*stage['S1']:.4f};  N/log^2 N = {N/math.log(N)**2:.4f};  ratio {lam[0]*stage['S1']/(N/math.log(N)**2):.4f}")
sel = [r for r in rows if r["p"] in (2, 3, 5, 7, 11, 13, 31, 97, 199, 499, 997, 1999, 3989) or r is rows[-1]]
print("     (prime: K_Q @ J) " + "  ".join(f"{r['p']}:{r['KQ']:.3f}@{r['JQ']}" for r in sel))
# mode populations against the random benchmark and the trivial bound
Js = np.array([1, 2, 3, 5, 10, 20, 50, 100, 200, 500, 1000])
Js = Js[Js <= N]
rand = []
for t in range(40):
    cr_ = c * rng.choice([-1.0, 1.0], size=N); br = V.T @ cr_; rand.append(np.cumsum(br**2))
rand = np.array(rand); rmean = rand.mean(axis=0)
print("   mode populations S(J) N^2/J^2: Moebius, random-sign mean (benchmark 3 pi/2 = 4.712), trivial bound J log^2 J / N  *  N^2/J^2")
for J in Js:
    print(f"     J = {J:5d}: Moebius {S[J-1]*N**2/J**2:8.4f}   random {rmean[J-1]*N**2/J**2:8.4f}   trivial {J*max(1, math.log(J))**2/N * N**2/J**2:10.1f}")
Rr = rand @ np.diff(np.concatenate([[0], np.cumsum(lam)]))  # placeholder, replaced below
Rrand = [float(np.sum(lam * np.diff(np.concatenate([[0], r_])))) for r_ in rand]
print(f"   R(N): Moebius {R:.4f};  random-sign mean {np.mean(Rrand):.4f} (3/pi log N = {3/math.pi*math.log(N):.4f});  K_N random mean {np.mean([float(np.max(lam * r_)) for r_ in rand]):.4f}")
# dyadic mode bands against the dyadic historical integral
print("   dyadic mode bands j in [2^k, 2^(k+1)) : band energy sum lambda_j b_j^2   against   sum_{N/2^(k+1) < x <= N/2^k} M(x)^2/(x(x+1))")
bands = []
for k in range(0, int(math.log2(N)) + 1):
    lo, hi = 2**k, min(2**(k+1), N + 1)
    be = float(np.sum(lam[lo-1:hi-1] * b[lo-1:hi-1]**2))
    xlo, xhi = N // 2**(k+1), N // 2**k
    hi_ = float(np.sum(M[xlo+1:xhi+1].astype(float)**2 / (np.arange(xlo+1, xhi+1) * (np.arange(xlo+1, xhi+1) + 1.0)))) if xhi > xlo else 0.0
    bands.append((k, be, hi_)); print(f"     k = {k:2d}: band {be:.5f}   historical {hi_:.5f}")
json.dump({"N": N, "R": R, "KN": KN, "Jstar": Jstar, "S": S[:1000].tolist(), "rand_S": rmean[:1000].tolist(),
           "induction": [{"p": r["p"], "KQ": r["KQ"], "JQ": r["JQ"], "S1": r["S1"], "J50": r[50]} for r in rows], "bands": bands}, open(f"data/subpower_closure_{N}.json", "w"))
# figure
BLUE, ORANGE, INK, INK2, SURF, GRID = "#2a78d6", "#eb6834", "#0b0b0b", "#52514e", "#fcfcfb", "#e6e5e1"
plt.rcParams.update({"font.size": 10, "axes.titlesize": 11, "axes.edgecolor": GRID, "axes.labelcolor": INK2, "xtick.color": INK2, "ytick.color": INK2,
                     "axes.spines.top": False, "axes.spines.right": False, "figure.facecolor": SURF, "axes.facecolor": SURF, "text.color": INK})
fig, ax = plt.subplots(1, 3, figsize=(17, 5.2))
for a in ax: a.grid(True, color=GRID, lw=0.6); a.set_axisbelow(True)
A = ax[0]; J = np.arange(1, N + 1)
A.plot(J, lam * S, color=BLUE, lw=1.4, label="Möbius source: λ_J S_N(J)")
A.plot(J, lam * rmean, color=ORANGE, lw=1.4, label="random-sign source, mean of 40: λ_J S_N(J)")
A.axhline(R, color=INK, lw=0.8, ls="--", label=f"R(N) = {R:.3f}")
A.scatter([Jstar], [KN], s=40, color=BLUE, zorder=4, edgecolor=SURF)
A.annotate(f"K_N = {KN:.3f} at J* = {Jstar}", (Jstar, KN), xytext=(-110, 12), textcoords="offset points", fontsize=8.5, color=INK2)
A.set_xscale("log"); A.set_xlabel("J"); A.set_ylabel("λ_J S_N(J)")
A.set_title(f"A.  The target along the modes, N = {N}: K_N = max_J λ_J S_N(J)", loc="left")
A.legend(frameon=False, fontsize=8.5, loc="center right")
Bx = ax[1]; ps = [r["p"] for r in rows]
Bx.plot(ps, [r["KQ"] for r in rows], color=BLUE, lw=1.6, label="K_Q(N) after adding the prime p, primes in increasing order")
Bx.axhline(KN, color=BLUE, lw=0.8, ls=":", label=f"final K_N = {KN:.3f}")
Bx.axhline(N / math.log(N)**2, color=INK2, lw=0.8, ls="--", label=f"N / log² N = {N/math.log(N)**2:.2f}")
Bx.axvline(N / 2, color=INK2, lw=0.8, ls=":"); Bx.text(N/2 * 1.05, KN * 1.5, "p = N/2", fontsize=8.5, color=INK2)
Bx.set_xscale("log"); Bx.set_yscale("log"); Bx.set_xlabel("prime p added"); Bx.set_ylabel("K_Q(N)")
Bx.set_title("B.  The prime induction is not an invariant region: K_Q(N) along the way", loc="left")
Bx.legend(frameon=False, fontsize=8.5, loc="center left")
C = ax[2]; Jt = 50
C.plot(ps, np.cumsum([r[Jt][1] for r in rows]), color=ORANGE, lw=1.4, label="Σ squared terms  p^{−2}‖Π_J T b_{Q,K}‖²")
C.plot(ps, np.cumsum([r[Jt][2] for r in rows]), color=INK2, lw=1.4, label="Σ cross terms  (2/p) Re⟨Π_J b_{Q,N}, Π_J T b_{Q,K}⟩")
C.plot(ps, [r[Jt][4] for r in rows], color=BLUE, lw=1.6, label=f"S_Q(J) after each prime, J = {Jt}")
C.axhline(Jt**2 / N**2, color=BLUE, lw=0.8, ls=":", label=f"J²/N² = {Jt**2/N**2:.1e}")
C.set_xscale("log"); C.set_xlabel("prime p added"); C.set_ylabel(f"population of the first {Jt} modes")
C.set_title("C.  The update (D) at J = 50: squared terms, cross terms, and what they leave", loc="left")
C.legend(frameon=False, fontsize=8.5, loc="upper left")
fig.suptitle("The closure attempt through the horizon modes: the target K_N, and the prime induction that builds the Möbius source", fontsize=11)
fig.tight_layout(); fig.savefig(out, dpi=160, bbox_inches="tight"); print("wrote", out)
