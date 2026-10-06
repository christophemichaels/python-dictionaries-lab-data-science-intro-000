"""The map through the horizon (THE_HORIZON_MAP.md): the light plane (u = log x, the dilation flow u -> u + t, the field
Hamiltonian H0 = H_dil^2 + 1/4, the propagator 1/(k^2 + 1/4)) and its arithmetic image, the lattice log N with the weights
mu(n) n^{-1/2}, which admits two evolutions: by integers (J_N, the Dirichlet polynomial) and by primes (J^(P) =
prod_{p<=P} (1 - p^{-1/2} S_p) delta_0, the finite Euler product). Computes: (1) the spectral colours of both, the integer
side showing lines at the zeros, the prime side the envelope exp(-2 Re S_P(k)) of the pole term; (2) Plancherel: the mean
colour density of J_N equals the delivered energy sum mu^2/n; (3) the prime-cutoff energy E_P from its spectrum (checked
against an exact enumeration at P = 7) and the prime update E' = (1 + 1/q) E - 2 q^{-1/2} C_q, which is the plane's
conservation under the shift plus one interference; the growth of E_P in sqrt(P)/log P; (4) the integer side, delivered
sum mu^2/n (slope 6/pi^2 per e-fold) against retained R(N).  Usage: python3 rh_horizon_map.py THE_HORIZON_MAP.png"""
import sys, math, json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
out = sys.argv[1]
X = 10**7
s = np.ones(X + 1, dtype=bool); s[:2] = False
for p in range(2, int(X**0.5) + 1):
    if s[p]: s[p*p::p] = False
primes = np.nonzero(s)[0]
mu = np.ones(X + 1, dtype=np.int8)
for p in primes: mu[p::p] *= -1; mu[p*p::p*p] = 0
mu[0] = 0
n = np.arange(X + 1, dtype=float); M = np.cumsum(mu.astype(np.int64))
def R(Nn): return float(np.sum(M[1:Nn].astype(float)**2 * (1/n[1:Nn] - 1/n[2:Nn+1])) + M[Nn]**2 / Nn)
zeros = np.loadtxt("data/zeros_6000.txt")[:, 1]
res = {}
# ---------- (1) spectral colours ----------
kk = np.arange(0, 60, 0.005)
def spec_integer(Nn, k):
    """|J_N(k)|^2, J_N(k) = sum_{n<=N} mu(n) n^{-1/2-ik}, by direct summation."""
    idx = np.nonzero(mu[1:Nn+1])[0] + 1; w = mu[idx] / np.sqrt(idx); ln = np.log(idx); acc = np.zeros(len(k), dtype=complex)
    for b in range(0, len(k), 20000):
        kb = k[b:b+20000]
        for a in range(0, len(idx), 1000):
            acc[b:b+20000] += (w[a:a+1000][None, :] * np.exp(-1j * kb[:, None] * ln[a:a+1000][None, :])).sum(axis=1)
    return np.abs(acc)**2
def spec_prime(P, k):
    """log of prod_{p<=P} |1 - p^{-1/2-ik}|^2 = sum log(1 + 1/p - 2 p^{-1/2} cos(k log p)), and -2 Re S_P(k), S_P = sum p^{-1/2-ik}."""
    ps = primes[primes <= P]; lp = np.log(ps); lg = np.zeros(len(k)); S = np.zeros(len(k))
    for a in range(0, len(ps), 500):
        c = np.cos(k[:, None] * lp[a:a+500][None, :]) / np.sqrt(ps[a:a+500])[None, :]
        lg += np.sum(np.log(1 + 1/ps[a:a+500][None, :] - 2*c), axis=1); S += np.sum(c, axis=1)
    return lg, -2*S
S3 = spec_integer(10**3, kk); S5 = spec_integer(10**5, kk)
lgE2, poleE2 = spec_prime(100, kk); lgE4, poleE4 = spec_prime(10**4, kk)
print("spectral colours on 0 <= k < 60   (zeros at 14.135, 21.022, 25.011)")
for name, S in (("integer cutoff N = 10^3", S3), ("integer cutoff N = 10^5", S5)):
    near = np.min(np.abs(kk[:, None] - zeros[None, :3]), axis=1) < 0.25
    print(f"  {name}: mean over colours = {np.mean(S):.3f};  mean within 0.25 of gamma_1,2,3 = {np.mean(S[near]):.3f} (ratio {np.mean(S[near])/np.mean(S):.1f});  peak colour {kk[np.argmax(S)]:.3f};  at k = 0: {S[0]:.4f}")
for name, lg, pole in (("prime cutoff P = 100 ", lgE2, poleE2), ("prime cutoff P = 10^4", lgE4, poleE4)):
    off = kk > 2
    print(f"  {name}: log10 of the product ranges over [{np.min(lg)/math.log(10):.1f}, {np.max(lg)/math.log(10):.1f}];  at gamma_1: {math.exp(lg[np.argmin(np.abs(kk-zeros[0]))]):.4g};  "
          f"log(product) - (-2 Re S_P): mean {np.mean((lg-pole)[off]):+.3f}, sd {np.std((lg-pole)[off]):.3f} for k > 2")
# ---------- (2) Plancherel: the mean colour density is the delivered energy ----------
kw = np.arange(0, 10**4, 0.01); SW = spec_integer(10**3, kw)
deliv3 = float(np.sum(mu[1:1001].astype(float)**2 / n[1:1001]))
print(f"Plancherel at N = 10^3: mean of |J_N(k)|^2 over 0 <= k < 10^4 is {np.mean(SW):.4f} against sum mu^2/n = {deliv3:.4f}")
del SW, kw
# ---------- (3) the prime-cutoff energy and the prime update ----------
def E_prime(P, K=4000.0, dk=0.001):
    """E_P = (1/2 pi) int |J^(P)(k)|^2 / (k^2 + 1/4) dk, the tail beyond K from the mean of the product, prod (1 + 1/p)."""
    k = np.arange(0, K, dk); lg, _ = spec_prime(P, k); ps = primes[primes <= P]
    return np.trapezoid(np.exp(lg) / (k**2 + 0.25), k) / math.pi + float(np.prod(1 + 1/ps)) / K / math.pi
sm = [1]
for p in (2, 3, 5, 7): sm = sm + [p*m for m in sm]
sm = np.array(sorted(sm)); musm = np.array([mu[m] for m in sm])
E7_exact = float(sum(musm[i]*musm[j]/max(sm[i], sm[j]) for i in range(16) for j in range(16)))
print(f"prime cutoff P = 7: energy from the spectrum {E_prime(7):.6f}; exact sum over the 16 squarefree 7-smooth numbers {E7_exact:.6f}")
print("prime update  E_(P,q) = (1 + 1/q) E_P - 2 q^(-1/2) C_q :")
print("   q    free prod(1+1/p)      E_P      interference C_q   E_P / free   log E_P / (sqrt P / log P)")
Plist = [int(q) for q in primes[primes <= 100]] + [int(q) for q in primes[(primes > 100) & (primes <= 400)][::4]]
rows = []; prev = 1.0
for q in Plist:
    Eq = E_prime(q); free = float(np.prod(1 + 1/primes[primes <= q])); Cq = ((1 + 1/q) * prev - Eq) * math.sqrt(q) / 2
    rows.append((q, free, Eq, Cq)); prev = Eq
    print(f"  {q:3d}   {free:10.4f}   {Eq:12.5f}   {Cq:+12.5f}   {Eq/free:9.4f}   {math.log(Eq)/(math.sqrt(q)/math.log(q)):8.4f}")
res["prime_update"] = rows
# ---------- (4) integer side: delivered against retained ----------
Ns = np.unique(np.concatenate([np.logspace(2, 7, 60).astype(int), [10**3, 10**4, 10**5, 10**6, 10**7]]))
self_e = np.array([float(np.sum(mu[1:Nn+1].astype(float)**2 / n[1:Nn+1])) for Nn in Ns]); Rs = np.array([R(int(Nn)) for Nn in Ns])
print("integer side: delivered sum mu^2/n (slope 6/pi^2 = 0.6079 per e-fold) against retained R(N)")
for Nn in (10**3, 10**4, 10**5, 10**6, 10**7):
    i = np.searchsorted(Ns, Nn); print(f"  N = {Nn:8d}: delivered {self_e[i]:.4f}   retained {Rs[i]:.4f}   fraction {Rs[i]/self_e[i]:.4f}")
sl = np.polyfit(np.log(Ns[Ns >= 10**4]), Rs[Ns >= 10**4], 1)[0]
print(f"  least-squares slope of R(N) in log N over 10^4..10^7: {sl:.4f} per e-fold, against 6/pi^2 = {6/math.pi**2:.4f} delivered")
res["integer"] = {"N": Ns.tolist(), "delivered": self_e.tolist(), "R": Rs.tolist(), "slope": sl}
json.dump(res, open("data/horizon_map.json", "w"), indent=1)
# ---------- figure ----------
BLUE, ORANGE, INK, INK2, SURF, GRID = "#2a78d6", "#eb6834", "#0b0b0b", "#52514e", "#fcfcfb", "#e6e5e1"
plt.rcParams.update({"font.size": 10, "axes.titlesize": 11, "axes.edgecolor": GRID, "axes.labelcolor": INK2, "xtick.color": INK2, "ytick.color": INK2,
                     "axes.spines.top": False, "axes.spines.right": False, "figure.facecolor": SURF, "axes.facecolor": SURF, "text.color": INK})
fig, ax = plt.subplots(2, 2, figsize=(14, 9))
for a in ax.flat: a.grid(True, color=GRID, lw=0.6); a.set_axisbelow(True)
A = ax[0, 0]
A.plot(kk, S3, color=BLUE, lw=1.0, label="N = 10³")
A.plot(kk, S5, color=ORANGE, lw=1.0, alpha=0.9, label="N = 10⁵")
for g in zeros[zeros < 60]: A.axvline(g, color=INK2, lw=0.6, ls=":")
A.set_title("A.  By integers: |Σ_{n≤N} μ(n) n^{−1/2−ik}|², lines at the zeros (dotted)", loc="left")
A.set_xlabel("colour k"); A.set_ylabel("|Ĵ_N(k)|²"); A.legend(frameon=False, fontsize=9, loc="upper left")
B = ax[0, 1]
B.plot(kk, lgE2 / math.log(10), color=BLUE, lw=1.0, label="P = 100")
B.plot(kk, lgE4 / math.log(10), color=ORANGE, lw=1.0, alpha=0.9, label="P = 10⁴")
B.plot(kk, poleE4 / math.log(10), color=INK, lw=0.8, ls="--", label="pole term exp(−2 Re Σ_{p≤10⁴} p^{−1/2−ik})")
for g in zeros[zeros < 60]: B.axvline(g, color=INK2, lw=0.6, ls=":")
B.set_title("B.  By primes: Π_{p≤P} |1 − p^{−1/2−ik}|², the pole's envelope; lines only above k ≈ √P/log P", loc="left")
B.set_xlabel("colour k"); B.set_ylabel("log₁₀ |Ĵ^{(P)}(k)|²"); B.legend(frameon=False, fontsize=9, loc="upper right")
C = ax[1, 0]
xs = [q for q, _, _, _ in rows]
C.plot(xs, [E for _, _, E, _ in rows], color=BLUE, lw=1.4, marker="o", ms=4, label="E_P, the energy of the prime-cutoff source")
C.plot(xs, [f for _, f, _, _ in rows], color=ORANGE, lw=1.4, marker="o", ms=4, label="Π_{p≤P}(1 + 1/p), the conserved copies alone")
C.set_xscale("log"); C.set_yscale("log")
C.set_title("C.  By primes the energy runs away: log E_P grows like 2√P/log P, the pole", loc="left")
C.set_xlabel("prime cutoff P"); C.set_ylabel("energy"); C.legend(frameon=False, fontsize=9, loc="upper left")
D = ax[1, 1]
D.plot(Ns, self_e, color=ORANGE, lw=1.6, label="delivered: Σ_{n≤N} μ(n)²/n, slope 6/π² per e-fold")
D.plot(Ns, Rs, color=BLUE, lw=1.6, label=f"retained: R(N), slope ≈ {sl:.3f} per e-fold over 10⁴–10⁷")
D.set_xscale("log")
D.set_title("D.  By integers the energy is tempered: delivered against retained", loc="left")
D.set_xlabel("integer cutoff N"); D.set_ylabel("energy"); D.legend(frameon=False, fontsize=9, loc="upper left")
fig.suptitle("The map through the horizon: the Möbius source seen in its colours, by integers (pole-free, tempered) and by primes (the pole passes)", fontsize=11)
fig.tight_layout(); fig.savefig(out, dpi=160, bbox_inches="tight"); print("wrote", out)
