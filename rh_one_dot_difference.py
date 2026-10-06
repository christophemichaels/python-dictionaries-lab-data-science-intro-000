"""Forward from the dot against backward from the stones (THE_OTHER_END.md, Section 4).
Forward: every zero placed exactly on the point (beta = 1/2), the first K of them, rebuilt into psi_K(x) by the explicit formula.
Backward: the actual psi(x) from the primes. The difference D_K(x) = psi(x) - psi_K(x), measured in the dot's own unit sqrt x and
averaged in log x, is predicted (Parseval for the almost periodic normalized error, under RH) to have mean square
sum_{n>K} 2/|rho_n|^2. Usage: python3 rh_one_dot_difference.py data/zeros_6000.txt THE_DIFFERENCE.png"""
import sys, math
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

zf, out = sys.argv[1], sys.argv[2]
g = np.array([float(l.split()[1]) for l in open(zf)]); Kmax = len(g)
X = 10**6
# backward: psi(x) from the primes
s = np.ones(X + 1, dtype=bool); s[:2] = False
for p in range(2, int(X**0.5) + 1):
    if s[p]: s[p*p::p] = False
Lam = np.zeros(X + 1)
for p in np.nonzero(s)[0]:
    q = p
    while q <= X: Lam[q] = math.log(p); q *= p
psi = np.cumsum(Lam)
# log grid of non-integer x in [10^2, 10^6]
x = np.floor(np.geomspace(100, X - 1, 20000)) + 0.5
u = np.log(x)
psi_x = psi[np.floor(x).astype(int)]
stones = (psi_x - x) / np.sqrt(x)                       # the normalized error, backward
# forward: the dot's rebuild with K zeros, all placed at beta = 1/2
rho = 0.5 + 1j * g
def rebuild(K):
    tot = np.zeros_like(x)
    for a in range(0, K, 500):
        b = min(K, a + 500)
        tot += 2 * np.real(np.exp(np.outer(u, 1j * g[a:b])) / rho[a:b]).sum(axis=1)   # sum of 2 Re x^{i gamma}/rho, times sqrt x below
    return -tot - (math.log(2 * math.pi) + 0.5 * np.log(1 - x**-2)) / np.sqrt(x)
rows = []
inv2 = 2 / (0.25 + g**2)
T = g[-1]; tail_beyond = (math.log(T / (2 * math.pi)) + 1) / (math.pi * T)        # integral of 2/t^2 dN beyond the last zero
for K in (10, 30, 100, 300, 1000, 3000, 6000):
    dot = rebuild(K); D = stones - dot
    pred = math.sqrt(inv2[K:].sum() + tail_beyond)
    rows.append((K, g[K-1], D.mean(), math.sqrt((D**2).mean()), pred, np.abs(D).max()))
    print(f"K = {K:5d} (height {g[K-1]:8.2f}): difference stones - dot, in units of sqrt x: mean {D.mean():+.4f}, rms {math.sqrt((D**2).mean()):.4f}, predicted rms {pred:.4f}, max {np.abs(D).max():.3f}")
# the stones alone and the dot alone
print(f"stones (psi - x)/sqrt x over the grid: mean {stones.mean():+.4f}, rms {math.sqrt((stones**2).mean()):.4f}; predicted rms if all zeros are on the point sqrt(sum 2/|rho|^2) = {math.sqrt(inv2.sum() + tail_beyond):.4f}")
dot6000 = rebuild(6000); D6000 = stones - dot6000

BLUE, ORANGE, INK, INK2, SURF, GRID = "#2a78d6", "#eb6834", "#0b0b0b", "#52514e", "#fcfcfb", "#e6e5e1"
plt.rcParams.update({"font.size": 10, "axes.titlesize": 11, "axes.edgecolor": GRID, "axes.labelcolor": INK2, "xtick.color": INK2, "ytick.color": INK2,
                     "axes.spines.top": False, "axes.spines.right": False, "figure.facecolor": SURF, "axes.facecolor": SURF, "text.color": INK})
fig, ax = plt.subplots(1, 3, figsize=(15.5, 4.8))
for a in ax: a.grid(True, color=GRID, lw=0.6); a.set_axisbelow(True)
w = (x > 2000) & (x < 20000)
A = ax[0]
A.plot(x[w], stones[w], color=BLUE, lw=1.0, label="backward: the stones, (ψ(x) − x)/√x from the primes")
A.plot(x[w], dot6000[w], color=ORANGE, lw=1.0, alpha=0.9, label="forward: the dot, rebuilt from 6,000 zeros placed at β = ½")
A.set_xscale("log"); A.set_xlabel("x (log scale)"); A.set_ylabel("(ψ − x) / √x"); A.set_title("A.  Forward and backward, 2,000 < x < 20,000", loc="left"); A.legend(frameon=False, fontsize=8.5, loc="upper left")
B = ax[1]
B.plot(x, D6000, color=BLUE, lw=0.6, label="stones − dot, 6,000 zeros")
B.axhline(0, color=INK2, lw=0.8); B.set_xscale("log"); B.set_ylim(-0.6, 0.6)
B.set_xlabel("x (log scale)"); B.set_ylabel("difference, in units of √x"); B.set_title(f"B.  The difference: rms {math.sqrt((D6000**2).mean()):.3f}, predicted {rows[-1][4]:.3f}", loc="left"); B.legend(frameon=False, fontsize=8.5, loc="upper left")
C = ax[2]
Ks = [r[0] for r in rows]; C.plot(Ks, [r[3] for r in rows], "o-", color=BLUE, lw=1.8, ms=6, label="measured rms of the difference")
C.plot(Ks, [r[4] for r in rows], "s--", color=ORANGE, lw=1.6, ms=5, label="predicted: √(Σ_{n>K} 2/|ρₙ|²), the rings beyond K")
C.set_xscale("log"); C.set_yscale("log"); C.set_xlabel("K zeros placed on the dot"); C.set_ylabel("rms, in units of √x")
C.set_title("C.  The difference is the rings not yet placed", loc="left"); C.legend(frameon=False, fontsize=8.5)
fig.suptitle("All zeros on one dot, forward, against the primes, backward: the difference in the dot's unit √x, averaged in log x", fontsize=11)
fig.tight_layout(); fig.savefig(out, dpi=160, bbox_inches="tight"); print("wrote", out)
