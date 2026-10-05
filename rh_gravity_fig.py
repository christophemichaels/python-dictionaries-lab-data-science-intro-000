"""Figure for THE_GRAVITY_PLANE.md: the field energies of the prime star and the Moebius star against log N, with the slopes the zeros
predict. Usage: python3 rh_gravity_fig.py THE_GRAVITY_PLANE.png [X=1e7]"""
import sys, math
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
out = sys.argv[1]; X = int(float(sys.argv[2])) if len(sys.argv) > 2 else 10**7
EG = 0.5772156649015329; CRAMER = 2 + EG - math.log(4 * math.pi); MOEB = 0.0288   # sum 2/|rho zeta'(rho)|^2 over the first 400 pairs
s = np.ones(X + 1, dtype=bool); s[:2] = False
for p in range(2, int(X**0.5) + 1):
    if s[p]: s[p*p::p] = False
primes = np.nonzero(s)[0]
Lam = np.zeros(X + 1)
for p in primes:
    q = p
    while q <= X: Lam[q] = math.log(p); q *= p
mu = np.ones(X + 1, dtype=np.int8)
for p in primes: mu[p::p] *= -1; mu[p*p::p*p] = 0
mu[0] = 0
psi = np.cumsum(Lam); M = np.cumsum(mu.astype(np.int64)); t = np.arange(X + 1, dtype=float)
w = 1/t[1:-1] - 1/t[2:]
IM = np.concatenate([[0.0], np.cumsum(M[1:-1].astype(float)**2 * w)])
fl = np.concatenate([[0.0], np.cumsum((psi[1:-1] - t[1:-1])**2 * w)])
selfE = np.cumsum(mu[1:].astype(float)**2 / t[1:])
Ns = np.array([int(round(10**e)) for e in np.arange(2, 7.0001, 0.02)]); L = np.log(Ns)
BLUE, ORANGE, INK, INK2, SURF, GRID = "#2a78d6", "#eb6834", "#0b0b0b", "#52514e", "#fcfcfb", "#e6e5e1"
plt.rcParams.update({"font.size": 10, "axes.titlesize": 11, "axes.edgecolor": GRID, "axes.labelcolor": INK2, "xtick.color": INK2, "ytick.color": INK2,
                     "axes.spines.top": False, "axes.spines.right": False, "figure.facecolor": SURF, "axes.facecolor": SURF, "text.color": INK})
fig, ax = plt.subplots(1, 2, figsize=(13, 4.8))
for a in ax: a.grid(True, color=GRID, lw=0.6); a.set_axisbelow(True)
A = ax[0]
y = fl[Ns - 1]; A.plot(L, y, color=BLUE, lw=2, label="prime star: fluctuation field energy ∫(ψ − t)²/t² dt")
A.plot(L, y[-1] + CRAMER * (L - L[-1]), color=ORANGE, lw=1.6, ls="--", label=f"slope set by the zeros: Σ 1/|ρ|² = 2 + γ − log 4π = {CRAMER:.4f}")
A.set_xlabel("log N"); A.set_ylabel("field energy"); A.set_title("A.  The prime star: its fluctuation energy grows at Cramér's rate", loc="left"); A.legend(frameon=False, fontsize=8.5, loc="upper left")
B = ax[1]
B.plot(L, selfE[Ns - 1], color=INK2, lw=1.2, label="self-energy of the shells Σ μ(n)²/n ~ (6/π²) log N")
B.plot(L, IM[Ns - 1] - selfE[Ns - 1], color=INK2, lw=1.2, ls=":", label="interaction between the shells (the pull): negative")
B.plot(L, IM[Ns - 1], color=BLUE, lw=2, label="Möbius star: net field energy I_M(N)")
B.plot(L, IM[Ns - 1][-1] + MOEB * (L - L[-1]), color=ORANGE, lw=1.6, ls="--", label=f"slope set by the zeros: Σ 2/|ρζ′(ρ)|² ≈ {MOEB:.4f}")
B.axhline(0, color=GRID, lw=0.8)
B.set_xlabel("log N"); B.set_ylabel("field energy"); B.set_title("B.  The Möbius star: bound, and growing at the zeros' rate", loc="left"); B.legend(frameon=False, fontsize=8.5, loc="upper left")
fig.suptitle("The plane as a gravitational field: masses on shells of radius n, potential 1/max(r, R), energy ∫ M(t)²/t² dt, modes = the zeros", fontsize=11)
fig.tight_layout(); fig.savefig(out, dpi=160, bbox_inches="tight"); print("wrote", out)
