"""The zero hole (THE_ZERO_HOLE.md). Invert the gravity plane, r -> 1/r (the Kelvin transform): the shells of radius n become the
layers of an onion at radius 1/n, the point at infinity becomes the centre, and the potential of all the shells at the centre is
sum m(n)/n, the Dirichlet series at s = 1. For the Moebius star that potential is h(N) = sum_{n<=N} mu(n)/n -> 0 (the prime number
theorem), at the square-root rate iff RH; for the prime star it is sum Lambda(n)/n = log N - gamma + o(1) (Mertens).
Usage: python3 rh_zero_hole.py THE_ZERO_HOLE.png [X=1e7]"""
import sys, math
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
out = sys.argv[1]; X = int(float(sys.argv[2])) if len(sys.argv) > 2 else 10**7
EG = 0.5772156649015329
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
n = np.arange(X + 1, dtype=float)
h = np.concatenate([[0.0], np.cumsum(mu[1:].astype(float) / n[1:])])          # potential of the Moebius star at the centre
mert = np.concatenate([[0.0], np.cumsum(Lam[1:] / n[1:])])                      # potential of the prime star at the centre
Ns = np.array([int(round(10**e)) for e in np.arange(2, 7.0001, 0.01)]); L = np.log(Ns)
print("potential at the zero hole (the centre of the onion), sum_{n<=N} m(n)/n:")
for N in (10**3, 10**4, 10**5, 10**6, X):
    print(f"  N = {N:9d}: Moebius star h(N) = {h[N]:+.6f}, sqrt N * h(N) = {math.sqrt(N)*h[N]:+.4f};  prime star {mert[N]:.6f} = log N - {math.log(N)-mert[N]:.6f}  (gamma = {EG:.6f}); scaled Mertens error (log N - gamma - sum) sqrt N / log^2 N = {(math.log(N)-EG-mert[N])*math.sqrt(N)/math.log(N)**2:+.4f}")
sh = np.sqrt(Ns) * h[Ns]
print(f"sqrt N * h(N) over 10^2..10^7: min {sh.min():+.3f}, max {sh.max():+.3f}, rms {math.sqrt((sh**2).mean()):.3f}")
BLUE, ORANGE, INK, INK2, SURF, GRID, MAG = "#2a78d6", "#eb6834", "#0b0b0b", "#52514e", "#fcfcfb", "#e6e5e1", "#c2178f"
plt.rcParams.update({"font.size": 10, "axes.titlesize": 11, "axes.edgecolor": GRID, "axes.labelcolor": INK2, "xtick.color": INK2, "ytick.color": INK2,
                     "axes.spines.top": False, "axes.spines.right": False, "figure.facecolor": SURF, "axes.facecolor": SURF, "text.color": INK})
fig = plt.figure(figsize=(16, 5.2))
A = fig.add_subplot(1, 3, 1); A.set_aspect("equal"); A.set_xticks([]); A.set_yticks([]); A.spines["left"].set_visible(False); A.spines["bottom"].set_visible(False)
th = np.linspace(0, 2 * math.pi, 400)
for k in range(1, 61):
    r = 1 / k; c = BLUE if mu[k] > 0 else (ORANGE if mu[k] < 0 else GRID)
    A.plot(r * np.cos(th), r * np.sin(th), color=c, lw=1.4 if mu[k] else 0.6, alpha=0.9 if mu[k] else 0.6)
A.plot([0], [0], "o", color=MAG, ms=9); A.annotate("the zero hole: r = ∞ folded to the centre;\npotential here = Σ m(n)/n", xy=(0, 0), xytext=(0.12, 0.78), fontsize=8.5, color=MAG, arrowprops=dict(arrowstyle="-", color=MAG, lw=0.8))
A.plot([], [], color=BLUE, lw=1.4, label="layer n with μ(n) = +1"); A.plot([], [], color=ORANGE, lw=1.4, label="layer n with μ(n) = −1"); A.plot([], [], color=GRID, lw=0.6, label="μ(n) = 0: no mass")
A.legend(frameon=False, fontsize=8.5, loc="lower left"); A.set_xlim(-1.05, 1.05); A.set_ylim(-1.05, 1.05)
A.set_title("A.  The Möbius onion: shell n at radius 1/n, n ≤ 60", loc="left")
B = fig.add_subplot(1, 3, 2); B.grid(True, color=GRID, lw=0.6); B.set_axisbelow(True)
B.plot(L, sh, color=BLUE, lw=1.2, label="√N · Σ_{n≤N} μ(n)/n  (Möbius star)")
B.axhline(0, color=INK2, lw=0.8); B.set_ylim(-1.2, 1.2)
B.set_xlabel("log N"); B.set_ylabel("potential × √N"); B.set_title("B.  The potential at the hole vanishes, at the square-root rate", loc="left"); B.legend(frameon=False, fontsize=8.5, loc="upper left")
C = fig.add_subplot(1, 3, 3); C.grid(True, color=GRID, lw=0.6); C.set_axisbelow(True)
C.plot(L, mert[Ns] - L, color=BLUE, lw=1.4, label="Σ_{n≤N} Λ(n)/n − log N  (prime star)")
C.axhline(-EG, color=ORANGE, lw=1.6, ls="--", label=f"−γ = {-EG:.4f}  (Mertens)")
C.set_ylim(-0.7, -0.45); C.set_xlabel("log N"); C.set_ylabel("potential − log N")
C.set_title("C.  The prime star's potential at the hole is log N − γ", loc="left"); C.legend(frameon=False, fontsize=8.5, loc="lower right")
fig.suptitle("The zero hole: the plane inverted into an onion, every shell's point at infinity folded to one centre, and the potential of the masses there", fontsize=11)
fig.tight_layout(); fig.savefig(out, dpi=160, bbox_inches="tight"); print("wrote", out)
