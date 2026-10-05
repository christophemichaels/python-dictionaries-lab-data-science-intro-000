"""The colours of the primes and the spectrum of light (THE_LIGHT_PLANE.md). Each prime power n is a colour: the wave
-2 Lambda(n) n^{-1/2} cos(t log n) in the height variable t. Their superposition over all n <= X, smoothly cut off, has its peaks at
the zeros: the primes compose the spectrum, and the zeros are its lines. This is the explicit formula read from the prime side.
Usage: python3 rh_spectrum_of_light.py data/zeros_6000.txt THE_SPECTRUM_OF_LIGHT.png [X=1e6]"""
import sys, math
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
zf, out = sys.argv[1], sys.argv[2]; X = int(float(sys.argv[3])) if len(sys.argv) > 3 else 10**6
g = np.array([float(l.split()[1]) for l in open(zf)])
s = np.ones(X + 1, dtype=bool); s[:2] = False
for p in range(2, int(X**0.5) + 1):
    if s[p]: s[p*p::p] = False
primes = np.nonzero(s)[0]
ns, lam = [], []
for p in primes:
    q = p
    while q <= X: ns.append(q); lam.append(math.log(p)); q *= p
ns = np.array(ns, dtype=float); lam = np.array(lam); ln = np.log(ns)
t = np.linspace(0, 60, 6001)
wt = lam / np.sqrt(ns) * (1 - ln / math.log(X))            # smooth (Fejer-type) cut-off at X
spec = np.zeros_like(t)
for a in range(0, len(ns), 2000):
    b = min(len(ns), a + 2000); spec += -2 * (np.cos(np.outer(t, ln[a:b])) * wt[a:b]).sum(axis=1)
# local maxima of the superposition and the nearest zeros
pk = [i for i in range(1, len(t)-1) if spec[i] > spec[i-1] and spec[i] > spec[i+1] and spec[i] > 0.4 * spec[t > 10].max()]
print("peaks of the prime superposition (t > 10) against the zeros:")
for i in pk:
    if t[i] > 10: j = int(np.argmin(np.abs(g - t[i]))); print(f"  peak at t = {t[i]:6.2f}  nearest zero gamma_{j+1} = {g[j]:7.3f}  (difference {t[i]-g[j]:+.2f})")
BLUE, ORANGE, AQUA, YELLOW, MAG, INK, INK2, SURF, GRID = "#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#0b0b0b", "#52514e", "#fcfcfb", "#e6e5e1"
plt.rcParams.update({"font.size": 10, "axes.titlesize": 11, "axes.edgecolor": GRID, "axes.labelcolor": INK2, "xtick.color": INK2, "ytick.color": INK2,
                     "axes.spines.top": False, "axes.spines.right": False, "figure.facecolor": SURF, "axes.facecolor": SURF, "text.color": INK})
fig, ax = plt.subplots(2, 1, figsize=(13, 7.8))
for a in ax: a.grid(True, color=GRID, lw=0.6); a.set_axisbelow(True)
A = ax[0]
for p, c in zip((2, 3, 5, 7), (BLUE, ORANGE, AQUA, YELLOW)):
    A.plot(t, -2 * math.log(p) / math.sqrt(p) * np.cos(t * math.log(p)), color=c, lw=1.4, label=f"colour {p}: −2 log {p} · {p}^(−½) · cos(t log {p})")
A.set_ylabel("amplitude"); A.set_title("A.  Four colours: the waves of the primes 2, 3, 5, 7 in the height variable t", loc="left"); A.legend(frameon=False, fontsize=8, ncol=4, loc="lower center", bbox_to_anchor=(0.5, 1.0)); A.set_ylim(-1.7, 1.7)
B = ax[1]
B.plot(t, spec, color=INK, lw=1.2, label=f"all colours to {X:.0e}, smoothly cut off: −2 Σ Λ(n) n^(−½) (1 − log n/log X) cos(t log n)")
for k, gm in enumerate(g[g < 60]):
    B.axvline(gm, color=MAG, lw=1.0, alpha=0.7, label="the zeros γₙ" if k == 0 else None)
B.set_xlim(5, 60); B.set_ylim(-12, 22); B.set_xlabel("height t"); B.set_ylabel("superposition"); B.set_title("B.  The spectrum, 5 ≤ t ≤ 60: the colours superposed peak at the zeros (the pole dominates below t = 5 and is cut from view)", loc="left"); B.legend(frameon=False, fontsize=8.5, loc="upper right")
fig.suptitle("The colours of the primes compose the spectrum of light, and the zeros are its lines (the explicit formula read from the prime side)", fontsize=11)
fig.tight_layout(); fig.savefig(out, dpi=160, bbox_inches="tight"); print("wrote", out)
