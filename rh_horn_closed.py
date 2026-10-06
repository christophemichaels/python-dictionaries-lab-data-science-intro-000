"""Closing the horn by unfolding (THE_HORN.md, Section 6). The Riemann-Siegel phase theta(t)/pi + 1 sends the n-th zero to
phi_n = n - S(gamma_n): in that coordinate every ring has unit radius and sits at integer height up to the wobble S, so the horn
becomes a cylinder, and the spacings of the rings follow the GUE law. Usage: python3 rh_horn_closed.py data/zeros_6000.txt THE_HORN_CLOSED.png"""
import sys, math
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

zf, out = sys.argv[1], sys.argv[2]
g = np.array([float(l.split()[1]) for l in open(zf)])
n = np.arange(1, len(g) + 1); r = g / (2 * math.pi)
theta = lambda t: (t / 2) * np.log(t / (2 * math.pi)) - t / 2 - math.pi / 8 + 1 / (48 * t)
phi = theta(g) / math.pi + 1                 # unfolded height, phi_n = n - S(gamma_n)
jit = phi - n; sp = np.diff(phi)
x = np.linspace(0, 3, 400)
gue = (32 / math.pi**2) * x**2 * np.exp(-4 * x**2 / math.pi); poi = np.exp(-x)
edges = np.linspace(0, 3, 31); h, _ = np.histogram(sp, bins=edges, density=True); c = (edges[1:] + edges[:-1]) / 2
gue_c = (32 / math.pi**2) * c**2 * np.exp(-4 * c**2 / math.pi); poi_c = np.exp(-c); w = edges[1] - edges[0]
print(f"unfolded rings: jitter mean {jit.mean():.4f}, std {jit.std():.4f}, range [{jit.min():.3f}, {jit.max():.3f}]")
print(f"spacings: mean {sp.mean():.5f}, std {sp.std():.4f}, min {sp.min():.4f}, max {sp.max():.3f}")
small = (sp < 0.25).mean(); gue_small = float(np.trapezoid(gue[x <= 0.25], x[x <= 0.25])); poi_small = 1 - math.exp(-0.25)
print(f"spacings below 0.25: measured {small:.4f}; GUE {gue_small:.4f}; Poisson {poi_small:.4f}")
print(f"histogram distance to GUE {math.sqrt(((h-gue_c)**2).sum()*w):.4f}, to Poisson {math.sqrt(((h-poi_c)**2).sum()*w):.4f}")

BLUE, ORANGE, INK, INK2, SURF, GRID = "#2a78d6", "#eb6834", "#0b0b0b", "#52514e", "#fcfcfb", "#e6e5e1"
plt.rcParams.update({"font.size": 10, "axes.titlesize": 11, "axes.edgecolor": GRID, "axes.labelcolor": INK2, "xtick.color": INK2, "ytick.color": INK2,
                     "axes.spines.top": False, "axes.spines.right": False, "figure.facecolor": SURF, "axes.facecolor": SURF, "text.color": INK})
fig, ax = plt.subplots(1, 3, figsize=(15, 4.8))
for a in ax: a.grid(True, color=GRID, lw=0.6); a.set_axisbelow(True)
A = ax[0]
A.plot(n[:300], r[:300] / r[299], color=BLUE, lw=2, label="the horn: ring radius γₙ/2π (scaled to ring 300)")
A.plot(n[:300], np.ones(300), color=ORANGE, lw=2, label="the closed horn: unit rings at the unfolded height φₙ")
A.set_xlabel("ring number n"); A.set_ylabel("radius"); A.set_title("A.  Unfolding closes the flare: the horn becomes a cylinder", loc="left")
A.legend(frameon=False, fontsize=8.5, loc="upper left"); A.set_ylim(0, 1.15)
B = ax[1]
B.plot(n[:200], jit[:200] + 0.5, color=BLUE, lw=1.2, label="ring n sits at height n − S(γₙ): the jitter S (first 200 rings, centred)")
B.axhline(0, color=INK2, lw=0.8); B.set_ylim(-1.2, 1.2)
B.set_xlabel("ring number n"); B.set_ylabel("displacement from integer height"); B.set_title("B.  What is left after closing: the wobble S, never past one ring", loc="left")
B.legend(frameon=False, fontsize=8.5, loc="upper left")
C = ax[2]
C.bar(c, h, width=w * 0.9, color=BLUE, alpha=0.35, label="measured spacings of the unfolded rings (5,999)")
C.plot(x, gue, color=BLUE, lw=2, label="GUE (Wigner surmise): rings repel")
C.plot(x, poi, color=ORANGE, lw=2, label="Poisson: rings independent")
C.set_xlabel("spacing between consecutive rings, in units of the mean"); C.set_ylabel("density")
C.set_title("C.  The law of the rings' spacing: the quantum signature", loc="left"); C.legend(frameon=False, fontsize=8.5)
fig.suptitle("Forcing the horn to close: in the phase coordinate θ/π the rings are unit rings at integer heights, and what remains is S and the GUE law", fontsize=11)
fig.tight_layout(); fig.savefig(out, dpi=160, bbox_inches="tight"); print("wrote", out)
