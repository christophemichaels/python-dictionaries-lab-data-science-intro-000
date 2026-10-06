"""The magenta horn, drawn: the first 1,000 zeros as rings of circumference gamma_n at height n, with the stones on its axis
(prime power m at the height where the ring radius is m, Crossing 1 of ARITHMOPHYSICS). Usage: python3 rh_magenta_horn.py data/zeros_6000.txt MAGENTA_HORN.png"""
import sys, math
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

zf, out = sys.argv[1], sys.argv[2]
g = np.array([float(l.split()[1]) for l in open(zf)])[:1000]
n = np.arange(1, len(g) + 1); r = g / (2 * math.pi)
th = np.linspace(0, 2 * math.pi, 181)
X = np.outer(r, np.cos(th)); Y = np.outer(r, np.sin(th)); Z = np.outer(n, np.ones_like(th))
cmap = LinearSegmentedColormap.from_list("magenta", ["#fbe4ef", "#e87ba4", "#c2178f", "#6e0a52"])
SURF, INK, INK2 = "#fcfcfb", "#0b0b0b", "#52514e"
fig = plt.figure(figsize=(9.5, 8.2), facecolor=SURF)
ax = fig.add_subplot(111, projection="3d", facecolor=SURF)
ax.plot_surface(X, Y, Z, facecolors=cmap((Z - 1) / (len(g) - 1)), rstride=4, cstride=4, linewidth=0, antialiased=True, alpha=0.92, shade=False)
for k in range(0, len(g), 40):                       # every 40th ring drawn as a line
    ax.plot(X[k], Y[k], Z[k], color="#6e0a52", lw=0.5, alpha=0.7)
# the stones on the axis: prime power m sits at the height where the ring radius equals m (horizon 2 pi m)
def is_pp(m):
    for p in range(2, m + 1):
        if m % p == 0:
            while m % p == 0: m //= p
            return m == 1
    return False
stones = [m for m in range(2, 230) if is_pp(m)]
hs = [int((g <= 2 * math.pi * m).sum()) for m in stones]
ax.plot([0] * len(stones), [0] * len(stones), hs, "o", color=INK, ms=3, zorder=10)
ax.plot([0, 0], [0, 0], [0, len(g)], color=INK2, lw=0.8, alpha=0.8)
for m, h in zip(stones, hs):
    if m in (11, 17, 31, 61, 127, 227):
        ax.text(0, 0, h, f"  {m}", color=INK, fontsize=8, zorder=11)
ax.set_xlim(-230, 230); ax.set_ylim(-230, 230); ax.set_zlim(0, 1000)
ax.set_xlabel("ring radius = γₙ / 2π", color=INK2, fontsize=9); ax.set_ylabel("", color=INK2); ax.set_zlabel("zero number n", color=INK2, fontsize=9)
ax.view_init(elev=22, azim=-58)
ax.xaxis.pane.fill = ax.yaxis.pane.fill = ax.zaxis.pane.fill = False
for a in (ax.xaxis, ax.yaxis, ax.zaxis): a.pane.set_edgecolor("#e6e5e1"); a._axinfo["grid"]["color"] = "#e6e5e1"
ax.tick_params(colors=INK2, labelsize=8)
ax.set_title("The magenta horn: the first 1,000 zeros as rings of circumference γₙ at height n.\nOn the axis, the stones: prime power m at the height where the ring radius is m (the horizon 2πm).", fontsize=10, color=INK, loc="left")
fig.tight_layout(); fig.savefig(out, dpi=170, bbox_inches="tight", facecolor=SURF); print("wrote", out, "stones", list(zip(stones, hs))[:8])
