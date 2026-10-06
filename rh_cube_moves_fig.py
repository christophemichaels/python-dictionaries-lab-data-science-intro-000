"""Figure for CUBE_MOVES.md from data/cube_moves.json: (A) the fixed-ratio increment Delta_6 R(N) across N <= 1e7, binned
min/mean/max, against the controls' maxima on the same bins; (B) the signed components of the fixed-packet difference (13)."""
import json, math, matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt, numpy as np
d = json.load(open("data/cube_moves.json"))
BLUE, ORANGE, INK, MUTED, GRID, SURF = "#2a78d6", "#eb6834", "#0b0b0b", "#52514e", "#e6e5e1", "#fcfcfb"
plt.rcParams.update({"font.size": 10, "axes.edgecolor": MUTED, "axes.labelcolor": INK, "xtick.color": MUTED, "ytick.color": MUTED, "text.color": INK, "axes.titlecolor": INK})
fig, (ax, bx) = plt.subplots(1, 2, figsize=(12.6, 4.6), facecolor=SURF)
for a in (ax, bx): a.set_facecolor(SURF); a.grid(True, color=GRID, lw=0.8); a.set_axisbelow(True); [a.spines[s].set_visible(False) for s in ("top", "right")]
env = np.array(d["envelope"]); x = np.sqrt(env[:, 0]*env[:, 1])
ax.fill_between(x, env[:, 2], env[:, 4], color=BLUE, alpha=0.18, lw=0, label="Möbius: bin min to max")
ax.plot(x, env[:, 3], color=BLUE, lw=1.6, label="Möbius: bin mean")
ctr = {c["name"]: np.array(c["envelope"]) for c in d["controls"]}
ax.plot(x, ctr["independent signs, squarefree support"][:, 4], color=ORANGE, lw=1.2, ls="-", label="random signs, bin max")
ax.plot(x, ctr["block-conditioned shuffle (hexadic)"][:, 4], color=ORANGE, lw=1.2, ls=":", label="hexadic shuffle, bin max")
ax.axhline(0, color=MUTED, lw=0.8)
for N, v in d["increment"]["peaks"]:
    if N >= 1000: ax.plot(N, v, "o", ms=4.5, color=INK, zorder=5)
ax.annotate("largest increments above 10³:\nnew maxima of |M|/√N", xy=(42968, 0.2276), xytext=(1.5e3, 0.62), fontsize=8.5, color=INK, arrowprops=dict(arrowstyle="-", color=MUTED, lw=0.8))
ax.set_xscale("log"); ax.set_xlim(6, 1.2e7); ax.set_ylim(-0.5, 1.6)
ax.set_xlabel("horizon N"); ax.set_ylabel("Δ₆R(N) = R(N) − R(⌊N/6⌋)")
ax.set_title("A. The fixed-ratio increment, every N ≤ 10⁷ (234 log bins)", loc="left", fontsize=10.5)
ax.legend(frameon=False, fontsize=8.5, loc="upper right")
rows = d["packet23"]; Ns = np.array([r["N"] for r in rows]); keys = [("D", "ΔD₂₃ (diagonal, bounded by 1.862)", ORANGE, "-"), ("kern", "Δ overlap (newly completed parents)", ORANGE, "--"), ("cross", "Δ cross 2⟨g,t⟩", MUTED, "-"), ("unf", "Δ unfinished ‖t‖²", MUTED, "--"), ("R", "Δ₆R = sum", BLUE, "-")]
for k, lab, col, ls in keys:
    bx.plot(Ns, [r["diff"][k] for r in rows], color=col, ls=ls, lw=1.8 if k == "R" else 1.3, marker="o" if k == "R" else None, ms=3.5, label=lab)
bx.axhline(0, color=MUTED, lw=0.8); bx.set_xscale("log"); bx.set_xlim(150, 1.5e7); bx.set_ylim(-0.45, 0.45)
bx.set_xlabel("horizon N  (L = ⌊N/6⌋)"); bx.set_ylabel("increment from L to N")
bx.set_title("B. The fixed-packet difference (13): diagonal and overlap cancel", loc="left", fontsize=10.5)
bx.legend(frameon=False, fontsize=8.5, loc="lower left", ncol=1)
fig.tight_layout(); fig.savefig("CUBE_MOVES.png", dpi=160, facecolor=SURF); print("wrote CUBE_MOVES.png")
