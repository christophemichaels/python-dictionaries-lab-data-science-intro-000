"""Figure for BALL_SHAPE.md from data/ball_shape.json: (A) the certified decade maxima of the normalized coefficient c_N and
the per-step values along the chain from 1e7, against the constant c* and the envelope curve 0.2453/B(N^(1/6))^5;
(B) the two-mode signed readout along the chain from 1e7: nu+, nu- and their trace c_N; (C) the Hardy head H_K(N) against
its proved allowance (9.10) along the same chain, with the tail increment against delta_N."""
import json, math, re, matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt, numpy as np
BLUE, ORANGE, INK, MUTED, GRID, SURF = "#2a78d6", "#eb6834", "#0b0b0b", "#52514e", "#e6e5e1", "#fcfcfb"
plt.rcParams.update({"font.size": 10, "axes.edgecolor": MUTED, "axes.labelcolor": INK, "xtick.color": MUTED, "ytick.color": MUTED, "text.color": INK, "axes.titlecolor": INK})
d = json.load(open("data/ball_shape.json"))
def mid(s):   # midpoint of a ball string "[m +/- r]" or a plain number
    m = re.match(r"\[?\s*(-?[0-9.eE+-]+)", str(s)); return float(m.group(1))
fig, axs = plt.subplots(1, 3, figsize=(14.5, 4.4), facecolor=SURF)
for a in axs: a.set_facecolor(SURF); a.grid(True, color=GRID, lw=0.8); a.set_axisbelow(True); [a.spines[s].set_visible(False) for s in ("top", "right")]
ax, bx, cx = axs
dec = d["decades"]; ax.plot([r["argmax"] for r in dec], [mid(r["c"]) for r in dec], "o", ms=6, color=ORANGE, label="certified decade maximum of c_N")
ch = [r for r in d["chains"] if r["N"] <= 10**7 and r["N"] in {10000000, 1666666, 277777, 46296, 7716, 1286, 214, 35}]
ax.plot([r["N"] for r in ch], [max(mid(r["c"]), 1e-5) for r in ch], "s", ms=4, color=BLUE, label="c_N on the chain from 10⁷ (clipped at 10⁻⁵)")
cstar = d["cstar"][0]/d["cstar"][1]; ax.axhline(cstar, color=INK, lw=1.0, ls="--", label=f"c* = 18307/480480 = {cstar:.4f} (N = 13)")
Ns = np.logspace(2, 7, 400); b5 = []
Sq = {1: 1, 2: 1, 3: 1, 4: 1, 5: 43/30, 6: 43/30, 7: 43/30, 8: 43/30, 9: 43/30, 10: 43/30, 11: 43/30, 12: 43/30, 13: 51629/30030, 14: 51629/30030}
for N in Ns:
    x = int(math.floor(N**(1/6) + 1e-9)); b5.append((1 + Sq[min(max(x, 1), 14)])**5)
ax.plot(Ns, 0.24525/np.array(b5), color=MUTED, lw=1.2, label="0.2453 / B(N^{1/6})⁵ (largest increment above 100 over the envelope)")
ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlim(8, 2e7); ax.set_ylim(8e-6, 0.08)
ax.set_xlabel("horizon N"); ax.set_ylabel("c_N = [δ_N]₊ / B(N^{1/6})⁵")
ax.set_title("A. The coefficient of Hypothesis 1.1, certified", loc="left", fontsize=10.5); ax.legend(frameon=False, fontsize=7.8, loc="lower left")
chain = [r for r in d["chains"]][:8]
x = np.arange(len(chain)); w = 0.38
bx.bar(x - w/2, [mid(r["nu_plus"]) for r in chain], w, color=BLUE, label="ν₊")
bx.bar(x + w/2, [mid(r["nu_minus"]) for r in chain], w, color=ORANGE, label="ν₋")
bx.plot(x, [mid(r["c"]) for r in chain], "o-", ms=4, color=INK, lw=1.0, label="trace ν₊ + ν₋ = δ_N / B⁵")
bx.axhline(0, color=MUTED, lw=0.8); bx.set_xticks(x); bx.set_xticklabels([f"{r['N']:.0e}".replace("e+0", "e") if r["N"] >= 1000 else str(r["N"]) for r in chain], fontsize=8)
bx.set_xlabel("chain N_j from 10⁷"); bx.set_ylabel("eigenvalues of the rank-two readout D_N")
bx.set_title("B. Two-mode signed readout (10.10) along the chain", loc="left", fontsize=10.5); bx.legend(frameon=False, fontsize=8.5, loc="upper right")
H = [r for r in d["hardy"] if r["N"] in {10000000, 1666666, 277777, 46296, 7716, 1286, 214}]
xn = [r["N"] for r in H]
cx.plot(xn, [mid(r["HK_N"]) for r in H], "o-", color=BLUE, ms=4, label="head H_K(N), K = ⌊(log N)^{3/4}⌋ (certified)")
cx.plot(xn, [mid(r["allow910"]) for r in H], "s--", color=ORANGE, ms=4, label="proved head allowance (9.10)")
cx.plot(xn, [abs(mid(r["tail_inc"])) for r in H], "^-", color=MUTED, ms=4, label="|T_K(N) − T_K(L)| (tail increment)")
cx.plot(xn, [abs(mid(r["delta"])) for r in H], "v:", color=INK, ms=4, label="|δ_N|")
cx.set_xscale("log"); cx.set_yscale("log"); cx.set_xlabel("horizon N on the chain from 10⁷"); cx.set_ylabel("energy")
cx.set_title("C. Hardy head, its allowance, and the tail increment", loc="left", fontsize=10.5); cx.legend(frameon=False, fontsize=8, loc="center left")
fig.tight_layout(); fig.savefig("BALL_SHAPE.png", dpi=160, facecolor=SURF); print("wrote BALL_SHAPE.png")
