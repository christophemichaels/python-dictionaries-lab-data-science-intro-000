"""Figure FUNCTION_FIELD_CONE.png: the light cone of a curve over a finite field beside the number-field cone (FUNCTION_FIELD_CONE.md).
Usage: python3 rh_ff_cone_fig.py . FUNCTION_FIELD_CONE.png   (reads data/ff_cone_*.json, floor_grid_K40.csv, floor_grid_arb.csv)"""
import json, math, csv, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

R = sys.argv[1]          # repo root
OUT = sys.argv[2]        # output png
MAG = "#c2178f"; GREY = "#555555"
plt.rcParams.update({"font.size": 10, "axes.titlesize": 11, "axes.spines.top": False, "axes.spines.right": False})

fig, (axA, axB, axC) = plt.subplots(1, 3, figsize=(16, 5.2), gridspec_kw=dict(width_ratios=[1.15, 1, 1]))

# ---------- Panel A: the cone over the genus-3 curve over F_3 ----------
d3 = json.load(open(f"{R}/data/ff_cone_g3_q3.json"))
q, g = d3["q"], d3["g"]; lq = math.log(q)
th = np.array(d3["thetas"])
def N_of(m): return q**m + 1 - (q**(m/2) * np.sum(np.exp(1j*m*th))).real
b = {}
for m in range(1, 2*g + 3):
    s = sum(d*b[d] for d in b if m % d == 0)
    b[m] = int(round((N_of(m) - s)/m))
xmax = (2*g + 2.6)*lq; ymax = (g + 1.35)*lq
# cone interior: log n < 2a
xs = np.linspace(0, xmax, 200)
axA.fill_between(xs, xs/2, ymax, color=MAG, alpha=0.07, lw=0)
axA.plot(xs, xs/2, color=MAG, lw=2, label="edge  log n = 2a")
# closed points: a vertical line from its entry a = d log q / 2 upward, count written at the base
for d in range(1, 2*g + 3):
    x = d*lq
    axA.plot([x, x], [x/2, ymax], color=GREY, lw=0.9, alpha=0.8)
    axA.plot([x], [x/2], marker="o", color=GREY, ms=4)
    axA.text(x, x/2 - 0.13*lq, f"{b[d]}", ha="center", va="top", fontsize=8.5, color=GREY)
    axA.text(x, ymax + 0.03*lq, f"deg {d}", ha="center", va="bottom", fontsize=8, color=GREY)
# horizon and determination
axA.axhline(g*lq, color=MAG, lw=2.2, ls="-")
axA.axhline(g*lq/2, color=MAG, lw=1.2, ls="--")
axA.fill_between([0, xmax], g*lq, ymax, color="#dddddd", alpha=0.55, lw=0)
axA.text(0.08*lq, g*lq + 0.07*lq, "horizon  a = g log q :  floor = 0,  $\\mathcal{P}_X(a)=\\{\\nu_X\\}$  (the cone has closed)", fontsize=8.8, color=MAG, va="bottom")
axA.text(0.08*lq, g*lq/2 + 0.07*lq, "a = g log q / 2 :  N₁ … N_g inside; with the genus known the zeros are fixed", fontsize=8.3, color=MAG, va="bottom")
# floors written along the left axis at a = D log q
axA.set_xlim(0, xmax); axA.set_ylim(0, ymax + 0.35*lq)
print("b_d:", b)
axA2 = axA.twinx()
axA2.set_ylim(axA.get_ylim())
levels = [r for r in d3["floors"] if r["D"] <= g + 1]
axA2.set_yticks([r["a"] for r in levels])
axA2.set_yticklabels([f"D = {r['D']}:  λ = {max(r['lam'], 0.0):.2f}" for r in levels], fontsize=8, color="#1f4e79")
axA2.spines["top"].set_visible(False); axA2.tick_params(axis="y", length=0)
axA2.set_ylabel("floor of the Weil form at support D", color="#1f4e79", fontsize=9)
axA.set_xlabel("log N(P) = (deg P) · log q     (closed points P of the curve; counts b_d at their entry)")
axA.set_ylabel("half-width a")
axA.set_title(f"A.  The cone over the curve  y² = 2+x²+x³+2x⁵+x⁷  over 𝔽₃  (genus {g})", loc="left")
axA.legend(loc="lower right", frameon=False, fontsize=8.5)

# ---------- Panel B: function-field floors ----------
cols = {1: "#2a9d8f", 2: "#1f4e79", 3: "#8a3ffc"}
for gg in (1, 2, 3):
    d = json.load(open(f"{R}/data/ff_cone_g{gg}_q3.json"))
    a = [r["a"] for r in d["floors"]]; le = [r["lam_even"] for r in d["floors"]]; lo = [r["lam_odd"] for r in d["floors"]]
    axB.plot(a, le, "-o", color=cols[gg], ms=5, lw=1.6, label=f"g = {gg}, even sector")
    axB.plot(a[1:], lo[1:], "--s", color=cols[gg], ms=4, lw=1.1, alpha=0.75, label=f"g = {gg}, odd sector")
    axB.plot([gg*lq], [0], marker="o", ms=11, mfc="none", mec=MAG, mew=1.8)
    axB.axvline(gg*lq, color=MAG, lw=0.8, alpha=0.5)
    axB.text(gg*lq, -0.45, f"closes\ng = {gg}", ha="center", va="top", fontsize=7.5, color=MAG)
axB.set_xlim(-0.2, 5.0); axB.set_ylim(-1.2, 9.4)
axB.set_xticks([D*lq for D in range(5)]); axB.set_xticklabels([f"D={D}\n{D*lq:.2f}" for D in range(5)], fontsize=8)
axB.set_xlabel("half-width a = D · log q   (q = 3)")
axB.set_ylabel("floor of the Weil form at support D")
axB.set_title("B.  Function field: the floor reaches zero exactly at a = g log q", loc="left")
axB.legend(frameon=False, fontsize=8, ncol=2, loc="upper right")

# ---------- Panel C: the number-field floor ----------
a40, l40 = [], []
for row in csv.DictReader(open(f"{R}/floor_grid_K40.csv")):
    a40.append(float(row["a"])); l40.append(float(row["lambda0"]))
aarb, larb = [], []
for row in csv.DictReader(open(f"{R}/floor_grid_arb.csv")):
    aarb.append(float(row["a"])); larb.append(float(row["lambda0"]))
axC.semilogy(a40, l40, color="#1f4e79", lw=1.6, label="odd Weil form, K = 40 modes (floor_grid_K40.csv)")
axC.semilogy(aarb, larb, "o", color=MAG, ms=5, label="certified enclosures (floor_grid_arb.csv)")
for n in (2, 3, 4, 5, 7, 8, 9, 11):
    an = math.log(n)/2
    axC.axvline(an, color=GREY, lw=0.6, ls=":", alpha=0.8)
    axC.text(an, 2.5e-100, f"{n}", ha="center", va="bottom", fontsize=8, color=GREY)
axC.set_xlim(0.25, 1.6); axC.set_ylim(1e-100, 5)
axC.set_xlabel("half-width a   (prime powers n enter at a = ½ log n, marked above)")
axC.set_ylabel("floor λ(a)   (log scale)")
axC.set_title("C.  Number field: the floor decays and is positive at every support computed", loc="left")
axC.text(0.27, 1e-74, "no finite horizon: the zeros are infinitely many,\nso the form has no null vector at any finite a.\nPositivity for all a is Weil's criterion, i.e. RH.", fontsize=8.5, color=GREY, va="bottom")
axC.legend(frameon=False, fontsize=8, loc="upper right")

fig.suptitle("The paperwork of the curves projected onto the cone: where the genus is finite the cone closes; for ζ it never does", fontsize=12, y=1.0)
fig.tight_layout()
fig.savefig(OUT, dpi=160, bbox_inches="tight")
print("wrote", OUT)
