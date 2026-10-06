"""Figure for PARITY_DECAY.md: A. N delta_N across horizons (the relative parity defect times N) with its two terms of (10);
B. the completed-packet remainder I_{2,3}(N) = R(N) - (2/3) sum_{m <= N/6, (m,6)=1} mu(m)^2/m for N <= 1e7, with the energy records.
The parity rows are the output of rh_parity_decay.py, Section 2.  Usage: python3 rh_parity_decay_fig.py PARITY_DECAY.png"""
import sys, math, numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
out = sys.argv[1]
rows = [(1024,3.857582,0.098345,1.346349),(1448,3.891684,0.016703,1.439313),(2048,3.985138,0.006928,1.474843),(2896,4.271878,0.169143,1.414384),(4096,4.247554,0.127990,1.447702),
        (5792,4.054816,0.004696,1.496186),(8192,4.222803,0.002770,1.561214),(11585,4.156579,0.000035,1.538062),(16384,4.292974,0.071891,1.516452),(23170,4.325999,0.058979,1.541114),
        (32768,4.258502,0.001721,1.572913),(46340,4.295435,0.002713,1.585781),(65536,4.271713,0.000672,1.578621),(92681,4.297644,0.009376,1.579393),(131072,4.326663,0.004707,1.594676),
        (185363,4.384628,0.002400,1.618369),(262144,4.365506,0.000013,1.613446),(370727,4.449707,0.013830,1.630842),(524288,4.512763,0.019716,1.648134),(741455,4.470483,0.000427,1.651800),
        (1000000,4.606942,0.015244,1.687405),(1048576,4.661279,0.024604,1.698138),(1482910,4.538943,0.008437,1.669061),(2097152,4.713314,0.035588,1.706363),(2965820,4.570259,0.000004,1.689061),
        (4194304,4.627553,0.004699,1.705543),(5931641,4.643700,0.001359,1.714841),(8388608,4.659775,0.000114,1.722025)]
X = 10**7
s = np.ones(X + 1, dtype=bool); s[:2] = False
for p in range(2, int(X**0.5) + 1):
    if s[p]: s[p*p::p] = False
MU = np.ones(X + 1, dtype=np.int8)
for p in np.nonzero(s)[0]: MU[p::p] *= -1; MU[p*p::p*p] = 0
MU[0] = 0
n_ = np.arange(X + 1, dtype=float); M = np.cumsum(MU.astype(np.int64))
inc = np.zeros(X + 1); inc[1:] = (MU[1:].astype(float)**2 + 2*MU[1:]*np.concatenate([[0], M[1:-1]]))/n_[1:]; Rn = np.cumsum(inc)
cop = (MU != 0) & (np.gcd(np.arange(X + 1), 6) == 1); Dcop = np.cumsum(np.where(cop, 1.0, 0.0)/np.maximum(n_, 1))
IC = Rn[1:] - (2/3)*Dcop[(np.arange(1, X + 1))//6]
rec = [1]; best = Rn[1]
for n in range(2, X + 1):
    if Rn[n] > best + 1e-12: rec.append(n); best = Rn[n]
BLUE, ORANGE, INK, INK2, SURF, GRID = "#2a78d6", "#eb6834", "#0b0b0b", "#52514e", "#fcfcfb", "#e6e5e1"
plt.rcParams.update({"font.size": 10, "axes.titlesize": 11, "axes.edgecolor": GRID, "axes.labelcolor": INK2, "xtick.color": INK2, "ytick.color": INK2,
                     "axes.spines.top": False, "axes.spines.right": False, "figure.facecolor": SURF, "axes.facecolor": SURF, "text.color": INK})
fig, ax = plt.subplots(1, 2, figsize=(14, 5.2))
for a in ax: a.grid(True, color=GRID, lw=0.6); a.set_axisbelow(True)
A = ax[0]; Ns = [r[0] for r in rows]
A.plot(Ns, [r[1] for r in rows], color=BLUE, lw=1.6, marker="o", ms=4, label="N·δ_N = N·E_P/A_N")
A.plot(Ns, [r[0]*r[3]/(r[1]/r[0]*r[0]) * 0 + r[3]/r[1]*r[1]*1.0 if False else r[3]*r[1]/(r[2]+r[3]) for r in rows], color=ORANGE, lw=1.2, marker="o", ms=3, label="directional part, 2ab(1−cos θ)·N/A_N")
A.plot(Ns, [r[2]*r[1]/(r[2]+r[3]) for r in rows], color=INK2, lw=1.0, marker="o", ms=3, label="amplitude part, (a−b)²·N/A_N")
A.set_xscale("log"); A.set_xlabel("N"); A.set_ylabel("N·δ_N"); A.set_ylim(0, 5.2)
A.set_title("A.  The relative parity defect times N: no power, a slow drift (3.86 → 4.66)", loc="left"); A.legend(frameon=False, fontsize=8.5, loc="center right")
Bx = ax[1]; xs = np.unique(np.logspace(1, 7, 3000).astype(int))
Bx.plot(xs, IC[xs-1], color=BLUE, lw=1.0, label="I_{2,3}(N) = R(N) − (2/3)·Σ_{m≤N/6,(m,6)=1} μ(m)²/m")
Bx.plot(xs, 1.45 - 0.19*np.log(xs), color=INK2, lw=0.8, ls="--", label="1.45 − 0.19 log N")
Bx.scatter(rec[4:], IC[np.array(rec[4:])-1], s=10, color=ORANGE, zorder=4, label="energy records (first passages of |M|)")
Bx.axhline(0, color=INK, lw=0.8); Bx.axvline(2837, color=INK2, lw=0.8, ls=":"); Bx.text(3200, 0.9, "last N with I > 0:\n2837", fontsize=8.5, color=INK2)
Bx.set_xscale("log"); Bx.set_xlabel("N"); Bx.set_ylabel("I_{2,3}(N)"); Bx.set_ylim(-1.8, 1.6)
Bx.set_title("B.  The completed-packet remainder is negative for every N from 2838 to 10⁷", loc="left"); Bx.legend(frameon=False, fontsize=8.5, loc="upper right")
fig.suptitle("Decay inside the growth: the parity defect, and the sign of the {2,3} packet remainder", fontsize=11)
fig.tight_layout(); fig.savefig(out, dpi=160, bbox_inches="tight"); print("wrote", out)
