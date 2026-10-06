"""The horn of rings (Primes, Folds, and the One Dot, Section 10, Figure 26), reopened: THE_HORN.md.
Ring n = the n-th zero 1/2 + i gamma_n drawn as a circle of circumference gamma_n at height n, radius r_n = gamma_n / 2 pi.
Usage: python3 rh_horn.py data/zeros_6000.txt THE_HORN.png   -> prints the checks, writes data/horn_rings.csv, data/horn_stones.json, the figure.
"""
import sys, math, json, csv
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

zf, out = sys.argv[1], sys.argv[2]
g = np.array([float(l.split()[1]) for l in open(zf)])
n = np.arange(1, len(g) + 1)
r = g / (2 * math.pi)
law_r = lambda rr: rr * (np.log(rr) - 1) + 7 / 8          # Riemann-von Mangoldt smooth count in the radius variable
theta = lambda t: (t / 2) * np.log(t / (2 * math.pi)) - t / 2 - math.pi / 8 + 1 / (48 * t)
wob = n - law_r(r)                                        # rings counted (just above the n-th zero) minus the law
dr = np.diff(r); mid = (r[1:] + r[:-1]) / 2
print(f"rings {len(g)}: circumference {g[0]:.4f} .. {g[-1]:.2f}; radius {r[0]:.4f} .. {r[-1]:.2f}")
print(f"growth per ring: first 1000 mean {dr[:999].mean():.4f} vs 1/log r {(1/np.log(mid[:999])).mean():.4f}; all {dr.mean():.4f} vs {(1/np.log(mid)).mean():.4f}")
print(f"rings minus law: first 1000 in [{wob[:1000].min():.3f}, {wob[:1000].max():.3f}]; all in [{wob.min():.3f}, {wob.max():.3f}], rms {np.sqrt((wob**2).mean()):.3f}")
print(f"Riemann-Siegel phase advance per ring, in units of pi: {np.diff(theta(g)/math.pi).mean():.5f}")
# Gram's law
def gram_point(k):
    lo, hi = 7.0, 1e5
    for _ in range(80):
        m = (lo + hi) / 2
        if theta(m) < k * math.pi: lo = m
        else: hi = m
    return (lo + hi) / 2
K = len(g) - 10
gram = np.array([gram_point(k) for k in range(-1, K)])
counts = np.histogram(g, bins=gram)[0]
gramlaw = dict(one=float((counts == 1).mean()), none=float((counts == 0).mean()), two=float((counts == 2).mean()))
print(f"Gram intervals k=-1..{K-1}: exactly one ring {gramlaw['one']:.3f}, none {gramlaw['none']:.3f}, two {gramlaw['two']:.3f}")
# the stones inside the horn: prime power m enters the cone at a = log(m)/2, where the horizon 2 pi e^{2a} = 2 pi m, i.e. ring radius m
def is_pp(m):
    for p in range(2, m + 1):
        if m % p == 0:
            while m % p == 0: m //= p
            return m == 1
    return False
stones = [m for m in range(2, 1001) if is_pp(m)]
rows = [(m, math.log(m) / 2, 2 * math.pi * m, int((g <= 2 * math.pi * m).sum()), float(law_r(m))) for m in stones]
print("stone m, entry a, horizon 2 pi m, rings resolved, law:")
for row in rows[:10]: print("  %4d  %.4f  %8.2f  %5d  %8.2f" % row)
print("  ...  %4d  %.4f  %8.2f  %5d  %8.2f" % rows[-1])
# curve comparison: zeros of a curve of genus g over F_q are periodic in height, rings at radius (theta_i/2pi + k)/log q: an exact cone of slope 1/(2 g log q)
gq = (3, 3); cone_slope = 1 / (2 * gq[0] * math.log(gq[1]))
print(f"curve genus {gq[0]} over F_{gq[1]}: exact cone, slope 1/(2g log q) = {cone_slope:.4f}; zeta's local slope at ring 1000 = {1/math.log(r[999]):.4f}, at ring 6000 = {1/math.log(r[-1]):.4f}")
with open("data/horn_rings.csv", "w", newline="") as f:
    w = csv.writer(f); w.writerow(["n", "gamma_n", "radius", "rings_minus_law"])
    for i in range(len(g)): w.writerow([int(n[i]), f"{g[i]:.10f}", f"{r[i]:.10f}", f"{wob[i]:.6f}"])
json.dump(dict(stones=[dict(m=m, entry_a=a, horizon=T, rings=N, law=L) for m, a, T, N, L in rows], gram=gramlaw,
               growth=dict(first1000=float(dr[:999].mean()), law_first1000=float((1/np.log(mid[:999])).mean())),
               wobble=dict(min=float(wob.min()), max=float(wob.max()), rms=float(np.sqrt((wob**2).mean()))),
               curve=dict(g=gq[0], q=gq[1], cone_slope=cone_slope)), open("data/horn_stones.json", "w"), indent=1)

# ---------------- figure ----------------
BLUE, ORANGE, INK, INK2, SURF, GRID = "#2a78d6", "#eb6834", "#0b0b0b", "#52514e", "#fcfcfb", "#e6e5e1"
plt.rcParams.update({"font.size": 10, "axes.titlesize": 11, "axes.edgecolor": GRID, "axes.labelcolor": INK2, "xtick.color": INK2, "ytick.color": INK2,
                     "axes.spines.top": False, "axes.spines.right": False, "figure.facecolor": SURF, "axes.facecolor": SURF, "text.color": INK})
fig, ax = plt.subplots(2, 2, figsize=(13.5, 9.2))
for a in ax.flat: a.grid(True, color=GRID, lw=0.6); a.set_axisbelow(True)
# A: profile
A = ax[0, 0]
A.plot(n, r, color=BLUE, lw=2, label="ζ: ring radius γₙ / 2π  (6,000 rings)")
A.plot(n, n * cone_slope, color=ORANGE, lw=2, label=f"curve of genus 3 over 𝔽₃: exact cone, slope 1/(2g log q) = {cone_slope:.3f}")
A.plot(n, n * (r[999] / 1000), color=INK2, lw=1, alpha=0.6, label="the cone through ring 1,000 (what the eye reads)")
A.set_xlabel("ring number n  (zero number)"); A.set_ylabel("ring radius r = circumference / 2π")
A.set_title("A.  The profile: a horn, not a cone", loc="left"); A.legend(frameon=False, fontsize=8.5, loc="upper left")
A.annotate("ring 1,000: r = 225.9", xy=(1000, r[999]), xytext=(1500, 150), fontsize=8.5, color=INK2, arrowprops=dict(arrowstyle="-", color=INK2, lw=0.8))
# B: slope
B = ax[0, 1]
edges = np.unique(np.round(np.geomspace(1, 5999, 40)).astype(int))
bm = [dr[edges[i]-1:edges[i+1]-1].mean() for i in range(len(edges)-1)]; bx = [(edges[i]+edges[i+1])/2 for i in range(len(edges)-1)]
B.plot(n[1:], 1 / np.log(mid), color=BLUE, lw=2, label="ζ: law 1 / log r  (mean spacing in radius units)")
B.plot(bx, bm, "o", color=BLUE, ms=5, mfc=SURF, mew=1.5, label="ζ: measured growth per ring (geometric bins)")
B.axhline(cone_slope, color=ORANGE, lw=2, label="curve of genus 3 over 𝔽₃: constant slope")
B.set_xscale("log"); B.set_ylim(0, 0.8); B.set_xlabel("ring number n  (log scale)"); B.set_ylabel("growth per ring  dr/dn")
B.set_title("B.  The slope: ζ's horn keeps steepening, the curve's cone does not", loc="left"); B.legend(frameon=False, fontsize=8.5)
# C: stones inside the horn
C = ax[1, 0]
ms = np.array([row[0] for row in rows]); Ns = np.array([row[3] for row in rows]); Ls = np.array([row[4] for row in rows])
mm = np.geomspace(2, 1000, 300)
C.plot(mm, law_r(mm), color=INK2, lw=1.2, label="law  m (log m − 1) + 7/8")
C.plot(ms, Ns, "o", color=BLUE, ms=4, label="stone m: rings resolved when m enters, N(2πm)")
for m_ in (2, 3, 5, 7, 11, 31, 97, 997):
    i = list(ms).index(m_); C.annotate(f"{m_}", xy=(ms[i], Ns[i]), xytext=(0, 7), textcoords="offset points", ha="center", fontsize=8, color=INK2)
C.set_xscale("log"); C.set_xlabel("prime power m  (log scale); it enters the cone at a = ½ log m, where the horizon is ring radius m")
C.set_ylabel("rings resolved N(2π m)"); C.set_title("C.  The stones inside the horn: stone m sits at ring radius m", loc="left"); C.legend(frameon=False, fontsize=8.5, loc="upper left")
# D: wobble
D = ax[1, 1]
D.plot(n, wob, color=BLUE, lw=0.7, label="ζ: rings counted minus the law  (= S(γₙ) + ½, the wobble)")
D.axhline(0, color=INK2, lw=0.8); D.set_ylim(-1.5, 2.2)
D.set_xlabel("ring number n"); D.set_ylabel("rings − law"); D.set_title(f"D.  The wobble stays within [{wob.min():.2f}, {wob.max():.2f}] over 6,000 rings", loc="left"); D.legend(frameon=False, fontsize=8.5, loc="upper left")
fig.suptitle("The horn of rings reopened: radius = height / 2π, profile = inverse counting function, and where π enters", fontsize=12)
fig.tight_layout(); fig.savefig(out, dpi=160, bbox_inches="tight"); print("wrote", out)
