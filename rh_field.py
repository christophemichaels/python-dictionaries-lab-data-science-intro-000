"""The arithmetic field (THE_FIELD.md). In u = log x the kernel 1/max(a,b) = (ab)^{-1/2} e^{-|log a - log b|/2} is the Green
function of H0 = -d^2/du^2 + 1/4 (a Yukawa field of mass 1/2). With the Moebius source J_N = sum mu(n) n^{-1/2} delta(u - log n)
the field is phi_N = H0^{-1} J_N and its energy int (phi'^2 + phi^2/4) du equals R(N) = I_M(N) + M(N)^2/N. Checked here, with the
decomposition of the field into its inward and outward solutions and the h-form R(N) = h(N)^2 + sum_{k<N} (h(N) - h(k))^2.
Usage: python3 rh_field.py THE_FIELD.png [N=1e4]"""
import sys, math
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
out = sys.argv[1]; N = int(float(sys.argv[2])) if len(sys.argv) > 2 else 10**4
X = max(N, 10**6)
s = np.ones(X + 1, dtype=bool); s[:2] = False
for p in range(2, int(X**0.5) + 1):
    if s[p]: s[p*p::p] = False
mu = np.ones(X + 1, dtype=np.int8)
for p in np.nonzero(s)[0]: mu[p::p] *= -1; mu[p*p::p*p] = 0
mu[0] = 0
n = np.arange(X + 1, dtype=float); M = np.cumsum(mu.astype(np.int64)); h = np.concatenate([[0.0], np.cumsum(mu[1:] / n[1:])])
def R(Nn): return float(np.sum(M[1:Nn].astype(float)**2 * (1/n[1:Nn] - 1/n[2:Nn+1])) + M[Nn]**2 / Nn)
for Nn in (10**3, 10**4, 10**5, 10**6):
    hform = h[Nn]**2 + float(np.sum((h[Nn] - h[1:Nn])**2))
    inward = R(Nn); outward = float(np.sum((h[Nn] - h[1:Nn])**2)); field = 0.5 * inward + 0.5 * outward + 0.5 * h[Nn]**2
    print(f"N = {Nn:7d}: R(N) = {inward:.6f};  h-form h(N)^2 + sum (h(N)-h(k))^2 = {hform:.6f};  field energy (1/2 inward + 1/2 outward + 1/2 h(N)^2) = {field:.6f}")
# the field for N = 10^4 on a grid in u, piecewise exact: phi = e^{-u/2} M(x) + e^{u/2} (h(N) - h(x)), x = e^u
u = np.linspace(-1, math.log(N) + 3, 20000); x = np.exp(u); k = np.minimum(np.floor(x).astype(int), N); k = np.maximum(k, 0)
Mx = M[k]; hx = h[k]
inw = np.exp(-u/2) * Mx; outw = np.exp(u/2) * (h[N] - hx); phi = inw + outw
dens = 0.5 * Mx.astype(float)**2 * np.exp(-u) + 0.5 * (h[N] - hx)**2 * np.exp(u)        # energy density, exact between pulses
print(f"N = {N}: total energy from the density on the grid {np.trapezoid(dens, u):.6f} against R(N) = {R(N):.6f}")
BLUE, ORANGE, INK, INK2, SURF, GRID = "#2a78d6", "#eb6834", "#0b0b0b", "#52514e", "#fcfcfb", "#e6e5e1"
plt.rcParams.update({"font.size": 10, "axes.titlesize": 11, "axes.edgecolor": GRID, "axes.labelcolor": INK2, "xtick.color": INK2, "ytick.color": INK2,
                     "axes.spines.top": False, "axes.spines.right": False, "figure.facecolor": SURF, "axes.facecolor": SURF, "text.color": INK})
fig, ax = plt.subplots(2, 1, figsize=(13, 7.4), sharex=True)
for a in ax: a.grid(True, color=GRID, lw=0.6); a.set_axisbelow(True)
A = ax[0]
A.plot(u, phi, color=INK, lw=1.0, label="the field φ_N(u) = H₀⁻¹ J_N, N = 10⁴")
A.plot(u, inw, color=BLUE, lw=1.0, alpha=0.9, label="inward solution e^{−u/2} M(e^u): the Möbius state")
A.plot(u, outw, color=ORANGE, lw=1.0, alpha=0.9, label="outward solution e^{u/2} (h(N) − h(e^u)): the mean obligation")
A.axvline(math.log(N), color=INK2, lw=0.8, ls=":"); A.text(math.log(N) + 0.1, A.get_ylim()[1]*0.85 if False else 0.9, "last pulse\nu = log N", fontsize=8, color=INK2)
A.set_ylabel("field"); A.set_title("A.  The arithmetic field of the Möbius source, and its two fundamental solutions", loc="left"); A.legend(frameon=False, fontsize=8.5, loc="upper left")
B = ax[1]
B.plot(u, dens, color=INK, lw=1.0, label="energy density φ′² + φ²/4  (exact between pulses)")
B.plot(u, 0.5 * Mx.astype(float)**2 * np.exp(-u), color=BLUE, lw=1.0, alpha=0.9, label="inward part ½ M(x)²/x")
B.plot(u, 0.5 * (h[N] - hx)**2 * np.exp(u), color=ORANGE, lw=1.0, alpha=0.9, label="outward part ½ x (h(N) − h(x))²")
B.set_xlabel("u = log x"); B.set_ylabel("energy density"); B.set_ylim(0, None)
B.set_title(f"B.  Energy density; total {R(N):.4f} = R(N), shared equally between the two solutions", loc="left"); B.legend(frameon=False, fontsize=8.5, loc="upper left")
fig.suptitle("The field equation H₀ φ = J with H₀ = −d²/du² + ¼: the Green energy of the record as the energy of a Yukawa field of mass ½", fontsize=11)
fig.tight_layout(); fig.savefig(out, dpi=160, bbox_inches="tight"); print("wrote", out)
