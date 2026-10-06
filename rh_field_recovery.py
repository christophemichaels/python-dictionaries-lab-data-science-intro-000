"""Exact recoverability of the arithmetic field (THE_FIELD.md, Sections 5-7). The field phi_N = H0^{-1} J_N of rh_field.py
is evaluated here directly as the sum of its pulses, mu(n) n^{-1/2} e^{-|u - log n|/2}, never through the closed form, and
the closed form's consequences are then checked against it: the recovery law M(e^u) = e^{u/2}(phi/2 - phi'),
h(N) - h(e^u) = e^{-u/2}(phi/2 + phi'); the endpoint values phi_N(0) = h(N), phi_N(log N) = M(N)/sqrt N; the two channels
A = e^{-u/2} M, B = e^{u/2}(h(N) - h) each integrating to R(N) over the whole line; the successor identity
R(n) - R(n-1) = (mu(n)^2 + 2 mu(n) M(n-1))/n with the interaction read off the field already present, phi_{n-1}(log n) =
M(n-1)/sqrt n; and the damped energies int_0^{log N} e^{-2 delta u} |m(u)|^2 du beside the undamped I_M(N).
Usage: python3 rh_field_recovery.py THE_FIELD_RECOVERY.png"""
import sys, math
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
out = sys.argv[1]
X = 10**7
s = np.ones(X + 1, dtype=bool); s[:2] = False
for p in range(2, int(X**0.5) + 1):
    if s[p]: s[p*p::p] = False
mu = np.ones(X + 1, dtype=np.int8)
for p in np.nonzero(s)[0]: mu[p::p] *= -1; mu[p*p::p*p] = 0
mu[0] = 0
n = np.arange(X + 1, dtype=float); M = np.cumsum(mu.astype(np.int64)); h = np.concatenate([[0.0], np.cumsum(mu[1:] / n[1:])])
def R(Nn): return float(np.sum(M[1:Nn].astype(float)**2 * (1/n[1:Nn] - 1/n[2:Nn+1])) + M[Nn]**2 / Nn)
# --- the field and its slope, directly from the pulses (independent of the closed form) ---
def field(Nn, u):
    """phi_N(u) and phi_N'(u) as the sum over n <= N of mu(n) n^{-1/2} e^{-|u - log n|/2}; u an array of points off the pulses."""
    nn = np.arange(1, Nn + 1); w = mu[1:Nn+1] / np.sqrt(nn); ln = np.log(nn)
    phi = np.empty(len(u)); dphi = np.empty(len(u))
    for i, ui in enumerate(u):
        d = ui - ln; g = np.exp(-np.abs(d) / 2)
        phi[i] = np.sum(w * g); dphi[i] = np.sum(w * g * (-0.5 * np.sign(d)))
    return phi, dphi
N = 10**4
k = np.arange(1, N); um = np.log(k + 0.5)                       # midpoints between consecutive pulses
phi, dphi = field(N, um)
M_rec = np.exp(um / 2) * (phi / 2 - dphi); T_rec = np.exp(-um / 2) * (phi / 2 + dphi)
print(f"recovery at N = {N}, {len(um)} midpoints: max |M recovered - M| = {np.max(np.abs(M_rec - M[k])):.2e};  "
      f"max |tail recovered - (h(N) - h)| = {np.max(np.abs(T_rec - (h[N] - h[k]))):.2e}")
for Nn in (10**3, 10**4):
    p0, _ = field(Nn, np.array([-1e-9])); pL, _ = field(Nn, np.array([math.log(Nn) + 1e-9]))
    print(f"endpoints at N = {Nn}: phi(0) = {p0[0]:.12f} against h(N) = {h[Nn]:.12f};  phi(log N) = {pL[0]:.12f} against M(N)/sqrt N = {M[Nn]/math.sqrt(Nn):.12f}")
# --- the two channels over the whole line, and the successor identity ---
print("channels over the whole line (truncated sums): int A^2 = I_M(N) + M(N)^2/N, int B^2 = h(N)^2 + sum (h(N)-h(k))^2")
for Nn in (10**3, 10**4, 10**5, 10**6, 10**7):
    IA = R(Nn); IB = h[Nn]**2 + float(np.sum((h[Nn] - h[1:Nn])**2))
    print(f"  N = {Nn:8d}: int A^2 = {IA:.6f}   int B^2 = {IB:.6f}   R(N) = {R(Nn):.6f}")
Rn = np.concatenate([[0.0], np.cumsum(((mu[1:].astype(float))**2 + 2 * mu[1:] * np.concatenate([[0], M[1:-1]])) / n[1:])])
err = np.max(np.abs(Rn[[10**3, 10**4, 10**5, 10**6, 10**7]] - np.array([R(v) for v in (10**3, 10**4, 10**5, 10**6, 10**7)])))
print(f"successor identity summed from n = 1: max |sum (mu^2 + 2 mu M(n-1))/n - R(N)| over N = 10^3..10^7 is {err:.2e}")
for nn in (30, 31, 1009, 9973):
    p, _ = field(nn - 1, np.array([math.log(nn)]))
    print(f"  interaction at n = {nn}: phi_(n-1)(log n) = {p[0]:+.9f}  against M(n-1)/sqrt n = {M[nn-1]/math.sqrt(nn):+.9f};  "
          f"R(n) - R(n-1) = {R(nn) - R(nn-1):+.9f} = (mu^2 + 2 mu M(n-1))/n = {(mu[nn]**2 + 2*mu[nn]*M[nn-1])/nn:+.9f}")
# --- damped energies ---
print("damped energy  int_0^{log N} e^{-2 delta u} |m(u)|^2 du,  m(u) = e^{-u/2} M(e^u)   (delta = 0 is the undamped I_M(N))")
print("  N        delta=0      delta=0.05   delta=0.1    delta=0.25")
for Nn in (10**3, 10**4, 10**5, 10**6, 10**7):
    row = []
    for d in (0.0, 0.05, 0.1, 0.25):
        a = 1 + 2 * d; row.append(float(np.sum(M[1:Nn].astype(float)**2 * (n[1:Nn]**(-a) - n[2:Nn+1]**(-a))) / a))
    print(f"  {Nn:8d}  " + "  ".join(f"{v:10.6f}" for v in row))
# --- figure ---
BLUE, ORANGE, INK, INK2, SURF, GRID = "#2a78d6", "#eb6834", "#0b0b0b", "#52514e", "#fcfcfb", "#e6e5e1"
plt.rcParams.update({"font.size": 10, "axes.titlesize": 11, "axes.edgecolor": GRID, "axes.labelcolor": INK2, "xtick.color": INK2, "ytick.color": INK2,
                     "axes.spines.top": False, "axes.spines.right": False, "figure.facecolor": SURF, "axes.facecolor": SURF, "text.color": INK})
fig, ax = plt.subplots(3, 1, figsize=(13, 11))
for a in ax: a.grid(True, color=GRID, lw=0.6); a.set_axisbelow(True)
K = 120; kk = np.arange(1, K + 1)
A = ax[0]
A.step(np.log(np.concatenate([kk, [K + 1]])), np.concatenate([M[1:K+1], [M[K]]]), where="post", color=INK, lw=1.2, label="M(x), from the sieve")
sel = k <= K
A.scatter(um[sel], M_rec[sel], s=22, color=BLUE, zorder=3, label="e^{u/2}(φ/2 − φ′) at the midpoints, from the pulses only")
A.set_ylabel("M"); A.set_title("A.  Recovery of the enclosed mass M(e^u) = e^{u/2}(φ_N/2 − φ_N′), N = 10⁴, shown to x = 120", loc="left")
A.legend(frameon=False, fontsize=8.5, loc="lower left")
B = ax[1]
B.step(np.log(np.concatenate([kk, [K + 1]])), np.concatenate([h[N] - h[1:K+1], [h[N] - h[K]]]), where="post", color=INK, lw=1.2, label="h(N) − h(x), from the sieve")
B.scatter(um[sel], T_rec[sel], s=22, color=ORANGE, zorder=3, label="e^{−u/2}(φ/2 + φ′) at the midpoints, from the pulses only")
B.axhline(0, color=INK2, lw=0.6)
B.set_ylabel("remaining tail"); B.set_title("B.  Recovery of the remaining obligation h(N) − h(e^u) = e^{−u/2}(φ_N/2 + φ_N′)", loc="left")
B.legend(frameon=False, fontsize=8.5, loc="lower right")
for a in (A, B): a.set_xlabel("u = log x")
C = ax[2]
nn = 31; ug = np.linspace(0, 5.2, 1200)
p_old, _ = field(nn - 1, ug); pulse = mu[nn] / math.sqrt(nn) * np.exp(-np.abs(ug - math.log(nn)) / 2)
C.plot(ug, p_old, color=BLUE, lw=1.6, label="φ₃₀, the field already present")
C.plot(ug, pulse, color=ORANGE, lw=1.6, label="the arriving pulse μ(31)·31^{−1/2} e^{−|u − log 31|/2}")
C.plot(ug, p_old + pulse, color=INK, lw=1.0, ls="--", label="φ₃₁ = φ₃₀ + pulse")
C.scatter([math.log(nn)], [M[nn-1] / math.sqrt(nn)], s=40, color=BLUE, zorder=4, edgecolor=SURF, linewidth=1.5)
C.annotate(f"φ₃₀(log 31) = M(30)/√31 = {M[nn-1]/math.sqrt(nn):+.4f}\nthe whole earlier source, read where the pulse lands",
           (math.log(nn), M[nn-1] / math.sqrt(nn)), xytext=(0.75, -0.36), fontsize=8.5, color=INK2, arrowprops=dict(arrowstyle="-", color=INK2, lw=0.7))
C.axvline(math.log(nn), color=INK2, lw=0.8, ls=":")
C.set_xlabel("u = log x"); C.set_ylabel("field")
C.set_title(f"C.  The successor as a field event: R(31) − R(30) = (μ(31)² + 2 μ(31) M(30))/31 = {(mu[nn]**2 + 2*mu[nn]*M[nn-1])/nn:+.5f}", loc="left")
C.legend(frameon=False, fontsize=8.5, loc="lower right")
fig.suptitle("Exact recoverability: the arithmetic builds the field, and the field's value and slope give the arithmetic back", fontsize=11)
fig.tight_layout(); fig.savefig(out, dpi=160, bbox_inches="tight"); print("wrote", out)
