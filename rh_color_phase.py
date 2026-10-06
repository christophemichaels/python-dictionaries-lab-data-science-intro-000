"""The colour phase (COLORS_AND_DESCENT.md, Sections 2-3). E(t) = ||sum_k e^{ikt} p_k||^2_B, the energy of the colour source at phase t
(t = 0 the positive control, t = pi the Moebius source), compared with the Selberg-Delange prediction
   E_SD(t) = |C(z)/(z Gamma(z))|^2 int_0^N |(log N)^z - (log x)^z|^2 dx,  z = e^{it},  C(z) = prod_p (1 + z/p)(1 - 1/p)^z,
which vanishes at z = -1 where 1/Gamma has its zero; and the energies E(f_{C_y}) of the Moebius source with the primes <= y removed,
along y, with the per-prime ratio E((I - D_p) v)/E(v).  Usage: python3 rh_color_phase.py COLOR_PHASE.png"""
import sys, math, cmath, numpy as np

def gamma_c(z):
    """complex Gamma by Lanczos (g = 7, n = 9)."""
    g = 7; coef = [0.99999999999980993, 676.5203681218851, -1259.1392167224028, 771.32342877765313, -176.61502916214059,
                   12.507343278686905, -0.13857109526572012, 9.9843695780195716e-6, 1.5056327351493116e-7]
    if z.real < 0.5: return cmath.pi/(cmath.sin(cmath.pi*z)*gamma_c(1 - z))
    z -= 1; x = coef[0]
    for i in range(1, g + 2): x += coef[i]/(z + i)
    t = z + g + 0.5; return cmath.sqrt(2*cmath.pi)*t**(z + 0.5)*cmath.exp(-t)*x
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
out = sys.argv[1]
X = 10**6
s = np.ones(X + 1, dtype=bool); s[:2] = False
for p in range(2, int(X**0.5) + 1):
    if s[p]: s[p*p::p] = False
primes = np.nonzero(s)[0]
MU = np.ones(X + 1, dtype=np.int8)
for p in primes: MU[p::p] *= -1; MU[p*p::p*p] = 0
MU[0] = 0
OM = np.zeros(X + 1, dtype=np.int8)
for p in primes: OM[p::p] += 1
def tails(c): return c[::-1].cumsum()[::-1]
def C_of(z, P=10**6):
    ps = primes[primes <= P].astype(float); return complex(np.prod((1 + z/ps) * np.exp(z*np.log(1 - 1/ps))))
def E_SD(t, N):
    z = cmath.exp(1j*t); x = np.arange(2, N + 1, dtype=float); llx = np.log(np.log(x)); llN = math.log(math.log(N))
    integrand = np.abs(np.exp(z*llN) - np.exp(z*llx))**2        # |(log N)^z - (log x)^z|^2 on [2, N]; x in [0, 2) adds at most 2 |log N|^(2 Re z)
    integrand = np.concatenate([[abs(np.exp(z*llN))**2]*1, integrand])
    pref = abs(C_of(z)/(z*gamma_c(z)))**2 if abs(z + 1) > 1e-12 else 0.0
    return pref*float(np.sum(integrand))
ts = np.linspace(0, math.pi, 181); prof = {}
print("== the colour phase: E(t) against Selberg-Delange ==")
for N in (16384, 10**5, 10**6):
    nn = np.arange(1, N + 1); sq = MU[1:N+1] != 0; om = OM[1:N+1]; K = int(om[sq].max())
    T = [tails(np.where(sq & (om == k), 1.0/nn, 0.0)) for k in range(K + 1)]
    G = np.array([[float(np.dot(T[k], T[l])) for l in range(K + 1)] for k in range(K + 1)])
    E = np.array([sum(G[k, l]*math.cos((k - l)*t) for k in range(K + 1) for l in range(K + 1)) for t in ts]); prof[N] = E
    sd = np.array([E_SD(t, N) for t in ts])
    print(f"   N = {N:7d}, colours 0..{K}: E(0) = {E[0]:.1f} (prediction 2N(6/pi^2)^2 = {2*N*(6/math.pi**2)**2:.1f}; SD {sd[0]:.1f});  mean over t = {np.mean(E):.1f} = diag {np.trace(G):.1f};  E(pi) = {E[-1]:.4f} = R(N) = {float(np.dot(tails(MU[1:N+1]/nn), tails(MU[1:N+1]/nn))):.4f}")
    for t in (math.pi/4, math.pi/2, 3*math.pi/4, 0.9*math.pi, 0.95*math.pi, 0.99*math.pi):
        i = int(round(t/math.pi*180)); print(f"      t = {t/math.pi:.2f} pi: E(t) = {E[i]:12.4f}   Selberg-Delange main term {sd[i]:12.4f}   ratio {E[i]/sd[i] if sd[i] > 0 else float('nan'):.3f}")
    prof[(N, "sd")] = sd
print("== prime packets: E(f_{C_y}), the Moebius source with the primes <= y removed, N = 16384 ==")
N = 16384; nn = np.arange(1, N + 1); c = MU[1:N+1]/nn
def energy(v): t = tails(v); return float(np.dot(t, t))
ys = [1, 2, 3, 5, 7, 11, 13, 20, 50, 100, 300, 1000, 3000, 8192, 16384]
prev = None
for y in ys:
    q = np.ones(N, dtype=bool)
    for p in primes[primes <= y]: q[p-1::p] = False
    f = np.where(q, c, 0.0); e = energy(f)
    print(f"   y = {y:5d}: E(f_C) = {e:10.4f}" + (f"   ratio to previous {e/prev:.3f}" if prev else "") + f"   (primes removed: {int(np.sum(primes <= y))})")
    prev = e
# figure
BLUE, ORANGE, INK, INK2, SURF, GRID = "#2a78d6", "#eb6834", "#0b0b0b", "#52514e", "#fcfcfb", "#e6e5e1"
plt.rcParams.update({"font.size": 10, "axes.titlesize": 11, "axes.edgecolor": GRID, "axes.labelcolor": INK2, "xtick.color": INK2, "ytick.color": INK2,
                     "axes.spines.top": False, "axes.spines.right": False, "figure.facecolor": SURF, "axes.facecolor": SURF, "text.color": INK})
fig, ax = plt.subplots(1, 1, figsize=(9, 5.2)); ax.grid(True, color=GRID, lw=0.6); ax.set_axisbelow(True)
for N, col in ((16384, BLUE), (10**6, ORANGE)):
    ax.plot(ts/math.pi, prof[N]/N, color=col, lw=1.6, label=f"E(t)/N, N = {N}")
    ax.plot(ts/math.pi, prof[(N, "sd")]/N, color=col, lw=1.0, ls="--", label=f"Selberg–Delange main term, N = {N}")
ax.set_yscale("log"); ax.set_xlabel("colour phase t/π   (t = 0: all colours positive;  t = π: the Möbius parity)"); ax.set_ylabel("energy / N")
ax.set_title("The energy of the colour source as a function of its phase: the parity phase is the zero of 1/Γ(e^{it})", loc="left")
ax.legend(frameon=False, fontsize=8.5, loc="lower left")
fig.tight_layout(); fig.savefig(out, dpi=160, bbox_inches="tight"); print("wrote", out)
