"""Figure for GLOBAL_CLOSURE.md: (A) the Mobius potential f_N(u) = V_N(e^u) at N = 39 (M(39) = 0); (B) its fold onto the
Connes-Consani orbit C_5 = R/(log 5)Z over one period, with the kinks at the classes of the unfinished 5-packets and the
smooth classes of the complete ones; (C) the chain-refined ledger (37)/(42) from N = 1e7 (data/global_closure.json)."""
import json, math, matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt, numpy as np
BLUE, ORANGE, INK, MUTED, GRID, SURF = "#2a78d6", "#eb6834", "#0b0b0b", "#52514e", "#e6e5e1", "#fcfcfb"
plt.rcParams.update({"font.size": 10, "axes.edgecolor": MUTED, "axes.labelcolor": INK, "xtick.color": MUTED, "ytick.color": MUTED, "text.color": INK, "axes.titlecolor": INK})
N, p = 39, 5
mu = [0, 1]
for n in range(2, N + 1):
    m, f, k = n, [], 2
    while k*k <= m:
        if m % k == 0:
            m //= k
            if m % k == 0: f = None; break
            f.append(k)
        else: k += 1
    if f is None: mu.append(0); continue
    if m > 1: f.append(m)
    mu.append((-1)**len(f))
M = np.cumsum(mu); h = np.cumsum([0] + [mu[n]/n for n in range(1, N + 1)])
def V(x):
    x = np.asarray(x, float); k = np.minimum(np.floor(x), N).astype(int); return (M[k] + x*(h[N] - h[k]))*(x > 0)
def fold(x): return sum(V(p**k*x) for k in range(-60, 8))
fig, axs = plt.subplots(1, 3, figsize=(14.5, 4.4), facecolor=SURF)
for a in axs: a.set_facecolor(SURF); a.grid(True, color=GRID, lw=0.8); a.set_axisbelow(True); [a.spines[s].set_visible(False) for s in ("top", "right")]
ax, bx, cx = axs
u = np.linspace(-1.5, math.log(N) + 1.2, 3000); ax.plot(u, V(np.exp(u)), color=BLUE, lw=1.6)
sq = [n for n in range(1, N + 1) if mu[n]]
ax.plot([math.log(n) for n in sq], V(np.array(sq, float)), "o", ms=3.2, color=INK)
ax.axhline(0, color=MUTED, lw=0.8); ax.set_xlabel("u = log x"); ax.set_ylabel("f₃₉(u) = V₃₉(eᵘ)")
ax.set_title("A. The Möbius potential at N = 39, M(39) = 0", loc="left", fontsize=10.5)
ax.annotate("slope h(N) at the source end", xy=(-1.2, V(np.exp(-1.2))), xytext=(-1.45, -0.45), fontsize=8.5, arrowprops=dict(arrowstyle="-", color=MUTED, lw=0.8))
ax.annotate("value M(N) = 0 beyond log N", xy=(math.log(N) + 0.8, 0.0), xytext=(1.6, -1.25), fontsize=8.5, arrowprops=dict(arrowstyle="-", color=MUTED, lw=0.8)); ax.set_ylim(-1.75, 0.3)
uu = np.linspace(0, math.log(p), 4000); F = fold(np.exp(uu)); bx.plot(uu, F, color=BLUE, lw=1.6)
unf = [m for m in sq if m % p and p*m > N]; comp = [m for m in sq if m % p and p*m <= N]
rep = lambda m: math.log(m) - math.floor(math.log(m)/math.log(p))*math.log(p)
bx.plot([rep(m) for m in unf], fold(np.exp(np.array([rep(m) for m in unf]))), "o", ms=4, color=ORANGE, label="unfinished parents (kinks, order −μ(m))")
bx.plot([rep(m) for m in comp], fold(np.exp(np.array([rep(m) for m in comp]))), "o", ms=5, mfc="none", color=INK, label="complete packets (no kink, constant μ(m))")
bx.set_xlim(0, math.log(p)); bx.set_xlabel("u on C₅ = ℝ/(log 5)ℤ"); bx.set_ylabel("fold Σₖ V₃₉(5ᵏ eᵘ)")
bx.set_title("B. Fold onto C₅: degree 0, kinks only at unfinished packets", loc="left", fontsize=10.5); bx.set_ylim(-2.58, -1.88)
bx.legend(frameon=False, fontsize=8.2, loc="upper right")
d = json.load(open("data/global_closure.json")); L = [r for r in d["ledger"] if r["N"] == 10**7][0]
J = len(L["dhat"]); xj = np.arange(J); w = 0.38
cx.bar(xj - w/2, L["dhat"], w, color=BLUE, label="Δ̂ⱼ (projected increment)")
cx.bar(xj + w/2, L["d23"], w, color=ORANGE, label="d₂₃(Nⱼ) (diagonal, ≤ 4/5)")
cx.plot(xj, L["eps"], "o", ms=4, color=INK, label="εⱼ ≥ 0 (storage), Σ ≤ 2^(−4/3)")
cx.axhline(0, color=MUTED, lw=0.8); cx.set_xticks(xj); cx.set_xticklabels([f"{n:.0e}".replace("e+0", "e") if n >= 1000 else str(n) for n in L["chain"][:-1]], fontsize=8)
cx.set_xlabel("chain Nⱼ from 10⁷ (Nⱼ₊₁ = ⌊Nⱼ/6⌋)"); cx.set_ylabel("increment")
cx.set_title("C. Ledger (37) on the chain from 10⁷: B̂₂₃ = 0", loc="left", fontsize=10.5)
cx.legend(frameon=False, fontsize=8.2, loc="upper left")
fig.tight_layout(); fig.savefig("GLOBAL_CLOSURE.png", dpi=160, facecolor=SURF); print("wrote GLOBAL_CLOSURE.png")
