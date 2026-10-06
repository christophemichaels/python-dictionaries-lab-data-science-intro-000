"""The horizon round trip (HORIZON_ROUND_TRIP.md). Verifies the supplied lemmas: the two-sided bound (2) with the constant
1 + 2 log(3N); the localization (6) ||x||_B <= A_L ||D_N x||_B on the tail n >= L; the diagonal bounds (8), (9); the kernel (13);
the high-mode bound. Computes the exact cost accounting of the localization at L = 21 (head energy, cross term 2 M(L-1)(h(N) - h(L-1)),
tail energy) to 1e7; the diagnostics N^2/J^2 (D, O, total) for the unweighted and log weights by FFT for N = 1e3..1e7 and all
J < J0(N) = ceil(N / log^2(eN)), with the random-sign and all-positive controls; the dyadic Gram blocks (11) at N = 2000.
Usage: python3 rh_round_trip.py HORIZON_ROUND_TRIP.png"""
import sys, math, json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
out = sys.argv[1]
rng = np.random.default_rng(11)
def sieve(N):
    s = np.ones(N + 1, dtype=bool); s[:2] = False
    for p in range(2, int(N**0.5) + 1):
        if s[p]: s[p*p::p] = False
    primes = np.nonzero(s)[0]
    mu = np.ones(N + 1, dtype=np.int8)
    for p in primes: mu[p::p] *= -1; mu[p*p::p*p] = 0
    mu[0] = 0
    return primes, mu
# ---------- lemmas at N = 1000, 2000, 4000 ----------
print("== lemmas ==")
for N in (1000, 2000, 4000):
    primes, mu = sieve(N); nn = np.arange(1, N + 1); jj = nn; D = 2*N + 1
    theta = (2*jj - 1) * math.pi / D; lam = 1 / (4 * np.sin(theta/2)**2)
    V = (2 / math.sqrt(D)) * np.sin(np.outer(nn, theta)); c = mu[1:N+1] / nn; b = V.T @ c
    S = np.cumsum(b**2); K = float(np.max(lam * S)); R = float(np.sum(lam * b**2))
    const2 = 1 + 2*math.log(3*N); actual = 1 + math.log(lam[0]/lam[-1])
    # (6) localization on the tail n >= L
    L = 21; g = 1 / (np.log(nn) - 1); Bt = np.minimum.outer(nn[L-1:], nn[L-1:]).astype(float); Lb = np.linalg.cholesky(Bt)
    A = np.linalg.solve(Lb, (Lb.T @ np.diag(g[L-1:])).T).T   # Lb^T G Lb^{-T}
    normL = float(np.linalg.norm(Lb.T @ np.diag(g[L-1:]) @ np.linalg.inv(Lb.T), 2)); AL = 1/(math.log(L)-1) + 2/(math.log(L)-1)**2
    # (8), (9) diagonal bounds, (13) kernel
    for J in (10, 100):
        Kdiag = np.sum(V[:, :J]**2, axis=1); D1 = float(np.sum((mu[1:N+1]**2 / nn**2) * Kdiag)); bound8 = 2*math.pi**2*J**2 / D**2
        w = np.where(nn >= 2, np.log(nn) - 1, 0.0); Dw = float(np.sum((mu[1:N+1]**2 * w**2 / nn**2) * Kdiag)); ell = math.log(math.e*N); bound9 = ell**2 * bound8
        a, bb = 17, 40; H = lambda t: 2*J if abs(math.sin(t)) < 1e-15 else math.sin(2*J*t)/math.sin(t)
        k13 = (H(math.pi*(a-bb)/D) - H(math.pi*(a+bb)/D)) / D; kdir = float(np.sum(V[a-1, :J] * V[bb-1, :J]))
        if J == 10: print(f"N = {N}: (2) multiplier actual 1 + log(l1/lN) = {actual:.3f} <= 1 + 2 log(3N) = {const2:.3f};  K <= R <= mult K: {K:.4f} <= {R:.4f} <= {const2*K:.3f};  (6) tail norm at L = 21: {normL:.4f} <= A_21 = {AL:.4f}  (full inverse norm 25.45)")
        print(f"     J = {J:3d}: (8) D^(1) = {D1:.3e} <= {bound8:.3e} (ratio {D1/bound8:.3f});  (9) D^(log) = {Dw:.3e} <= {bound9:.3e} (ratio {Dw/bound9:.3f});  (13) kernel at (17,40): {k13:.6e} vs direct {kdir:.6e}")
    # high-mode bound for F = -z^(log): sum_{j>J0} lambda_j F_j^2 <= lambda_{J0+1} * Clog
    w = np.where(nn >= 2, np.log(nn) - 1, 0.0); F = V.T @ (c * w); ell = math.log(math.e*N); J0 = math.ceil(N / ell**2)
    Clog = float(np.sum((np.log(np.arange(2, 10**6)) - 1)**2 / np.arange(2, 10**6)**2.0))
    print(f"     high modes j > J0 = {J0}: sum lambda_j F_j^2 = {float(np.sum(lam[J0:] * F[J0:]**2)):.4f} <= (9/4) Clog ell^4 = {2.25*Clog*ell**4:.1f};  Clog = {Clog:.4f};  R^log(N) = {float(np.sum(lam*F**2)):.4f}")
# ---------- exact cost accounting of the localization at L = 21, to 1e7 ----------
print("== localization accounting, L = 21: R(N) = R(20) + 2 M(20) (h(N) - h(20)) + R_tail(N) ==")
X = 10**7; primes, mu = sieve(X); n = np.arange(X + 1, dtype=float); M = np.cumsum(mu.astype(np.int64)); h = np.concatenate([[0.0], np.cumsum(mu[1:]/n[1:])])
def R(Nn): return float(np.sum(M[1:Nn].astype(float)**2 * (1/n[1:Nn] - 1/n[2:Nn+1])) + M[Nn]**2 / Nn)
L = 21; head = R(L-1); print(f"   head R(20) = {head:.6f}, M(20) = {M[20]}, h(20) = {h[20]:.6f}")
for Nn in (10**3, 10**4, 10**5, 10**6, 10**7):
    cross = 2 * M[L-1] * (h[Nn] - h[L-1]); tail = R(Nn) - head - cross
    print(f"   N = {Nn:8d}: R = {R(Nn):.6f} = head {head:.6f} + cross {cross:+.6f} + tail {tail:.6f};  tail/R = {tail/R(Nn):.4f}")
# ---------- diagnostics by FFT ----------
print("== diagnostics: N^2/J^2 x (D, O, total) for J < J0(N); weights w = 1 and w = log n - 1 ==")
def modes(cvec, N):
    """b_j, j = 1..N, from c (length N) by one FFT: b_j = (2/sqrt D) sum_n c_n sin(n (2j-1) pi / D) = (2/sqrt D) Im FFT_{4N+2}(c)[2j-1]."""
    D = 2*N + 1; x = np.zeros(2*D); x[1:N+1] = cvec
    f = np.fft.rfft(x); return (2/math.sqrt(D)) * (-np.imag(f[1:2*N+1:2]))   # e^{-2 pi i k n/(2D)} -> Im gives -sin
def diag_terms(c2, N, J):
    """D_{N,J} = (4/D) sum_{j<=J} sum_n c_n^2 sin^2(n theta_j) = (2/D) sum_{j<=J} [sum c^2 - sum c^2 cos(2 n theta_j)], by one FFT of c^2."""
    D = 2*N + 1; x = np.zeros(D); x[1:N+1] = c2; f = np.fft.rfft(x); cosT = np.real(f[1:2*J+1:2])      # frequencies (2j-1)/D
    return (2/D) * np.cumsum(np.sum(c2) - cosT)
res = {}
for Nn in (10**3, 10**4, 10**5, 10**6, 10**7):
    nn = np.arange(1, Nn + 1); c1 = mu[1:Nn+1] / nn; w = np.where(nn >= 2, np.log(nn) - 1, 0.0); cl = c1 * w
    ell = math.log(math.e*Nn); J0 = math.ceil(Nn / ell**2)
    b1 = modes(c1, Nn); bl = modes(cl, Nn); S1 = np.cumsum(b1**2); Sl = np.cumsum(bl**2)
    D1 = diag_terms(c1**2, Nn, J0); Dl = diag_terms(cl**2, Nn, J0)
    Js = [J for J in (1, 10, 100, 1000, 10000) if J < J0] + [J0 - 1]
    print(f"   N = {Nn:8d}, J0 = {J0:6d}:")
    for J in Js:
        f = Nn**2 / J**2
        print(f"      J = {J:6d}: w=1  D {f*D1[J-1]:9.4f}  O {f*(S1[J-1]-D1[J-1]):+10.4f}  total {f*S1[J-1]:9.4f}   |  w=log  D {f*Dl[J-1]:9.3f}  O {f*(Sl[J-1]-Dl[J-1]):+10.3f}  total {f*Sl[J-1]:9.3f}   (log^2(eN) = {ell**2:.1f})")
    res[Nn] = {"J0": J0, "J": Js, "D1": [float(D1[J-1]) for J in Js], "S1": [float(S1[J-1]) for J in Js], "Dl": [float(Dl[J-1]) for J in Js], "Sl": [float(Sl[J-1]) for J in Js]}
    if Nn in (10**4, 10**5):
        sq = (mu[1:Nn+1] != 0)
        rs = [modes(c1 * rng.choice([-1.0, 1.0], size=Nn), Nn) for _ in range(5)]; Sr = np.mean([np.cumsum(r**2) for r in rs], axis=0)
        bp = modes(np.abs(c1), Nn); Sp = np.cumsum(bp**2)
        print(f"      controls (w = 1), N^2/J^2 x total: " + "  ".join(f"J={J}: random {Nn**2/J**2*Sr[J-1]:.3f}, positive {Nn**2/J**2*Sp[J-1]:.3e}, Moebius {Nn**2/J**2*S1[J-1]:.3f}" for J in Js))
        res[Nn]["random"] = [float(Sr[J-1]) for J in Js]; res[Nn]["positive"] = [float(Sp[J-1]) for J in Js]
# ---------- dyadic Gram blocks (11) at N = 2000 ----------
print("== dyadic Gram blocks of F_N at N = 2000: G_rs = <Pi_J q_r, Pi_J q_s>, q_r = sum_{2^r <= d < 2^(r+1)} ((Lambda(d)-1)/d) T_d b_{N/d} ==")
N = 2000; primes, mu2 = sieve(N); nn = np.arange(1, N + 1); D = 2*N + 1; theta = (2*nn - 1)*math.pi/D; V = (2/math.sqrt(D))*np.sin(np.outer(nn, theta))
Lam = np.zeros(N + 1)
for p in primes:
    q = p
    while q <= N: Lam[q] = math.log(p); q *= p
c = mu2[1:N+1] / nn; blocks = []
for r in range(1, int(math.log2(N)) + 1):
    v = np.zeros(N)
    for d in range(2**r, min(2**(r+1), N + 1)):
        K = N // d; v[d*np.arange(1, K+1) - 1] += (Lam[d] - 1) / d * c[:K]
    blocks.append(V.T @ v)
Q = np.array(blocks); F = Q.sum(axis=0); Fdir = V.T @ (-np.where(nn >= 2, (np.log(nn) - 1) * c, 0.0))
print(f"   sum of blocks = F_N: error {np.max(np.abs(F - Fdir)):.1e};  blocks r = 1..{len(blocks)}")
gram = {}
for J in (10, 50):
    G = Q[:, :J] @ Q[:, :J].T; tot = float(np.sum(G)); gram[J] = G.tolist()
    print(f"   J = {J}: ||Pi_J F||^2 = {tot:.4e} = sum of Gram;  diagonal blocks G_rr: " + " ".join(f"{G[r,r]:.1e}" for r in range(len(blocks))))
    print(f"           row sums sum_s G_rs (block r's net share): " + " ".join(f"{np.sum(G[r]):+.1e}" for r in range(len(blocks))))
    off = G - np.diag(np.diag(G)); print(f"           sum of diagonal blocks {np.trace(G):.4e}, sum of off-diagonal {np.sum(off):+.4e};  largest off-diagonal pairs: " + ", ".join(f"({i+1},{j+1}) {off[i,j]:+.1e}" for i, j in zip(*np.unravel_index(np.argsort(np.abs(off).ravel())[-4:][::-1], off.shape)) if i < j or True)[:140])
json.dump({"diag": res, "gram": gram}, open("data/round_trip.json", "w"))
# ---------- figure ----------
BLUE, ORANGE, INK, INK2, SURF, GRID = "#2a78d6", "#eb6834", "#0b0b0b", "#52514e", "#fcfcfb", "#e6e5e1"
plt.rcParams.update({"font.size": 10, "axes.titlesize": 11, "axes.edgecolor": GRID, "axes.labelcolor": INK2, "xtick.color": INK2, "ytick.color": INK2,
                     "axes.spines.top": False, "axes.spines.right": False, "figure.facecolor": SURF, "axes.facecolor": SURF, "text.color": INK})
fig, ax = plt.subplots(1, 2, figsize=(14, 5.2))
for a in ax: a.grid(True, color=GRID, lw=0.6); a.set_axisbelow(True)
A = ax[0]
for Nn, col, lw in ((10**4, BLUE, 1.0), (10**6, ORANGE, 1.2)):
    nn = np.arange(1, Nn + 1); c1 = mu[1:Nn+1] / nn; ell = math.log(math.e*Nn); J0 = math.ceil(Nn / ell**2)
    b1 = modes(c1, Nn); S1 = np.cumsum(b1**2); Js = np.arange(1, J0); A.plot(Js, Nn**2/Js**2 * S1[:J0-1], color=col, lw=lw, label=f"Möbius, N = 10^{int(math.log10(Nn))}, J < J₀ = {J0}")
Nn = 10**4; nn = np.arange(1, Nn + 1); c1 = mu[1:Nn+1] / nn; J0 = math.ceil(Nn / math.log(math.e*Nn)**2); Js = np.arange(1, J0)
Sr = np.mean([np.cumsum(modes(c1 * rng.choice([-1.0, 1.0], size=Nn), Nn)**2) for _ in range(5)], axis=0)
A.plot(Js, Nn**2/Js**2 * Sr[:J0-1], color=INK2, lw=1.0, ls="--", label="random signs, N = 10⁴, mean of 5")
A.axhline(2*math.pi**2/4, color=INK, lw=0.8, ls=":", label="diagonal bound (8): π²/2")
A.set_xscale("log"); A.set_yscale("log"); A.set_xlabel("J"); A.set_ylabel("N²/J² · ‖Π_J b_N‖²")
A.set_title("A.  The low-mode population against its target scale J²/N²: Möbius and the random control", loc="left")
A.legend(frameon=False, fontsize=8.5, loc="lower right")
Bx = ax[1]; G = np.array(gram[50]); r = len(G)
im = Bx.imshow(np.sign(G) * np.log10(1 + np.abs(G) * 1e6), cmap="RdBu_r", vmin=-4, vmax=4)
Bx.set_xticks(range(r)); Bx.set_xticklabels([f"2^{k+1}" for k in range(r)], fontsize=8); Bx.set_yticks(range(r)); Bx.set_yticklabels([f"2^{k+1}" for k in range(r)], fontsize=8)
Bx.set_xlabel("dilation block s (d ∈ [2^s, 2^{s+1}))"); Bx.set_ylabel("dilation block r"); Bx.grid(False)
Bx.set_title("B.  Gram blocks of F_N at N = 2000, J = 50: sign·log₁₀(1 + 10⁶|G_rs|)", loc="left")
cb = fig.colorbar(im, ax=Bx, fraction=0.046, pad=0.03); cb.ax.tick_params(labelsize=8)
fig.suptitle("The round trip: the signed correlations that remain after the diagonal is bounded", fontsize=11)
fig.tight_layout(); fig.savefig(out, dpi=160, bbox_inches="tight"); print("wrote", out)
