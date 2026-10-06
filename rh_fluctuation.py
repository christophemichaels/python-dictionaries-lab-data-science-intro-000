"""The coupled prime-fluctuation equation (SUBPOWER_CLOSURE.md, Section 9). With L_N = V^T diag(log n) V, u_N = V^T e_1 and
F_N = sum_{d=2}^N ((Lambda(d)-1)/d) T_d b_{N/d}, the exact equation (L_N - I) b_N = -u_N - F_N (the d = 1 term of the fluctuation
identity separated). Computes: the completion formula M_{Q_P}(N) = M(N) + sum_{P<p<=N} M(N/p) against the naive one; the equation's
error; the inversion norm ||(L_N - I)^{-1}|| in the energy norm (as the 2-norm of C^T diag(1/(log n - 1)) C^{-T}); the chain
sqrt R <= norm (1 + ||F||); and the exact collapse of F_N: its coefficient vector is f_n = -(log n - 1) mu(n)/n for n >= 2, f_1 = 0,
so ||F_N||^2 = R^log(N) = sum_{a,b>=2} mu(a)mu(b)(log a - 1)(log b - 1)/max(a,b), computed to 1e7 by partial sums.
Usage: python3 rh_fluctuation.py"""
import math, numpy as np
def sieve(N):
    s = np.ones(N + 1, dtype=bool); s[:2] = False
    for p in range(2, int(N**0.5) + 1):
        if s[p]: s[p*p::p] = False
    primes = np.nonzero(s)[0]
    mu = np.ones(N + 1, dtype=np.int8)
    for p in primes: mu[p::p] *= -1; mu[p*p::p*p] = 0
    mu[0] = 0
    Lam = np.zeros(N + 1)
    for p in primes:
        q = p
        while q <= N: Lam[q] = math.log(p); q *= p
    Pplus = np.zeros(N + 1, dtype=np.int64)
    for q in primes: Pplus[q::q] = q
    return primes, mu, Lam, Pplus
# completion formula at N = 4000, P = 1000
N = 4000; primes, mu, Lam, Pplus = sieve(N); M = np.cumsum(mu.astype(np.int64))
P = 1000; MQ = int(np.sum(mu[1:N+1][Pplus[1:N+1] <= P])); pi = lambda x: int(np.sum(primes <= x))
print(f"N = {N}, P = {P}: M_Q(N) = {MQ};  M(N) + sum_(P<p<=N) M(N/p) = {M[N] + sum(int(M[N//p]) for p in primes[primes > P])};  M(N) + pi(N) - pi(P) = {M[N] + pi(N) - pi(P)}")
for N in (1000, 2000, 4000):
    primes, mu, Lam, Pplus = sieve(N); M = np.cumsum(mu.astype(np.int64))
    nn = np.arange(1, N + 1); jj = np.arange(1, N + 1); theta = (2*jj - 1) * math.pi / (2*N + 1); lam = 1 / (4 * np.sin(theta/2)**2)
    V = (2 / math.sqrt(2*N + 1)) * np.sin(np.outer(nn, theta)); c = mu[1:N+1] / nn; b = V.T @ c
    logn = np.log(nn); e1 = np.zeros(N); e1[0] = 1.0; u = V.T @ e1
    fcoef = np.zeros(N)
    for dd in range(2, N + 1):
        K = N // dd; fcoef[dd*np.arange(1, K+1) - 1] += (Lam[dd] - 1) / dd * c[:K]
    F = V.T @ fcoef
    lhs = V.T @ ((logn - 1) * c); err_eq = np.max(np.abs(lhs + u + F))
    f_pred = np.where(nn >= 2, -(logn - 1) * c, 0.0); err_collapse = np.max(np.abs(fcoef - f_pred))
    g = 1 / (logn - 1)
    # operator norm of (L - I)^{-1} in the energy norm: ||C^T diag(g) C^{-T}||_2, C^T upper-triangular ones, (C^T)^{-1} = I - superdiagonal
    CTinv = np.eye(N) - np.eye(N, k=1); A = np.cumsum((g[:, None] * CTinv)[::-1], axis=0)[::-1]      # C^T (g * C^{-T}): column sums from row n downward
    nrm = float(np.linalg.norm(A, 2))
    Rn = float(np.sum(lam * b**2)); nF = float(np.sum(lam * F**2)); nu = float(np.sum(lam * u**2))
    print(f"N = {N}: equation (L-I)b = -u - F: error {err_eq:.1e};  collapse f_n = -(log n - 1) mu(n)/n (n>=2): error {err_collapse:.1e};  ||u||^2 = {nu:.6f}")
    print(f"        ||(L-I)^(-1)||_(energy) = {nrm:.4f} (claimed <= 216; multiplier values -1, {g[1]:.3f}, {g[2]:.3f}, {g[3]:.3f}, ...);  sqrt R(N) = {math.sqrt(Rn):.5f};  ||F_N|| = {math.sqrt(nF):.5f};  chain sqrt R <= norm (1 + ||F||) = {nrm*(1+math.sqrt(nF)):.3f}")
    print(f"        ||F_N||^2 = {nF:.6f};  R(N) = {Rn:.6f};  R(N) log^2 N = {Rn*math.log(N)**2:.4f};  ratio ||F||^2 / R = {nF/Rn:.3f}")
# the log-weighted energy to 1e7
X = 10**7; primes, mu, Lam, Pplus = sieve(X); n = np.arange(X + 1, dtype=float)
g = mu.astype(float) * (np.log(np.maximum(n, 1)) - 1); g[0] = 0; g[1] = 0
G = np.cumsum(g); M = np.cumsum(mu.astype(np.int64))
def Rlog(Nn): return float(np.sum(G[2:Nn].astype(float)**2 * (1/n[2:Nn] - 1/n[3:Nn+1])) + G[Nn]**2 / Nn)
def R(Nn): return float(np.sum(M[1:Nn].astype(float)**2 * (1/n[1:Nn] - 1/n[2:Nn+1])) + M[Nn]**2 / Nn)
print("the collapsed coupled energy R^log(N) = sum_{a,b>=2} mu(a)mu(b)(log a - 1)(log b - 1)/max(a,b), beside R(N):")
for Nn in (10**3, 10**4, 10**5, 10**6, 10**7):
    print(f"   N = {Nn:8d}: R^log = {Rlog(Nn):9.4f}   R = {R(Nn):.4f}   R^log / (R log^2 N) = {Rlog(Nn)/(R(Nn)*math.log(Nn)**2):.4f}   G(N) = {G[Nn]:.2f}   G(N)^2/N = {G[Nn]**2/Nn:.4f}")
