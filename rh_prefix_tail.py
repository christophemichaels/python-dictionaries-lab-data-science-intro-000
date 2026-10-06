"""The prefix-tail mechanism (HORIZON_ROUND_TRIP.md, Section 6). Reproduces the table of Theorem 2.1 of
received/Michaels_Arithmophysics_Growing_Mode_Bounds.pdf; tests the sharper lowest-mode estimate for the completed prefix,
   lambda_J ||Pi_J h||^2 <= (8 pi^2/3) (J/N) M*(L)^2 (1 + log+(pi J L / N))^2,   M*(L) = max_{x<=L} |M(x)|,
against the actual value and against the prefix energy R(L) used in the paper; and quantifies the Cauchy-Schwarz loss on
the prefix-tail interaction Gamma.  Usage: python3 rh_prefix_tail.py"""
import math, numpy as np
def sieve(N):
    s = np.ones(N + 1, dtype=bool); s[:2] = False
    for p in range(2, int(N**0.5) + 1):
        if s[p]: s[p*p::p] = False
    mu = np.ones(N + 1, dtype=np.int8)
    for p in np.nonzero(s)[0]: mu[p::p] *= -1; mu[p*p::p*p] = 0
    mu[0] = 0; return mu
def setup(N):
    mu = sieve(N); nn = np.arange(1, N + 1); D = 2*N + 1; th = (2*nn - 1)*math.pi/D; lam = 1/(4*np.sin(th/2)**2)
    V = (2/math.sqrt(D))*np.sin(np.outer(nn, th)); c = mu[1:N+1]/nn; M = np.cumsum(mu.astype(np.int64))
    return mu, nn, D, th, lam, V, c, M
def R_of(M, L):
    if L == 0: return 0.0
    k = np.arange(1, L, dtype=float); return float(np.sum(M[1:L].astype(float)**2/(k*(k+1))) + M[L]**2/L)
print("== Theorem 2.1 table (paper: 0.169033/0.749076/2.343828; 0.093853/1.193002/4.288237; 0.090570/1.183036/4.595958; 0.035052/1.396416/9.851200) ==")
for N, J, L in ((64, 8, 16), (257, 8, 65), (1024, 32, 64), (4096, 32, 256)):
    mu, nn, D, th, lam, V, c, M = setup(N); b = V.T @ c
    S = float(np.sum(b[:J]**2)); E = float(np.sum(lam[J-1:]*b[J-1:]**2)); T = float(np.sum(c[L:]**2)); allow = (math.sqrt(R_of(M, L)) + math.sqrt(lam[J-1]*T))**2
    print(f"   N = {N:4d}, J = {J:2d}, L = {L:3d}: lambda_J S = {lam[J-1]*S:.6f}   E>=J = {E:.6f}   allowance {allow:.6f}")
print("== the prefix's lowest modes: actual lambda_J ||Pi_J h||^2, the new bound (8 pi^2/3)(J/N) M*(L)^2 (1 + log+(pi J L/N))^2, and R(L) ==")
for N in (1024, 4096, 16384):
    mu, nn, D, th, lam, V, c, M = setup(N)
    Mstar = np.maximum.accumulate(np.abs(M))
    for J in (4, 32, 256):
        for L in (N//(2*J), N//J, 4*N//J, min(N, 16*N//J)):
            if L < 2: continue
            h = V.T @ np.where(nn <= L, c, 0.0); t = V.T @ np.where(nn > L, c, 0.0)
            P = lam[J-1]*float(np.sum(h[:J]**2)); newb = (8*math.pi**2/3)*(J/N)*Mstar[L]**2*(1 + max(0.0, math.log(math.pi*J*L/N)))**2; RL = R_of(M, L)
            G = float(np.dot(h[:J], t[:J])); cs = math.sqrt(float(np.sum(h[:J]**2))*float(np.sum(t[:J]**2)))
            Tt = lam[J-1]*float(np.sum(t[:J]**2)); Tpar = lam[J-1]*float(np.sum(c[L:]**2))
            print(f"   N = {N:5d} J = {J:3d} L = {L:5d} (JL/N = {J*L/N:5.1f}): actual {P:9.5f}   new bound {newb:10.4f}   R(L) {RL:8.4f}  | tail: actual {Tt:9.5f} Parseval {Tpar:9.4f}  | Gamma {lam[J-1]*G:+9.5f} vs CS {lam[J-1]*cs:9.5f} (ratio {G/cs if cs > 0 else 0:+.3f})")
print("== how far the truth sits below the envelope: max_J lambda_J S_N(J) (J/N) at N = 4096, 16384 ==")
for N in (4096, 16384):
    mu, nn, D, th, lam, V, c, M = setup(N); b = V.T @ c; S = np.cumsum(b**2); J = np.arange(1, N + 1); q = lam*S*J/N
    print(f"   N = {N}: max over J of lambda_J S(J) J/N = {q.max():.4f} at J = {int(np.argmax(q))+1};  at J = N^(1/2) = {int(math.sqrt(N))}: {q[int(math.sqrt(N))-1]:.4f};  at J = N^(3/4) = {int(N**0.75)}: {q[int(N**0.75)-1]:.4f};  lambda_J S(J) at J = N^(1/2): {lam[int(math.sqrt(N))-1]*S[int(math.sqrt(N))-1]:.4f}")
print("== the completed prefix's lowest modes are M(L): h_j = (2 theta_j / sqrt D) (M(L) + e_j), |e_j| <= theta_j^2 L^2 M*(L)/3 in the monotone regime ==")
for N, J, L in ((4096, 32, 64), (4096, 4, 512), (16384, 256, 32), (16384, 32, 256)):
    mu, nn, D, th, lam, V, c, M = setup(N); Mstar = np.maximum.accumulate(np.abs(M))
    h = V.T @ np.where(nn <= L, c, 0.0)
    e = h[:J]*math.sqrt(D)/(2*th[:J]) - M[L]; ebound = th[:J]**2 * L**2 * Mstar[L]/3
    sharp = lam[J-1]*(4/D)*Mstar[L]**2*float(np.sum(th[:J]**2)); actual = lam[J-1]*float(np.sum(h[:J]**2)); lead = lam[J-1]*(4/D)*M[L]**2*float(np.sum(th[:J]**2))
    print(f"   N = {N:5d} J = {J:3d} L = {L:4d} (theta_J L = {th[J-1]*L:.3f}): M(L) = {M[L]:4d}, M*(L) = {Mstar[L]:3d};  max|e_j| = {np.max(np.abs(e)):.4f} <= {np.max(ebound):.4f};  lambda_J ||Pi_J h||^2: actual {actual:.5f}, leading (4/D) M(L)^2 sum theta^2 {lead:.5f}, monotone bound with M* {sharp:.5f}, R(L) {R_of(M, L):.4f}")
