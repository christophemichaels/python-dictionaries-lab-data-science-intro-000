"""The plane as a gravitational field (THE_GRAVITY_PLANE.md). Masses on concentric shells of radius n: the Newtonian potential of a
shell is 1/max(r, R) (the shell theorem), so the Green energy of the record, sum m_i m_j / max(i, j), is the gravitational self-energy
of the shells, and its field-energy form is int M(t)^2 / t^2 dt + M(N)^2 / N with M(t) the enclosed mass. Two stars: the prime star,
masses Lambda(n); the Moebius star, masses mu(n). Usage: python3 rh_gravity.py [X=1e7]"""
import sys, math
import numpy as np
X = int(float(sys.argv[1])) if len(sys.argv) > 1 else 10**7
EG = 0.5772156649015329
s = np.ones(X + 1, dtype=bool); s[:2] = False
for p in range(2, int(X**0.5) + 1):
    if s[p]: s[p*p::p] = False
primes = np.nonzero(s)[0]
Lam = np.zeros(X + 1)
for p in primes:
    q = p
    while q <= X: Lam[q] = math.log(p); q *= p
mu = np.ones(X + 1, dtype=np.int8)
for p in primes: mu[p::p] *= -1; mu[p*p::p*p] = 0
mu[0] = 0
psi = np.cumsum(Lam); M = np.cumsum(mu.astype(np.int64)); t = np.arange(X + 1, dtype=float)
# 1. the shell theorem identity, checked directly at N = 2000: sum m_i m_j / max(i,j) = int_1^N M(t)^2 dt / t^2 + M(N)^2 / N
N0 = 2000; m = mu[1:N0+1].astype(float); idx = np.arange(1, N0+1, dtype=float)
G = 1.0 / np.maximum.outer(idx, idx); direct = m @ G @ m
field = sum(M[k]**2 * (1/k - 1/(k+1)) for k in range(1, N0)) + M[N0]**2 / N0
print(f"shell theorem, Moebius star, N = {N0}: double sum of masses over 1/max(i,j) = {direct:.6f}; field energy int M^2/t^2 + M(N)^2/N = {field:.6f}")
# 2. the Moebius star: field energy I_M(N), its diagonal (self-energy of the shells) and its interaction (the pull between shells)
Ns = [10**3, 10**4, 10**5, 10**6, X]
w = 1/t[1:-1] - 1/t[2:]          # weight of the interval [k, k+1)
IM = np.concatenate([[0.0], np.cumsum(M[1:-1].astype(float)**2 * w)])   # I_M(k) for k = 1..X-1  (index k-1)
self_energy = np.cumsum(mu[1:].astype(float)**2 / t[1:])                # sum_{n<=k} mu(n)^2 / n
print("Moebius star (masses mu(n)):")
for N in Ns:
    R = IM[N-1] + M[N]**2 / N; diag = self_energy[N-1]; inter = R - diag
    print(f"  N = {N:9d}: field energy R(N) = {R:8.4f} = {R/math.log(N):.4f} log N;  self-energy of the shells {diag:8.4f} = {diag/math.log(N):.4f} log N (6/pi^2 = {6/math.pi**2:.4f});  interaction {inter:+8.4f} = {inter/math.log(N):+.4f} log N;  M(N) = {M[N]}")
# 3. the prime star: masses Lambda(n); enclosed mass psi(t) ~ t (a uniform star); fluctuation field (psi - t)/t^2 and its energy
print("prime star (masses Lambda(n)):")
fluct = np.concatenate([[0.0], np.cumsum((psi[1:-1] - t[1:-1])**2 * w)])
for N in Ns:
    E = np.sum(psi[1:N]**2 * w[:N-1]) + psi[N]**2 / N
    print(f"  N = {N:9d}: field energy {E:12.1f} = {E/N:.4f} N (a uniform star has energy ~ N);  fluctuation energy int (psi - t)^2/t^2 = {fluct[N-1]:.4f} = {fluct[N-1]/math.log(N):.5f} log N;  Cramer's constant under RH: sum 1/|rho|^2 = 2 + gamma - log 4 pi = {2 + EG - math.log(4*math.pi):.5f}")
