"""Shared-checkpoint checks (HORIZON_ROUND_TRIP.md, Section 5): the table ||F_N||_E^2 at N = 64, 256, 1024 against R^log(N);
the first-mode residual-to-diagonal ratio S_N(1)/D_{N,1} against M(N)^2/N and against its exact form (pi^2/(2N^3)) Mtilde(N)^2 / D_{N,1}
with Mtilde(N) = sum mu(n) sinc(n theta_1); and the max of |M(x)|/sqrt x to 1e7.  Usage: python3 rh_checkpoint_check.py"""
import math, numpy as np
X = 10**7
s = np.ones(X + 1, dtype=bool); s[:2] = False
for p in range(2, int(X**0.5) + 1):
    if s[p]: s[p*p::p] = False
mu = np.ones(X + 1, dtype=np.int8)
for p in np.nonzero(s)[0]: mu[p::p] *= -1; mu[p*p::p*p] = 0
mu[0] = 0
n = np.arange(X + 1, dtype=float); M = np.cumsum(mu.astype(np.int64))
g = mu.astype(float) * (np.log(np.maximum(n, 1)) - 1); g[:2] = 0; G = np.cumsum(g)
def Rlog(N): return float(np.sum(G[2:N].astype(float)**2 * (1/n[2:N] - 1/n[3:N+1])) + G[N]**2 / N)
print("checkpoint table ||F_N||_E^2 (reported 0.430376710, 1.078227825, 2.161216072) against R^log(N):")
for N in (64, 256, 1024): print(f"   N = {N:5d}: R^log(N) = {Rlog(N):.9f}")
print("first mode: S_N(1)/D_{N,1}, M(N)^2/N, and the exact (pi^2/(2 N^3)) Mtilde^2 / D form")
for N in (10**3, 10**4, 10**5, 10**6, 10**7):
    nn = np.arange(1, N + 1); D = 2*N + 1; th = math.pi / D; c = mu[1:N+1] / nn
    b1 = (2/math.sqrt(D)) * float(np.sum(c * np.sin(nn*th))); D1 = (4/D) * float(np.sum(c**2 * np.sin(nn*th)**2))
    Mt = float(np.sum(mu[1:N+1] * np.sin(nn*th)/(nn*th)))
    print(f"   N = {N:8d}: S/D = {b1**2/D1:.4f};  M(N)^2/N = {M[N]**2/N:.4f};  Mtilde(N)^2/N = {Mt**2/N:.4f};  (4 theta^2/D) Mtilde^2 / D1 = {(4*th**2/D)*Mt**2/D1:.4f};  M(N) = {M[N]}, Mtilde = {Mt:.1f};  D1 N^2 = {D1*N**2:.4f}")
q = np.abs(M[1:]) / np.sqrt(n[1:]); i = int(np.argmax(q)) + 1
print(f"max |M(x)|/sqrt x for x <= 1e7: {q.max():.4f} at x = {i};  the sign O_(N,1) <= 0 is |Mtilde(N)| <= {math.sqrt(2.3211/(math.pi**2/2)):.3f} sqrt N (from D_1 N^2 = 2.3211)")
