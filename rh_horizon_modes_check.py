"""Checks of the Exact Horizon Map (received/Michaels_Light_Plane_Horizon_Map.pdf) at N beyond its own script's N = 64:
the sine modes of min(a,b) with gains 1/(4 sin^2(theta_j/2)), the mode resolution R(N) = sum lambda_j |b_j|^2, the mode-population
allowance under (G), the local sample inversion mu(n) = n(2 w_n - w_{n-1} - w_{n+1}), the phase readout R_N(tau) and its transport
bound with q(tau) = sqrt(1 + tau^2) + |tau|, and the long phase average = sum mu^2/n.  Usage: python3 rh_horizon_modes_check.py"""
import math, numpy as np
X = 10**4
s = np.ones(X + 1, dtype=bool); s[:2] = False
for p in range(2, int(X**0.5) + 1):
    if s[p]: s[p*p::p] = False
mu = np.ones(X + 1, dtype=np.int8)
for p in np.nonzero(s)[0]: mu[p::p] *= -1; mu[p*p::p*p] = 0
mu[0] = 0
M = np.cumsum(mu.astype(np.int64)); n = np.arange(X + 1, dtype=float)
def R(N): return float(np.sum(M[1:N].astype(float)**2 * (1/n[1:N] - 1/n[2:N+1])) + M[N]**2 / N)
for N in (10**3, 10**4):
    nn = np.arange(1, N + 1); j = np.arange(1, N + 1); theta = (2*j - 1) * math.pi / (2*N + 1); lam = 1 / (4 * np.sin(theta/2)**2)
    V = (2 / math.sqrt(2*N + 1)) * np.sin(np.outer(nn, theta))                  # columns v_j(n)
    c = mu[1:N+1] / nn
    b = V.T @ c
    orth = np.max(np.abs(V[:, :50].T @ V[:, :50] - np.eye(50)))
    Bc = np.minimum.outer(nn, nn) @ c if N <= 10**3 else None
    print(f"N = {N}: sum lambda_j |b_j|^2 = {float(np.sum(lam * b**2)):.6f} against R(N) = {R(N):.6f};  orthonormality of the first 50 modes to {orth:.1e}"
          + (f";  |B v_1 - lambda_1 v_1| = {np.max(np.abs(np.minimum.outer(nn, nn) @ V[:, 0] - lam[0] * V[:, 0])):.1e}" if N <= 10**3 else ""))
    for J in (1, 10, 100):
        pop = float(np.sum(b[:J]**2)); allow = 4 * R(N) * math.sin((2*J - 1) * math.pi / (4*N + 2))**2
        print(f"   modes j <= {J:3d}: population sum |b_j|^2 = {pop:.3e};  allowance 4 R(N) sin^2((2J-1)pi/(4N+2)) = {allow:.3e};  ratio {pop/allow:.3f}")
    # local sample inversion from the field at the lattice
    w = np.concatenate([[0.0], np.minimum.outer(nn, nn) @ c]) if N <= 10**3 else None
    if w is not None:
        rec = np.array([k * (2*w[k] - w[k-1] - w[k+1]) for k in range(1, N)] + [N * (w[N] - w[N-1])])
        print(f"   sample inversion mu(n) = n(2 w_n - w_(n-1) - w_(n+1)): max error {np.max(np.abs(rec - mu[1:N+1])):.1e}")
# phase readout and the transport bound at N = 10^3
N = 10**3; nn = np.arange(1, N + 1); idx = np.nonzero(mu[1:N+1])[0] + 1; w = mu[idx] / np.sqrt(idx); ln = np.log(idx)
def RN(tau):
    z = np.sum(w * np.exp(-1j * tau * ln) * np.sqrt(idx))        # not used
    a = w * np.exp(-1j * tau * ln)
    G = np.sqrt(np.minimum.outer(idx, idx) / np.maximum.outer(idx, idx))
    return float(np.real(np.conj(a) @ G @ a))
R0 = R(N); deliv = float(np.sum(mu[1:N+1].astype(float)**2 / nn))
print(f"phase readout at N = {N}: R_N(0) = {RN(0.0):.6f} = R(N) = {R0:.6f}")
for tau in (0.5, 1.0, 3.0, 14.134725):
    q = math.sqrt(1 + tau**2) + abs(tau); v = RN(tau)
    print(f"   tau = {tau:9.6f}: R_N(tau) = {v:9.4f};  bounds [{R0/q**2:.4f}, {R0*q**2:.1f}];  inside: {R0/q**2 <= v <= R0*q**2}")
taus = np.linspace(0, 2000, 4001); avg = np.mean([RN(t) for t in taus[::8]])
print(f"   phase average of R_N(tau) over 0..2000 (sampled): {avg:.4f} against the delivered sum mu^2/n = {deliv:.4f}")
