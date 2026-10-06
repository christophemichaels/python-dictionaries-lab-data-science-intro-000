"""The boundary carry (BOUNDARY_CARRY.md). Verifies the boundary replacement of received/Michaels_Arithmophysics_Boundary_Carry_
Continuation_Prompt.pdf (q = c_{<=L} - (M(L)/L) e_L, w = c_{>L} + (M(L)/L) e_L: q^T B q = I_L, q^T B w = 0, r_j = sum_{k<L} M(k) D_j(k),
the table R(L)/lambda_J and I_L Xi, the bound (B) at L = floor((N/J)^(5/6))); carries out the first calculation (the retained
population W split exactly into boundary square, tail square and mixed term; the max-|M| comparison and its loss) and the second
(interval identities (I); the partition with lengths (N/J)^(5/3)/a; per-interval errors, their sum, the exact total and its
coherence; the Gram form (P) of the increments; coherence diagnostics for Moebius, random signs and all-positive).
Usage: python3 rh_boundary_carry.py"""
import math, numpy as np
rng = np.random.default_rng(3)
def sieve(N):
    s = np.ones(N + 1, dtype=bool); s[:2] = False
    for p in range(2, int(N**0.5) + 1):
        if s[p]: s[p*p::p] = False
    mu = np.ones(N + 1, dtype=np.int8)
    for p in np.nonzero(s)[0]: mu[p::p] *= -1; mu[p*p::p*p] = 0
    mu[0] = 0; return mu
MU = sieve(10**6)
def modes_of(coef, N, J):
    """(V_N^T coef)_j for j <= J, coef indexed 1..N (array length N)."""
    nn = np.arange(1, N + 1); D = 2*N + 1; th = (2*np.arange(1, J + 1) - 1)*math.pi/D
    out = np.empty(J)
    for a in range(0, N, 100000):
        blk = (2/math.sqrt(D))*np.sin(np.outer(th, nn[a:a+100000]))
        out = out + blk @ coef[a:a+100000] if a else blk @ coef[a:a+100000]
    return out
def lam_of(N, j): return 1/(4*math.sin((2*j - 1)*math.pi/(2*(2*N + 1)))**2)
def Dj(N, j, k):
    """D_j(k) = v_j(k)/k - v_j(k+1)/(k+1), vectorised in k."""
    D = 2*N + 1; th = (2*j - 1)*math.pi/D; return (2/math.sqrt(D))*(np.sin(k*th)/k - np.sin((k+1)*th)/(k+1))
print("== 1. the boundary replacement: identities and the supplied table ==")
for N, J, L in ((257, 1, 101), (1024, 1, 322), (4096, 1, 1024), (4096, 4, 322)):
    mu = MU[:N+1]; nn = np.arange(1, N + 1); c = mu[1:N+1]/nn; M = np.cumsum(mu.astype(np.int64))
    k = np.arange(1, L, dtype=float); IL = float(np.sum(M[1:L]**2/(k*(k+1)))); RL = IL + M[L]**2/L
    q = np.where(nn <= L, c, 0.0); q[L-1] -= M[L]/L; w = c - q
    Bq = np.minimum.outer(nn, nn) @ q if N <= 4096 else None
    Xi = sum(float(np.sum(k*(k+1)*Dj(N, j, k)**2)) for j in range(1, J + 1))
    r = modes_of(q, N, J); rabel = np.array([float(np.sum(M[1:L]*Dj(N, j, k))) for j in range(1, J + 1)])
    print(f"   (N,J,L) = ({N},{J},{L}): R(L)/lambda_J = {RL/lam_of(N, J):.6e}  I_L Xi = {IL*Xi:.6e};  q^T B q = {float(q @ Bq):.9f} vs I_L = {IL:.9f};  q^T B w = {float(w @ Bq):.1e};  "
          f"r_j Abel error {np.max(np.abs(r - rabel)):.1e};  ||Pi_J r||^2 = {float(np.sum(r**2)):.4e} <= I_L min(Xi, 1/lambda_J) = {IL*min(Xi, 1/lam_of(N, J)):.4e}")
    if (N, J, L) == (4096, 1, 1024):
        h = modes_of(np.where(nn <= L, c, 0.0), N, 1); t = modes_of(np.where(nn > L, c, 0.0), N, 1); e = np.zeros(N); e[L-1] = M[L]/L; bnd = modes_of(e, N, 1)
        print(f"      inner products at J = 1: <h, t> = {float(h[0]*t[0]):.6e} (reported 2.509743e-9);  boundary part <(M(L)/L) e_L, t> = {float(bnd[0]*t[0]):.6e} (reported 2.491264e-9)")
print("== 2. the bound (B) and the first calculation: W = boundary^2 + tail^2 + mixed, at L = floor((N/J)^(5/6)) ==")
Cd = 2*math.pi**6/315
for N in (4096, 16384, 65536):
    mu = MU[:N+1]; nn = np.arange(1, N + 1); c = mu[1:N+1]/nn; M = np.cumsum(mu.astype(np.int64)); D = 2*N + 1
    for J in (1, 4, 32):
        L = int((N/J)**(5/6)); th = (2*np.arange(1, J + 1) - 1)*math.pi/D
        q = np.where(nn <= L, c, 0.0); q[L-1] -= M[L]/L; r = modes_of(q, N, J)
        bnd = (2/math.sqrt(D))*(M[L]/L)*np.sin(L*th); tail = modes_of(np.where(nn > L, c, 0.0), N, J); b = modes_of(c, N, J)
        W = float(np.sum((bnd + tail)**2)); Bs = float(np.sum(bnd**2)); Ts = float(np.sum(tail**2)); mix = 2*float(np.dot(bnd, tail)); S = float(np.sum(b**2))
        Mstar = float(np.max(np.abs(M[L:N+1]))); kk = np.arange(L, N, dtype=float)
        maxcmp = sum((Mstar*float(np.sum(np.abs(Dj(N, j, kk)))) + abs(M[N])/N*abs(math.sin(N*th[j-1]))*(2/math.sqrt(D)))**2 for j in range(1, J + 1))
        print(f"   N = {N:5d} J = {J:2d} L = {L:5d}: ||Pi_J r||^2 = {float(np.sum(r**2)):.3e} <= C_d J^2/N^2 = {Cd*J**2/N**2:.3e} (ratio {float(np.sum(r**2))/(Cd*J**2/N**2):.4f});  "
              f"W = {W:.4e} = bnd {Bs:.3e} + tail {Ts:.3e} + mixed {mix:+.3e};  S = {S:.4e};  |sqrt S - sqrt W| = {abs(math.sqrt(S)-math.sqrt(W)):.2e} <= {math.sqrt(Cd)*J/N:.2e};  "
              f"W N^2/J^2 = {W*N**2/J**2:.4f};  max-|M| comparison bound {maxcmp:.3e} (loss x{maxcmp/W:.0f}; M*_[L,N] = {Mstar:.0f})")
print("== 3. interval identities (I) and B-orthogonality at N = 1024, dyadic partition ==")
N = 1024; mu = MU[:N+1]; nn = np.arange(1, N + 1); c = mu[1:N+1]/nn; M = np.concatenate([[0], np.cumsum(mu[1:].astype(np.int64))]); B = np.minimum.outer(nn, nn).astype(float)
xs = [0] + [2**r for r in range(0, 11)]; qs = []; err = 0.0; errI = 0.0
for a, b in zip(xs[:-1], xs[1:]):
    cr = np.where((nn > a) & (nn <= b), c, 0.0); Q = M[b] - M[a]; qr = cr.copy(); qr[b-1] -= Q/b; qs.append(qr)
    k = np.arange(a + 1, b); modesI = np.array([float(np.sum((M[k] - M[a])*Dj(N, j, k.astype(float)))) if len(k) else 0.0 for j in range(1, 4)])
    err = max(err, float(np.max(np.abs(modes_of(qr, N, 3) - modesI)))); Ir = float(np.sum((M[k]-M[a])**2/(k*(k+1.0)))) if len(k) else 0.0; errI = max(errI, abs(float(qr @ B @ qr) - Ir))
G = np.array([[float(qi @ B @ qj) for qj in qs] for qi in qs]); off = np.max(np.abs(G - np.diag(np.diag(G))))
print(f"   (V^T q^(r))_j vs sum (M(k)-M(a)) D_j(k): max error {err:.1e};  I_r formula error {errI:.1e};  B-orthogonality of disjoint remainders: max off-diagonal {off:.1e};  sum of remainders + carried charges = c: {np.max(np.abs(sum(qs) + sum(np.where(nn == b, (M[b]-M[a])/b, 0.0) for a, b in zip(xs[:-1], xs[1:])) - c)):.1e}")
print("== 4. the partition above L with lengths (N/J)^(5/3)/a: errors, coherence, the Gram form (P), controls ==")
def partition(N, J, L):
    xs = [L]; a = L
    while a < N:
        ell = max(1, int((N/J)**(5/3)/a)); a = min(N, a + ell); xs.append(a)
    return xs
for N, J in ((65536, 1), (65536, 4), (10**6, 1), (10**6, 4)):
    mu = MU[:N+1]; nn = np.arange(1, N + 1); M = np.concatenate([[0], np.cumsum(mu[1:].astype(np.int64))]); D = 2*N + 1; L = int((N/J)**(5/6)); xs = partition(N, J, L); m = len(xs) - 1
    th = (2*np.arange(1, J + 1) - 1)*math.pi/D
    sources = {"Moebius": mu[1:N+1].astype(float), "random": (mu[1:N+1] != 0)*rng.choice([-1.0, 1.0], size=N), "positive": (mu[1:N+1] != 0).astype(float)}
    line = f"   N = {N:7d} J = {J}: L = {L}, intervals m = {m} (lengths {xs[1]-xs[0]} ... {xs[-1]-xs[-2]});"
    for name, s in sources.items():
        Ms = np.concatenate([[0], np.cumsum(s)]); Q = np.array([Ms[b] - Ms[a] for a, b in zip(xs[:-1], xs[1:])]); x = np.array(xs[1:], dtype=float)
        K = (4/D)*np.sin(np.outer(x, th)) @ np.sin(np.outer(x, th)).T                         # K_{N,J}(x_r, x_s)
        WP = float((Q/x) @ K @ (Q/x)); WP0 = (4/D)*float(np.sum(np.sin(L*th)**2))*(Ms[L]/L)**2
        bnd_and_tail = (2/math.sqrt(D))*((Ms[L]/L)*np.sin(L*th) + (Q/x) @ np.sin(np.outer(x, th)))
        coh = Q.sum()**2/float(np.sum(Q**2))
        if N <= 65536:
            cvec = s/nn; tailc = np.where(nn > L, cvec, 0.0); Wtrue = float(np.sum(((2/math.sqrt(D))*(Ms[L]/L)*np.sin(L*th) + modes_of(tailc, N, J))**2))
            # per-interval remainders and their projected errors
            errs = []; tot = np.zeros(J); mom = []
            for a, b in zip(xs[:-1], xs[1:]):
                qr = np.where((nn > a) & (nn <= b), cvec, 0.0); qr[b-1] -= (Ms[b]-Ms[a])/b; rr = modes_of(qr, N, J); errs.append(float(np.sum(rr**2))); tot += rr
                k = np.arange(a + 1, b); mom.append(float(np.sum(k*(Ms[k] - Ms[a]))) if len(k) else 0.0)
            mom = np.array(mom); line += f"\n      {name:8s}: W (bnd+tail, exact) = {Wtrue:.4e};  W^P (increments at right endpoints) = {float(np.sum(bnd_and_tail**2)):.4e};  sum of per-interval errors {sum(errs):.3e};  exact total error {float(np.sum(tot**2)):.3e} (C_d J^2/N^2 = {Cd*J**2/N**2:.3e});  first-moment coherence (sum m_r)^2/sum m_r^2 = {mom.sum()**2/float(np.sum(mom**2)):.2f};  increment coherence (sum Q)^2/sum Q^2 = {coh:.2f}"
        else:
            line += f"\n      {name:8s}: W^P (increments at right endpoints, with the boundary) = {float(np.sum(bnd_and_tail**2)):.4e} (N^2/J^2 x = {float(np.sum(bnd_and_tail**2))*N**2/J**2:.4f});  increment coherence (sum Q)^2/sum Q^2 = {coh:.2f};  sum Q^2 = {float(np.sum(Q**2)):.0f}, (N - L) 6/pi^2 = {(N-L)*6/math.pi**2:.0f}"
    print(line)
print("== 5. Theorem P: intervals of length floor(N^(3/2) J^(-5/2) / a) above L: total projected error <= C_2 J^2/N^2 for every bounded charge, C_2 = 64 pi^6/378 ==")
C2 = 64*math.pi**6/378
def partition2(N, J, L):
    xs = [L]; a = L; Lam = N**1.5 * J**-2.5
    while a < N:
        ell = max(1, int(Lam/a)); a = min(N, a + ell); xs.append(a)
    return xs
for N, J in ((4096, 1), (16384, 1), (65536, 1), (65536, 4), (16384, 8)):
    mu = MU[:N+1]; nn = np.arange(1, N + 1); D = 2*N + 1; L = int((N/J)**(5/6)); xs = partition2(N, J, L); m = len(xs) - 1
    th = (2*np.arange(1, J + 1) - 1)*math.pi/D
    sources = {"Moebius": mu[1:N+1].astype(float), "random": (mu[1:N+1] != 0)*rng.choice([-1.0, 1.0], size=N), "positive": (mu[1:N+1] != 0).astype(float)}
    line = f"   N = {N:5d} J = {J}: L = {L}, m = {m} (<= sqrt N J^(5/2) = {int(math.sqrt(N)*J**2.5)}), lengths {xs[1]-xs[0]} ... {xs[-1]-xs[-2]};  budget C_2 J^2/N^2 = {C2*J**2/N**2:.3e}"
    for name, s in sources.items():
        Ms = np.concatenate([[0], np.cumsum(s)]); cvec = s/nn; tot = np.zeros(J); esum = 0.0
        for a, b in zip(xs[:-1], xs[1:]):
            qr = np.where((nn > a) & (nn <= b), cvec, 0.0); qr[b-1] -= (Ms[b]-Ms[a])/b; rr = modes_of(qr, N, J); tot += rr; esum += float(np.sum(rr**2))
        Q = np.array([Ms[b] - Ms[a] for a, b in zip(xs[:-1], xs[1:])]); x = np.array(xs[1:], dtype=float)
        WP = float(np.sum(((2/math.sqrt(D))*((Ms[L]/L)*np.sin(L*th) + (Q/x) @ np.sin(np.outer(x, th))))**2))
        Wtrue = float(np.sum(((2/math.sqrt(D))*(Ms[L]/L)*np.sin(L*th) + modes_of(np.where(nn > L, cvec, 0.0), N, J))**2))
        line += f"\n      {name:8s}: total error {float(np.sum(tot**2)):.3e} (sum of per-interval {esum:.3e}), ratio to budget {float(np.sum(tot**2))/(C2*J**2/N**2):.2e};  W = {Wtrue:.4e}, W^P = {WP:.4e}, |sqrt W - sqrt W^P| = {abs(math.sqrt(Wtrue)-math.sqrt(WP)):.2e} <= {math.sqrt(C2)*J/N:.2e};  increment coherence {Q.sum()**2/float(np.sum(Q**2)):.2f}"
    print(line)
