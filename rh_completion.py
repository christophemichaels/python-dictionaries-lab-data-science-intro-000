"""The completion of the prime construction (SUBPOWER_CLOSURE.md, Section 8). For P >= sqrt N every n <= N that is not P-smooth
has exactly one prime factor above P, so c_N = c_{Q_P,N} - sum_{P<p<=N} p^{-1} E_p c_{N/p} with the complete Moebius vectors at the
small horizons. Computes, at N = 1000, 2000, 4000: the first-mode energy before and after the last packet (P = N/2), the packet's
own first-mode energy against 0.30338 N / log^2 N; for P = sqrt N, N/4, N/2 the populations of the unfinished state, of the completion
and of the survivor in the first J modes with their cosine; and the identity that governs the completion, the mode form of
mu log = -(Lambda * mu) with the uniform completion sum_d d^{-1} E_d c_{N/d} = e_1 removed:
   b^log_j := sum mu(n) log n v_j(n)/n  =  -v_j(1) - sum_{d<=N} ((Lambda(d) - 1)/d) (V^T E_d c_{N/d})_j.
Usage: python3 rh_completion.py"""
import math, numpy as np
def run(N):
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
    nn = np.arange(1, N + 1); jj = np.arange(1, N + 1); theta = (2*jj - 1) * math.pi / (2*N + 1); lam = 1 / (4 * np.sin(theta/2)**2)
    V = (2 / math.sqrt(2*N + 1)) * np.sin(np.outer(nn, theta)); c = mu[1:N+1] / nn; b = V.T @ c
    I = 0.0; xs = np.linspace(0.5, 1, 20001); I = np.trapezoid(np.sin(math.pi*xs/2)/xs, xs)
    print(f"N = {N}:  integral int_(1/2)^1 sin(pi x/2)/x dx = {I:.8f};  8 I^2/pi^2 = {8*I**2/math.pi**2:.8f}")
    for P in (int(math.ceil(math.sqrt(N))), N // 4, N // 2):
        cQ = np.where(Pplus[1:N+1] <= P, c, 0.0)                         # the unfinished source: P-smooth part
        comp = np.zeros(N)
        for p in primes[primes > P]:
            K = N // p; comp[p*np.arange(1, K+1) - 1] += c[:K] / p         # sum_{P<p<=N} p^{-1} E_p c_{N/p}, complete small sources
        assert np.max(np.abs(cQ - comp - c)) < 1e-15
        a = V.T @ cQ; d = V.T @ comp
        line = f"   P = {P:5d}: "
        for J in (1, 10, 50):
            na, nd, ns = float(np.sum(a[:J]**2)), float(np.sum(d[:J]**2)), float(np.sum((a[:J]-d[:J])**2))
            cos = float(np.dot(a[:J], d[:J]) / math.sqrt(na*nd))
            line += f"J={J:2d}: unfinished {lam[J-1]*na:9.4f}  completion {lam[J-1]*nd:9.4f}  survivor {lam[J-1]*ns:8.5f}  cos {cos:.6f} | "
        print(line + "(entries are lambda_J times the population of the first J modes)")
        if P == N // 2:
            print(f"      first mode: before {lam[0]*a[0]**2:.6f}  after {lam[0]*b[0]**2:.6f};  packet lambda_1 d_1^2 = {lam[0]*d[0]**2:.6f} = {lam[0]*d[0]**2/(N/math.log(N)**2):.5f} N/log^2 N  (asymptotic 8I^2/pi^2 = {8*I**2/math.pi**2:.5f})")
    # the governing identity
    unif = np.zeros(N); fluct = np.zeros(N)
    for dd in range(1, N + 1):
        K = N // dd; idx = dd*np.arange(1, K+1) - 1
        unif[idx] += c[:K] / dd; fluct[idx] += (Lam[dd] - 1) / dd * c[:K]
    e1 = np.zeros(N); e1[0] = 1.0
    print(f"   uniform completion sum_d d^(-1) E_d c_(N/d) = e_1: max error {np.max(np.abs(unif - e1)):.1e}")
    blog = V.T @ (c * np.log(nn)); rhs = -(V.T @ e1) - (V.T @ fluct)
    print(f"   identity b^log = -V^T e_1 - sum_d ((Lambda(d)-1)/d) V^T E_d c_(N/d): max error {np.max(np.abs(blog - rhs)):.1e}")
    for J in (1, 10, 50):
        print(f"      J = {J:2d}: lambda_J * population of  b^log {lam[J-1]*float(np.sum(blog[:J]**2)):9.5f}   unit term {lam[J-1]*float(np.sum((V.T @ e1)[:J]**2)):9.5f}   fluctuation term {lam[J-1]*float(np.sum((V.T @ fluct)[:J]**2)):9.5f}   log N * b: {lam[J-1]*float(np.sum((math.log(N)*b)[:J]**2)):9.5f}")
for N in (1000, 2000, 4000): run(N)
