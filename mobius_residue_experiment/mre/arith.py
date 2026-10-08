"""Section 9A: the signed arithmetic remainder.  Exact rational arithmetic for fixtures and small cutoffs, float64 with
identity checks for the sweep, ball arithmetic (python-flint arb) for certified values at selected cutoffs.

Definitions (manuscript, 7 October revision, Sections 1, 12, 13, 18, as restated in the master prompt):
  R(N) = sum_{k<N} M(k)^2/(k(k+1)) + M(N)^2/N,   L = floor(N/6),   K_red = floor(N^(1/6)),   H = h(N) - h(L),
  V_N = sum_{k=L}^{N-1} (h(N) - h(k))^2,   kappa_N = 2 M(L) H + L H^2,   delta_N = R(N) - R(L) = kappa_N + V_N,
  A_N = sum_{L<a,b<=N, g=(a,b), a/g<=K_red, b/g<=K_red} mu(a) mu(b) (min(a,b) - L)/(ab)   (selected family),
  C_N = V_N - A_N (complementary aggregate),   Gamma_N = delta_N - A_N = kappa_N + C_N,
  S(x) = max{1, R(n): n <= x},  B = 1 + S,  c_N = [delta_N]_+ / B(K_red)^5,  Gamma_N / B(K_red)^5 (signed diagnostic).
"""
import math, json
from fractions import Fraction as Fr
import numpy as np

def sieve_mu(X):
    """Mobius function mu(0..X) as int8 by a prime sieve."""
    s = np.ones(X + 1, dtype=bool); s[:2] = False
    for p in range(2, int(X**0.5) + 1):
        if s[p]: s[p*p::p] = False
    primes = np.nonzero(s)[0]
    mu = np.ones(X + 1, dtype=np.int8)
    for p in primes:
        mu[p::p] *= -1; mu[p*p::p*p] = 0
    mu[0] = 0
    return mu

MU_FIRST_30 = [1, -1, -1, 0, -1, 1, -1, 0, 0, 1, -1, 0, -1, 1, 1, 0, -1, 0, -1, 0, 1, 1, -1, 0, 0, 1, 0, 0, -1, -1]

def check_mu(mu):
    """Initial values test for the sieve (n = 1..30)."""
    return [int(v) for v in mu[1:31]] == MU_FIRST_30

def isqrt6(N):
    """Exact integer sixth root floor(N^(1/6)), correct at perfect sixth powers."""
    x = int(round(N**(1/6)))
    while x**6 > N: x -= 1
    while (x + 1)**6 <= N: x += 1
    return x

# ------------------------------------------------------------------ exact rational pieces (small N)
def R_pairwise_exact(a):
    """R_a(N) = sum_{m,n<=N} a(m)a(n)/max(m,n), exact; a is a list of ints/Fractions indexed from 1 (a[0] unused)."""
    N = len(a) - 1
    return sum(Fr(a[m])*Fr(a[n])/max(m, n) for m in range(1, N + 1) for n in range(1, N + 1))

def R_cumulative_exact(a):
    """R_a(N) by the cumulative formula, exact."""
    N = len(a) - 1; A = 0; tot = Fr(0)
    for k in range(1, N + 1):
        A += a[k]
        tot += Fr(A*A, k*(k + 1)) if k < N else Fr(A*A, N)
    return tot

def h_exact(mu, N):
    """h(k) = sum_{n<=k} mu(n)/n for k = 0..N as Fractions (prefix list)."""
    h = [Fr(0)]
    for n in range(1, N + 1): h.append(h[-1] + Fr(int(mu[n]), n))
    return h

def V_direct_exact(mu, N):
    L = N//6; h = h_exact(mu, N)
    return sum((h[N] - h[k])**2 for k in range(L, N))

def A_pairsum_exact(mu, N):
    """Selected family by its defining pair sum (exact, O(N^2))."""
    L = N//6; K = isqrt6(N); tot = Fr(0)
    for a in range(L + 1, N + 1):
        if mu[a] == 0: continue
        for b in range(L + 1, N + 1):
            if mu[b] == 0: continue
            g = math.gcd(a, b)
            if a//g <= K and b//g <= K:
                tot += Fr(int(mu[a])*int(mu[b])*(min(a, b) - L), a*b)
    return tot

def kappa_exact(mu, N):
    L = N//6; h = h_exact(mu, N); M = int(np.sum(mu[1:L+1])); H = h[N] - h[L]
    return 2*M*H + L*H*H

def exact_case(mu, N):
    """All Section-9A quantities at one N as exact rationals (N <= ~600 practical)."""
    L = N//6; K = isqrt6(N)
    a = [0] + [int(v) for v in mu[1:N+1]]; aL = [0] + [int(v) for v in mu[1:L+1]]
    RN = R_cumulative_exact(a); RL = R_cumulative_exact(aL) if L >= 1 else Fr(0)
    V = V_direct_exact(mu, N); A = A_pairsum_exact(mu, N); kap = kappa_exact(mu, N)
    return dict(N=N, L=L, K_red=K, R=RN, R_L=RL, V=V, A=A, C=V - A, kappa=kap, delta=RN - RL, delta_via=kap + V, Gamma=kap + V - A)

# ------------------------------------------------------------------ float sweep with prefix arrays
class ArithSweep:
    def __init__(self, Nmax, mu=None):
        self.Nmax = Nmax
        self.mu = sieve_mu(Nmax) if mu is None else mu
        assert check_mu(self.mu), "Mobius sieve failed its initial-value test"
        mu = self.mu; n = np.arange(Nmax + 1, dtype=float); n[0] = 1.0
        self.M = np.cumsum(mu.astype(np.int64))
        self.h = np.cumsum(mu.astype(float)/n); self.h[0] = 0.0
        inner = np.cumsum(self.M.astype(float)**2/(n*(n + 1)))
        self.R = np.zeros(Nmax + 1); self.R[1:] = inner[:-1] + self.M[1:].astype(float)**2/n[1:]
        self.S = np.maximum.accumulate(np.maximum(self.R, 1.0)); self.S[0] = 1.0        # S(x) = max{1, R(n): n <= x}
        self.hsum = np.cumsum(self.h); self.h2sum = np.cumsum(self.h*self.h)                 # prefix sums for the O(1) V formula
        self.sqf = (mu != 0)
        self._F = {}
    def F(self, q):
        """prefix arrays F_{q,1}, F_{q,2}: sum_{g<=z, (g,q)=1} mu(g)^2/g^j."""
        if q not in self._F:
            g = np.arange(self.Nmax + 1); mask = self.sqf & (np.gcd(g, q) == 1); mask[0] = False
            gf = g.astype(float); gf[0] = 1.0
            self._F[q] = (np.cumsum(np.where(mask, 1.0/gf, 0.0)), np.cumsum(np.where(mask, 1.0/(gf*gf), 0.0)))
        return self._F[q]
    def pairs(self, K):
        """ordered coprime squarefree pairs (u, v) with u, v <= K."""
        us = [u for u in range(1, K + 1) if self.mu[u] != 0]
        return [(u, v) for u in us for v in us if math.gcd(u, v) == 1]
    def A_prefix(self, N):
        L = N//6; K = isqrt6(N); tot = 0.0
        for u, v in self.pairs(K):
            gmin = L//min(u, v) + 1; gmax = N//max(u, v)
            if gmin > gmax: continue
            F1, F2 = self.F(u*v)
            dF1 = F1[gmax] - F1[gmin - 1]; dF2 = F2[gmax] - F2[gmin - 1]
            tot += int(self.mu[u])*int(self.mu[v])/(u*v)*(min(u, v)*dF1 - L*dF2)
        return tot
    def V_direct(self, N):
        L = N//6; return float(np.sum((self.h[N] - self.h[L:N])**2))
    def V_prefix(self, N):
        L = N//6; hN = self.h[N]
        sh = self.hsum[N-1] - (self.hsum[L-1] if L >= 1 else 0.0); sh2 = self.h2sum[N-1] - (self.h2sum[L-1] if L >= 1 else 0.0)
        return (N - L)*hN*hN - 2*hN*sh + sh2
    def row(self, N):
        L = N//6; K = isqrt6(N); H = self.h[N] - self.h[L]; ML = int(self.M[L])
        V = self.V_direct(N); Vp = self.V_prefix(N); A = self.A_prefix(N); kap = 2*ML*H + L*H*H
        delta = self.R[N] - self.R[L]; C = V - A; Gam = kap + C
        B = 1.0 + float(self.S[K]); B5 = B**5
        return dict(N=N, L=L, K_red=K, R=float(self.R[N]), S=float(self.S[N]), M=int(self.M[N]), h=float(self.h[N]), V=V, V_prefix_formula=Vp,
                    A_selected=A, C_complement=C, carry=kap, delta=delta, delta_via_kappa_V=kap + V, Gamma=Gam, B_Kred=B,
                    c_normalized=max(delta, 0.0)/B5, Gamma_over_B5=Gam/B5, err_delta=abs(delta - (kap + V)), err_V=abs(V - Vp))
    def sweep(self, Nmin=6, Nmax=None):
        Nmax = self.Nmax if Nmax is None else Nmax
        return [self.row(N) for N in range(Nmin, Nmax + 1)]

# ------------------------------------------------------------------ ball arithmetic at selected N
def ball_case(mu, N, prec=200):
    """R, V, A, kappa, delta, C, Gamma at one N as 200-bit balls (python-flint); V by the direct positive sum."""
    from flint import arb, ctx
    ctx.prec = prec
    L = N//6; K = isqrt6(N)
    h = [arb(0)]; M = [0]
    for n in range(1, N + 1):
        h.append(h[-1] + (arb(int(mu[n]))/n if mu[n] else arb(0))); M.append(M[-1] + int(mu[n]))
    def Rb(n):
        tot = arb(0)
        for k in range(1, n): tot += arb(M[k]*M[k])/(k*(k + 1))
        return tot + (arb(M[n]*M[n])/n if n >= 1 else arb(0))
    RN = Rb(N); RL = Rb(L) if L >= 1 else arb(0)
    V = arb(0)
    for k in range(L, N):
        dk = h[N] - h[k]; V += dk*dk                                   # products, not **2: arb pow of a negative ball is nan
    H = h[N] - h[L]; kap = 2*M[L]*H + L*H*H
    # selected family by the prefix formula in balls
    A = arb(0)
    us = [u for u in range(1, K + 1) if mu[u] != 0]
    for u in us:
        for v in us:
            if math.gcd(u, v) != 1: continue
            gmin = L//min(u, v) + 1; gmax = N//max(u, v)
            if gmin > gmax: continue
            q = u*v; dF1 = arb(0); dF2 = arb(0)
            for g in range(gmin, gmax + 1):
                if mu[g] != 0 and math.gcd(g, q) == 1: dF1 += arb(1)/g; dF2 += arb(1)/(g*g)
            A += arb(int(mu[u])*int(mu[v]))/(u*v)*(min(u, v)*dF1 - L*dF2)
    return dict(N=N, L=L, K_red=K, R=RN, R_L=RL, V=V, A=A, C=V - A, kappa=kap, delta=RN - RL, delta_via=kap + V, Gamma=kap + V - A)

# ------------------------------------------------------------------ block decomposition (manuscript Section 13)
def dyadic_blocks(lo, hi):
    """dyadic intervals (U, 2U] with U a power of two (U = 1/2 for u = 1) covering the integers in (lo, hi]."""
    blocks = []; U = 0.5
    while U < hi:
        a = max(int(math.floor(U)), int(math.floor(lo))); b = min(int(math.floor(2*U)), int(math.floor(hi)))
        idx = [u for u in range(a + 1, b + 1) if u > lo and u > U and u <= 2*U and u <= hi]
        if idx: blocks.append((U, idx))
        U *= 2
    return blocks

def block_decomposition(mu, N, g_budget=None, entries_budget=2*10**8):
    """C_N = sum_g mu(g)^2 sum_{(U,V) ordered} S_{g,U,V}, S = m_U^T B m_V, with the SVD split E+ - E-.
    Returns (rows, aggregate, complete_flag, stats)."""
    L = N//6; K = isqrt6(N); rows = []; agg = 0.0; entries = 0; complete = True
    gmax = N//(K + 1)
    for g in range(1, gmax + 1):
        if mu[g] == 0: continue
        if g_budget is not None and g > g_budget: complete = False; break
        lo = L/g; hi = N/g
        idx_all = [u for u in range(int(math.floor(lo)) + 1, int(math.floor(hi)) + 1) if mu[u] != 0 and math.gcd(u, g) == 1 and u > lo]
        if not idx_all: continue
        blocks = dyadic_blocks(lo, hi)
        blocks = [(U, [u for u in idx if mu[u] != 0 and math.gcd(u, g) == 1]) for U, idx in blocks]
        blocks = [(U, idx) for U, idx in blocks if idx]
        for U, iu in blocks:
            for Vb, iv in blocks:
                ua = np.array(iu); va = np.array(iv)
                entries += len(ua)*len(va)
                if entries > entries_budget: complete = False; break
                mask = (np.gcd(ua[:, None], va[None, :]) == 1) & (np.maximum(ua[:, None], va[None, :]) > K)
                W = (np.minimum(ua[:, None], va[None, :]) - L/g)/(g*ua[:, None].astype(float)*va[None, :].astype(float))
                B = np.where(mask, W, 0.0)
                if not np.any(mask): continue
                mU = mu[ua].astype(float); mV = mu[va].astype(float)
                Sval = float(mU @ B @ mV)
                Ul, sig, Vt = np.linalg.svd(B, full_matrices=False)
                al = Ul.T @ mU; be = Vt @ mV
                Ep = 0.25*float(np.sum(sig*(al + be)**2)); Em = 0.25*float(np.sum(sig*(al - be)**2))
                lead = int(np.argmax(np.abs(sig*al*be))) if len(sig) else -1
                rows.append(dict(N=N, g=g, U=U, V=Vb, dim_u=len(iu), dim_v=len(iv), E_plus=Ep, E_minus=Em, signed_block=Sval, svd_error=abs(Sval - (Ep - Em)),
                                 normalized_block=g*max(U, Vb)/math.sqrt(U*Vb)*Sval, sigma_max=float(sig[0]) if len(sig) else 0.0,
                                 lead_mode=lead, lead_overlap=float(sig[lead]*al[lead]*be[lead]) if lead >= 0 else 0.0))
                agg += Sval
            if not complete: break
        if not complete: break
    return rows, agg, complete, dict(g_max=gmax, entries=entries)

# ------------------------------------------------------------------ records
def running_records(Ns, vals):
    """indices where vals attains a new running maximum (strict)."""
    rec = []; best = -math.inf
    for N, v in zip(Ns, vals):
        if v > best: best = v; rec.append((int(N), float(v)))
    return rec

def ball_str(x, digits=15):
    return x.str(digits, radius=True)

if __name__ == "__main__":
    import time
    t = time.time(); sw = ArithSweep(3000)
    # exact identities at small N
    for N in (6, 7, 12, 36, 64, 100, 128, 256):
        ex = exact_case(sw.mu, N); row = sw.row(N)
        assert ex['delta'] == ex['delta_via'], N
        assert abs(float(ex['A']) - row['A_selected']) < 1e-12 and abs(float(ex['V']) - row['V']) < 1e-12 and abs(float(ex['C']) - row['C_complement']) < 1e-12
        print(f"N={N}: exact delta = {ex['delta']} = kappa+V (True); A exact {float(ex['A']):.12f} prefix {row['A_selected']:.12f}; C = {float(ex['C']):.12f}; Gamma = {float(ex['Gamma']):.12f}")
    rows = sw.sweep(6, 3000); print("sweep to 3000:", len(rows), "rows; max err_delta", max(r['err_delta'] for r in rows), "max err_V", max(r['err_V'] for r in rows), "time", round(time.time() - t, 1))
    # block decomposition check at N = 64, 128
    for N in (64, 128):
        rws, agg, comp, st = block_decomposition(sw.mu, N); row = sw.row(N)
        print(f"blocks N={N}: {len(rws)} blocks, aggregate {agg:.12f} vs C = {row['C_complement']:.12f}, complete {comp}, max svd err {max(r['svd_error'] for r in rws):.1e}")
    b = ball_case(sw.mu, 256); print("ball N=256: delta", ball_str(b['delta']), " via", ball_str(b['delta_via']), " C", ball_str(b['C']), " Gamma", ball_str(b['Gamma']))
