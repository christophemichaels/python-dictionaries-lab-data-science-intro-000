"""
Minimal high-precision engine for the odd-sector Weil form Q_a on [-a,a].

Q_a(f) = W_inf(g) - sum_n Lambda(n) n^{-1/2} (g(log n)+g(-log n)) + polar,  g = f * f~,
W_inf(g) = -(gamma+log pi) g(0) + 2 int_0^inf [g(0) e^{-2x} - g(x) e^{-x/2}] / (1-e^{-2x}) dx,
polar (odd f) = -2 (int f sinh(x/2))^2.

Basis: f(y) = a^{-1/2} sum_i c_i N_i P_{2i+1}(y/a), N_i = sqrt((4i+3)/2)  (orthonormal).
"""
import mpmath as mp
import numpy as np

mp.mp.dps = 50

def gl_nodes(M):
    """Gauss-Legendre nodes/weights on [-1,1] at mp precision (numpy seed + Newton)."""
    x0, _ = np.polynomial.legendre.leggauss(M)
    xs, ws = [], []
    for x in x0:
        x = mp.mpf(x)
        for _ in range(4):
            p0, p1 = mp.mpf(1), x
            for k in range(2, M+1):
                p0, p1 = p1, ((2*k-1)*x*p1 - (k-1)*p0)/k
            dp = M*(x*p1 - p0)/(x*x-1)
            x -= p1/dp
        p0, p1 = mp.mpf(1), x
        for k in range(2, M+1):
            p0, p1 = p1, ((2*k-1)*x*p1 - (k-1)*p0)/k
        dp = M*(x*p1 - p0)/(x*x-1)
        xs.append(x); ws.append(2/((1-x*x)*dp*dp))
    return xs, ws

def legendre_all(u, nmax):
    P = [mp.mpf(1), u]
    for k in range(2, nmax+1):
        P.append(((2*k-1)*u*P[-1] - (k-1)*P[-2])/k)
    return P

class OddWeil:
    def __init__(self, K=24, M=64, Mx=80):
        self.K = K
        self.N = [mp.sqrt(mp.mpf(4*i+3)/2) for i in range(K)]
        self.xs, self.ws = gl_nodes(M)       # inner (u) quadrature
        self.xx, self.wx = gl_nodes(Mx)      # outer (x) quadrature
        self.deg = 2*K-1

    def phi(self, u):
        P = legendre_all(u, self.deg)
        return [self.N[i]*P[2*i+1] for i in range(self.K)]

    def gmat(self, s):
        """g_ij(s) = int_{s-1}^{1} phi_i(u) phi_j(u-s) du, symmetrized; s in [0,2]."""
        K = self.K
        G = mp.zeros(K, K)
        if s >= 2:
            return G
        lo, hi = s-1, mp.mpf(1)
        half = (hi-lo)/2; mid = (hi+lo)/2
        for x, w in zip(self.xs, self.ws):
            u = mid + half*x
            A = self.phi(u); B = self.phi(u-s)
            ww = w*half
            for i in range(K):
                Ai = A[i]*ww
                for j in range(K):
                    G[i,j] += Ai*B[j]
        return (G + G.T)/2

    def matrix(self, a, primes=True, delete=None, shift=None):
        a = mp.mpf(a); K = self.K
        # archimedean
        A = mp.zeros(K, K)
        half = a; mid = a                     # x in [0,2a]
        for x, w in zip(self.xx, self.wx):
            xv = mid + half*x
            G = self.gmat(xv/a)
            den = 1 - mp.e**(-2*xv)
            e2 = mp.e**(-2*xv); eh = mp.e**(-xv/2)
            for i in range(K):
                for j in range(K):
                    A[i,j] += 2*w*half*((e2 if i==j else 0) - G[i,j]*eh)/den
        tail = -mp.log(1 - mp.e**(-4*a))/2       # int_{2a}^inf e^{-2x}/(1-e^{-2x}) dx
        for i in range(K):
            A[i,i] += -(mp.euler + mp.log(mp.pi)) + 2*tail
        Q = A
        # primes
        if primes:
            n = 2
            while mp.log(n) < 2*a:
                lam = self.vonmangoldt(n)
                if lam and n != delete:
                    ln = mp.log(n)
                    if shift:
                        for p, eps in shift.items():
                            m, k = n, 0
                            while m % p == 0: m //= p; k += 1
                            if m == 1 and k > 0: ln = k*(mp.log(p) + mp.mpf(eps))
                    Q -= 2*lam/mp.sqrt(n)*self.gmat(ln/a)
                n += 1
        # polar
        s = mp.zeros(K, 1)
        for x, w in zip(self.xs, self.ws):
            ph = self.phi(x); sh = mp.sinh(a*x/2)
            for i in range(K):
                s[i] += w*ph[i]*sh
        s = s*mp.sqrt(a)
        Q -= 2*s*s.T
        return Q

    @staticmethod
    def vonmangoldt(n):
        for p in range(2, n+1):
            if n % p == 0:
                m = n
                while m % p == 0: m //= p
                return mp.log(p) if m == 1 else 0
        return 0

    def floor(self, a, **kw):
        Q = self.matrix(a, **kw)
        E, V = mp.eigsy(Q)
        k = min(range(self.K), key=lambda i: E[i])
        c = V[:, k]
        fa = sum(c[i]*self.N[i] for i in range(self.K))/mp.sqrt(a)   # f(a) = a^{-1/2} sum c_i N_i
        return E[k], fa, E

if __name__ == "__main__":
    import sys, time
    W = OddWeil(K=24)
    t0 = time.time()
    for a in [0.4, 0.5, 0.6, 0.7, 0.8]:
        lam, fa, E = W.floor(a)
        print(f"a={a}: lambda={mp.nstr(lam,6)}  f(a)={mp.nstr(fa,6)}  neg-count={sum(1 for e in E if e<0)}  ({time.time()-t0:.0f}s)")
