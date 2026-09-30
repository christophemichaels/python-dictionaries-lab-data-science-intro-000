"""
Ball-arithmetic engine for the odd-sector Weil form Q_a on [-a,a]  (python-flint / arb).

Same form, basis and normalization as rh_weil_odd.py:
  Q_a(f) = W_inf(g) - sum_n Lambda(n) n^{-1/2} (g(log n)+g(-log n)) + polar,   g = f * f~,
  W_inf(g) = -(gamma+log pi) g(0) + 2 int_0^inf [g(0) e^{-2x} - g(x) e^{-x/2}] / (1-e^{-2x}) dx,
  polar (odd f) = -2 (int f sinh(x/2))^2,
  f(y) = a^{-1/2} sum_i c_i N_i P_{2i+1}(y/a),  N_i = sqrt((4i+3)/2).

Differences from the mpmath engine: every number is an arb ball (rounding is tracked rigorously), the
inner Gauss-Legendre rule has M = 2K+8 nodes (exact for the polynomial products), the outer rule for the
archimedean integral has Mx = 2K+100 nodes (its truncation error is not enclosed; it is far below the
working precision for a <= 1.5), and the lowest eigenvalue is returned as a rigorous enclosure of the
lowest eigenvalue of the assembled ball matrix (acb_mat.eig).  The Legendre recurrence is evaluated with a midpoint
reset at each degree (its ball version doubles the radius per degree), so the enclosure covers the assembly and the
eigenvalue computation but not the polynomial evaluations, which are 500-bit floating point.  Hundreds of times
faster than mpmath.

Usage: python3 rh_weil_arb.py a K [prec_bits]
"""
import sys, time, math
from flint import arb, arb_mat, acb_mat, ctx

def gl_rule(M):
    """Gauss-Legendre nodes and weights on [-1,1] as arb balls."""
    xs, ws = [], []
    for k in range(M):
        x, w = arb.legendre_p_root(M, k, weight=True); xs.append(x); ws.append(w)
    return xs, ws

class OddWeilArb:
    def __init__(self, K=40, prec=500, Mx=None):
        ctx.prec = prec
        self.K = K; self.prec = prec
        self.M = 2*K + 8
        self.Mx = Mx if Mx else 2*K + 100
        self.N = [(arb(4*i + 3)/2).sqrt() for i in range(K)]
        self.xs, self.ws = gl_rule(self.M)
        self.xx, self.wx = gl_rule(self.Mx)
        self.deg = 2*K - 1

    def phi_rows(self, us, scale=None):
        """List of rows [N_i P_{2i+1}(u)]_i for u in us, optionally each row multiplied by scale[m]."""
        K, deg = self.K, self.deg
        rows = []
        for m, u in enumerate(us):
            p0, p1 = arb(1), u
            vals = [p1]
            for k in range(2, deg + 1):
                p0, p1 = p1, arb((((2*k - 1)*u*p1 - (k - 1)*p0)/k).mid())   # midpoint reset: the ball
                if k % 2 == 1: vals.append(p1)                                # recurrence doubles radii per degree
            sc = scale[m] if scale is not None else None
            rows.append([(self.N[i]*vals[i]*sc if sc is not None else self.N[i]*vals[i]) for i in range(K)])
        return rows

    def gmat(self, s):
        """g_ij(s) = int_{s-1}^{1} phi_i(u) phi_j(u-s) du, symmetrized; s in [0,2] (arb)."""
        K = self.K
        if s >= 2: return arb_mat(K, K)
        lo, hi = s - 1, arb(1)
        half = (hi - lo)/2; mid = (hi + lo)/2
        us = [mid + half*x for x in self.xs]
        A = arb_mat(self.phi_rows(us))                                   # M x K
        B = arb_mat(self.phi_rows([u - s for u in us], scale=[w*half for w in self.ws]))
        G = A.transpose()*B
        return (G + G.transpose())/2

    def matrix(self, a, primes=True):
        a = arb(a) if not isinstance(a, arb) else a
        K = self.K
        Q = arb_mat(K, K)
        for x, w in zip(self.xx, self.wx):
            xv = a + a*x                                                  # x in [0, 2a]
            G = self.gmat(xv/a)
            e2 = (-2*xv).exp(); eh = (-xv/2).exp(); den = 1 - e2
            c = 2*w*a/den
            Q = Q - G*(c*eh)
            d = c*e2
            for i in range(K): Q[i, i] = Q[i, i] + d
        tail = -(1 - (-4*a).exp()).log()/2
        diag = -(arb.const_euler() + arb.pi().log()) + 2*tail
        for i in range(K): Q[i, i] = Q[i, i] + diag
        if primes:
            n = 2
            while math.log(n) < 2*float(a.mid()):
                lam = vonmangoldt(n)
                if lam:
                    ln = arb(n).log()
                    Q = Q - self.gmat(ln/a)*(2*lam/arb(n).sqrt())
                n += 1
        # polar
        sh = [(a*x/2).sinh() for x in self.xs]
        rows = self.phi_rows(self.xs)
        s = [sum(w*rows[m][i]*sh[m] for m, w in enumerate(self.ws))*a.sqrt() for i in range(K)]
        S = arb_mat(K, 1, s)
        Q = Q - S*S.transpose()*2
        return Q

    def floor(self, a, eigvec=True):
        a = arb(a); Q = self.matrix(a)
        E = acb_mat(Q).eig(nonstop=True)
        E = sorted(E, key=lambda z: z.real.mid())
        lam0, lam1 = E[0].real, E[1].real
        fa = None
        if eigvec:                                  # inverse iteration on the midpoint matrix (not rigorous)
            K = self.K
            Qm = arb_mat(K, K, [arb(Q[i, j].mid()) for i in range(K) for j in range(K)])
            for shift in (1 - arb(2)**(-20), arb(1)/2):
                sigma = lam0.mid()*shift
                Qs = arb_mat(Qm)
                for i in range(K): Qs[i, i] = Qs[i, i] - sigma
                try:
                    v = Qs.solve(arb_mat(K, 1, [arb(1) for _ in range(K)]), nonstop=True)
                    if any(v[i, 0].is_nan() for i in range(K)): continue
                    nrm = sum(v[i, 0]*v[i, 0] for i in range(K)).sqrt()
                    c = [v[i, 0]/nrm for i in range(K)]
                    if c[K//2] < 0: c = [-x for x in c]
                    fa = sum(c[i]*self.N[i] for i in range(K))/a.sqrt()
                    break
                except ZeroDivisionError:
                    continue
        return lam0, lam1, fa

def vonmangoldt(n):
    for p in range(2, n + 1):
        if n % p == 0:
            m = n
            while m % p == 0: m //= p
            return arb(p).log() if m == 1 else 0
    return 0

if __name__ == "__main__":
    a = float(sys.argv[1]); K = int(sys.argv[2]); prec = int(sys.argv[3]) if len(sys.argv) > 3 else 500
    t0 = time.time()
    W = OddWeilArb(K=K, prec=prec)
    t1 = time.time()
    lam0, lam1, fa = W.floor(a)
    print(f"a={a} K={K} prec={prec}: lambda0 = {lam0.str(25)}  radius {lam0.rad().str(3)}  lambda1 = {lam1.str(12)}  f(a) = {fa.str(12) if fa is not None else None}  (rules {t1-t0:.0f}s, total {time.time()-t0:.0f}s)", flush=True)
