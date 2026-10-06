"""
The even sector of the Weil form on the window (paper, verification notes of Section 9.2 and the external review of 2026-10-01).

For even f on [-a,a] the polar term of the explicit formula has the opposite sign to the odd sector's: with c(y) = cosh(y/2),
    Q^even(f) = <f, K_a f> + 2 <f, c>^2,
K_a the window operator with the symbol Psi_a (archimedean part plus the prime shifts), the same operator as in the odd sector.
Off the diagonal the kernel of K_a is negative (-J_Gamma(x-y), and -c_n at the shifts), so on even functions, folded to the half
window, K_a is an M-operator: its ground state is positive (Perron-Frobenius), and it overlaps the positive function c.  The Weil
criterion on even windows, Q^even >= 0 for all a, is then: K_a^even has at most one negative eigenvalue mu_1, and the rank-one
positive perturbation closes it, 1 + 2 <c, K_a^{-1} c> <= 0.  This script computes, on the even Legendre basis N_i P_{2i}(x/a)
with the arb engine of rh_weil_arb.py: the two lowest eigenvalues mu_1 <= mu_2 of K_a^even, the sign pattern of the ground
state, the secular number 1 + 2 <c, K^{-1} c>, and the even floor lambda^even(a) = min spec Q^even, for a grid of supports.

Usage: python3 rh_weil_even.py [K] [prec] [a ...]      (default K = 30, prec = 200, a = 0.2 ... 2.0)
"""
import sys, math
from flint import arb, arb_mat, acb_mat, ctx
from rh_weil_arb import OddWeilArb, gl_rule, vonmangoldt

class EvenWeilArb(OddWeilArb):
    def __init__(self, K=30, prec=200, Mx=None):
        ctx.prec = prec
        self.K = K; self.prec = prec
        self.M = 2*K + 8
        self.Mx = Mx if Mx else 2*K + 100
        self.N = [(arb(4*i + 1)/2).sqrt() for i in range(K)]            # N_i P_{2i} orthonormal on [-1,1]
        self.xs, self.ws = gl_rule(self.M)
        self.xx, self.wx = gl_rule(self.Mx)
        self.deg = 2*K - 2

    def phi_rows(self, us, scale=None):
        K, deg = self.K, self.deg
        rows = []
        for m, u in enumerate(us):
            p0, p1 = arb(1), u
            vals = [p0]
            for k in range(2, deg + 1):
                p0, p1 = p1, arb((((2*k - 1)*u*p1 - (k - 1)*p0)/k).mid())
                if k % 2 == 0: vals.append(p1)
            sc = scale[m] if scale is not None else None
            rows.append([(self.N[i]*vals[i]*sc if sc is not None else self.N[i]*vals[i]) for i in range(K)])
        return rows

    def even_parts(self, a):
        """K^even = A + sum_n P_n (no polar term), and the polar vector c_i = sqrt(a) int phi_i(x) cosh(a x/2) dx."""
        a = arb(a)
        A, P, _ = self.parts(a)
        Kmat = arb_mat(A)
        for n in P: Kmat = Kmat + P[n]
        ch = [(a*x/2).cosh() for x in self.xs]
        rows = self.phi_rows(self.xs)
        c = [sum(w*rows[m][i]*ch[m] for m, w in enumerate(self.ws))*a.sqrt() for i in range(self.K)]
        return Kmat, arb_mat(self.K, 1, c)

def eigs(M):
    E = acb_mat(M).eig(algorithm="approx")
    return sorted([float(z.real.mid()) for z in E])

if __name__ == "__main__":
    args = sys.argv[1:]
    K = int(args[0]) if args else 30
    prec = int(args[1]) if len(args) > 1 else 200
    supports = [float(x) for x in args[2:]] or [0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0]
    W = EvenWeilArb(K=K, prec=prec)
    print(f"even sector, K = {K} modes (degree {2*K-2}), prec = {prec} bits")
    print("   a     mu_1(K^even)   mu_2(K^even)   #neg   ground state sign-changes   1 + 2<c,K^-1 c>   lambda^even = min spec(K + 2 c c^T)   lambda^even_2")
    for a in supports:
        Kmat, c = W.even_parts(a)
        ev = eigs(Kmat); mu1, mu2 = ev[0], ev[1]; nneg = sum(1 for e in ev if e < 0)
        # ground state by inverse iteration on the midpoint matrix
        n = W.K
        Km = arb_mat(n, n, [arb(Kmat[i, j].mid()) for i in range(n) for j in range(n)])
        Ks = arb_mat(Km)
        for i in range(n): Ks[i, i] = Ks[i, i] - arb(mu1)*(1 - arb(2)**(-20))
        v = Ks.solve(arb_mat(n, 1, [arb(1) for _ in range(n)]))
        # evaluate the ground state on the Gauss nodes and count sign changes
        rows = W.phi_rows(W.xs)
        vals = [sum(float(v[i, 0].mid())*float(rows[m][i].mid()) for i in range(n)) for m in range(len(W.xs))]
        sgn = [1 if x > 0 else -1 for x in vals]
        changes = sum(1 for i in range(1, len(sgn)) if sgn[i] != sgn[i - 1])
        # secular number
        sol = Km.solve(arb_mat(n, 1, [arb(c[i, 0].mid()) for i in range(n)]))
        sec = 1 + 2*sum(float(c[i, 0].mid())*float(sol[i, 0].mid()) for i in range(n))
        Q = Kmat + c*c.transpose()*2
        eq = eigs(Q)
        print(f"  {a:4.2f}   {mu1:12.6f}   {mu2:12.4e}   {nneg:3d}        {changes:3d}                 {sec:12.4e}        {eq[0]:14.4e}        {eq[1]:12.4e}", flush=True)
