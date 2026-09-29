"""
Edge-adapted P1 Galerkin solver for the odd-sector Weil form Q_a on [-a,a].

Same form as rh_weil_odd.py:
  Q(f) = c1 g(0) + 2 int_0^inf J(x) [g(0) - g(x)] dx - sum_n 2 Lambda(n) n^{-1/2} g(log n) - 2 (int f sinh(u/2) du)^2
  g(x) = int f(u) f(u-x) du,   J(x) = e^{-x/2}/(1-e^{-2x}),   c1 = -(gamma + log pi) - pi/2 - 3 log 2.

Basis: odd extensions of hat functions on a mesh of the half window that is uniform in the interior and
geometrically graded towards the endpoint (down to delta_min), so that the profile of the minimizer can be read
at log(1/delta) up to ~25 instead of the ~1/K^2 resolution of a Legendre truncation.

Coordinates: w = u - a (right edge at w = 0, centre at w = -a). Hats h_j live on [-a, 0] with nodes at
w = -delta_j; the odd partner is h~_j(w) = -h_j(-2a - w). For x >= 0 the cross-correlation of two odd-extended
hats is  Gamma_ij(x) = corr(h_i,h_j,x) + corr(h_j,h_i,x) - corr(h_i, h_j^R, x - 2a),  h^R(v) = h(-v),
which keeps every breakpoint near the edge exact in floating point.

All cross-correlations of piecewise-linear functions are computed exactly (Simpson on merged breakpoints);
the x-integration is Gauss-Legendre on the pieces on which Gamma_ij is a single cubic, closed-form elsewhere.
"""
import numpy as np, math, sys, json, time
from numpy.polynomial.legendre import leggauss

GAMMA = 0.57721566490153286
C1 = -(GAMMA + math.log(math.pi)) - math.pi/2 - 3*math.log(2)
GL_N = 16
_gx, _gw = leggauss(GL_N)

def J(x):
    return np.exp(-x/2)/(-np.expm1(-2*x))

def F(x):
    """int_x^inf J = artanh(e^{-x/2}) + arctan(e^{-x/2}), written to keep precision at small x."""
    w = np.exp(-x/2)
    return 0.5*(np.log1p(w) - np.log(-np.expm1(-x/2))) + np.arctan(w)

class PL:
    """Piecewise-linear function with compact support: breakpoints xp (sorted), values fp; zero outside."""
    __slots__ = ("xp", "fp")
    def __init__(self, xp, fp):
        self.xp = np.asarray(xp, float); self.fp = np.asarray(fp, float)
    def __call__(self, u):
        return np.interp(u, self.xp, self.fp, left=0.0, right=0.0)
    def refl(self):
        return PL(-self.xp[::-1], self.fp[::-1])

def corr(p, q, xs):
    """int p(u) q(u - x) du for each x in xs; exact for piecewise-linear p, q."""
    xs = np.atleast_1d(np.asarray(xs, float))
    n = len(xs)
    B = np.concatenate([np.broadcast_to(p.xp, (n, len(p.xp))), q.xp[None, :] + xs[:, None]], axis=1)
    lo = np.maximum(p.xp[0], q.xp[0] + xs); hi = np.minimum(p.xp[-1], q.xp[-1] + xs)
    B = np.minimum(np.maximum(B, lo[:, None]), hi[:, None])
    B.sort(axis=1)
    L = B[:, :-1]; R = B[:, 1:]; Mid = 0.5*(L + R); d = R - L
    def prod(U):
        return p(U.ravel()).reshape(U.shape) * q((U - xs[:, None]).ravel()).reshape(U.shape)
    return (d/6*(prod(L) + 4*prod(Mid) + prod(R))).sum(axis=1)

def corr_support(p, q):
    return (p.xp[0] - q.xp[-1], p.xp[-1] - q.xp[0])

def corr_breaks(p, q):
    return (p.xp[:, None] - q.xp[None, :]).ravel()

def mesh(a, n_int=200, rho=0.1, r=0.5, dmin=1e-10, r1=None, dswitch=None):
    """Distances delta_j = a - y_j from the right edge, decreasing from a (centre) to 0 (edge).
    Uniform on [rho a, a]; then geometric with ratio r1 down to dswitch*a (if given), then ratio r down to dmin."""
    d = [np.linspace(a, rho*a, n_int + 1)]                            # uniform in the interior, includes rho*a
    top = rho*a
    if r1 is not None and dswitch is not None:
        m1 = int(math.ceil(math.log(top/(dswitch*a))/math.log(1/r1)))
        d.append(top*r1**np.arange(1, m1 + 1)); top = d[-1][-1]
    m = int(math.ceil(math.log(top/dmin)/math.log(1/r)))
    d.append(top*r**np.arange(1, m + 1))
    d.append([0.0])
    return np.concatenate(d)

class EdgeFEM:
    def __init__(self, a, deltas, primes=True, verbose=False):
        self.a = float(a); self.d = np.asarray(deltas, float); self.primes = primes
        w = -self.d                                                # node coordinates, increasing
        self.w = w; self.n = len(w) - 1                            # unknowns j = 1..n (w_0 = -a is the centre, f = 0)
        hats = []
        for j in range(1, self.n + 1):
            if j < self.n: hats.append(PL([w[j-1], w[j], w[j+1]], [0, 1, 0]))
            else:          hats.append(PL([w[j-1], w[j]], [0, 1]))  # half hat at the edge
        self.h = hats; self.hR = [h.refl() for h in hats]
        self.verbose = verbose

    def gamma(self, i, j, xs, ss):
        """Gamma_ij at x = xs (>= 0), with ss = xs - 2a supplied separately for precision."""
        hi, hj = self.h[i], self.h[j]
        return corr(hi, hj, xs) + corr(hj, hi, xs) - corr(hi, self.hR[j], ss)

    def dij(self, i, j):
        """2 int_0^inf J(x) [Gamma_ij(0) - Gamma_ij(x)] dx  and  Gamma_ij(0)."""
        a2 = 2*self.a
        hi, hj = self.h[i], self.h[j]; hRj = self.hR[j]
        G0 = 2*corr(hi, hj, [0.0])[0]
        # breakpoints: (x, s) pairs; corr ones exact in x, conv ones exact in s
        bx = np.concatenate([corr_breaks(hi, hj), corr_breaks(hj, hi)])
        bx = bx[bx >= 0]
        bs = corr_breaks(hi, hRj)                                   # conv breakpoints in s = x - 2a
        bs = bs[bs >= -a2]
        pts = [(0.0, -a2)] + [(x, x - a2) for x in bx] + [(a2 + s, s) for s in bs]
        pts.sort()
        # dedupe
        P = [pts[0]]
        for t in pts[1:]:
            if t[0] - P[-1][0] > 1e-300: P.append(t)
        sup = [corr_support(hi, hj), corr_support(hj, hi)]
        supS = corr_support(hi, hRj)
        total = 0.0
        for (xl, sl), (xr, sr) in zip(P[:-1], P[1:]):
            xm = 0.5*(xl + xr)
            active = any(lo < xm < up for lo, up in sup) or (supS[0] < 0.5*(sl + sr) < supS[1])
            if active:
                hx = 0.5*(xr - xl); xn = 0.5*(xl + xr) + hx*_gx; sn = 0.5*(sl + sr) + hx*_gx
                total += hx*np.dot(_gw, J(xn)*(G0 - self.gamma(i, j, xn, sn)))
            elif G0 != 0.0:
                total += G0*(F(xl) - F(xr)) if xl > 0 else 0.0     # xl == 0 only happens on an active piece
        if G0 != 0.0:
            total += G0*F(P[-1][0])
        return 2*total, G0

    def assemble(self):
        n = self.n; a = self.a
        Q = np.zeros((n, n)); M = np.zeros((n, n))
        t0 = time.time()
        for i in range(n):
            for j in range(i, n):
                d, g0 = self.dij(i, j)
                Q[i, j] = Q[j, i] = d + C1*g0
                M[i, j] = M[j, i] = g0
            if self.verbose and i % 50 == 0:
                print(f"  row {i}/{n}  ({time.time()-t0:.0f}s)", flush=True)
        # primes: Lambda(n) n^{-1/2} for log n < 2a
        if self.primes:
            nn = 2
            while math.log(nn) < 2*a:
                lam = vonmangoldt(nn)
                if lam:
                    x = math.log(nn); s = x - 2*a
                    for i in range(n):
                        for j in range(i, n):
                            g = self.gamma(i, j, [x], [s])[0]
                            Q[i, j] -= 2*lam/math.sqrt(nn)*g
                            if j != i: Q[j, i] = Q[i, j]
                nn += 1
        # polar: s_i = 2 int h_i(w) sinh((w+a)/2) dw
        gx, gw = leggauss(12)
        sv = np.zeros(n)
        for i, h in enumerate(self.h):
            for l, r in zip(h.xp[:-1], h.xp[1:]):
                hx = 0.5*(r - l); u = 0.5*(l + r) + hx*gx
                sv[i] += hx*np.dot(gw, h(u)*np.sinh((u + a)/2))
        sv *= 2
        Q -= 2*np.outer(sv, sv)
        self.Q, self.M, self.sv = Q, M, sv
        return Q, M

    def solve(self):
        Q, M = self.Q, self.M
        S = 1/np.sqrt(np.diag(M))
        Qs = S[:, None]*Q*S[None, :]; Ms = S[:, None]*M*S[None, :]
        L = np.linalg.cholesky(Ms)
        X = np.linalg.solve(L, Qs); A = np.linalg.solve(L, X.T).T
        A = 0.5*(A + A.T)
        E, Z = np.linalg.eigh(A)
        z = Z[:, 0]
        c = S*np.linalg.solve(L.T, z)
        if c[len(c)//2] < 0: c = -c
        self.E, self.f = E, c                                        # f = nodal values at w_1..w_n
        return E[0], c

def vonmangoldt(n):
    for p in range(2, n + 1):
        if n % p == 0:
            m = n
            while m % p == 0: m //= p
            return math.log(p) if m == 1 else 0.0
    return 0.0

def profile(fem, fracs):
    """f at y/a = frac (i.e. delta = a(1-frac)) by linear interpolation on the mesh."""
    y = -fem.d[1:]          # w-coordinates of unknowns (increasing)
    out = {}
    for fr in fracs:
        wq = -fem.a*(1 - fr)
        out[str(fr)] = float(np.interp(wq, y, fem.f, left=0.0))
    return out

if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--a", type=float, default=0.3)
    ap.add_argument("--primes", type=int, default=0)
    ap.add_argument("--n_int", type=int, default=200)
    ap.add_argument("--rho", type=float, default=0.1)
    ap.add_argument("--r", type=float, default=0.5)
    ap.add_argument("--dmin", type=float, default=1e-10)
    ap.add_argument("--r1", type=float, default=None)
    ap.add_argument("--dswitch", type=float, default=None, help="switch from ratio r1 to r at delta = dswitch*a")
    ap.add_argument("--out", type=str, default="")
    args = ap.parse_args()
    d = mesh(args.a, args.n_int, args.rho, args.r, args.dmin, args.r1, args.dswitch)
    t0 = time.time()
    fem = EdgeFEM(args.a, d, primes=bool(args.primes), verbose=True)
    fem.assemble(); lam, f = fem.solve()
    print(f"a={args.a} primes={args.primes} nodes={fem.n} lambda={lam:.10e} f(a)={f[-1]:.8f} ({time.time()-t0:.0f}s)")
    fr = [0.5, 0.9, 0.95, 0.98, 0.99, 0.995, 0.999, 1.0]
    print("profile", json.dumps(profile(fem, fr)))
    # edge law: t = log(a/delta) against f
    dd = fem.d[1:]; msk = (dd > 0) & (dd < args.rho*args.a)
    t = np.log(args.a/dd[msk]); fv = f[msk]
    print(" t         delta        f          f^2*t      1/f^2")
    for tt, de, ff in zip(t, dd[msk], fv):
        print(f"{tt:7.3f}  {de:11.3e}  {ff:10.6f}  {ff*ff*tt:9.5f}  {1/(ff*ff):9.4f}")
    if args.out:
        json.dump({"a": args.a, "primes": args.primes, "lambda": lam, "deltas": dd.tolist(), "f": f.tolist(),
                   "profile": profile(fem, fr)}, open(args.out, "w"))
