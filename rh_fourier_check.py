"""
Independent Fourier-side check of the sliding-window identity (paper Computation 5.9, memo section 2.8) at K = 12, a = 0.6.

Matrix side: a*lambda_K' = a <c, Q_K'(a) c> from rh_weil_arb.parts (central differences at step 1e-30 in ball
arithmetic), and, independently of any identity, central differences of the K-mode floor lambda_K(a +- 1e-8).
Fourier side, on the same coefficient vector c: D_inf = (1/2pi) int |F|^2 t d/dt Re psi(1/4+it/2) dt, with F(t) in
closed form through spherical Bessel functions (int_{-a}^{a} P_n(x/a) e^{itx} dx = 2a i^n j_n(ta)), the t-integral by
composite Gauss-Legendre on [0,T] and the Parseval tail (t d/dt Re psi -> 1, so (1/pi) int_0^inf |F|^2 = ||f||^2 is
used exactly and int_T^inf |F|^2 (w-1) = O(T^-3) is dropped); D_P = -sum 2 Lambda(n) n^{-1/2} log n g'(log n) from
central differences of the exact autocorrelation g (a polynomial on [0,2a]); P, M by Gauss-Legendre.  The two sides
agree to twelve digits.  The weight t d/dt Re psi(1/4+it/2) has poles at t = i(2k+1/2), so the panels near t = 0 are
short (width 0.2 on [0,4]); wide panels there cost 1e-9 in D_inf, which is 6e-5 of a*lambda'.

Usage: python3 rh_fourier_check.py [T]      (T = cutoff of the t-integral, default 800; ~1 minute)
"""
import sys, os, time, mpmath as mp
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from flint import arb, arb_mat, acb_mat, ctx
from rh_weil_arb import OddWeilArb
import rh_weil_odd as eng
K = 12; prec = 500; ctx.prec = prec; a_arb = arb("0.6")
W = OddWeilArb(K=K, prec=prec)
t0 = time.time()
def floor_vec(aa):
    Q = W.matrix(aa)
    E = acb_mat(Q).eig(nonstop=True); E = sorted(E, key=lambda z: z.real.mid()); lam0 = E[0].real
    Qm = arb_mat(K, K, [arb(Q[i, j].mid()) for i in range(K) for j in range(K)])
    sigma = lam0.mid()*(1 - arb(2)**(-20))
    for i in range(K): Qm[i, i] = Qm[i, i] - sigma
    v = Qm.solve(arb_mat(K, 1, [arb(1) for _ in range(K)]), nonstop=True)
    nrm = sum(v[i, 0]*v[i, 0] for i in range(K)).sqrt()
    return lam0, arb_mat(K, 1, [v[i, 0]/nrm for i in range(K)])
lam0, c = floor_vec(a_arb)
def quad(X): return (c.transpose()*X*c)[0, 0]
h = arb(10)**(-30)
Ap, Pp, sp = W.parts(a_arb + h); Am, Pm, sm = W.parts(a_arb - h)
dA = (Ap - Am)/(2*h); dP = {n: (Pp[n] - Pm[n])/(2*h) for n in Pp}
dS = (sp*sp.transpose() - sm*sm.transpose())/(2*h)
lam_prime = quad(dA) + sum(quad(dP[n]) for n in dP) - 2*quad(dS)
print(f"arb engine K={K} a=0.6: lambda_K = {lam0.str(15)}; a lambda_K' (parts, step 1e-30) = {(a_arb*lam_prime).str(15)}", flush=True)
# central differences of the floor itself (independent of parts)
hh = arb(10)**(-8)
lp, _ = floor_vec(a_arb + hh); lm, _ = floor_vec(a_arb - hh)
print(f"  a lambda_K' (central difference of the floor, h = 1e-8) = {(a_arb*(lp - lm)/(2*hh)).str(12)}   ({time.time()-t0:.0f}s)", flush=True)
# --- Fourier side on the same c ---
mp.mp.dps = 30
a = mp.mpf("0.6")
Nn = [mp.sqrt(mp.mpf(4*i+3)/2) for i in range(K)]
coef = [mp.mpf(c[i, 0].mid().str(40, radius=False))*Nn[i] for i in range(K)]
def f(x):
    P = eng.legendre_all(x/a, 2*K-1); return sum(coef[i]*P[2*i+1] for i in range(K))/mp.sqrt(a)
def gl(n):
    xs, ws = [], []
    for k in range(1, n+1):
        x = mp.cos(mp.pi*(k - mp.mpf(1)/4)/(n + mp.mpf(1)/2))
        for _ in range(60):
            p0, p1 = mp.mpf(1), x
            for j in range(2, n+1): p0, p1 = p1, ((2*j-1)*x*p1 - (j-1)*p0)/j
            dp = n*(x*p1 - p0)/(x*x - 1); dx = p1/dp; x -= dx
            if abs(dx) < mp.mpf(10)**(-mp.mp.dps+3): break
        xs.append(x); ws.append(2/((1 - x*x)*dp*dp))
    return xs, ws
def integrate(fun, lo, hi, xs, ws):
    m, hlen = (lo+hi)/2, (hi-lo)/2
    return hlen*sum(w*fun(m + hlen*x) for x, w in zip(xs, ws))
X40, W40 = gl(40); X12, W12 = gl(12)
mp.mp.dps = 20
def Fabs2(t):
    z = t*a
    if z == 0: return mp.mpf(0)
    s = sum(coef[i]*(-1)**i*mp.sqrt(mp.pi/(2*z))*mp.besselj(2*i + mp.mpf(3)/2, z) for i in range(K))
    return 4*a*s*s
def w(t): return -t*mp.im(mp.polygamma(1, mp.mpf(1)/4 + 1j*t/2))/2
# psi'(1/4+it/2) has poles at t = i(2k+1/2): the weight w is analytic only in the strip |Im t| < 1/2, so the panels
# near t = 0 must be short (width 0.2 on [0,4], then width 2): with width 2 throughout the first panels are only
# accurate to ~1e-9, which was the source of a 6e-5 discrepancy in a first run.
T_END = float(sys.argv[1]) if len(sys.argv) > 1 else 800; acc = mp.mpf(0); lo = mp.mpf(0)
pars = mp.mpf(0)
while lo < T_END:
    step = mp.mpf("0.2") if lo < 4 else mp.mpf(2); hi = lo + step
    acc += integrate(lambda t: Fabs2(t)*(w(t) - 1), lo, hi, X12, W12)
    pars += integrate(Fabs2, lo, hi, X12, W12); lo = hi
fa = f(a*(1 - mp.mpf(10)**-25))
print(f"  Parseval check: (1/pi) int_0^T |F|^2 + 2 f(a)^2/(pi T) = {mp.nstr(pars/mp.pi + 2*fa*fa/(mp.pi*T_END), 12)}  (should be ||f||^2 = 1 up to O(T^-2))")
norm2 = integrate(lambda x: f(x)**2, -a, a, X40, W40)
Dinf = norm2 + acc/mp.pi
def g(u):
    if abs(u) >= 2*a: return mp.mpf(0)
    return integrate(lambda x: f(x)*f(x-u), u-a, a, X40, W40)
def gp(u, h=mp.mpf("1e-12")): return (g(u+h) - g(u-h))/(2*h)
prime = [n for n in range(2, 40) if mp.log(n) < 2*a and eng.OddWeil.vonmangoldt(n)]
mp.mp.dps = 30      # g is a polynomial on [0,2a] (40-node rule exact); central differences at h = 1e-12 need 30 digits
DP = -sum(2*eng.OddWeil.vonmangoldt(n)/mp.sqrt(n)*mp.log(n)*gp(mp.log(n)) for n in prime)
mp.mp.dps = 20
P = integrate(lambda x: f(x)*mp.sinh(x/2), -a, a, X40, W40); M = integrate(lambda x: f(x)*x*mp.cosh(x/2), -a, a, X40, W40)
rhs = -(Dinf + DP + 2*P*P + 2*P*M)
print(f"  ||f||^2 = {mp.nstr(norm2, 12)}  primes {prime}")
print(f"  Fourier side: D_inf = {mp.nstr(Dinf, 15)}  D_P = {mp.nstr(DP, 15)}  2P^2+2PM = {mp.nstr(2*P*P+2*P*M, 15)}")
print(f"  -[D_inf + D_P + 2P^2 + 2PM] = {mp.nstr(rhs, 15)}")
alp = mp.mpf((a_arb*lam_prime).mid().str(40, radius=False))
print(f"  relative difference vs arb parts: {mp.nstr(abs(rhs/alp - 1), 3)}   ({time.time()-t0:.0f}s)")
