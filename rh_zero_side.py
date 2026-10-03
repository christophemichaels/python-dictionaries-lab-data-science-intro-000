"""
The floor and the sliding-window identity on the zero side (paper Section 7.1).

For the K-mode minimizer f = f_K at support a, the explicit formula gives Q(f) = sum_rho ghat(rho) with g = f * f~,
and every zero used here lies on the line, so ghat(1/2 + i gamma) = |F(gamma)|^2:  lambda_K = 2 sum_{gamma>0} |F(gamma)|^2
(the tail above the last zero is estimated with the Riemann-von Mangoldt density), and the dilation derivative reads
a lambda_K' = lambda_K + 2 sum_{gamma>0} gamma (|F|^2)'(gamma).  F is evaluated in closed form through spherical Bessel
functions, int_{-a}^{a} P_n(x/a) e^{itx} dx = 2a i^n j_n(ta), at 35 digits; the zeros come from mpmath.zetazero (21 digits).
Reported: the partial sums against lambda_K and a lambda_K' (the latter from rh_weil_arb.parts), the share of the floor
carried by the zeros below the horizon T* = 2 pi e^{2a}, and the zeros carrying most of it.

Usage: python3 rh_zero_side.py a K zeros.txt [bits] [max_zeros]        (zeros.txt: lines "n gamma_n")
"""
import sys, os, time, mpmath as mp
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from flint import arb, arb_mat, acb_mat, ctx
from rh_weil_arb import OddWeilArb

a_f = float(sys.argv[1]); K = int(sys.argv[2]); zfile = sys.argv[3]
prec = int(sys.argv[4]) if len(sys.argv) > 4 else 500; NMAX = int(sys.argv[5]) if len(sys.argv) > 5 else 10**9
ctx.prec = prec; a_arb = arb(str(a_f))
t0 = time.time()
W = OddWeilArb(K=K, prec=prec); Q = W.matrix(a_arb)
E = acb_mat(Q).eig(nonstop=True); E = sorted(E, key=lambda z: z.real.mid()); lam0 = E[0].real
Qm = arb_mat(K, K, [arb(Q[i, j].mid()) for i in range(K) for j in range(K)])
sigma = lam0.mid()*(1 - arb(2)**(-20))
for i in range(K): Qm[i, i] = Qm[i, i] - sigma
v = Qm.solve(arb_mat(K, 1, [arb(1) for _ in range(K)]), nonstop=True); v = Qm.solve(v, nonstop=True)
nrm = sum(v[i, 0]*v[i, 0] for i in range(K)).sqrt(); c = arb_mat(K, 1, [v[i, 0]/nrm for i in range(K)])
def quad(X): return (c.transpose()*X*c)[0, 0]
h = arb(10)**(-30)
Ap, Pp, sp = W.parts(a_arb + h); Am, Pm, sm = W.parts(a_arb - h)
lam_prime = quad((Ap - Am)/(2*h)) + sum(quad((Pp[n] - Pm[n])/(2*h)) for n in Pp) - 2*quad((sp*sp.transpose() - sm*sm.transpose())/(2*h))
alp = mp.mpf((a_arb*lam_prime).mid().str(30, radius=False)); lamK = mp.mpf(lam0.mid().str(30, radius=False))
mp.mp.dps = 35
a = mp.mpf(str(a_f)); Nn = [mp.sqrt(mp.mpf(4*i+3)/2) for i in range(K)]
coef = [mp.mpf(c[i, 0].mid().str(40, radius=False))*Nn[i] for i in range(K)]
fa = sum(coef)/mp.sqrt(a)                                   # f_K(a)
def G(t):                                                   # |F(t)|^2 = G(t)^2
    z = t*a
    return 2*mp.sqrt(a)*sum(coef[i]*(-1)**i*mp.sqrt(mp.pi/(2*z))*mp.besselj(2*i + mp.mpf(3)/2, z) for i in range(K))
hd = mp.mpf("1e-12")
zeros = []
for line in open(zfile):
    p = line.split()
    if len(p) == 2: zeros.append((int(p[0]), mp.mpf(p[1])))
zeros.sort(); zeros = zeros[:NMAX]
Ts = 2*mp.pi*mp.exp(2*a)
print(f"a = {a_f}  K = {K}  lambda_K = {mp.nstr(lamK, 10)}  a lambda_K' = {mp.nstr(alp, 10)}  f_K(a) = {mp.nstr(fa, 4)}  T* = {mp.nstr(Ts, 6)}  zeros: {len(zeros)} up to {mp.nstr(zeros[-1][1], 8)}  ({time.time()-t0:.0f}s)", flush=True)
S0 = mp.mpf(0); S1 = mp.mpf(0); rows = []; below = mp.mpf(0)
checkpoints = {10, 30, 100, 300, 1000, 2000, 3000, 4000, 5000, 6000, len(zeros)}
for n, g in zeros:
    Gp = G(g + hd); Gm = G(g - hd); G0 = (Gp + Gm)/2; dG = (Gp - Gm)/(2*hd)
    v0 = 2*G0*G0; v1 = 4*g*G0*dG
    S0 += v0; S1 += v1; rows.append((v0, n, g))
    if g < Ts: below += v0
    if n in checkpoints:
        # smooth tails from the Riemann-von Mangoldt density rho(t) = log(t/2pi)/2pi, |F|^2 averaging to 2 f(a)^2/t^2
        T = g; rho = lambda t: mp.log(t/(2*mp.pi))/(2*mp.pi)
        tail0 = 2*mp.quad(lambda t: G(t)**2*rho(t), mp.linspace(T, T + 60, 31)) + 4*fa*fa*(mp.log((T+60)/(2*mp.pi)) + 1)/(2*mp.pi*(T+60))
        tail1 = -2*tail0                                      # t d/dt (2 f(a)^2/t^2) = -4 f(a)^2/t^2
        print(f"  N={n:5d} (T={mp.nstr(g,7)}): 2 sum |F|^2 = {mp.nstr(S0, 8)}  + tail {mp.nstr(tail0, 3)} -> {mp.nstr(S0 + tail0, 8)}  vs lambda_K {mp.nstr(lamK, 8)}  ratio {mp.nstr((S0+tail0)/lamK, 7)} |"
              f"  lambda + 2 sum gamma(|F|^2)' = {mp.nstr(lamK + S1, 8)} + tail {mp.nstr(tail1, 3)} -> {mp.nstr(lamK + S1 + tail1, 8)}  vs a lambda' {mp.nstr(alp, 8)}  ratio {mp.nstr((lamK + S1 + tail1)/alp, 7)}", flush=True)
nb = sum(1 for n, g in zeros if g < Ts)
print(f"  zeros below the horizon T* = {mp.nstr(Ts, 6)}: {nb}; their share of 2 sum |F|^2: {mp.nstr(below/S0, 5)}")
rows.sort(reverse=True)
print("  largest terms 2|F(gamma)|^2 / lambda_K:  " + "  ".join(f"gamma_{n}={mp.nstr(g, 6)}: {mp.nstr(v/lamK, 4)}" for v, n, g in rows[:8]))
print("  first zeros 2|F(gamma)|^2 / lambda_K:     " + "  ".join(f"{mp.nstr(v/lamK, 3)}" for v, n, g in sorted(rows, key=lambda r: r[1])[:12]))
