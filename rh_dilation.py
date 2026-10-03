"""
The sliding-window identity on the minimizer (memo section 2.8).

At support a, with the K-mode minimizer c (normalized), Hellmann-Feynman gives lambda'(a) = <c, Q'(a) c>, and Q'
splits into the archimedean, the prime and the polar parts.  Each part is differentiated by a central difference
with step h = 1e-30 in ball arithmetic.  Reported: a*<c, X' c> for each part X (so that -a<c,A'c> = D_inf, the
archimedean dilation form, expected close to ||f||^2 = 1, and -a<c,P_n' c> the prime dilation forms), their sum
a*lambda', and Phi' = -lambda'/lambda against the grid value.

Usage: python3 rh_dilation.py a K [prec]
"""
import sys, time, math
from flint import arb, arb_mat, acb_mat, ctx
from rh_weil_arb import OddWeilArb

K = int(sys.argv[2]); prec = int(sys.argv[3]) if len(sys.argv) > 3 else 500
PRIMES = not (len(sys.argv) > 4 and sys.argv[4] == 'noprimes')
ctx.prec = prec
a = arb(sys.argv[1])                       # after setting the precision: an exact-to-prec ball
t0 = time.time()
W = OddWeilArb(K=K, prec=prec)
Q = W.matrix(a, primes=PRIMES)
E = acb_mat(Q).eig(nonstop=True); E = sorted(E, key=lambda z: z.real.mid()); lam0 = E[0].real
# eigenvector by inverse iteration (midpoint matrix)
Qm = arb_mat(K, K, [arb(Q[i, j].mid()) for i in range(K) for j in range(K)])
sigma = lam0.mid()*(1 - arb(2)**(-20))
for i in range(K): Qm[i, i] = Qm[i, i] - sigma
v = Qm.solve(arb_mat(K, 1, [arb(1) for _ in range(K)]), nonstop=True)
nrm = sum(v[i, 0]*v[i, 0] for i in range(K)).sqrt()
c = arb_mat(K, 1, [v[i, 0]/nrm for i in range(K)])
def quad(X): return (c.transpose()*X*c)[0, 0]
h = arb(10)**(-30)
Ap, Pp, sp = W.parts(a + h); Am, Pm, sm = W.parts(a - h); A0, P0, s0 = W.parts(a)
if not PRIMES: Pp, Pm, P0 = {}, {}, {}
dA = (Ap - Am)/(2*h); dP = {n: (Pp[n] - Pm[n])/(2*h) for n in P0}
dS = (sp*sp.transpose() - sm*sm.transpose())/(2*h)
qA, qP, qS = quad(dA), {n: quad(dP[n]) for n in dP}, quad(dS)*(-2)
lam_check = quad(A0) + sum(quad(P0[n]) for n in P0) - 2*quad(s0*s0.transpose())
lam_prime = qA + sum(qP.values()) + qS
Ts = 2*math.pi*math.exp(2*float(a.mid()))
print(f"a = {a.str(5)}  K = {K}  lambda = {lam0.str(12)}  (Rayleigh check {lam_check.str(12)})")
print(f"  archimedean:  a<c,A'c>   = {(a*qA).str(15)}      -> D_inf = {(-a*qA).str(12)}")
for n in sorted(qP): print(f"  prime {n:2d}:     a<c,P_n'c> = {(a*qP[n]).str(15)}")
print(f"  primes total: a<c,P'c>   = {(a*sum(qP.values())).str(15)}   -> D_P = {(-a*sum(qP.values())).str(12)}")
print(f"  polar:        a<c,S'c>   = {(a*qS).str(15)}")
print(f"  sum = a lambda' = {(a*lam_prime).str(15)}")
print(f"  Phi' = -lambda'/lambda = {(-lam_prime/lam0.mid()).str(10)}   Phi'/T* = {(-lam_prime/lam0.mid()/Ts).str(8)}   ({time.time()-t0:.0f}s)")
# the pencil (-Q', Q): largest generalized eigenvalue nu_max = max_f (-Q'(f))/Q(f); strong form of Conjecture A
# would be nu_max <= c T*.  Q is positive definite in the K-mode space here; eigenvalues of Q^{-1}(-Q') are real.
Qp = dA + sum(dP.values()) - dS*2
Qmid = arb_mat(K, K, [arb(Q[i, j].mid()) for i in range(K) for j in range(K)])
Qpm = arb_mat(K, K, [arb(-Qp[i, j].mid()) for i in range(K) for j in range(K)])
Mp = Qmid.solve(Qpm, nonstop=True)
ev = acb_mat(Mp).eig(algorithm="approx")
re = sorted([z.real.mid() for z in ev], key=lambda x: float(x))
print(f"  pencil (-Q', Q): eigenvalues/T*: max = {(re[-1]/Ts).str(8)}, next = {(re[-2]/Ts).str(6)}, {(re[-3]/Ts).str(6)};  min = {(re[0]/Ts).str(6)};  minimizer direction gives {(-lam_prime/lam0.mid()/Ts).str(6)}")
