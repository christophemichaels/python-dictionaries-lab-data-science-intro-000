"""
The cutoff dilation identity (PROGRAM.md item 2, exact finite-height form of the tail law), tested on the prime-free
edge-FEM critical point.  With D = -a lambda' = (1/2pi) int |F|^2 t Psi' + 2P^2 + 2PM, cutting the integral at |t| = T,
integrating by parts and using the eigen-equation Psi F = lambda F + 2P Sigma + G on |t| < T:
   D = tail(T) + (1/pi) T |F(T)|^2 (Psi(T) - lambda) + tailD(T) - 4P <s, ((xf)')_{>T}> + 2 B_T,
   tail(T) = (1/2pi) int_{|t|>T} |F|^2 (Psi - lambda),  tailD(T) = (1/2pi) int_{|t|>T} |F|^2 t Psi',
   B_T = <g_out, ((xf)')_T>  (outside force against the low-pass of the dilation generator);  D = 2 B_inf, and with
   the boundary law lambda' = -2C^2, B_inf = a C^2.  The tail law is B_inf - B_T = B_inf/(pi a T)(1 + o(1)) up to the
   explicit oscillating terms.
Every term is computed independently (B_T in x-space from the archimedean kernel).  Usage:
   python3 rh_cutoff_identity.py data/fem_a0.5_primefree.json 3,5,10,20
"""
import sys, json, math, numpy as np, mpmath as mp
from numpy.polynomial.legendre import leggauss
fn = sys.argv[1]; ks = [float(k) for k in sys.argv[2].split(",")] if len(sys.argv) > 2 else [3, 5, 10, 20]
d = json.load(open(fn)); a = d["a"]; lam = d["lambda"]
dl = np.array(d["deltas"]); f = np.array(d["f"]); x = a - dl; o = np.argsort(x); x, f = x[o], f[o]
if x[0] > 0: x = np.concatenate([[0.0], x]); f = np.concatenate([[0.0], f])
B = np.diff(f)/np.diff(x); A = f[:-1] - B*x[:-1]; fa = f[-1]
Ts = 2*math.pi*math.exp(2*a)
def fval(y):
    y = np.atleast_1d(y); s = np.sign(y); u = np.abs(y); out = np.zeros_like(u)
    m = u < a; idx = np.clip(np.searchsorted(x, u[m], side="right") - 1, 0, len(B)-1); out[m] = A[idx] + B[idx]*u[m]; return s*out
def S(t):     # F = 2i S, S(t) = int_0^a f sin(tx)
    t = np.atleast_1d(t)[:, None]
    return (-fa*np.cos(t[:, 0]*a)/t[:, 0]) + (B[None, :]*(np.sin(t*x[1:][None, :]) - np.sin(t*x[:-1][None, :]))).sum(1)/t[:, 0]**2
def Sp(t):    # S'(t)
    t = np.atleast_1d(t)[:, None]; t0 = t[:, 0]
    return fa*(a*np.sin(t0*a)/t0 + np.cos(t0*a)/t0**2) + (B[None, :]*((x[1:][None, :]*np.cos(t*x[1:][None, :]) - x[:-1][None, :]*np.cos(t*x[:-1][None, :]))/t**2 - 2*(np.sin(t*x[1:][None, :]) - np.sin(t*x[:-1][None, :]))/t**3)).sum(1)
tg = np.concatenate([np.linspace(0.02, 60, 3000), np.linspace(60, 40000, 8000)[1:]])
psi_g = np.array([float(mp.re(mp.psi(0, mp.mpf(1)/4 + 1j*mp.mpf(float(t))/2))) for t in tg])
dpsi_g = np.array([float(-mp.im(mp.polygamma(1, mp.mpf(1)/4 + 1j*mp.mpf(float(t))/2))/2) for t in tg])
Psi = lambda t: np.interp(t, tg, psi_g) - math.log(math.pi); dPsi = lambda t: np.interp(t, tg, dpsi_g)
gx, gw = leggauss(400); ys = a*(gx+1)/2; ws = a*gw/2
P = 2*np.sum(ws*fval(ys)*np.sinh(ys/2)); M = 2*np.sum(ws*fval(ys)*ys*np.cosh(ys/2))
tt = np.concatenate([np.linspace(1e-4, 200, 400001), np.linspace(200, 40000, 800001)[1:]])
F2 = 4*S(tt)**2
Dinf = np.trapezoid(F2*tt*dPsi(tt), tt)/math.pi
D = Dinf + 2*P*P + 2*P*M
print(f"{fn}: a={a} lambda={lam:.6e} P={P:.6f} M={M:.6f} D_inf={Dinf:.6f}  D = -a lambda' = {D:.6f}  (arb K-mode a lambda' from the dilation log for comparison)")
def J(s): return np.exp(-s/2)/(-np.expm1(-2*s))
def Jsmooth(s): return J(s) - 1/(2*s)
gq, gwq = leggauss(24)
def gout(y):   # prime-free: -int J(y-z) f(z) dz for y > a
    y = np.atleast_1d(y); out = np.zeros_like(y)
    for i, yy in enumerate(y):
        tot = 0.0
        for j in range(len(B)):
            z0, z1 = x[j], x[j+1]; Aj, Bj = A[j], B[j]
            tot += 0.5*((Aj + Bj*yy)*math.log((yy - z0)/(yy - z1)) - Bj*(z1 - z0))
            zz = z0 + (z1 - z0)*(gq + 1)/2; wz = (z1 - z0)*gwq/2; tot += np.sum(wz*(Aj + Bj*zz)*Jsmooth(yy - zz))
        for j in range(len(B)):
            z0, z1 = -x[j+1], -x[j]; zz = z0 + (z1 - z0)*(gq + 1)/2; wz = (z1 - z0)*gwq/2; tot += np.sum(wz*fval(zz)*J(yy - zz))
        out[i] = -tot
    return out
def lowpass_xfprime(y, T):   # ((xf)')_T(y) = (2/pi) int_0^T (-t S'(t)) sin(t y) dt
    n = int(max(40000, 60*T*(np.max(y) + a))); t = np.linspace(1e-6, T, n); h = -t*Sp(t)
    out = np.empty(len(y))
    for i0 in range(0, len(y), 200):
        yy = y[i0:i0+200]; out[i0:i0+200] = (2/math.pi)*np.trapezoid(h[None, :]*np.sin(t[None, :]*yy[:, None]), t, axis=1)
    return out
Binf = D/2
for k in ks:
    T = k*Ts
    m = tt > T
    tail = np.trapezoid(F2[m]*(Psi(tt[m]) - lam), tt[m])/math.pi; tailD = np.trapezoid(F2[m]*tt[m]*dPsi(tt[m]), tt[m])/math.pi
    bdry = T*4*S(np.array([T]))[0]**2*(Psi(np.array([T]))[0] - lam)/math.pi
    yin = ys; lp_in = lowpass_xfprime(yin, T); s_lp = 2*np.sum(ws*np.sinh(yin/2)*lp_in)          # <s, ((xf)')_T>
    polar_hp = 4*P*(-M/2 - s_lp)                                                                # 4P <s, ((xf)')_{>T}> = 4P(<s,(xf)'> - <s,((xf)')_T>)
    yo = np.concatenate([a + np.logspace(-7, -2, 50), a + np.linspace(0.011, 8, int(8/min(0.01, 0.15*2*math.pi/T)) + 1)])
    BT = 2*np.trapezoid(gout(yo)*lowpass_xfprime(yo, T), yo)
    rhs = tail + bdry + tailD - polar_hp + 2*BT
    law = D/(math.pi*a*T)
    print(f"  T/T*={k:5.1f}: D = {D:.6f}  vs  rhs = {rhs:.6f} (ratio {rhs/D:.4f});  tail = {tail:.5e} = {tail/law:.4f} of the law D/(pi a T);  boundary {bdry:+.3e}  tailD {tailD:+.3e}  4P<s,((xf)')_(>T)> {polar_hp:+.3e}  B_T = {BT:.5f} (B_inf = D/2 = {Binf:.5f})", flush=True)
