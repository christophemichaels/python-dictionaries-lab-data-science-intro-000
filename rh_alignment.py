"""
The zeros' rotation and the entries (paper Computation 7.18): the tail of the floor carries the average of cos(2 gamma a)
over the zeros, weighted by the tail envelope 1/(gamma^2 log(gamma/2pi)).  R(a)/W = <cos(2 gamma a)> from the first
6000 zeros, on a grid of a and at the entries a_n = (1/2) log n; the zeros anti-align (R/W = -0.5 to -0.6) at the
primes 2, 3, 5, 7, 11 and nowhere else, and around an entry R(a) is a V-shaped cusp: the kink of the floor at the
entries seen from the zero side.  Usage: python3 rh_alignment.py [data/zeros_6000.txt]
"""
import sys, math, numpy as np
zf = sys.argv[1] if len(sys.argv) > 1 else "data/zeros_6000.txt"
g = np.array([float(l.split()[1]) for l in open(zf)])
w = 1/(g**2*np.log(g/(2*np.pi))); W = w.sum()
a = np.arange(0.30, 1.30, 0.0005); R = np.array([(np.cos(2*g*x)*w).sum() for x in a])/W
print(f"<cos(2 gamma a)> over the tail-weighted zeros ({len(g)} zeros): mean |R/W| = {np.abs(R).mean():.4f}, 95th percentile {np.percentile(np.abs(R), 95):.4f}, min {R.min():.4f} at a = {a[R.argmin()]:.4f}")
for n in (2, 3, 4, 5, 7, 8, 9, 11, 13):
    an = 0.5*math.log(n); i = np.argmin(abs(a - an)); j = slice(max(0, i-40), i+41); k = j.start + R[j].argmin()
    q = min(p for p in range(2, n+1) if n % p == 0); m = n
    while m % q == 0: m //= q
    lam = math.log(q) if m == 1 else 0.0
    print(f"  n={n:2d}  a_n={an:.4f}  R/W(a_n)={R[i]:+.4f}  min within 0.02: {R[k]:+.4f} at a={a[k]:.4f}  Lambda(n)/sqrt(n)={lam/math.sqrt(n):.3f}")
a3 = 0.5*math.log(3); ds = np.arange(-0.004, 0.0041, 0.001)
print("around the entry of 3, offsets " + ", ".join(f"{d:+.3f}" for d in ds) + ":")
print("  all zeros, weight 1/(gamma^2 log):      " + "  ".join(f"{(np.cos(2*g*(a3+d))*w).sum()/W:+.4f}" for d in ds))
Ts = 2*math.pi*math.exp(2*a3); m = g > 3*Ts; w3 = np.where(m, 1/g**2, 0.0); W3 = w3.sum()
print("  zeros above 3T* only, weight 1/gamma^2: " + "  ".join(f"{(np.cos(2*g*(a3+d))*w3).sum()/W3:+.4f}" for d in ds))
