"""What comes out the other end when every zero is forced onto the point (THE_OTHER_END.md).
Four consequences of RH on the integers, checked numerically: Schoenfeld's bounds on pi(x) and psi(x), the Moebius sum at the
square-root scale, Robin's and Lagarias's divisor-sum inequalities. Usage: python3 rh_other_end.py [X=10^7] [N=10^6]"""
import sys, math
import numpy as np
X = int(float(sys.argv[1])) if len(sys.argv) > 1 else 10**7
N = int(float(sys.argv[2])) if len(sys.argv) > 2 else 10**6
EG = math.exp(0.5772156649015329)
# primes, Lambda, mu up to X
s = np.ones(X + 1, dtype=bool); s[:2] = False
for p in range(2, int(X**0.5) + 1):
    if s[p]: s[p*p::p] = False
primes = np.nonzero(s)[0]
pi = np.cumsum(s)
# Lambda(n): log p at prime powers
Lam = np.zeros(X + 1)
for p in primes:
    q = p
    while q <= X: Lam[q] = math.log(p); q *= p
psi = np.cumsum(Lam)
# mu(n) by sieve
mu = np.ones(X + 1, dtype=np.int8); 
for p in primes:
    mu[p::p] *= -1; mu[p*p::p*p] = 0
mu[0] = 0; M = np.cumsum(mu.astype(np.int64))
x = np.arange(1, X + 1, dtype=float)
# li(x) by the series li(x) = gamma + log log x + sum (log x)^k/(k k!) (Ramanujan's form) -- use mpmath-free numeric integration via cumulative trapezoid of 1/log t from 2 plus li(2)
li2 = 1.045163780117492784844588889194613136522615578151
t = np.arange(2, X + 1, dtype=float); li = np.concatenate([[0.0], np.cumsum(np.concatenate([[0.0], (1/np.log(t[1:]) + 1/np.log(t[:-1]))/2]))]) + li2
li_full = np.concatenate([[0.0], li])   # index n -> li(n), with li(1) set to 0 (unused)
n = np.arange(0, X + 1, dtype=float)
with np.errstate(divide='ignore', invalid='ignore'):
    sch_pi = np.sqrt(n) * np.log(n) / (8 * math.pi); sch_psi = np.sqrt(n) * np.log(n)**2 / (8 * math.pi)
err_pi = np.abs(pi - li_full); m = np.arange(X + 1) >= 2657
print(f"pi(x) vs li(x), x <= {X:.0e}: max |pi - li| / (sqrt x log x / 8 pi) over x >= 2657 = {(err_pi[m]/sch_pi[m]).max():.4f}  (RH says < 1, Schoenfeld)")
err_psi = np.abs(psi - n); m2 = np.arange(X + 1) >= 74
print(f"psi(x) vs x:            max |psi - x| / (sqrt x log^2 x / 8 pi) over x >= 74 = {(err_psi[m2]/sch_psi[m2]).max():.4f}  (RH says < 1, Schoenfeld)")
ratio = np.abs(M[1:]) / np.sqrt(x); ratio[:99] = 0; k = int(ratio.argmax())
print(f"Moebius sum:            max |M(x)| / sqrt x over 100 <= x <= {X:.0e} = {ratio.max():.4f} at x = {k+1}  (RH says M(x) = O(x^(1/2+eps)); the trivial bound is x)")
m5 = np.arange(X + 1) >= 10**5
print(f"pi(x) vs li(x), x >= 1e5: max ratio to Schoenfeld's bound = {(err_pi[m5]/sch_pi[m5]).max():.4f}; the bound at x = {X:.0e} allows |pi - li| <= {sch_pi[X]:.0f} and the actual error there is {err_pi[X]:.1f}")
# divisor sums to N
sig = np.zeros(N + 1, dtype=np.int64)
for d in range(1, N + 1): sig[d::d] += d
nn = np.arange(N + 1, dtype=float)
with np.errstate(divide='ignore', invalid='ignore'):
    robin = sig / (EG * nn * np.log(np.log(nn)))
m3 = np.arange(N + 1) > 5040; j = int(np.argmax(np.where(m3, robin, -1)))
print(f"Robin:                  max sigma(n) / (e^gamma n log log n) over 5040 < n <= {N:.0e} = {robin[j]:.5f} at n = {j}  (RH <=> this is < 1 for every n > 5040)")
H = np.concatenate([[0.0], np.cumsum(1 / nn[1:])])
lag = sig[1:] - (H[1:] + np.exp(H[1:]) * np.log(H[1:]))
print(f"Lagarias:               max of sigma(n) - (H_n + e^(H_n) log H_n) over 1 <= n <= {N:.0e} = {lag.max():.3f} (at n = {int(lag.argmax())+1}); RH <=> this is <= 0 for every n, equality only at n = 1")
# prime gaps vs Cramer's conditional bound
gaps = np.diff(primes); pg = primes[:-1]
cr = gaps / (np.sqrt(pg) * np.log(pg)); cr[pg < 100] = 0; kk = int(cr.argmax())
print(f"prime gaps:             max (p' - p) / (sqrt p log p) over p >= 100 = {cr.max():.4f} at p = {pg[kk]} (gap {gaps[kk]}); RH gives O(1) (Cramer 1920)")
