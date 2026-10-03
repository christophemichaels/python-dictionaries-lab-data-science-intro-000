"""Sandbox (side project, 2026-10-02): sequence versus value.

Paper statement tested: none. This is the 'strange math' sandbox of SANDBOX_STRANGE_MATH.md. The question: which parts of the
window picture depend only on the SEQUENCE of the primes (their ranks, their density) and which depend on their VALUES?
The K-mode odd Weil form is computed with four prime systems in place of the primes (positions d, weights c in the term
-2 c g_f(d) of the form; for the true primes d = k log p, c = log p / p^{k/2}):
  T   true primes and prime powers (the baseline of rh_weil_odd.py);
  R1  'two is one, three is two': every integer n >= 2 is a prime of weight log n (the rank system, Beurling style);
  R2  true positions log p, rank weights: Lambda(p_n) := log(n+1), p^{-1/2} := (n+1)^{-1/2} ('the energy is the rank');
  R3  rank positions log(n+1) for p_n, true weights ('the place is the rank, the energy is the value').
Output: the floor lambda(a) for each system at a = 0.5, ..., 1.0, the prime energy E(a) = sum c^2 (Mertens: 2a^2 for T),
and the sign of the floor. Usage: python3 rh_sandbox_rank_primes.py [K] [dps]. Log: data/sandbox_rank_primes.log.
"""
import sys, time
import mpmath as mp
from rh_weil_odd import OddWeil

def is_prime(n):
    if n < 2: return False
    p = 2
    while p*p <= n:
        if n % p == 0: return False
        p += 1
    return True

def primes_below(x):
    return [p for p in range(2, int(x)+2) if is_prime(p) and p < x]

def system(name, a):
    """List of (position, weight) with position < 2a."""
    a = mp.mpf(a); terms = []
    if name == 'T':
        n = 2
        while mp.log(n) < 2*a:
            lam = OddWeil.vonmangoldt(n)
            if lam: terms.append((mp.log(n), lam/mp.sqrt(n)))
            n += 1
    elif name == 'R1':
        n = 2
        while mp.log(n) < 2*a:
            terms.append((mp.log(n), mp.log(n)/mp.sqrt(n)))   # every integer a prime of weight log n
            n += 1
    elif name in ('R2', 'R3'):
        ps = primes_below(mp.e**(2*a) * 4)   # generous; filtered by position below
        for idx, p in enumerate(ps):
            r = idx + 1                     # rank: 2 -> 1, 3 -> 2, 5 -> 3
            k = 1
            while True:
                if name == 'R2':
                    pos = k*mp.log(p); w = mp.log(r+1)/mp.mpf(r+1)**(mp.mpf(k)/2)
                else:
                    pos = k*mp.log(r+1); w = mp.log(p)/mp.mpf(p)**(mp.mpf(k)/2)
                if pos >= 2*a: break
                terms.append((pos, w)); k += 1
    return sorted(terms)

class SandboxWeil(OddWeil):
    def matrix_system(self, a, terms):
        Q = self.matrix(a, primes=False)
        a = mp.mpf(a)
        for pos, w in terms:
            Q -= 2*w*self.gmat(pos/a)
        return Q
    def floor_system(self, a, terms):
        Q = self.matrix_system(a, terms)
        E, V = mp.eigsy(Q)
        k = min(range(self.K), key=lambda i: E[i])
        return E[k], sum(1 for e in E if e < 0)

if __name__ == '__main__':
    K = int(sys.argv[1]) if len(sys.argv) > 1 else 30
    mp.mp.dps = int(sys.argv[2]) if len(sys.argv) > 2 else 40
    W = SandboxWeil(K=K)
    print(f"# sandbox rank primes: K={K} dps={mp.mp.dps}", flush=True)
    t0 = time.time()
    for a in [0.5, 0.6, 0.7, 0.8, 0.9, 1.0]:
        row = []
        for name in ['T', 'R1', 'R2', 'R3']:
            terms = system(name, a)
            energy = sum(w*w for _, w in terms)
            lam, neg = W.floor_system(a, terms)
            row.append((name, len(terms), energy, lam, neg))
        print(f"a={a:.2f}  2a^2={2*a*a:.3f}", flush=True)
        for name, nt, energy, lam, neg in row:
            print(f"   {name:3s} terms={nt:2d}  E={mp.nstr(energy,5):>9s}  lambda={mp.nstr(lam,8):>14s}  sign={'+' if lam>0 else '-'}  neg={neg}  ({time.time()-t0:.0f}s)", flush=True)
