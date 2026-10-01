"""
The second and fourth moments of the edge's spectral measure, exactly, to large supports (paper Proposition 9.22(v), Computation 9.25).

m_2 |b|^2 = |C b|^2 = sum over the second-generation points q of (sum of the weights of the two-step walks from the edge landing at q)^2,
a two-step walk being (m, m') with a step -log m then -+ log m' (the first-generation point a - log m, then inward or outward), landing
at a - log(m m') or a - log(m/m') (inside the window, not at the edge).  Its landing point is the rational m m' or m/m' (> 1), so the
walks are grouped exactly by that rational, with no floating-point coincidence test.  Likewise m_4 |b|^2 = |C^2 b|^2 groups the
three-step walks by their rational.  Both are finite sums over the entries (prime powers below e^{2a}), exact, and are compared with
the Chebyshev bounds on the mass of mu_b above the horizon's level 2a and above five horizons, and with the asymptotic bound
m_2 <= (5/12)(2a)^2 (1 + O(1/a)) of Proposition 9.22(v).

Usage: python3 rh_edge_moments.py [a ...]      (default 1.0 1.5 2.0 2.5 3.0 3.5 4.0 4.5; m_4 only where the entry count allows)
"""
import math, sys
from fractions import Fraction

def entries(a):
    """Prime powers m < e^{2a} with their weights c_m = log p / sqrt m."""
    X = math.exp(2*a); out = []
    N = int(X) + 1; sieve = bytearray([1])*(N + 1); sieve[0] = sieve[1] = 0
    for i in range(2, int(N**0.5) + 1):
        if sieve[i]: sieve[i*i::i] = bytearray(len(sieve[i*i::i]))
    for p in range(2, N + 1):
        if not sieve[p]: continue
        m = p
        while m < X: out.append((m, math.log(p)/math.sqrt(m), p)); m *= p
    return out

def moment2(a, E):
    """Exact m_2 |b|^2 = sum_q (Cb)_q^2 with q keyed by the rational landing point (numerator, denominator)."""
    X = math.exp(2*a); land = {}
    for m, cm, _ in E:
        for mp, cmp_, _ in E:
            if m*mp < X:                                   # inward-inward: a - log(m m') > -a
                key = (m*mp, 1); land[key] = land.get(key, 0.0) + cm*cmp_
            if mp < m:                                     # outward: a - log(m/m') < a, inside since m/m' < m < e^{2a}
                g = math.gcd(m, mp); key = (m//g, mp//g); land[key] = land.get(key, 0.0) + cm*cmp_
    return sum(v*v for v in land.values()), land

def moment4(a, E, land2):
    """Exact m_4 |b|^2 = sum_q (C^2 b)_q^2: three-step walks grouped by rational; (C^2 b)_q = sum over second-generation
    points r adjacent to q of (Cb)_r c_{|r-q|}, with the step from r = (num, den) by m: inward num*m/den, outward num/(den*m)."""
    X = math.exp(2*a); land = {}
    for (num, den), w in land2.items():
        for m, cm, _ in E:
            n1, d1 = num*m, den                           # inward
            if n1 < X*d1:
                g = math.gcd(n1, d1); key = (n1//g, d1//g); land[key] = land.get(key, 0.0) + w*cm
            n2, d2 = num, den*m                           # outward: position a - log(num/(den m)); inside if num/(den m) > 1 (strictly, not the edge)
            if n2 > d2:
                g = math.gcd(n2, d2); key = (n2//g, d2//g); land[key] = land.get(key, 0.0) + w*cm
    return sum(v*v for v in land.values())

if __name__ == "__main__":
    supports = [float(x) for x in sys.argv[1:]] or [1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5]
    print("   a    entries   |b|^2    2a^2    m_2      m_2/(2a)^2   (5/12)   m_2/(2a+log5)^2   m_4      m_4^(1/4)/|b|   m_4/(2a+log5)^4")
    for a in supports:
        E = entries(a); b2 = sum(c*c for _, c, _ in E); s5 = 2*a + math.log(5)
        M2, land2 = moment2(a, E); m2 = M2/b2
        row = f"{a:5.2f}  {len(E):7d}   {b2:7.3f}  {2*a*a:6.2f}  {m2:7.3f}   {m2/(2*a)**2:8.4f}    0.4167   {m2/s5**2:10.4f}"
        if len(E)*len(land2) < 3e7:
            m4 = moment4(a, E, land2)/b2
            row += f"   {m4:10.3f}   {m4**0.25/math.sqrt(b2):8.4f}   {m4/s5**4:10.4f}"
        print(row, flush=True)
