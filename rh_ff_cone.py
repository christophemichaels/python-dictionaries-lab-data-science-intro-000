"""The light cone of a curve over a finite field (FUNCTION_FIELD_CONE.md).
Usage: python3 rh_ff_cone.py q "c0,c1,...,c_{2g+1}"   (f = sum c_k x^k, monic, squarefree, odd degree)  -> JSON on stdout
Data: data/ff_cone_g{g}_q{q}.json
Hyperelliptic curve y^2 = f(x) over F_q (q odd prime), deg f = 2g+1 (one point at infinity).
Counts N_m = #X(F_{q^m}); zeta polynomial P(T) via log Z = sum N_m T^m/m; zeros alpha_i, |alpha_i| = sqrt q (Weil).
Weil form at degree cutoff D: Q(h) = sum_i |h^(theta_i)|^2 = count side (explicit formula), h supported on [-D, D].
Floor lambda(D) = min eigenvalue of the Toeplitz matrix [nu^(j-k)] of the zero measure; potential set collapses when singular.
"""
import itertools, math, sys, json
import numpy as np

q = int(sys.argv[1]) if len(sys.argv) > 1 else 3
fcoef = [int(c) for c in sys.argv[2].split(",")] if len(sys.argv) > 2 else [1, 2, 0, 1, 0, 0, 0, 1]  # f(x)=1+2x+x^3+x^7 over F_3
g = (len(fcoef) - 1 - 1) // 2
assert len(fcoef) - 1 == 2 * g + 1

def poly_gcd_deg(a, b):
    # degree of gcd over F_q (lists low->high)
    def trim(p):
        while p and p[-1] == 0: p.pop()
        return p
    a, b = trim(a[:]), trim(b[:])
    while b:
        inv = pow(b[-1], q - 2, q)
        while len(a) >= len(b):
            if a[-1]:
                c = a[-1] * inv % q
                for k in range(len(b)):
                    a[len(a) - len(b) + k] = (a[len(a) - len(b) + k] - c * b[k]) % q
            a.pop()
            trim(a)
            if not a: break
        a, b = b, a
    return len(a) - 1
fder = [(i * c) % q for i, c in enumerate(fcoef)][1:]
squarefree = poly_gcd_deg(fcoef, fder) == 0
assert squarefree, "f is not squarefree: the curve is singular"

# ---------- F_{q^m} arithmetic (polynomials mod an irreducible of degree m) ----------
def poly_mulmod(a, b, mod, m):
    # a, b: lists of length m (coeffs), mod: monic irreducible of degree m (list length m+1)
    res = [0] * (2 * m - 1)
    for i, ai in enumerate(a):
        if ai:
            for j, bj in enumerate(b):
                if bj:
                    res[i + j] = (res[i + j] + ai * bj) % q
    for d in range(2 * m - 2, m - 1, -1):
        c = res[d]
        if c:
            for k in range(m + 1):
                res[d - m + k] = (res[d - m + k] - c * mod[k]) % q
    return res[:m]

def is_irreducible(mod, m):
    # brute force: no monic factor of degree 1..m//2
    def polymod(a, b):
        a = a[:]
        while len(a) >= len(b):
            if a[-1]:
                c = a[-1]
                for k in range(len(b)):
                    a[len(a) - len(b) + k] = (a[len(a) - len(b) + k] - c * b[k]) % q
            a.pop()
        return a
    for d in range(1, m // 2 + 1):
        for tail in itertools.product(range(q), repeat=d):
            b = list(tail) + [1]
            r = polymod(mod, b)
            if all(x == 0 for x in r):
                return False
    return True

def find_irreducible(m):
    if m == 1:
        return [0, 1]
    for tail in itertools.product(range(q), repeat=m):
        mod = list(tail) + [1]
        if mod[0] != 0 and is_irreducible(mod, m):
            return mod
    raise RuntimeError

def count_points(m):
    mod = find_irreducible(m)
    e = (q ** m - 1) // 2
    N = 1  # point at infinity
    for coeffs in itertools.product(range(q), repeat=m):
        x = list(coeffs)
        # evaluate f(x) by Horner in F_{q^m}
        val = [0] * m
        for c in reversed(fcoef):
            val = poly_mulmod(val, x, mod, m)
            val[0] = (val[0] + c) % q
        if all(v == 0 for v in val):
            N += 1
            continue
        # quadratic character: val^e
        r = [1] + [0] * (m - 1); base = val; ee = e
        while ee:
            if ee & 1: r = poly_mulmod(r, base, mod, m)
            base = poly_mulmod(base, base, mod, m); ee >>= 1
        one = [1] + [0] * (m - 1)
        N += 2 if r == one else 0
    return N

Ns = [count_points(m) for m in range(1, 2 * g + 1)]
# zeta polynomial: Z(T) = exp(sum N_m T^m/m) = P(T)/((1-T)(1-qT)); get P mod T^{2g+1}
from fractions import Fraction
L = [Fraction(0)] + [Fraction(Ns[m - 1], m) for m in range(1, 2 * g + 1)]  # log Z coefficients
# exp of power series
Z = [Fraction(1)] + [Fraction(0)] * (2 * g)
for n in range(1, 2 * g + 1):
    Z[n] = sum(k * L[k] * Z[n - k] for k in range(1, n + 1)) / n
# multiply by (1-T)(1-qT) = 1 - (q+1)T + qT^2
P = [Fraction(0)] * (2 * g + 3)
for i, z in enumerate(Z):
    P[i] += z; P[i + 1] -= (q + 1) * z; P[i + 2] += q * z
P = [int(c) for c in P[:2 * g + 1]]
assert all(Fraction(c).denominator == 1 for c in P[:2 * g + 1])
fe_ok = all(P[2 * g - k] == q ** (g - k) * P[k] for k in range(g + 1))  # functional equation
roots = np.roots(P[::-1])  # P(T) = sum P[k] T^k ; roots in T; alpha = 1/T
alphas = 1 / roots
mods = np.abs(alphas)
thetas = np.angle(alphas)
# closed points by degree: N_m = sum_{d|m} d b_d
b = {}
for m in range(1, 2 * g + 1):
    s = sum(d * b[d] for d in b if m % d == 0)
    b[m] = (Ns[m - 1] - s) // m
# Fourier coefficients of the zero measure nu = sum delta_{theta_i}: nu^(m) = sum_i e^{i m theta_i} = q^{m/2} + q^{-m/2} - N_m q^{-m/2}
def nuhat_count(m):
    m = abs(m)
    if m == 0: return 2 * g
    return q ** (m / 2) + q ** (-m / 2) - Ns[m - 1] * q ** (-m / 2)
def nuhat(m):
    return float(np.sum(np.exp(1j * m * thetas)).real)
# check the counted coefficients against the zeros (m <= 2g, where N_m was counted)
chk = max(abs(nuhat_count(m) - nuhat(m)) for m in range(0, 2 * g + 1))
# explicit-formula check for a random h supported on [-D, D]
rng = np.random.default_rng(1)
D = g
h = rng.standard_normal(2 * D + 1)
zero_side = sum(abs(sum(h[n + D] * np.exp(1j * n * th) for n in range(-D, D + 1))) ** 2 for th in thetas)
count_side = sum(h[j] * h[k] * nuhat(j - k) for j in range(2 * D + 1) for k in range(2 * D + 1))
# floors
rows = []
for Dc in range(0, 2 * g + 2):
    T = np.array([[nuhat(j - k) for k in range(2 * Dc + 1)] for j in range(2 * Dc + 1)])
    ev = np.linalg.eigvalsh(T)
    # even/odd sectors
    n = 2 * Dc + 1
    Jrev = np.eye(n)[::-1]
    Pe = (np.eye(n) + Jrev) / 2; Po = (np.eye(n) - Jrev) / 2
    def sector_min(Pr):
        U, s, _ = np.linalg.svd(Pr); B = U[:, s > 0.5]
        return np.linalg.eigvalsh(B.T @ T @ B).min() if B.shape[1] else float('nan')
    rank = int(np.linalg.matrix_rank(T, tol=1e-7))
    rows.append(dict(D=Dc, a=Dc * math.log(q), lam=float(ev.min()), lam_even=float(sector_min(Pe)), lam_odd=float(sector_min(Po)), rank=rank, size=n))
# reconstruction at the horizon D = g: the null vector of the singular Toeplitz matrix is a polynomial whose roots are the zeros
Tg = np.array([[nuhat(j - k) for k in range(2 * g + 1)] for j in range(2 * g + 1)])
w, V = np.linalg.eigh(Tg)
c = V[:, 0]                       # eigenvector of the zero eigenvalue
zr = np.roots(c[::-1])            # polynomial sum c_n z^n, n = 0..2g
rec_thetas = np.sort(np.angle(zr)); true_thetas = np.sort(thetas)
reconstruction = dict(null_eigenvalue=float(w[0]), root_moduli=[float(abs(z)) for z in zr],
                      max_angle_error=float(np.max(np.abs(rec_thetas - true_thetas))))
# a fake count sequence (one zero pair off the circle, |alpha| = q^{1/2} r with r > 1): where does the cone catch it?
def fake_catch(r):
    al = [math.sqrt(q) * r, math.sqrt(q) / r] + [math.sqrt(q) * np.exp(1j * t) for t in thetas[2:]]  # keep the other zeros
    def nh(m):
        return float(np.sum([(a / math.sqrt(q)) ** abs(m) for a in al]).real)
    for Dc in range(0, 12):
        T = np.array([[nh(j - k) for k in range(2 * Dc + 1)] for j in range(2 * Dc + 1)])
        if np.linalg.eigvalsh(T).min() < -1e-9:
            return Dc
    return None
fakes = {str(r): fake_catch(r) for r in (1.05, 1.1, 1.2, 1.5, 2.0)}
out = dict(q=q, f=fcoef, g=g, squarefree=squarefree, reconstruction=reconstruction, fake_catch_D=fakes, N=Ns, P=P, functional_equation=fe_ok, moduli=[float(x) for x in mods], sqrtq=math.sqrt(q),
           thetas=[float(t) for t in thetas], closed_points=b, nuhat_check=float(chk), explicit_formula_check=float(abs(zero_side - count_side)),
           floors=rows, horizon_a=g * math.log(q), determination_a=g * math.log(q) / 2)
print(json.dumps(out, indent=1))
