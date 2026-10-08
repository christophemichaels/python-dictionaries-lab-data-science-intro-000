"""Section 9B: direct prime-zero comparison.  Observable: e(x) = psi(x) - x (von Mangoldt staircase with prime powers)
against the finite reconstruction e_Z(x) = -sum_{rho in Z} x^rho/rho - log(2 pi) - (1/2) log(1 - x^-2) from a listed
conjugation-symmetric multiset Z of computed zeros, in the window norm <f,g>_X = X^-2 int_X^{aX} conj(f) g dx.
P = ||e||^2 (exact in the staircase weights), Z = ||e_Z||^2, r = e - e_Z, D = ||r||^2, I = 2 Re <e_Z, r>, with P = Z + I + D.
Finite mode Gram bound: E_Z = c^* G c <= lambda_Z W_Z.  Zero provenance is recorded by zeros_cached()."""
import math, json, hashlib, os, time
import numpy as np

# ------------------------------------------------------------------ zeros
def zeros_cached(n_pos, path, dps=30):
    """First n_pos zeros of zeta on the critical line from mpmath.zetazero (cached as strings), with provenance.
    mpmath.zetazero(n, info=True) returns the zero, the bracketing Gram points, the index inside the Rosser block and the
    block pattern; it locates zeros as sign changes of the Riemann-Siegel Z function.  It does not certify that no zero
    off the critical line exists below the reported height; the list is 'the listed critical-line zeros'."""
    import mpmath
    cache = {}
    if os.path.exists(path):
        cache = json.load(open(path))
    have = len(cache.get('ordinates', []))
    if have < n_pos or cache.get('dps') != dps:
        mpmath.mp.dps = dps
        ords = cache.get('ordinates', []) if cache.get('dps') == dps else []
        infos = cache.get('gram_info', []) if cache.get('dps') == dps else []
        t0 = time.time()
        for n in range(len(ords) + 1, n_pos + 1):
            z, gram, idx, pattern = mpmath.zetazero(n, info=True)
            ords.append(mpmath.nstr(z.imag, dps)); infos.append([int(gram[0]), int(gram[1]), int(idx), str(pattern)])
        cache = dict(source="mpmath.zetazero(n, info=True)", mpmath_version=mpmath.__version__, dps=dps, ordinates=ords, gram_info=infos,
                     real_part="1/2 by construction (critical-line zero finder)", multiplicity="each listed once; conjugates added at use",
                     completeness="not certified by this computation: the finder reports Gram-point brackets and Rosser-block patterns for the zeros it indexes; no census of possible off-line zeros is produced",
                     elapsed_s=round(time.time() - t0, 2))
        json.dump(cache, open(path, 'w'), indent=1)
    ords = np.array([float(s) for s in cache['ordinates'][:n_pos]])
    h = hashlib.sha256(("\n".join(cache['ordinates'][:n_pos])).encode()).hexdigest()[:16]
    return ords, h, cache

def zero_multiset(ords):
    """conjugation-symmetric multiset: rho = 1/2 + i gamma and its conjugate."""
    g = np.asarray(ords, dtype=float)
    return np.concatenate([0.5 + 1j*g, 0.5 - 1j*g])

# ------------------------------------------------------------------ primes
def mangoldt(limit):
    """Lambda(n) for n <= limit (prime powers included) and the list of jump locations."""
    s = np.ones(limit + 1, dtype=bool); s[:2] = False
    for p in range(2, int(limit**0.5) + 1):
        if s[p]: s[p*p::p] = False
    primes = np.nonzero(s)[0]
    Lam = np.zeros(limit + 1)
    for p in primes:
        q = int(p)
        while q <= limit:
            Lam[q] = math.log(p); q *= int(p)
    return Lam

class Staircase:
    def __init__(self, limit):
        self.limit = limit; self.Lam = mangoldt(limit); self.psi = np.cumsum(self.Lam)
        self.jumps = np.nonzero(self.Lam)[0]
    def intervals(self, X, aX):
        """intervals (b, c) between prime-power jumps covering [X, aX], with the constant value s = psi on each."""
        pts = [X] + [float(j) for j in self.jumps if X < j < aX] + [aX]
        out = []
        for b, c in zip(pts[:-1], pts[1:]):
            s = float(self.psi[int(math.floor(b))]) if b == X else float(self.psi[int(b)])   # value just right of b
            if b == X:
                s = float(self.psi[int(math.floor(X))])
            out.append((b, c, s))
        return out
    def P_exact_weights(self, X, a):
        """P(X) = X^-2 sum over intervals of Delta [ (s - (b+c)/2)^2 + Delta^2/12 ]."""
        tot = 0.0
        for b, c, s in self.intervals(X, a*X):
            D = c - b; tot += D*((s - (b + c)/2)**2 + D*D/12)
        return tot/(X*X)
    def psi0(self, x):
        """midpoint convention at jumps (for pointwise plots only)."""
        x = np.asarray(x, dtype=float); k = np.floor(x).astype(int); v = self.psi[k].copy()
        at = (x == k) & (self.Lam[k] > 0); v[at] -= 0.5*self.Lam[k[at]]
        return v

# ------------------------------------------------------------------ reconstruction
LOG2PI = math.log(2*math.pi)
def e_Z(x, rhos):
    """-sum_rho x^rho/rho - log 2pi - (1/2) log(1 - x^-2), real for a conjugation-symmetric multiset."""
    x = np.asarray(x, dtype=float); lx = np.log(x)
    z = np.exp(np.outer(lx, rhos))/rhos[None, :]
    return -np.real(z.sum(axis=1)) - LOG2PI - 0.5*np.log1p(-1.0/(x*x))

def gauss_nodes(n=12):
    return np.polynomial.legendre.leggauss(n)

def window_integrals(stair, X, a, rhos, nodes=12, max_phase=0.5):
    """Z, cross <e_Z, e>, D, I and the quadrature check of P, by Gauss-Legendre on each jump interval, subdivided in
    log x so that gamma_max * (log span per piece) <= max_phase.  Returns dict with all pieces."""
    gx, gw = gauss_nodes(nodes); gmax = float(np.max(np.abs(np.imag(rhos)))) if len(rhos) else 1.0
    Zq = 0.0; XE = 0.0; Dq = 0.0; Iq = 0.0; Pq = 0.0
    for b, c, s in stair.intervals(X, a*X):
        span = math.log(c/b); m = max(1, int(math.ceil(gmax*span/max_phase)))
        edges = b*np.exp(np.linspace(0.0, span, m + 1))
        for lo, hi in zip(edges[:-1], edges[1:]):
            x = (hi - lo)/2*gx + (hi + lo)/2; w = (hi - lo)/2*gw
            ez = e_Z(x, rhos); e = s - x; r = e - ez
            Zq += float(np.sum(w*ez*ez)); XE += float(np.sum(w*ez*e)); Dq += float(np.sum(w*r*r)); Iq += float(np.sum(w*2*ez*r)); Pq += float(np.sum(w*e*e))
    X2 = X*X
    return dict(Z=Zq/X2, cross_eZ_e=XE/X2, D=Dq/X2, I=Iq/X2, P_quadrature=Pq/X2)

def mode_gram(X, a, rhos):
    """E_Z = c^* G c with c_rho = X^(rho-1/2)/rho, G = (exp((1+conj rho+sigma) l) - 1)/(1+conj rho+sigma); the weighted bound."""
    ell = math.log(a); c = np.exp((rhos - 0.5)*math.log(X))/rhos
    S = 1.0 + np.conj(rhos)[:, None] + rhos[None, :]
    G = np.expm1(S*ell)/S
    E = float(np.real(np.conj(c) @ G @ c))
    diag = float(np.real(np.sum(np.abs(c)**2*np.diag(G))))
    Dw = np.log(2 + np.abs(np.imag(rhos)))
    Gt = G/np.sqrt(np.outer(Dw, Dw)); lam = float(np.max(np.linalg.eigvalsh((Gt + Gt.conj().T)/2)))
    W = float(np.sum(np.abs(c)**2*Dw))
    return dict(E_modes=E, E_diag=diag, E_offdiag=E - diag, lambda_finite=lam, W_weighted=W, bound=lam*W, utilization=E/(lam*W) if lam*W > 0 else float('nan'))

def zero_part_quadrature(X, a, rhos, nodes=12, max_phase=0.5):
    """X^-2 int |sum x^rho/rho|^2 dx by quadrature, to cross-check the Gram expression."""
    gx, gw = gauss_nodes(nodes); gmax = float(np.max(np.abs(np.imag(rhos)))); span = math.log(a); m = max(1, int(math.ceil(gmax*span/max_phase)))
    edges = X*np.exp(np.linspace(0.0, span, m + 1)); tot = 0.0
    for lo, hi in zip(edges[:-1], edges[1:]):
        x = (hi - lo)/2*gx + (hi + lo)/2; w = (hi - lo)/2*gw; lx = np.log(x)
        z = np.real((np.exp(np.outer(lx, rhos))/rhos[None, :]).sum(axis=1)); tot += float(np.sum(w*z*z))
    return tot/(X*X)

def campaign(stair, ords_all, Xs, counts, a=6, nodes=12, conv_nodes=24, zero_hash="", max_height_all=None):
    rows = []
    for X in Xs:
        P = stair.P_exact_weights(X, a)
        for n in counts:
            rhos = zero_multiset(ords_all[:n]); t0 = time.time()
            q = window_integrals(stair, X, a, rhos, nodes=nodes); q2 = window_integrals(stair, X, a, rhos, nodes=conv_nodes, max_phase=0.25)
            g = mode_gram(X, a, rhos)
            ident = P - (q['Z'] + q['I'] + q['D'])
            rows.append(dict(X=X, window_ratio=a, zero_list_hash=zero_hash, positive_zero_count=n, actual_max_height=float(ords_all[n-1]),
                             completeness_status="reconstruction from the listed critical-line zeros (mpmath.zetazero); completeness not certified here",
                             P_prime=P, P_quadrature=q['P_quadrature'], Z_reconstruction=q['Z'], I_cross=q['I'], D_residual=q['D'], D_over_P=q['D']/P if P > 0 else float('nan'),
                             identity_residual=ident, bound1_ok=bool(abs(math.sqrt(P) - math.sqrt(q['Z'])) <= math.sqrt(q['D']) + 1e-12), bound2_ok=bool(abs(q['I']) <= 2*math.sqrt(q['Z']*q['D']) + 1e-12),
                             E_modes=g['E_modes'], E_diag=g['E_diag'], E_offdiag=g['E_offdiag'], W_weighted_diagonal=g['W_weighted'], lambda_finite=g['lambda_finite'], bound_utilization=g['utilization'], gram_bound_ok=bool(g['E_modes'] <= g['bound']*(1 + 1e-12)),
                             precision="float64", quadrature_Z_change=abs(q2['Z'] - q['Z']), quadrature_D_change=abs(q2['D'] - q['D']), quadrature_P_change=abs(q['P_quadrature'] - P), elapsed_s=round(time.time() - t0, 3)))
    return rows

if __name__ == "__main__":
    import sys
    path = sys.argv[1] if len(sys.argv) > 1 else "zeros_test.json"
    ords, zh, cache = zeros_cached(50, path); print("zeros:", len(ords), "max height", ords[-1], "hash", zh, "elapsed", cache.get('elapsed_s'))
    st = Staircase(6*256 + 10)
    for X in (2, 16, 256):
        P = st.P_exact_weights(X, 6); rh = zero_multiset(ords[:50])
        q = window_integrals(st, X, 6, rh); g = mode_gram(X, 6, rh); zq = zero_part_quadrature(X, 6, rh)
        print(f"X={X}: P={P:.6f} (quad {q['P_quadrature']:.6f}) Z={q['Z']:.6f} I={q['I']:+.6f} D={q['D']:.6f} sum={q['Z']+q['I']+q['D']:.6f}; E_modes={g['E_modes']:.6f} quad {zq:.6f}; lambda={g['lambda_finite']:.4f} W={g['W_weighted']:.6f} util={g['utilization']:.4f}")
