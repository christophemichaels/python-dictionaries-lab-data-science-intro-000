"""
Parallel ball-arithmetic engine for the odd-sector Weil form Q_a on [-a,a] (python-flint / arb), for K up to ~800.

Same form, basis, quadrature and normalization as rh_weil_arb.py:
  Q_a(f) = W_inf(g) - sum_n Lambda(n) n^{-1/2} (g(log n)+g(-log n)) - 2 (int f sinh(x/2))^2,   g = f * f~,
  f(y) = a^{-1/2} sum_i c_i N_i P_{2i+1}(y/a),  N_i = sqrt((4i+3)/2),
  g_ij(s) = int_{s-1}^{1} phi_i(u) phi_j(u-s) du by a Gauss-Legendre rule with M = 2K+8 nodes (exact).

What is different from rh_weil_arb.py.
 * The archimedean integral is written in the shift variable s = x/a on [0,2],
     A(a) = int_0^2 2a [ I e^{-2as} - G(s) e^{-as/2} ] / (1-e^{-2as}) ds - log(1-e^{-4a}) I - (gamma + log pi) I,
   so the overlap matrices G(s_i) at the outer nodes do not depend on a: one pass gives A(a), A'(a) (derivative of
   the weights) and A(a +- delta).  The outer rule has Mx = 2K + 600 nodes (the polynomial G has degree 4K in s and
   the weight is analytic within |Im s| < pi/a; the truncation error is ~4^{-(2Mx-4K)}).
 * One Legendre table per outer node serves both factors of the overlap (the inner nodes of P(u-s) are the mirror
   images of those of P(u)); the square roots of the inner weights ride in the recurrence's starting values.
 * Every matrix is a weighted sum of overlap matrices G(s): the outer nodes (weights for Q, Q', Q(a +- delta)) and the
   prime shifts log n / a (and log n / a +- h for d/da).  The tasks are shared out among worker processes; the partial
   sums travel through binary files (mantissa, exponent).  G is symmetric: only the upper block triangle is multiplied.
 * The lowest eigenpair comes from one inverse-iteration step with shift 0 (floating point at the working precision;
   lambda_0 is hundreds of digits below lambda_1), the Rayleigh quotient, and lambda'(a) = <c, Q'(a) c> (Hellmann-Feynman);
   a deflated second step gives lambda_1.  The nested basis makes the K'-mode matrices leading blocks of the K-mode one.

Usage: python3 rh_weil_fast.py a K prec out.json [delta] [subK1,subK2,...] [nproc] [Mx]
"""
import sys, os, json, time, math, struct
from multiprocessing import Pool
from flint import arb, arb_mat, fmpz, ctx

def vonmangoldt(n):
    for p in range(2, n + 1):
        if n % p == 0:
            m = n
            while m % p == 0: m //= p
            return p if m == 1 else 0
    return 0

def gl_rule(M):
    xs, ws = [], []
    for k in range(M):
        x, w = arb.legendre_p_root(M, k, weight=True); xs.append(x); ws.append(w)
    return xs, ws

def legendre_table(us, kmax, scale=None, reset=64):
    """vals[k] = list over nodes of P_k(u_m) (times scale[m] if given), k = 0..kmax: three-term recurrence per node
    with precomputed coefficients, midpoint reset every `reset` degrees (the ball radius doubles per degree)."""
    c1 = [None, None] + [arb(2*k - 1)/k for k in range(2, kmax + 1)]
    c2 = [None, None] + [arb(k - 1)/k for k in range(2, kmax + 1)]
    cols = []
    for m, u in enumerate(us):
        p0 = scale[m] if scale is not None else arb(1)
        p1 = u*p0
        col = [p0, p1]
        for k in range(2, kmax + 1):
            p0, p1 = p1, u*p1*c1[k] - p0*c2[k]
            if k % reset == 0: p0 = arb(p0.mid()); p1 = arb(p1.mid())
            col.append(p1)
        cols.append(col)
    return [[col[k] for col in cols] for k in range(kmax + 1)]

class Engine:
    def __init__(self, K, prec):
        ctx.prec = prec
        self.K, self.prec = K, prec
        self.M = 2*K + 8
        self.xs, self.ws = gl_rule(self.M)
        assert all((self.xs[m] + self.xs[self.M - 1 - m]).contains(0) for m in range(self.M))
        self.sqw = [w.sqrt() for w in self.ws]
        self.split = 4 if K >= 64 else 1
        self.c1 = [None, None] + [arb(2*k - 1)/k for k in range(2, 2*K)]
        self.c2 = [None, None] + [arb(k - 1)/k for k in range(2, 2*K)]

    def gmat_raw(self, s):
        """G_raw(s)_{ij} = int_{s-1}^{1} P_{2i+1}(u) P_{2j+1}(u-s) du (without the N_i N_j), s in (0,2).
        One table serves both factors: with u_m = s/2 + (1-s/2) x_m and the symmetric rule x_{M-1-m} = -x_m,
        u_m - s = -u_{M-1-m}, so P_{2j+1}(u_m - s) = -P_{2j+1}(u_{M-1-m}); the square roots of the weights are
        folded into the recurrence's starting values."""
        K, M = self.K, self.M
        half = (2 - s)/2; mid = s/2
        us = [mid + half*x for x in self.xs]
        sh = half.sqrt()
        c1, c2 = self.c1, self.c2
        cols = []
        for m, u in enumerate(us):
            p0 = self.sqw[m]*sh; p1 = u*p0
            col = [p1]
            for k in range(2, 2*K):
                p0, p1 = p1, u*p1*c1[k] - p0*c2[k]
                if k & 1: col.append(p1)
                if k % 64 == 0: p0 = arb(p0.mid()); p1 = arb(p1.mid())
            cols.append(col)
        p = self.split; bounds = [(K*r)//p for r in range(p + 1)]
        blocks = [arb_mat(M, bounds[r + 1] - bounds[r], [v for col in cols for v in col[bounds[r]:bounds[r + 1]]]) for r in range(p)]
        rblocks = [arb_mat(M, bounds[r + 1] - bounds[r], [v for col in reversed(cols) for v in col[bounds[r]:bounds[r + 1]]]) for r in range(p)]
        G = arb_mat(K, K)
        for r in range(p):
            for q in range(r, p):
                Grq = -(blocks[r].transpose()*rblocks[q])                     # G is symmetric: only r <= q is computed
                for i in range(bounds[r + 1] - bounds[r]):
                    for j in range(bounds[q + 1] - bounds[q]):
                        G[bounds[r] + i, bounds[q] + j] = Grq[i, j]
                        if q != r: G[bounds[q] + j, bounds[r] + i] = Grq[i, j]
        return G

def weights(s, a):
    """(w_G, w_I): the a-dependent weights of G(s) and of I in the archimedean integral, and their a-derivatives."""
    e2 = (-2*a*s).exp(); eh = (-a*s/2).exp(); den = 1 - e2
    wG = 2*a*eh/den; wI = 2*a*e2/den
    corr = 2*s*e2/den
    return wG, wI, wG*(1/a - s/2 - corr), wI*(1/a - 2*s - corr)

# ---------------------------------------------------------------- binary transfer of arb matrices (midpoints)
def write_mats(path, mats):
    with open(path, "wb") as fh:
        for Mt in mats:
            r, c = Mt.nrows(), Mt.ncols()
            fh.write(struct.pack("<ii", r, c))
            for i in range(r):
                for j in range(c):
                    m, e = Mt[i, j].mid().man_exp()
                    mb = int(m).to_bytes((abs(int(m)).bit_length() + 8)//8, "little", signed=True)
                    fh.write(struct.pack("<qi", int(e), len(mb))); fh.write(mb)

def read_mats(path, count):
    out = []
    with open(path, "rb") as fh:
        for _ in range(count):
            r, c = struct.unpack("<ii", fh.read(8))
            flat = []
            two = arb(2)
            for _ in range(r*c):
                e, ln = struct.unpack("<qi", fh.read(12))
                m = int.from_bytes(fh.read(ln), "little", signed=True)
                flat.append(arb(fmpz(m))*two**e if m else arb(0))
            out.append(arb_mat(r, c, flat))
    return out

# ---------------------------------------------------------------- workers
_W = {}
def _init(K, prec):
    ctx.prec = prec
    _W["E"] = Engine(K, prec)

def _chunk(args):
    tasks, labels, path = args                     # tasks: list of (s_string, {label: coefficient_string})
    E = _W["E"]; K = E.K
    acc = {l: arb_mat(K, K) for l in labels}
    t0 = time.time()
    for n, (s_s, coefs) in enumerate(tasks):
        G = E.gmat_raw(arb(s_s))
        for l, c_s in coefs.items():
            acc[l] += G*arb(c_s)
        if n % 50 == 0: print(f"    worker {os.getpid()}: task {n+1}/{len(tasks)} ({time.time()-t0:.0f}s)", flush=True)
    write_mats(path, [acc[l] for l in labels])
    return path

# ---------------------------------------------------------------- assembly
def assemble(a_s, K, prec, delta=None, nproc=4, Mx=None, log=print, tmpdir="."):
    ctx.prec = prec
    a = arb(a_s)
    Mx = Mx if Mx else 2*K + 600
    digits = int(prec*0.30103) + 6
    def S(x): return x.mid().str(digits, radius=False)
    avals = {"Q": a, "dQ": a}
    if delta: avals["Qm"] = a - arb(delta); avals["Qp"] = a + arb(delta)
    labels = list(avals)
    t0 = time.time()
    xs, ws = gl_rule(Mx)
    tasks = []; diag = {l: arb(0) for l in labels}
    for x, w in zip(xs, ws):                                    # archimedean integral, s = 1 + x on [0,2]
        s = 1 + x; coefs = {}
        for l in labels:
            wG, wI, dwG, dwI = weights(s, avals[l])
            if l == "dQ": coefs[l] = S(-w*dwG); diag[l] += w*dwI
            else: coefs[l] = S(-w*wG); diag[l] += w*wI
        tasks.append((S(s), coefs))
    h = arb(2)**(-(prec//3))
    entries = []; n = 2
    while math.log(n) < 2*float(a.mid()) + 0.02:
        p = vonmangoldt(n)
        if p: entries.append((n, p))
        n += 1
    used = []
    for n, p in entries:                                        # prime matrices G(log n / a) and d/da of them
        wgt = -2*arb(p).log()/arb(n).sqrt()
        for l in labels:
            aa = avals[l]; s = arb(n).log()/aa
            if s >= 2: continue
            if l == "dQ":
                fac = wgt*(-arb(n).log()/(aa*aa))/(2*h)
                tasks.append((S(s + h), {"dQ": S(fac)})); tasks.append((S(s - h), {"dQ": S(-fac)}))
            else:
                if l == "Q": used.append(n)
                tasks.append((S(s), {l: S(wgt)}))
    log(f"outer rule Mx={Mx}, entries {used}: {len(tasks)} tasks ({time.time()-t0:.0f}s)")
    chunks = [(tasks[c::nproc], labels, os.path.join(tmpdir, f"partial_{c}.bin")) for c in range(nproc)]
    t0 = time.time()
    acc = None
    with Pool(nproc, initializer=_init, initargs=(K, prec)) as pool:
        for path in pool.imap_unordered(_chunk, chunks):
            mats = read_mats(path, len(labels)); os.remove(path)
            if acc is None: acc = dict(zip(labels, mats))
            else:
                for l, Mt in zip(labels, mats): acc[l] += Mt
            log(f"  chunk done ({time.time()-t0:.0f}s)")
    log(f"quadrature and primes: {time.time()-t0:.0f}s")
    E = Engine(K, prec)
    N = [(arb(4*i + 3)/2).sqrt() for i in range(K)]
    def c0(aa): return -(1 - (-4*aa).exp()).log() - (arb.const_euler() + arb.pi().log())
    scal = {"Q": c0(a), "dQ": (c0(a + h) - c0(a - h))/(2*h)}
    if delta: scal["Qm"] = c0(avals["Qm"]); scal["Qp"] = c0(avals["Qp"])
    Q = {}
    for l in labels:
        Mt = arb_mat(K, K, [N[i]*N[j]*acc[l][i, j] for i in range(K) for j in range(K)])
        for i in range(K): Mt[i, i] = Mt[i, i] + diag[l] + scal[l]
        Q[l] = (Mt + Mt.transpose())/2
    t0 = time.time()
    for l in labels:
        S_, Sd = polar(E, avals[l])
        if l == "dQ": Q["dQ"] -= (S_*Sd.transpose() + Sd*S_.transpose())*2
        else: Q[l] -= S_*S_.transpose()*2
    log(f"polar: {time.time()-t0:.0f}s")
    return Q, Mx, used

def polar(E, a):
    K = E.K
    sh = [(a*x/2).sinh() for x in E.xs]; ch = [(a*x/2).cosh()*(x/2) for x in E.xs]
    T = legendre_table(E.xs, 2*K - 1, scale=E.ws)
    N = [(arb(4*i + 3)/2).sqrt() for i in range(K)]
    ra = a.sqrt()
    s = [sum((v*sv for v, sv in zip(T[2*i + 1], sh)), arb(0))*N[i] for i in range(K)]
    sd = [sum((v*cv for v, cv in zip(T[2*i + 1], ch)), arb(0))*N[i] for i in range(K)]
    S = arb_mat(K, 1, [ra*x for x in s])
    Sd = arb_mat(K, 1, [ra*sd[i] + s[i]*ra/(2*a) for i in range(K)])
    return S, Sd

def lowest(Q, Qd=None, steps=3):
    """Inverse iteration (floating point at the working precision): one step with shift 0, then steps-1 steps with the
    shift sigma = lambda(1 - 2^-30); Rayleigh quotients for lambda_0 and lambda'(a) = <c, Q' c>; a deflated step for lambda_1."""
    K = Q.nrows()
    b = arb_mat(K, 1, [arb(1) + arb(i)/(7*K) for i in range(K)])
    Qm = arb_mat(K, K, [arb(Q[i, j].mid()) for i in range(K) for j in range(K)])
    def normalize(v):
        nrm = sum((v[i, 0]*v[i, 0] for i in range(K)), arb(0)).sqrt()
        c = arb_mat(K, 1, [v[i, 0]/nrm for i in range(K)])
        return -c if c[K//2, 0] < 0 else c
    c = normalize(Qm.solve(b, algorithm="approx"))
    lam0 = (c.transpose()*Q*c)[0, 0]
    for _ in range(steps - 1):
        sigma = lam0.mid()*(1 - arb(2)**(-30))
        Qs = arb_mat(Qm)
        for i in range(K): Qs[i, i] = Qs[i, i] - sigma
        c = normalize(Qs.solve(c, algorithm="approx"))
        lam0 = (c.transpose()*Q*c)[0, 0]
    res = Q*c - c*lam0
    resn = sum((res[i, 0]*res[i, 0] for i in range(K)), arb(0)).sqrt()
    bc = (c.transpose()*b)[0, 0]
    v2 = Qm.solve(b - c*bc, algorithm="approx")
    v2 = v2 - c*(c.transpose()*v2)[0, 0]
    c2 = normalize(v2)
    lam1 = (c2.transpose()*Q*c2)[0, 0]
    lamd = (c.transpose()*Qd*c)[0, 0] if Qd is not None else None
    return lam0, lam1, c, lamd, resn

def subblock(Q, Kp):
    return arb_mat(Kp, Kp, [Q[i, j] for i in range(Kp) for j in range(Kp)])

def kappa(lam0, lamd, a):
    return float(-lamd/lam0/(2*arb.pi()*arb.pi()*(2*a).exp()))

if __name__ == "__main__":
    a_s, K, prec, out = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), sys.argv[4]
    delta = sys.argv[5] if len(sys.argv) > 5 and sys.argv[5] not in ("0", "-") else None
    subK = [int(x) for x in sys.argv[6].split(",")] if len(sys.argv) > 6 and sys.argv[6] not in ("-", "") else []
    nproc = int(sys.argv[7]) if len(sys.argv) > 7 else 4
    Mx = int(sys.argv[8]) if len(sys.argv) > 8 else None
    T0 = time.time()
    def log(msg): print(msg, flush=True)
    Q, Mx, used = assemble(a_s, K, prec, delta, nproc, Mx, log, tmpdir=os.path.dirname(os.path.abspath(out)))
    t0 = time.time()
    lam0, lam1, c, lamd, resn = lowest(Q["Q"], Q["dQ"])
    a = arb(a_s)
    log(f"a={a_s} K={K} prec={prec}: lambda_K = {lam0.str(20)}  a lambda' = {(a*lamd).str(15)}  Phi'/(pi T*) = {kappa(lam0, lamd, a):.5f}  lambda_1 = {lam1.str(8)}  residual {resn.str(3)}  ({time.time()-t0:.0f}s)")
    res = {"a": a_s, "K": K, "prec": prec, "lambda": lam0.mid().str(60, radius=False), "lambda_rad": lam0.rad().str(3),
           "a_lambda_prime": (a*lamd).mid().str(60, radius=False), "lambda1": lam1.mid().str(20, radius=False),
           "residual": resn.str(3), "Mx": Mx, "entries": used,
           "c": [c[i, 0].mid().str(int(prec*0.30103) + 2, radius=False) for i in range(K)], "sub": {}, "engine": "rh_weil_fast.py"}
    if delta:
        lm = lowest(Q["Qm"])[0]; lp = lowest(Q["Qp"])[0]
        res["lambda_am"] = lm.mid().str(60, radius=False); res["lambda_ap"] = lp.mid().str(60, radius=False)
        d = arb(delta); res["delta"] = delta
        res["a_lambda_prime_cdiff"] = (a*(lp - lm)/(2*d)).mid().str(30, radius=False)
        res["Phi_prime_cdiff"] = (-(lp.log() - lm.log())/(2*d)).mid().str(20, radius=False)
        log(f"  a-{delta}: lambda = {lm.str(15)}   a+{delta}: lambda = {lp.str(15)}")
        log(f"  central difference: a lambda' = {(a*(lp-lm)/(2*d)).str(12)}  Phi' = {(-(lp.log()-lm.log())/(2*d)).str(12)}  (HF: {(-lamd/lam0).str(12)})")
    for Kp in subK:
        if Kp >= K: continue
        l0, l1, _, ld, _ = lowest(subblock(Q["Q"], Kp), subblock(Q["dQ"], Kp))
        res["sub"][str(Kp)] = {"lambda": l0.mid().str(40, radius=False), "a_lambda_prime": (a*ld).mid().str(40, radius=False)}
        log(f"  K'={Kp}: lambda = {l0.str(15)}  a lambda' = {(a*ld).str(12)}  Phi'/(pi T*) = {kappa(l0, ld, a):.5f}")
    res["seconds"] = round(time.time() - T0)
    json.dump(res, open(out, "w"))
    log(f"wrote {out} ({res['seconds']}s)")
