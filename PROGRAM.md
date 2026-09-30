# Program: from the tail law to the why

Agreed 2026-09-30. The rate of decay of the floor is pi times the effective height of the zeros that carry it,
Phi' = pi T_eff (tail law, paper Conjecture 7.7). Conjecture A (which implies RH) is "T_eff is a bounded multiple of
the horizon". Three items, in order; this file is updated at the end of every round so that the work does not drift.

## 1. Edge law => tail law, as a proposition  [DONE 2026-09-30: paper Proposition 7.8]
Under RH and the edge law with two derivatives, lim T sum_{|gamma|>T} |F(gamma)|^2 = 2C^2/pi. Proof: edge asymptotics
of the transform (splitting at delta = 1/t), then the sum over zeros through the Riemann-von Mangoldt formula with
Littlewood's bound S_1(t) = O(log t/(log log t)^2). With the boundary law lambda' = -2C^2 this is Conjecture 7.7.
Still open inside item 1: the edge law as an asymptotic equality with derivatives (only the upper bound is a theorem,
under (H_inf)); the boundary law itself (Conjecture 8.1).

## 2. The sum rule  [DONE 2026-09-30: Proposition 7.9, Lemma 7.11, Lemma 7.12, Corollary 7.13]
(a) DONE: the edge overlap is 2C^2/(pi T) J(L(1/T)) with J = 1 + O(1/L^2): the kernel sin(s+v)/(s+v) is symmetric in the
two distances and the log ratio antisymmetric, so the 1/L term vanishes; J = 1.020 ... 1.003 for L = 3 ... 8
(rh_sum_rule.py jl). This is why the tail law is sharp at finite T.
(b) DONE (Lemma 7.12): the regular part of K_inf f is continuous across the edge; inside, the wall's +(1/2)log(1/delta) f
and the log-Laplacian's -C L^{1/2} cancel, outside -C L^{1/2} stands alone, same regular remainder (the difference is
-C(beta + log 2a) L^{-1/2} -> 0). Hence g_reg(a) = (K f)(a-) = 2P sinh(a/2), h = 2P s + g_reg is continuous at the edge,
and polar + smooth overlap = <f, h_{>T}> = O(C^2/(T L)). Corollary 7.13: tail energy = 2C^2/(pi T)(1 + O(1/L)).
Loose end (not blocking): the coefficient of the O(1/L) term; the data say it is small.
(c) NEW 2026-09-30, reframes 2 and 3: the tail law is ANALYTIC. On the prime-free critical points (form W_inf - 2P^2,
no zeros, floors of either sign) T (1/2pi) int_{|t|>T} |F|^2 (Psi_inf - lambda) / (-lambda'/pi) = 1 +- 0.02 from three
horizons on (Computation 7.14). And the profile is NOT in its asymptotic edge regime at delta = 1/T for T = 5-20 T*
(f^2 L / C^2 = 1.5-3.2 at delta = 1/(e T*)); Lemma 7.11 + 7.12 explain only the far field. So the finite-T sum rule is
an exact analytic identity still to be found: a Pohozaev-type identity for the eigen-equation (the edge pairing of f
with the dilation generator x f' gives -lambda' by the boundary law; the same pairing cut at height T gives the tail).
NEXT for item 2: derive it. Start from tail(T) = 2P <f, s_{>T}> - <f_T, g_out> and the dilation identity
a lambda' = -(D_inf + 2P^2 + 2PM), D_inf = (1/2pi) int |F|^2 t dPsi_inf/dt, and look for the identity
T tail(T) = -lambda'/pi + (terms that vanish like the prime-free data say, i.e. by 3 horizons).
Status: the exact identity is paper Proposition 7.9 (tail identity); rh_sum_rule.py verifies it on the FEM minimizers
to 1-2% at 5-20 horizons (Computation 7.10; beyond, the y-grid under-resolves the leak's oscillation of period 2 pi/T). The edge overlap alone is the law to +-10%, and the smooth overlap plus
the polar piece cancel to +-10% of the law. Remaining: (a) the asymptotic expansion of the edge overlap in 1/L (leading
term log-free; the naive first-order estimate 0.79/L is too large, the data allow at most 0.2/L); (b) the cancellation of
the smooth overlap against the polar piece, analytically. The polar coefficient in the eigen-equation is 2P sinh(y/2),
not 4P (the gradient); with 4P the identity fails by 30%.
Route: the Euler-Lagrange equation in Fourier form. On the window, K f = lambda f + 4P sinh(x/2); outside it, K f =: g_out
is unconstrained (K the operator with symbol Psi; g_out contains the archimedean kernel acting on the edge, which
behaves like -2C L(delta)^{1/2} just outside the edge, and the reflected primes, smooth). Hence
  (Psi - lambda) F = G + 4P Sigma      (G = transform of g_out, Sigma = transform of sinh(x/2) 1_[-a,a]),
and, since f and g_out have disjoint supports, for every T the exact identity
  int_{|t|>T} |F|^2 (Psi - lambda) dt = 4P int_{|t|>T} conj(F) Sigma dt - 2 pi < f_T , g_out >,
where f_T is the low-pass of f at T. The leak of f_T outside the edge is ~ C L(1/T)^{-1/2} (pi/2 - Si(T delta))/pi and
g_out ~ -2C L(delta)^{1/2}: their overlap is -2C^2/(pi T) with L^{-1/2} L^{1/2} = 1, no logarithm, and
-2 pi times it is 4C^2/T, which is the tail law (2 pi int |F|^2 rho = 4C^2/T). The O(C L^{-1/2}/T) pieces (smooth part
of g_out against the leak, and the polar tail) must cancel each other; the first-order 1/L correction from the
overlap is not yet derived correctly (a naive estimate gives 1 + 0.79/L, five times the measured 0.07 at L = 2).
Next step: compute g_out for the FEM minimizer (archimedean convolution with J_Gamma plus reflected primes) and verify
the exact identity and the size of each piece at T = 5..100 horizons; then derive the boundary overlap with a smooth
cutoff to get the correction coefficient. If it comes out as measured, the sum rule is this identity.
The law holds at 1% from five horizons on, where the proof's O(1/log T) correction is 7% (density integral at a = 0.5:
1.068, 1.060, 1.032, 1.020, 1.012, 1.006 of the law at 5-100 horizons; zero sums 1.00 +- 0.02). The zero sum is
sharper than the density integral. Target: an exact identity for T int_{|t|>T} |F|^2 rho dt, or for the zero sum, valid
for the minimizer at finite T. Candidates: the Wiener-Hopf structure of the Euler-Lagrange equation ((Psi - lambda) F =
G + polar with G the transform of a function supported outside the window); the dilation identity applied to a
truncated form; the explicit formula for the high-pass part of f (the fluctuation that cancels the log correction is the
prime sum on frequencies above T). Whatever cancels the logarithms is prime-side, so this is likely the statement we
are looking for.

## 3. The onset at a single support  [REFORMULATED 2026-09-30: Lemma 7.16, Proposition 7.17; NEXT]
Exhaustion height: 2N(T) = 2aT/pi exactly at T = e T* (Lemma 7.16). Data: medians of the floor's mass 2.39 ... 2.74 T*
(rising to e), theta(e T*) = 0.43 ... 0.50. Proposition 7.17: under RH, (S_kappa) [zeros above kappa T* sample the tail
of F with constant M] + (E_kappa) [edge overlap holds at kappa T* within M'] + boundary law => Conjecture A with
c = pi kappa M M'. Measured: M ~ 1, M' ~ 1, kappa = e -> bound pi e = 8.5 against the observed 3.3-3.9 = pi e theta(e T*).
Remaining: prove (S_e) (a sampling inequality on the half-line above e T* for high-pass window functions with a log
edge; density (2a+1)/2pi vs Nyquist a/pi; needs no gaps wider than pi/a) and (E_e) (the minimizer's edge layer is no
narrower than 1/(e T*): it hides from no zero above the exhaustion height). To prove RH by A, both must eventually be
replaced by their prime-side forms.
Measured 2026-09-30 (Computation 7.14): (E_kappa) is not about the profile's asymptotic regime (the profile is 1.5-3x
above the edge law at delta = 1/(e T*)); it holds because the tail law is analytic (item 2c), from ~3 T* on, at 1%.
(S_kappa): M(kappa) = 1.01-1.09 for kappa >= e at a = 0.4, 0.45, 0.5 (8.9, 3.6, 1.8 at kappa = 1; 0.7-1.5 at 1.5-2).
So ALL the arithmetic of Conjecture A is in (S_e): the zeros above e T* sample the analytic tail of F at their mean
density, constant ~1. NEXT for item 3: the analytic form of M through the explicit formula for the high-pass part of
f (the fluctuation is the prime sum on frequencies above kappa T*, whose only structure between entries is the
alignment of Computation 7.18), and the gap condition it needs (no gap wider than pi/a above e T*).
A sampling-type lower bound: for the minimizer at one support a, sum_{gamma > kappa T*} |F(gamma)|^2 >= c C^2/T* with
explicit kappa, c. Inputs: density of zeros above the horizon exceeds a/pi; no gaps wider than pi/a above a few horizons
(true on average from T*, for the largest gaps from about 3T*); the minimizer's transform above kappa T* is
edge-dominated (observed, no argument yet). Under RH this is Conjecture A at that support. To prove RH by A, item 3
must eventually be replaced by its prime-side form: the dilation derivative of the ground state controlled by T* lambda
using only the Euler product and the archimedean symbol.

## Not to drift into
Refitting the decay constant; more grids; more supports for the tail law. The data are sufficient; the gap is analytic.
