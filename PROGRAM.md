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

## 2. The sum rule  [IDENTITY PROVED AND VERIFIED; edge-overlap expansion done (Lemma 7.11); one local statement remains]
(a) DONE: the edge overlap is 2C^2/(pi T) J(L(1/T)) with J = 1 + O(1/L^2): the kernel sin(s+v)/(s+v) is symmetric in the
two distances and the log ratio antisymmetric, so the 1/L term vanishes; J = 1.020 ... 1.003 for L = 3 ... 8
(rh_sum_rule.py jl). This is why the tail law is sharp at finite T.
(b) OPEN, local: the smooth overlap cancels the polar piece iff g_reg(a) = 2P sinh(a/2), i.e. the regular part of K f
is continuous across the edge (the singularity -/+ C L^{1/2} carried by the outside). Supported to ~20%. A property of
the log-Laplacian kernel on the (log)^{-1/2} edge; provable by a local computation.
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

## 3. The onset at a single support  [REFORMULATED 2026-09-30: Lemma 7.13, Proposition 7.14]
Exhaustion height: 2N(T) = 2aT/pi exactly at T = e T* (Lemma 7.13). Data: medians of the floor's mass 2.39 ... 2.74 T*
(rising to e), theta(e T*) = 0.43 ... 0.50. Proposition 7.14: under RH, (S_kappa) [zeros above kappa T* sample the tail
of F with constant M] + (E_kappa) [edge overlap holds at kappa T* within M'] + boundary law => Conjecture A with
c = pi kappa M M'. Measured: M ~ 1, M' ~ 1, kappa = e -> bound pi e = 8.5 against the observed 3.3-3.9 = pi e theta(e T*).
Remaining: prove (S_e) (a sampling inequality on the half-line above e T* for high-pass window functions with a log
edge; density (2a+1)/2pi vs Nyquist a/pi; needs no gaps wider than pi/a) and (E_e) (the minimizer's edge layer is no
narrower than 1/(e T*): it hides from no zero above the exhaustion height). To prove RH by A, both must eventually be
replaced by their prime-side forms.
A sampling-type lower bound: for the minimizer at one support a, sum_{gamma > kappa T*} |F(gamma)|^2 >= c C^2/T* with
explicit kappa, c. Inputs: density of zeros above the horizon exceeds a/pi; no gaps wider than pi/a above a few horizons
(true on average from T*, for the largest gaps from about 3T*); the minimizer's transform above kappa T* is
edge-dominated (observed, no argument yet). Under RH this is Conjecture A at that support. To prove RH by A, item 3
must eventually be replaced by its prime-side form: the dilation derivative of the ground state controlled by T* lambda
using only the Euler product and the archimedean symbol.

## Not to drift into
Refitting the decay constant; more grids; more supports for the tail law. The data are sufficient; the gap is analytic.
