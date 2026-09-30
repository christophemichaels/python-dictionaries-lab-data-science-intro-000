# Program: from the tail law to the why

Agreed 2026-09-30. The rate of decay of the floor is pi times the effective height of the zeros that carry it,
Phi' = pi T_eff (tail law, paper Conjecture 7.7). Conjecture A (which implies RH) is "T_eff is a bounded multiple of
the horizon". Three items, in order; this file is updated at the end of every round so that the work does not drift.

## 1. Edge law => tail law, as a proposition  [DONE 2026-09-30: paper Proposition 7.8]
Under RH and the edge law with two derivatives, lim T sum_{|gamma|>T} |F(gamma)|^2 = 2C^2/pi. Proof: edge asymptotics
of the transform (splitting at delta = 1/t), then the sum over zeros through the Riemann-von Mangoldt formula with
Littlewood's bound S_1(t) = O(log t/(log log t)^2). With the boundary law lambda' = -2C^2 this is Conjecture 7.7.
Still open inside item 1: the edge law as an asymptotic equality with derivatives (only the upper bound is a theorem,
under (H_inf)). The boundary law is now a theorem under that same hypothesis (paper Proposition 8.3, 2026-09-30): B_T ->
a C^2 by the overlap of the leak of the dilation generator with the edge force, J_1(L) = 1 + O(1/L^2) with no 1/L term (the
antisymmetry of Lemma 7.11 again), and Lemma 8.2 (the edge force to second order). So under the edge law with derivatives:
boundary law => tail law (Prop 7.8), and both are the flux identity 2|A_0|^2 = -lambda' of the Wiener-Hopf form.

## 2. The sum rule  [DONE 2026-09-30: Proposition 7.9, Lemma 7.11, Lemma 7.12, Corollary 7.13; sharp form Conjecture 7.19 verified, Computation 7.20]
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
DONE (Proposition 7.15, Corollary 7.16, Computation 7.17): the exact identity. Cutting the dilation integral at T,
integrating by parts and using the eigen-equation on |t| < T:
  D = tail(T) + (1/pi) T |F(T)|^2 (Psi(T)-lambda) + tailD(T) - 4P <s, ((xf)')_{>T}> + 2 B_T,
  B_T = <g_out, ((xf)')_T>,  D = 2 B_inf,  so  -a lambda' = 2 lim <g_out, ((xf)')_T>  and  B_inf = a C^2.
Verified on the prime-free FEM critical point at a = 0.5: D = 1.129084 (engine: 1.129075), B_T -> 0.551 vs D/2 = 0.5645,
identity to 1-2% (rh_cutoff_identity.py, data/cutoff_identity_a0.5_primefree.log). The tail law is the statement
B_inf - B_T = B_inf/(pi a T)(1 + o(1)) up to explicit oscillating terms: the Hadamard pairing converges at the rate set
by its own value.
CORRECTION: the polar high-pass term in the pairing identity is O(L^{-1/2}), five times the tail, and B_inf - B_T
cancels it at that order; the tail is the O(1/T) remainder. The pairing identity is the dilation identity in Fourier
form; it does not simplify the law.
THE SHARP FORM (Computation 7.18, Conjecture 7.19): the LOCAL envelope <t^2 (Psi - lambda) |F(t)|^2> over one
interference period equals -lambda' from ~5 horizons on: to 0.1% with no drift over two decades on the prime-free
critical points, +-3% with primes. Mechanism: Wiener-Hopf. Near an edge the eigen-equation is a half-line problem with
symbol sigma = Psi - lambda; its solution has transform e^{ita} Phi_+(t)/(t sigma_+(t)) with sigma = |sigma_+|^2, so
t^2 sigma |F|^2 = |Phi_+|^2 -> |A|^2 with NO logarithms (they all sit in sigma_+; the edge law C L^{-1/2} is the x-space
asymptotics of 1/sigma_+, which is why the profile can be 1.5-3x above its asymptotic law at scale 1/T while the tail
law holds there). |A_+|^2 + |A_-|^2 = -lambda' is the boundary law as a flux. The onset is NOT the last sign change of
the symbol (0.7, 2.2, 3.5, 33, >150 T* at a = 0.5, 0.6, 0.8, 1.0, 1.25): the oscillating part of the symbol is a bounded
perturbation whose effect averages out.
NEXT for item 2 (the last analytic piece): prove Conjecture 7.19. Wiener-Hopf factorization of sigma = Psi_inf - lambda
(order-zero symbol with logarithmic growth; index from its sign change near t_lambda), the two edges and the polar
rank-one term as regular perturbations, the prime cosines as a bounded perturbation averaging out; then the tail law is
exact beyond the decoupling height and the boundary law is the flux identity. Numerical first step: compute sigma_+ for
the prime-free symbol at a = 0.5 by the Cauchy integral of log sigma and compare t sigma_+(t) F(t) e^{-ita} with a
constant.
DONE 2026-09-30 (Computation 7.20, rh_wiener_hopf.py, data/wiener_hopf_a0.5_*.log): the numerical step, and it is decisive. With
sigma~ the symbol with its real zeros divided out and theta = (1/2) H[log sigma~] (Hilbert transform, odd; ~ -pi/(4 log(t/2pi))),
the nulls of S = int f sin(tx) on the prime-free critical point at a = 0.5 sit at t a = pi/2 + k pi - theta(t) + kappa/t with rms
residual 1.2e-4 rad over 319 nulls (2.8-120 T*), kappa = 4.34 (vs t_lambda = 4.56, the arctan of the zero's factor); the raw shift is
0.34 -> 0.14, the local law c/log leaves 3.5e-3. Pointwise, S = -|A| cos(ta + theta - kappa/t)/(t sigma~^{1/2}) with |A|^2 = -lambda'/2
and nothing else fitted: L^2 error 4.5e-3, 3.4e-4, 8e-5, 3e-5 at 3, 10, 20, 50 horizons, decaying like T^-2. With the prime 2 the
phase oscillates at frequency log 2 (amplitude ~0.98/2sigma~) and the nulls follow it (rms 0.010 vs 0.070 for the prime-free phase).
Conjecture 7.19 is now stated in this sharp form: A purely imaginary (no free phase), one real kappa. The paper carries the
half-line derivation: cutoff chi, (K - lambda)(chi f) = h with h = 2P chi sinh + [K, chi] f smooth up to the edge; sigma F~ = G~ + H~
with F~ analytic below (O(1/t)), G~ above (O(t^-1 log^1/2 t), Lemma 7.12), H~ = i h(a)/t + O(t^-2); factor, split, Liouville
=> F~ = i h(a)/(sigma_+(0) t sigma_-(t)) + O(t^-2); real zeros of the symbol = standing waves of the interior = simple poles of E;
the reality symmetry keeps A imaginary and makes the next term a real phase kappa/t; |A| is fixed by the far edge (quantization of
lambda(a)) and |A|^2 = -lambda'/2 is the boundary law as a flux (|A| = C).
DONE 2026-09-30, second pass (the derivation rewritten in the paper; rh_wiener_hopf.py part (iii)): the exact two-edge equation
needs no cutoff: sigma_a F~ = H~ + G~_+ with F~ = e^{-ita} F analytic below, G~_+ (outside force on the right) analytic above,
H~ = polar step - 2iP sinh(a/2)/t + e^{-2ita}(far edge) analytic below; the primes OUTSIDE the window cancel identically (their
copies of f lie in g_out), so the symbol is sigma_a = Psi_a - lambda with the entries log n < 2a only. Liouville + the reality
symmetry F~(-t) = conj F~(t) => A(t) = A_0 + A_1/t with A_0 imaginary, A_1 REAL: the first correction is a pure phase kappa/t
(proved at the level of the expansion). With primes: the entries inside the window copy the edge to a - log n (echo), the echo
feeds back on the edge as a real term -c^2/sigma_inf (c = Lambda(n)/sqrt n), so the edge's symbol is sigma_eff = sigma_inf -
c^2/sigma_inf and the echo piece is (c/sigma_inf) e^{-it log n} times the edge piece; longer chains give a finite continued
fraction; the window's symbol sigma_inf - 2c cos = the INFINITE chain (fixed point). Verified at a = 0.5 with the prime 2: null
residual 1.1e-3 rms on 5-80 T* (window's symbol 9.8e-3), pointwise 0.010, 0.0021, 0.0007 at 10, 20, 50 T*, ~T^-2 (window's
symbol 0.024, 0.016, 0.010, not decaying). Conjecture 7.19 restated: (i) prime-free sharp form, (ii) the echo form.
DONE 2026-09-30, third pass: THE PRIME-FREE FORM IS A THEOREM. Paper Lemma 7.21 (the factor of a logarithmic symbol: the
outer function of sigma~^{1/2}, bounds c <= |sigma_+| <= C (log)^{1/2}), Lemma 7.22 (the kernel of 1/sigma_+ beyond the origin
decays like e^{-eta xi}: sigma~ has no zeros in a strip below the line), Proposition 7.23 (exact Wiener-Hopf representation:
for EVERY L^2 odd critical point of W_inf - 2P^2, F~ = P_-[H^_-/sigma_+]/sigma_-, with the real zero of the symbol removed by
the rational factor and Liouville on an entire function of growth (log)^{1/2}), Corollary 7.24 (the form of the edge: F~ =
A(t)/(t sigma_-) + e^{-2ita} Phi_-/sigma~ exactly, A(t) = A_0 + O(log t/t), A_0 purely imaginary, A(-t) = -conj A(t); the
far-edge feedback through the kernel of 1/sigma_+ beyond 2a; the mirror: F = A e^{ita}/(t sigma_-) - conj(A) e^{-ita}/(t sigma_+)
+ O(t^-2 (log t)^{-1/2}), hence S = -|A_0| cos(ta + theta - kappa(t)/t)/(t sigma~^{1/2}) + ...). Consequence: for every prime-free
critical point the tail energy is 2|A_0|^2/(pi T)(1 + o(1)) at every height beyond the onset, so the tail law is the single
identity 2|A_0|^2 = -lambda' (the boundary law as a flux, |A_0| = C). What the theorem does not decide: |A_0|^2 = -lambda'/2;
kappa constant (bound O(log t), data 1e-4); the rate (O(log t/t) against the measured t^-2); the echo form with primes.
DONE 2026-09-30 (a): the flux identity B_inf = a C^2 is Proposition 8.3, proved in x-space (not from the Fourier tail, which
only gives the rate: B_inf - B_T ~ a C g_reg L_T^{-1/2}, the polar cross term, g_reg = 2P sinh(a/2)); it needs the edge law
with two derivatives (the hypothesis of Prop 7.8) and Lemma 8.2. Verified: C from the edge nodes 1.0615, C^2 = 1.127 vs
-lambda'/2 = 1.129; the gap a C^2 - B_T = 0.014-0.0125 at 5-20 T* against the leading term 0.024-0.019 (L_T only 2-4 there).
NEXT for item 2 (what now separates the tail law from a theorem for the prime-free critical points): the DERIVATIVE CONTROL of
the edge from the representation: show that the inverse transform of A(t)/(t sigma_-(t)) with A - A_0 = O(log t/t) (Corollary
7.24) is C(L + beta)^{-1/2} + psi with psi^{(k)} = o(delta^{-k} L^{-3/2}), k <= 2, and C = |A_0|. Then Prop 8.3 gives
2|A_0|^2 = -lambda', Corollary 7.24 gives the tail energy 2|A_0|^2/(pi T), and the tail law is a theorem without primes.
(b) the echo form with primes as a theorem (the exact representation holds with the window's symbol; the edge asymptotics
reorganize into the finite echo chain). Old list, still valid for the remaining parts: (a) the factorization of a symbol growing
like log t, factors (log t)^{1/2} e^{+-i theta}, theta = O(1/log t); (b) the a priori decay of the three transforms: F~ = O(1/t) is the
edge-law upper bound (Theorem 3.2 under (H_inf)), G~ from Lemma 7.12, H~ from the smoothness of the commutator [K, chi] f near the
edge; (c) Liouville with the poles at the real zeros, and the far edge at relative order t^-2; (d) with primes, the echo chains
in general (several entries, chains of several links, the mirror chains, coincidences at the entries a = a_n). Then the tail law
is exact beyond the decoupling height and the boundary law is the flux identity |A|^2 = -lambda'/2.
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

## 3. The onset at a single support  [REFORMULATED 2026-09-30: Lemma 7.26, Proposition 7.27; NEXT]
Exhaustion height: 2N(T) = 2aT/pi exactly at T = e T* (Lemma 7.26). Data: medians of the floor's mass 2.39 ... 2.74 T*
(rising to e), theta(e T*) = 0.43 ... 0.50. Proposition 7.27: under RH, (S_kappa) [zeros above kappa T* sample the tail
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
alignment of Computation 7.28), and the gap condition it needs (no gap wider than pi/a above e T*).
A sampling-type lower bound: for the minimizer at one support a, sum_{gamma > kappa T*} |F(gamma)|^2 >= c C^2/T* with
explicit kappa, c. Inputs: density of zeros above the horizon exceeds a/pi; no gaps wider than pi/a above a few horizons
(true on average from T*, for the largest gaps from about 3T*); the minimizer's transform above kappa T* is
edge-dominated (observed, no argument yet). Under RH this is Conjecture A at that support. To prove RH by A, item 3
must eventually be replaced by its prime-side form: the dilation derivative of the ground state controlled by T* lambda
using only the Euler product and the archimedean symbol.

## Not to drift into
Refitting the decay constant; more grids; more supports for the tail law. The data are sufficient; the gap is analytic.
