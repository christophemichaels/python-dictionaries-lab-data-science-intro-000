# Program: from the tail law to the why

Agreed 2026-09-30. The rate of decay of the floor is pi times the effective height of the zeros that carry it,
Phi' = pi T_eff (tail law, paper Conjecture 7.7). Conjecture A (which implies RH) is "T_eff is a bounded multiple of
the horizon". Three items, in order; this file is updated at the end of every round so that the work does not drift.

## 0. The unconditional statement  [2026-10-01: paper Section 9, Theorem 9.1 (the reduction), Computation 9.2; THE BOUNDARY LAW AT ALMOST EVERY SUPPORT 2026-10-01: Section 9.2, Theorem 9.8, Corollary 9.9]
Nothing proves RH. What is proved (under (H_inf), at every support below a_inf = 0.843): the boundary law lambda' = -2|A_0|^2 and
the analytic tail law, both with no zeros involved. So RH reduces to ONE INEQUALITY about the minimizer, 2|A_0(a)|^2 <= c T* lambda(a)
for all a >= 0.8 (Theorem 9.1; the boundary law that links the two is now a THEOREM at every Diophantine support, i.e. almost
everywhere, Theorem 9.8, so the reduction no longer depends on it: Corollary 9.9), equivalently:
the energy of the minimizer against the shifted symbol above kappa horizons is at most c lambda/(pi kappa). Under RH this is the
zero sum above kappa T*, bounded by lambda because its terms are squares (Props 7.27, 7.32) - that is the ONLY place RH is used.
Unconditionally the inequality CONTAINS RH: if RH fails the floor is negative from some support on while the left side is >= 0.
So no argument that treats lambda as a given number can work; a proof must use the minimality of f to exclude a negative floor,
as the Rellich-Pohozaev identity does for the Laplacian (where the remainder is a sum of squares); here the identity is the
dilation identity and the obstruction is the sign of the prime part D_P. A bootstrap in height from the verified zeros gives
positivity only to a ~ (1/3) log(H/2pi) ~ 9 and gains O(1/T*) per step. The prime-side form: at a finite-tree support
A_0 = P A_0^(1)(lambda, a), P^2 = 1/|phi_1|^2, lambda the lowest root of the secular equation <phi_1, s> = 1, so the inequality is
2|A_0^(1)|^2 <= c T* lambda |phi_1|^2 at that root, with only the archimedean symbol and the entries entering; beyond a_inf the
entries enter only through the spectral measure mu_b of the echo lattice at the edge, whose moments are sums over the closed
walks through the edge: MULTIPLICATIVE RELATIONS m_1^{e_1} ... m_{k+2}^{e_{k+2}} = 1 among prime powers below e^{2a} with partial
products in the window, weighted by prod Lambda(m_i)/sqrt m_i (Computation 9.2: m_1 = 0.1-0.25 since a ratio of prime powers is
rarely a prime power; m_8^{1/8} sits 2.1-2.3 below 2a + log 5 at every a to 2.5, a factor 9 in height; Chebyshev with m_8 puts
< 4% of mu_b above five horizons). The uniformity of the onset = these relations do not accumulate at the top of the spectrum:
m_k^{1/k} <= 2a + log kappa - delta for k growing with a. This is where the problem is arithmetic and where the paper stops.
THE PLANE (2026-10-01, paper 9.1: Theorem 9.3, Corollary 9.4, Computation 9.5; rh_plane.py). The third leg of the triangulation
taken unconditionally: lambda = sum over ALL zeros of q(gamma_rho), q(z) = 4 S(z)^2, the quadratic functional of the minimizer's
transform on the plane; on the line q >= 0, off it a quartet +-gamma_0 +- i delta contributes 16 Re S(gamma_0 - i delta)^2, and to
second order q(x + iy) has real part 4S^2 - 4y^2 S'^2, negative near the nulls of S, where the zeros are. Theorem 9.3: (i) RH <=>
lambda(a) >= 0 for all a (Weil on the odd sector); (ii) under (H_inf) + exponential bound on the variation, limsup log^+(-lambda)/2a
<= delta_max = sup(Re rho - 1/2) (Gallagher + Plancherel-Polya); (iii) if the sup is attained, lambda(a) <= -e^{2 delta_max a}/(2 delta_max)
(the test function sinh(delta u) cos(gamma_0 u) phi(u/a)): THE FLOOR'S NEGATIVE RATE IS THE ZEROS' DISTANCE FROM THE LINE.
Corollary 9.4 (with the boundary law): either the edge amplitude grows like e^{2 delta_max a} (RH false) or it decays faster than any
exponential (RH + Conj A): no third behaviour; THE EDGE READS THE PLANE. Every unconditional theorem of the paper holds on both
branches (none uses the sign of the floor): that is why they cannot decide, and why inequality (25) is exactly "the negative branch
does not occur"; a proof must find a quantity the branches force apart BEFORE the floor changes sign. Computation 9.5 (K-mode
minimizers, a = 0.6, 0.8, 1.0): the stiffness kappa_2 = sum 16 S'(gamma)^2 = 1-5e-3, entirely on the zeros below the horizon; the
critical displacement y_c = (lambda/kappa_2)^{1/2} = 1.1e-2, 2.8e-6, 3.5e-12: the floor certifies the zeros below the horizon to
be within y_c of the line, a resolution e^{-cT*/2}; the growth of |q| into the strip above the horizon is the edge's cosh(2ya) to
3-4 digits at every y <= 0.45 (the plane reads the edge), and departs from it in the mass region; Re q < 0 on 7-15% of the
heights at y = 0.2 and 19-27% at 0.45, independent of height and support.
THE BOUNDARY LAW AT EVERY SUPPORT (2026-10-01, paper 9.2: Lemma 9.6, Proposition 9.7, Theorem 9.8, Corollary 9.9, Computation 9.10;
UPGRADED 2026-10-01 evening: the Diophantine condition DROPPED, Theorem 9.8 holds at EVERY support that is not an entry, resonant
supports included (the chain landing exactly on the edge is a relative O(L^-2) correction inside the remainder class; near-misses
are words of length >= 3 for small eps at a non-resonant support, weight O(s^-3), by Lemma 9.6(i) applied to w concatenated with the
reverse of the resonant word): no exceptional set at all;
rh_mirror_chains.py, rh_band_phase.py). The one conditional link of Theorem 9.1 beyond a_inf removed. Lemma 9.6: the lattice
Lambda_a = {log(M/M')} of the entries; a walk of k steps cannot return to within (1/2) e^{-2ak} of its start (|log(M/M')| >= 1/(2 min)
for M != M' < e^{2ak}); a is Diophantine if |2a - u_w| >= c e^{-kappa |w|} for all words, which holds outside a null set
(Borel-Cantelli, kappa > log(2|D_a|)). Proposition 9.7: above the periodic band (s > X = sum 2c_d) the window's symbol has the
absolutely convergent expansion log sigma_s = log s + sum_u beta_u e^{i tau u}, beta_u <= 0 indexed by the walks, sum |beta_u| =
-log(1 - X/s), beta_0 = log G(s) - log s (the geometric mean of Lemma 8.18); the Wiener-Hopf factors are G^{1/2} exp(sum_{+-u>0}
beta_u e^{i tau u}); 1/sigma_{s,+} = G^{-1/2} sum_{u>=0} rho_u e^{i tau u} with rho_u >= 0, the words of total length k weighing at
most (X/s)^k, so the weight within eps of the edge is <= (2 eps)^{alpha(s)}, alpha = log(s/X)/2a, and within eps of 2a at a
Diophantine support <= (eps/c)^{log(s/X)/kappa}; (s/G)^{1/2} = 1 + |b|^2/(2s^2) + ... = the effective symbol's correction, and
rho_{+-d} = c_d/s + ...: the first-generation echoes. Theorem 9.8 (under (H_inf), a Diophantine, not an entry): (i) A_a(t) -> A_0,
A_a = A_0 + (mirror chains: sum over lattice points v > 2a of m_v e^{it(v-2a)}, |m_v| <= C rho_v(s), the walks from the far edge
across the window) + O(1/log^2 t); (ii) the edge law with two derivatives with the SAME constants C = |A_0|, beta = -gamma - log(2 pi a)
- lambda, remainder O(delta^{-k} L^{-5/2}), the echoes within delta of the edge weighing O(delta^{alpha(L)}); (iii) lambda' = -2|A_0|^2
at every point of differentiability and, under RH, the tail law. Hence lambda(a_2) - lambda(a_1) = -2 int |A_0|^2 for ALL a_1 < a_2:
THE BOUNDARY LAW HOLDS ALMOST EVERYWHERE AT EVERY SUPPORT. Corollary 9.9: (H_inf) + the inequality at a.e. a >= 0.8 => RH. Proof
method: Corollary 8.6 and Proposition 8.4 with the window's symbol, the lattice expansion replacing the kernel lemma; the factors on
the line by the Hilbert estimate H[g e^{itu}] = -i sign(u) g e^{itu} + O(|g|/(|u| t)) summed with the Diophantine split; the far
feedback through the kernel of 1/sigma_{a,+} (rho-weighted spikes at the lattice points), its singularities at v - 2a being the
mirror chains; the inversion with the cutoff t_1 = delta^{-1/2} and the one-sided copies. The theorem is about lambda', NOT about the
onset: its asymptotics begin at s > X e^{2a}, astronomically far; the onset is item 3. Computation 9.10 (FEM at a = 0.6, 114 nulls
on 5-34 T*): the window's symbol phase residual 0.0365 (tree 0.0009, archimedean 0.145) is removed to 0.0186 by the two-step
mirror chains at v - 2a = 0.186, 0.592, 0.997 (log 4, log 6, log 9 minus 1.2) with amplitudes 0.19, 0.51, 0.35 against the chain
weights 0.24, 0.62, 0.40 (common factor 0.8); the tree's residual shows nothing at those frequencies (0.002); by height band the
residual after the fit falls 0.033 -> 0.013 as (X/s)^3 falls 0.34 -> 0.13.
WHAT THE K-MODE MINIMIZERS RESOLVE (2026-10-01, Computation 9.10; corrects Computation 8.19's pointwise claims). The K-mode
minimizer is a polynomial of degree 2K-1: it resolves the edge to delta ~ a/(2K-1)^2 (heights ~ (2K-1)^2/a, 500-700 T*) but the
interior cusps only to heights ~ (2K-1)/a (3-8 T*). At a = 0.6 the K-mode transform agrees with the FEM to 0.1-0.6% at 2-5 T* and
departs at 7 T* (= (2K-1)/a = 7.6 T*) to 15-20% at 10-30 T*; the first-generation echo phase is in the FEM nulls with the predicted
amplitudes (regression 0.502, 0.635 vs c_2 = 0.490, c_3 = 0.634; rms 0.139 -> 0.022) and absent from the K-mode nulls (-0.003,
0.16); at a = 0.8 and 1.0 the K-mode nulls follow the archimedean phase alone to 0.01 rad from 10-30 T* to 400 T* with echo
coefficients < 0.001. So above (2K-1)/a (3.4, 2.9, 2.4 T* at a = 1.0, 1.25, 1.5) the K-mode data say nothing about the pointwise
structure of the minimizer; their envelope (a mean) is right to the few % of sigma_eff/sigma_inf - 1, their zero sums and tail laws
stand, but the claims "|F|^2 uncorrelated with 1/sigma_a, the minimizer does not follow the window's symbol pointwise" of
Computation 8.19 are WITHDRAWN (they describe the polynomial's edge). The true pointwise statement above the band top is Proposition
9.7: |F|^2 ~ 1/sigma_a to first order in c_d/s. Testing the band's pointwise structure at a >= 1 needs an interior-resolving solver
(FEM with spacing < 1/(50 T*)), not more K-modes.
THE GLOBAL AND THE LOCAL (2026-10-01, paper 9.3: Proposition 9.11, Theorem 9.12, Proposition 9.13). Two global objects in the
companion work: the Atlas's complete Weil form W(q,q) on all of C_c^inf (positivity on every support = RH; Suzuki defect E_N; local
support theorems, CP20: Q_a^-(q) > 5278/10^6 |q|^2 for a <= 3/8) and the Mobius Green energy R(N) = sum mu(m)mu(n)/max(m,n)
(R(N) = O(N^eps) <=> RH; a zero at beta_0 forces R(N) > N^{2 beta_0 - 1 - delta} i.o.). Proposition 9.11: the Atlas's odd operator
A_a^- = C_a - log 2pi - Pi_a - K_alpha - 2|s><s| IS the window form Q (Weil's kernel kappa, the shifts with b_n = c_n, the polar
term), so CP20 is a CERTIFIED LOWER BOUND lambda(a) >= 5.28e-3 on (0, 3/8], against the engine's upper bound lambda_40(3/8) = 3.14e-2.
Theorem 9.12 (two readings of delta_max): unconditionally limsup log^+ R(N)/(2 log N) = delta_max (Mobius paper + M(x) << x^{Theta+eps});
and limsup log^+(-lambda(a))/(2a) = delta_max (Theorem 9.3): in the common variable a = log N the global energy and the local floor
grow at the same rate 2 delta_max, both zero iff RH. Detection scales: the window's floor turns negative only from a_det =
O(delta_0^{-1} log(1/delta_0)) (the quartet must outweigh the on-line zeros' O(a^2)), the Mobius budget exceeds C_M log N only from
N^{delta_0} >> |rho_0 zeta'(rho_0)| sqrt(C_M): the same scale, exponentially late in 1/delta_0 - THE WINDOWS THAT TAKE INCREDIBLY LONG
TO CLOSE. Proposition 9.13 (late windows are allowed): lambda' >= -(c T* + e(a)) lambda a.e. with ANY locally integrable excess e
still gives lambda > 0 everywhere and RH (Gronwall); conversely RH <=> -lambda'/lambda in L^1_loc. So the content of Conjecture A
is the RATE c T*, not positivity; windows where the rate is exceeded are permitted as long as the excess integrates: the resonant
supports (null set) cost nothing, a late onset costs its excess over the interval of lateness; a NON-integrable excess is a sign
change = a zero off the line. The all-support comparison the Atlas leaves open IS inequality (25) in the form of Proposition 9.13.
NOT done: the transfer of the quantitative witness delta = -lambda(a)|f|^2 ~ e^{2 delta_max a} to the Suzuki defect E_N needs the
Suzuki norm S(f,f) of the minimizer; the Atlas's complete energy identity (1/tau_R) K_a(I - K_a) = eta_p W*W + N_a + sum b_n P_{a,log n}
is not yet related to the edge/echo structure of the minimizer.
TWO TRANSFERS AND THE BLIND SPOT (2026-10-01, paper 9.4: Proposition 9.14, Proposition 9.15, Proposition 9.16, Theorem 9.17, Computation 9.18;
rh_blind_spot.py). (1) The Suzuki transfer done: with S(f,f) <= M_a^2 |f'|_1^2/pi (Atlas (B.7)) and the Atlas's CP9 (D1) transfer
applied to a mollified minimizer, lambda(a) < 0 gives E_N >= (1 + pi(1-eps)|lambda(a)|/(3 M_a^2 V_a^2))^2 for N >= N_0, V_a = |f'|_1 =
4 sup|f| (sup|f| = 1.46, 1.42, 1.40, 1.38 at a = 0.6-1.25, at the fixed position u = 0.22-0.23); with Theorem 9.3(iii) the defect's
permanent floor grows like e^{2 delta_max a}/(M_a^2 V_a^2). Uncontrolled: Suzuki's feature constant M_a and the variation of the
negative-branch minimizer. (2) The paired-prime matrices done: <q, Pi_a q> = 2 sum c_n h_q(log n) (the prime part of Q = the sampling
defect, Prop 7.29); on the minimizer with a finite tree its tail part is (M_chi - 1) x density = <e^T C_a e / sigma~> (Cor 8.14(vi)),
C_a the tree's weighted adjacency with the edges; in the periodic model the mean of 2cos(tau d)|E|^2 is e^T C^{(d)} e: THE ATLAS'S
PAIRED-PRIME MATRIX OF THE ENTRY n, COMPRESSED TO THE EDGE AND ITS ECHOES, IS THE ADJACENCY OF THE STEP log n; P_a is C_a.
(3) THE BLIND SPOT IS A THEOREM (9.17): 16 Re S(gamma - i delta)^2 = 16 S^2 - 16 delta^2 (S'^2 + S S'') + O(delta^4): every
explicit-formula functional (floor, Mobius energy) is EVEN in the distance from the line - no first-order detector exists anywhere.
The floor's single-zero sensitivity delta_c(gamma) = (lambda/16 S'(gamma)^2)^{1/2} is e^{-cT*/2} below the horizon (3.5e-12 for
the first zero at a = 1, 2e-24 at a = 1.25: NO BLIND SPOT BELOW THE HORIZON) and exceeds 1/2 from about 1.5 T* on (10-100 from
five horizons: BLIND ABOVE, the slopes being the edge's cosine, delta_c^2 >= gamma^2 sigma~ lambda/(8a^2(-lambda'))). The Mobius
energy (weight 1/|rho|^2) is blind in the same place. Verified: the second-order formula to 1.0000 at delta <= 1e-2 and 0.4% at 0.1.
CONSEQUENCE FOR ITEM (a): the search for a branch-separating quantity BELOW the sign change is closed: at supports where lambda > 0
every functional differs between the branches at second order, with the floor itself carrying the sharpest coefficient; the branch
is decided only at the sign change. RH is not detection at any finite scale; it is the structure, i.e. inequality (25) in the form
of Proposition 9.13. Item (a) is withdrawn in favour of (b), (c), (d).

## 1. Edge law => tail law, as a proposition  [DONE 2026-09-30: paper Proposition 7.8]
Under RH and the edge law with two derivatives, lim T sum_{|gamma|>T} |F(gamma)|^2 = 2C^2/pi. Proof: edge asymptotics
of the transform (splitting at delta = 1/t), then the sum over zeros through the Riemann-von Mangoldt formula with
Littlewood's bound S_1(t) = O(log t/(log log t)^2). With the boundary law lambda' = -2C^2 this is Conjecture 7.7.
Still open inside item 1: the edge law as an asymptotic equality with derivatives (only the upper bound is a theorem,
under (H_inf)). The boundary law is now a theorem under that same hypothesis (paper Proposition 8.3, 2026-09-30): B_T ->
a C^2 by the overlap of the leak of the dilation generator with the edge force, J_1(L) = 1 + O(1/L^2) with no 1/L term (the
antisymmetry of Lemma 7.11 again), and Lemma 8.2 (the edge force to second order). So under the edge law with derivatives:
boundary law => tail law (Prop 7.8), and both are the flux identity 2|A_0|^2 = -lambda' of the Wiener-Hopf form.
DONE 2026-10-01: the edge law with derivatives is a THEOREM under (H_inf) for the prime-free form (Prop 8.4) and for the full
form at a_2 < a < a_3, the prime 2 inside the window (Cor 8.8(ii)), with the constants C = |A_0|, beta = -gamma - log(2 pi a) - lambda.
At those supports item 1 is closed end to end: edge law (theorem) => boundary law (Prop 8.3) => with RH the tail law (Prop 7.8).

## 2. The sum rule  [THEOREM AT ALMOST EVERY SUPPORT 2026-10-01: Theorem 9.8, the boundary law and the edge law at every Diophantine support; DONE 2026-09-30: Proposition 7.9, Lemma 7.11, Lemma 7.12, Corollary 7.13; sharp form Conjecture 7.19 verified, Computation 7.20; THEOREM for the prime-free form 2026-10-01: Propositions 7.23, 8.3, 8.4, Corollaries 7.24, 8.5; THEOREM for a_2 < a < a_3 (the prime 2 inside the window) 2026-10-01: Proposition 8.7, Corollary 8.8; THEOREM AT EVERY SUPPORT WITH A FINITE ECHO TREE, a < a_inf = 0.843, 2026-10-01: Lemma 8.12, Proposition 8.13, Corollary 8.14]
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
DONE 2026-10-01: THE TAIL LAW IS A THEOREM FOR THE PRIME-FREE FORM under (H_inf). Paper Proposition 8.4 (the edge law from the
representation): the inverse transform of A_0/(t sigma_-(t)) is the edge law with its constants, C = |A_0| and
beta = -gamma_E - log(2 pi a) - lambda, i.e. f(a - delta) = |A_0| (sigma_inf(1/delta) - gamma_E)^{-1/2} (1 + O(L^-2)): the
profile is the inverse square root of the symbol at frequency 1/delta shifted by Euler's constant. Proof: Mellin expansion of
the transform of (log 1/delta)^{-1/2} (the real L^-2 corrections cancel; the phase e^{-i pi/(4M)} appears with M = log t + gamma),
the expansion theta = -pi/(4 l) - pi^3/(24 l^3) + O(l^-5) + O(log log t/t) from the x-integral form of the Hilbert phase (even
powers vanish by x -> 1/x), matching to relative order l^-3, and a splitting lemma for the inverse transform of the remainder
with three derivatives. Corollary 8.5: lambda_inf' = -2|A_0|^2 and T (1/2pi) int_{|t|>T} |F|^2 sigma -> -lambda'/pi: the
analytic law of Computation 7.14 is a theorem, Conjecture 7.19(i) holds with |A|^2 = -lambda'/2. Computation 8.6 /
rh_edge_constants.py: with NOTHING fitted, f / [(-lambda'/2)^{1/2} (L + beta)^{-1/2}] = 1.001-1.003 on 1e-10 < delta < 1e-5 at
a = 0.3, 0.4, 0.5, 0.6 (1.0003-1.0009 on the finer mesh): the boundary law's C and the symbol's beta, both predicted.
DONE 2026-10-01 (b): THE FULL FORM IS A THEOREM FOR a_2 < a < a_3. Paper Proposition 8.7: with exactly one entry n inside the
window (a < d = log n < 2a, so n = 2) the entries outside the window move it off itself, and the eigen-equation on the whole
line reads K_inf f = lambda f + 2P s 1_[-a,a] + c f_L(. - d) + c f_R(. + d) + g_inf 1_{|y|>a}, c = Lambda(n)/sqrt n: the copies
of the two thirds L, R of the window are SOURCES for the archimedean operator, and Proposition 7.23 applies verbatim. The copy
on L is the echo, f(a+u) = c (k * f_R(a + . + d))(u) + smooth near u = -d, k the kernel of 1/sigma~_inf (eq. 18). The copy on
R near the edge is, by the echo at u - d, c^2 (k * f_R(a + .))(u) 1_{u<0} + a step + a far source, i.e. c^2 P_-[F~/sigma~_inf]
up to smooth sources; moving c^2 F~/sigma~_inf to the left gives sigma_eff F~ = H'_- + G'_+ with
sigma_eff = sigma~_inf - c^2/sigma~_inf and the representation F~ = P_-[H'_-/sigma_eff,+]/sigma_eff,- (eq. 19).
Corollary 8.8: (i) the edge form of Cor 7.24 with sigma_eff, A_0 purely imaginary; (ii) the edge law with derivatives with
the SAME constants C = |A_0|, beta = -gamma - log(2 pi a) - lambda (the return term changes the symbol at relative order
c^2/l^2 only), f smooth on the window except at the edges and at the echo points +-(a - d); (iii) the echo form of
Conjecture 7.19(ii) with rho = c/sigma~_inf; (iv) lambda' = -2|A_0|^2; (v) under RH, T sum_{|gamma|>T} |F|^2 -> -lambda'/pi:
CONJECTURE 7.7 IS A THEOREM FOR a_2 < a < a_3 under (H_inf) and RH. Computation 8.9 (rh_edge_constants.py,
data/edge_constants.log): on the full minimizers at a = 0.4, 0.45, 0.5 (lambda' = -0.466, -0.101, -0.0112 from the dilation
engine) the ratio f / [(-lambda'/2)^{1/2} (L + beta)^{-1/2}], beta = -1.514, -1.619, -1.722, nothing fitted, is 1.000-1.002 on
1e-12 < delta < 1e-6 and 1.004-1.008 at 1e-4; the alternative beta - 2c (a first-order return of the echo) is off by 2-4%:
the prime inside the window changes neither constant.
DONE 2026-10-01 (c): THE ECHO TREE. The entries d = log m < 2a generate walks (steps +-d) inside the window; the echo set S is
the set of points reached from +-a, and the partition of the window by S is invariant under the entries (Lemma 8.12: the entries
move the pieces onto pieces). Where the walks are confined (every orbit finite, (F_a)) the tree is finite, and Proposition 8.13
gives, with the weighted adjacency C on S° (C_pq = c_d for |p - q| = d) and the couplings b_p = c_{a-p}: (i) f smooth off S°;
(ii) the local echo relations f(p + v) = sum_q c_{|p-q|} (k * f(q + .))(v) + r_p, i.e. (sigma~ I - C) Phi_S° = b Phi_a + bbar Phi_{-a}
+ O(t^-2 log^-3/2); (iii) the exact representation at the edge with the EFFECTIVE SYMBOL = SCHUR COMPLEMENT
sigma_eff = sigma~ - b^T (sigma~ I - C)^{-1} b = det(sigma~ - C_a)/det(sigma~ - C), a rational function of the archimedean symbol
(a finite continued fraction along a chain), with finitely many real zeros (eigenvalues of C_a) and poles (eigenvalues of C: the
standing waves of the interior tree, giving rational corrections). Corollary 8.14: the edge law with the SAME constants; the echo
form S = -|A_0||E| cos(ta + theta_eff + arg E - kappa/t)/(t sigma_eff^1/2), E = 1 + sum rho_p e^{-it(a-p)}, rho = (sigma~ - C)^{-1} b;
THE SHARP FORM IS AN ALGEBRAIC IDENTITY OF THE TREE: <sigma_a |E|^2> = sigma_eff (mean over the oscillations), and the mirror
term's mean is -Rbar (zero when the trees of the two edges are disjoint, as at every support computed); the boundary law; under RH
the tail law (Conjecture 7.7) AT EVERY SUPPORT WITH A FINITE TREE; and the sampling constant M_chi - 1 =
<e^T C_a e/sigma~>/<sigma_eff/sigma~> averaged over the tail with the weight chi_T/t^2, = 2|b|^2 <sigma~^-2> to first order,
|b|^2 = sum_{log m < 2a} Lambda(m)^2/m. Computation 8.15 (rh_echo_tree.py, data/echo_tree.log): |S°| = 2 (one entry), 6 for
a_3 < a < a_4 (the echo of an echo at a - log 3 + log 2), 12, 22, 148 (a = 0.825), 426, 3434, 8854 near the threshold; the walks
DECONFINE at a_inf = 0.8430 (between 0.84297 and 0.84336), where the echo set becomes dense: the infinite tree has a band, which
is why the window's symbol is negative far above the horizon for a >= 1 (last sign change 33 T* at a = 1.0, > 150 T* at 1.25)
while the finite tree has only finitely many resonances, all below the horizon. Computation 8.16 (edge FEM at a = 0.6 with the
entries 2, 3, lambda = 5.946e-7, lambda' = -4.1545e-5 from the engine at K = 64): the constants unchanged (ratio 1.0003-1.003 on
1e-12 < delta < 1e-6, nothing fitted); the nulls follow the tree's phase with rms 0.0114 on 5-80 T* against 0.0162 for the first
generation alone and 0.036 for the window's symbol; pointwise S is the tree form to 0.038, 0.0088, 0.0034 at 5, 10, 20 T* against
0.053, 0.028, 0.020 for the first generation: THE DIFFERENCE IS THE SECOND ECHO c_2 c_3/(sigma~^2 - c_2^2) = 0.040, 0.026, 0.018.
The explicit-formula identity holds to 1e-5 of the zero sum at every height.
CORRECTION to Corollary 7.30 (found with the sharper cutoff 1 - exp(-(t/T)^8), which the quartic cutoff's leakage below the
horizon had masked): the sampling defect is 2c^2 <sigma~^-2>_T averaged over the tail with the edge's weight chi_T/t^2, not
2c^2/sigma(T)^2; the average is 0.7 times the value at T at ten horizons (int_1^inf du/(u^2 (L + log u)^2) = L^-2 (1 - 2/L + ...)).
DONE 2026-10-01 (d): THE BAND. Beyond a_inf the exact statement is the representation with the WINDOW'S SYMBOL sigma_a = Psi_a - lambda
(entries inside the window only), two sources (polar step, far force), valid at EVERY support (Prop 8.17; Lemma 7.21 extended to
symbols whose log-derivative is bounded but not L^2: theta, theta' = O(log log t)); its edge amplitude A_a(t) does not converge
(the kernel of 1/sigma_{a,+} is singular on the lattice: the mirror chains; and the rational corrections of the N zeros are
O(sum t_j^2/t^2)), so it is asymptotic only far above the last sign change. The finite tree is the resummation: for one entry the
factorization s - 2c cos tau = sigma_fp (1 - rho e^{i tau})(1 - rho e^{-i tau}) is exact with sigma_fp the fixed point of
sigma -> s - c^2/sigma (Lemma 8.18: the periodic limit; the geometric mean of the window's symbol over the torus is the log
potential of the law of sum 2c_d cos tau_d, = s - |b|^2/s + ..., = the fixed point outside the band [-2 sum c_d, 2 sum c_d] and
pinned at c inside for one entry). THE BAND THAT MATTERS IS THE CONFINED LATTICE'S (Prop 8.20): the walks are confined to the
window even when S is dense, the adjacency C on l^2(S°) has rho(C) <= D_max <= sum_{d<=a} 2c_d + sum_{a<d<2a} c_d < 2 sum c_d, and
rho(C) >= rhobar(a) = sum_d 2c_d (1 - d/2a) (the mean weighted degree, by equidistribution + Folner), rhobar ~ 4e^a/a. Numerically
(rh_lattice_band.py): rho(C) = 1.94, >= 2.91, >= 3.86 at a = 1.0, 1.25, 1.5 (rhobar = 1.74, 2.80, 4.10; periodic 5.85, 8.52, 13.0),
so the confined band ends at 0.9, 1.5, ~3 T* while the periodic band reaches 47, 412, 22000 T* and the last sign change 33, 259, 397 T*.
THE SHARP FORM HOLDS ABOVE THE CONFINED BAND, INSIDE THE PERIODIC ONE (rh_band.py on the K-mode minimizers at a = 1.0, 1.25, 1.5,
Computation 8.19): <t^2 (Psi - lambda)|F|^2>/(-lambda') = 1.07/1.22/1.08 at 5 T*, 1.05/1.05/1.03 at 7, 1.00/1.01/1.01 at 10,
within 3% to 50 T*, with the window's symbol negative on 2.6-6% of the heights at 5 T* carrying -0.02 of the envelope, and
|F|^2 UNCORRELATED with 1/|sigma_a| (|corr| < 0.04): the minimizer does not follow the symbol's sign changes; the law is
carried by the mean. The prime part (envelope against sigma_inf minus the full one) is 15-30% at 5 T*, a few % at 10.
THE LATTICE'S BAND AGAINST THE EDGE'S (the output of this round, in two steps). (1) rhobar(a) - 2a = -0.25, +0.30, +1.1 at
a = 1.0, 1.25, 1.5 but +3.5 at a = 2 and +7.6 at 2.5, so the confined lattice's band top e^{rho - 2a} T* is ~3 T* at 1.5, ~35 T*
at 2, ~2000 T* at 2.5: if the onset followed the SPECTRAL RADIUS, T_eff would outrun the horizon from a ~ 2 and the route of
Prop 7.32 (fixed kappa) would fail. (2) IT DOES NOT FOLLOW IT. The edge couples to the lattice through b alone (the first-generation
points), so what enters sigma_eff = s - R(s) is the SPECTRAL MEASURE mu_b OF C AT b: R(s) = int dmu_b/(s - nu), real above
supp mu_b, complex inside (the edge radiates) (Prop 8.22; the moments of mu_b are exact under truncation to depth g up to order
2g - 2, so Lanczos from b resolves it). Computed (rh_spectral_measure.py, Lanczos 80 steps on 150000 points, checked against the
dense eigenproblem on 6000): mu_b sits at the BOTTOM of the spectrum (mean 0.1-0.25, width 1-2.3 against rho(C) = 2-8), its
weight above the horizon's height 2a is <= 6% at every a, and its 1% quantile ends at 0.8, 1.3 T* for a = 1.0, 1.25 (where
rho(C) itself is resolved: 1.94, 2.93 < 2a + log 5); at a >= 1.5 the Lanczos measure's top node lies below 2a + log 5, so its
'zero weight above 5 T*' is the truncation's [CORRECTED 2026-10-01, Proposition 9.16]: the support of mu_b reaches rho(C) >=
rho_bar(a) (connected lattice), which exceeds 2a + log 5 from a = 1.75; what the EXACT moments bound is the MASS above 5 T*:
<= 2e-6, 1e-4, 1e-3, 0.8%, 2.2%, 4.7%, 5.6% at a = 1.0 ... 2.5 (Chebyshev with the highest exact even moment, orders 14,14,14,10,8,6,6). So THE PREDICTION OF A RECEDING ONSET IS WITHDRAWN: the
spectral radius is the band of the worst coupling, the edge's coupling sits far below it. What remains predicted for a = 2:
the onset of the law a little above five horizons (7-10 T*, where sigma_eff/s = 0.72-0.77; at 5 T* s is just above the resolved
top of supp mu_b and the return is large, sigma_eff/s = 0.5), and the sampling defect sigma_eff/s at ten horizons going from
0.89 (a = 1) to 0.77 (a = 2). Conjecture A's route through (T_kappa) with a fixed kappa ~ 10 survives this test; the a = 2
minimizer remains the check.
NEXT for item 2: (a) the a = 2 minimizer (the check of the edge's band: onset at 7-10 T*?); (b) the uniformity of the edge's band in
a: prove that the upper tail of mu_b stays within a bounded multiple of the horizon, i.e. that the weight of mu_b above 2a + log kappa
vanishes for a fixed kappa and all a (the first-generation points' walks that stay near the high-degree region are few: a counting
statement about the entries); this is the prime-side form of item 3; (c) kappa constant and the rate, not needed for the tail law.
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

## 3. The onset at a single support  [REFORMULATED 2026-09-30: Lemma 7.26, Proposition 7.27; THE SAMPLING DEFECT IS THE PRIME SUM AT THE ENTRIES 2026-10-01: Proposition 7.29, Corollary 7.30; ONE HYPOTHESIS ON THE MINIMIZER ALONE 2026-10-01: Proposition 7.32 (the sharp form); NEXT]
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
density, constant ~1. DONE 2026-10-01: THE ANALYTIC FORM OF M. Paper Proposition 7.29: under RH, with the entire high-pass chi_T(t) = 1 - exp(-(t/T)^4)
and h = |F|^2 chi_T, the explicit formula gives EXACTLY
  sum_gamma |F(gamma)|^2 chi_T(gamma) = (1/2pi) int |F|^2 chi_T Psi_inf - 2 sum_n Lambda(n) n^{-1/2} h^_T(log n) - 2P^2 (1 - e^{-1/(16 T^4)}),
h^_T(x) = (1/2pi) int |F|^2 chi_T e^{itx} dt the autocorrelation of f with its structure below the scale 1/T removed. So the
defect of the zeros against their density is the prime sum at the entries log n < 2a (the autocorrelation of f lives on
[-2a, 2a] and is singular exactly at the entries, where the minimizer echoes its edge), and NO GAP CONDITION ON THE ZEROS
ENTERS. Corollary 7.30: for one entry inside the window (a_2 < a < a_3), by the echo form, M_chi(T) = 1 + 2c^2/sigma_inf(T)^2
(1 + o(1)), c = Lambda(n)/sqrt n: the defect is twice the square of the echo's amplitude rho = c/sigma_inf and tends to zero
like 1/log^2 T. Computation 7.31 (rh_sampling_defect.py, data/sampling_defect.log): on the full minimizers at a = 0.4, 0.45,
0.5 with the first 6000 zeros, every term evaluated independently at T = kappa T*, kappa = 1..20: the identity holds to 3e-5
of the zero sum or better at every height; the prime sum is n = 2 to three digits from kappa = e on (92-97% at kappa = 1,
the rest the smearing to n = 3 just outside the window); M_chi against 1 + 2c^2/sigma^2 is 1.045/1.050, 1.048/1.047,
1.085/1.044 at kappa = 10 and 1.028/1.034, 1.027/1.032, 1.030/1.030 at kappa = 20. At kappa = e the entire cutoff still
reaches the heights where the minimizer hides from the zeros (weight 0.21 at 0.7 T), so M_chi is 1.27, 1.50, 3.33 there
while the sharp M(e) is 1.01, 1.09, 1.02.
THE ARITHMETIC OF CONJECTURE A IN ONE PLACE: (S_kappa) with bounded M <=> the prime sum at the entries,
2 sum_{log n < 2a} Lambda(n) n^{-1/2} h^_T(log n), is at most a fixed fraction of the tail at T = kappa T*, uniformly in a.
A statement about the autocorrelation of the minimizer at the entries (the echoes), no longer about the zeros. For a window
with many entries the echoes form a finite tree and sum_{log n < 2a} Lambda(n)^2/n ~ 2a^2 (Mertens), so the defect need not
tend to zero, but it is bounded as long as each echo keeps its size c_n/sigma_inf.
DONE 2026-10-01: ONE HYPOTHESIS. The prime sum at the entries is the oscillating part of the symbol, so the zero sum with the
cutoff is (1/2pi) int |F|^2 chi_T (Psi_inf - sum 2 c cos) = (1/2pi) int |F|^2 chi_T Psi minus the polar term, exactly. Paper
Proposition 7.32: under RH, if the SHARP tail law holds at kappa T* within a factor M'',
  (T_kappa)  (1/2pi) int |F|^2 chi_T^(m) (Psi - lambda) >= -lambda'/(pi M'' T),   T = kappa T*,  chi^(m) = 1 - exp(-(t/T)^{2m}), m ~ 1.5 T even,
then lambda' >= -pi kappa M'' T* lambda - M'' e^{-8 T*}: Conjecture A with c = pi kappa M'' (the cutoff with m ~ 1.5 T is bounded in
the strip of the explicit formula and makes the polar term 2P^2 (2T)^{-2m}, which Gronwall absorbs for c < 16). (S_kappa) and
(E_kappa) are gone: no sampling constant, no density, no gap condition; the zeros enter only through the exact identity (17).
At every support with a finite tree (T_kappa) holds for kappa large with M'' -> 1 (Corollary 8.14), and the sharp form is there an
algebraic identity of the echo amplitudes. The measured onset is a few horizons at every support computed, including a = 1.0, 1.25
where the tree is infinite; c = pi e M'' at kappa = e is again 8.5 against the observed 3.3-3.9.
NEXT for item 3: UNIFORMITY IN a of the onset of the sharp form, i.e. (T_kappa) with one kappa and one M'' for all a. The band
analysis makes this concrete: the onset is where the archimedean symbol 2a + log kappa clears the band the edge sees, the upper
tail of the spectral measure mu_b of the confined lattice at the edge's coupling. Measured, that tail ends within a few horizons
up to a = 2.5 (mass above 5 T* <= 2e-6 ... 5.6% by the exact moments; the SUPPORT reaches rho(C) ~ 4e^a/a, Proposition 9.16),
so a fixed kappa ~ 10 suffices as far as computed. The prime-side statement to prove: the weight of mu_b above 2a + log kappa
is below a fixed fraction eta < 1 for all a (it is NOT zero: the support outruns the horizon), i.e. m_k <= eta (2a + log kappa)^k
for ONE even k and all a, a statement about closed walks through the
edge of the window with steps log m and weights Lambda(m)/sqrt m, i.e. about the multiplicative structure of the integers below
e^{2a}. Together with the finite-tree theorems this would make (T_kappa) uniform. The a = 2 minimizer checks the picture.
A sampling-type lower bound: for the minimizer at one support a, sum_{gamma > kappa T*} |F(gamma)|^2 >= c C^2/T* with
explicit kappa, c. Inputs: density of zeros above the horizon exceeds a/pi; no gaps wider than pi/a above a few horizons
(true on average from T*, for the largest gaps from about 3T*); the minimizer's transform above kappa T* is
edge-dominated (observed, no argument yet). Under RH this is Conjecture A at that support. To prove RH by A, item 3
must eventually be replaced by its prime-side form: the dilation derivative of the ground state controlled by T* lambda
using only the Euler product and the archimedean symbol.

## Next (as of 2026-10-01, evening)
(a) [CLOSED by Theorem 9.17: no first-order quantity exists, the floor is the sharpest second-order one] The quantity the branches force apart BEFORE the floor changes sign, now sharpened by 9.3: both readings are blind below the
detection scale a_det ~ delta^{-1} log(1/delta), so the quantity must be visible at supports where lambda > 0 and the Mobius budget is
tame; candidates: the stiffness kappa_2 = sum 16 S'(gamma)^2 against y_c (Comp 9.5); the overlap of the minimizer with the Mobius test
function sum mu(n) n^{-1/2} phi(x - log n) (the global sequence evaluated in the local metric); the Suzuki norm S(f,f) of the minimizer.
(b) The long multiplicative relations m_k^{1/k} for k ~ a (the uniformity of the onset, item 3).
(c) The Atlas bridge [DONE in 9.4 up to Suzuki's constant M_a]: compute M_a from Suzuki's feature (B.2) and S(f,f) for the K-mode minimizer
at a <= 1 to make Proposition 9.14 numerical; the paired-prime/echo-tree identification is Proposition 9.15.
(d) An interior-resolving solver at a = 1 (FEM, spacing < 1/(50 T*)) before any further pointwise claim in the band; a K = 200
    minimizer at a = 1 (interior resolved to 8.6 T*) is being computed (data/kmode_K200/) to see the echo phase appear at 4-8 T*.
(e) The mass of mu_b above the horizon, rigorously: prove m_k <= eta (2a + log kappa)^k for one even k and all a (Proposition 9.16(iv));
    the moments are explicit sums over closed walks, m_1 is in closed form (Prop 9.16(iii)), m_2 = |Cb|^2/|b|^2 next.

## Not to drift into
Refitting the decay constant; more grids; more supports for the tail law. The data are sufficient; the gap is analytic.
