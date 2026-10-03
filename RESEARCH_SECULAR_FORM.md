# The prime-side comparison through the polar source

Research note, 2026-10-03 (rewritten the same day in the fixed notation of DICTIONARY.md). Status labels: DEFINITION,
PROVED (with source), HYPOTHESIS. Companion to RH_IF_THEN_TREE.md, item N1d.2, and to
received/Prime_Side_Completion_Addendum.pdf (Section 6 of which states the target this note rewrites).

## 1. Audit of the addendum against the paper

| addendum | paper | status |
|---|---|---|
| Lemma 1, ||C_a|| <= B(a) = 2 sum_{n<e^{2a}} Lambda(n)/sqrt(n) and the resolvent cutoff | the band top 2 sum c_d of Computation comp:band | PROVED, elementary; known |
| Lemma 2, bounded variation of f gives a Lipschitz autocorrelation, hence (H_fin) | the mechanism of Proposition 9.15 | PROVED; known |
| the Cantor caution: an a.e. derivative identity does not integrate without integration regularity | why Proposition 9.17 needs (H_fin) on the exceptional set | correct; Step 1's remainder |
| Theorem 3, the first-null argument with K integrable through every prospective null | Proposition prop:AimpliesRH; Proposition 9.30 | PROVED; known |
| the measure inequality nu(E) <= int_E K lambda as the combined target | new formulation of regularity plus comparison | adopted |
| Section 6, 2|A_0^(1)|^2 <= K lambda ||phi_1||^2 at the root of the secular equation | the paper's prime-side form after Theorem thm:reduction | the same target; Section 2 below |
| an absolute bound |A_0|^2 <= M(a) does not establish the comparison | K = 2M/lambda may diverge at a first null | correct |

## 2. The exact separation of the polar term (PROVED; the user's calculation of 2026-10-03, verified)

On the admissible space, Q_a(f) = <f, H_a f> - |<v_a, f>|^2 with H_a the prime-side operator and v_a = sqrt(2) sinh(x/2).
HYPOTHESIS (H+): H_a >= h_a I > 0 at the support in question. Under (H+):

(P1) Decomposition. With y = H_a^{1/2} f, w = H_a^{-1/2} v_a, m_a = ||w||^2 = <v_a, H_a^{-1} v_a>, delta_a = 1 - m_a,
        Q_a(f) = || y - (<w,y>/m_a) w ||^2 + (delta_a/m_a) |<v_a, f>|^2 .
     Hence Q_a >= 0 <=> m_a <= 1; and if delta_a > 0 then Q_a(f) >= delta_a <f, H_a f> >= h_a delta_a ||f||^2, so
        lambda(a) >= h_a delta_a .
     The full prime interaction stays inside H_a^{-1}; nothing is discarded.
(P2) Secular equation. If 0 < lambda(a) < h_a and <v_a, f_a> != 0 (not the polar branch), then with phi = (H_a - lambda)^{-1} v_a,
        f_a = phi/||phi||,   m(lambda(a), a) = <v_a, phi> = 1 ,
     and lambda(a) is the lowest root of m(., a) = 1 (m increases in lambda from 0 to +infinity on (-infinity, h_a)).
(P3) The margin through the resolvent identity.
        delta_a = m(lambda, a) - m(0, a) = lambda <v_a, H_a^{-1} (H_a - lambda)^{-1} v_a>   at lambda = lambda(a),
     and by spectral calculus (H_a and its resolvent commute, mu_j >= mu_j - lambda > 0)
        delta_a <= lambda(a) ||phi||^2 .
(P4) The sufficient boundary estimate. With |A_0|^2 = |alpha|^2/||phi||^2,
        2|alpha(a)|^2 <= K(a) delta_a   implies   L(a) = 2|A_0(a)|^2 <= K(a) delta_a/||phi||^2 <= K(a) lambda(a),
     the prime-side comparison. TARGET (HYPOTHESIS): 2|alpha(a)|^2 <= K(a) (1 - <v_a, H_a^{-1} v_a>) with K integrable.

Verification on the K-mode form (600 bits; K = 48, 56, 60, 64, 72; all finite-tree supports, a < a_inf = 0.843 except
the last, which is at the edge of that range):

      a      h_a          m(lambda*)         delta_a        ||phi||^2   h_a delta_a   lambda*       delta/(lambda ||phi||^2)
      0.6    3.07e-3      1.0000000000000    2.475e-5       41.50       7.6e-8        5.966e-7      0.99987
      0.7    6.54e-6      1.0000000000000    9.546e-9       37.34       6.2e-14       2.557e-10     0.99998
      0.75   9.18e-8      1.0000000000000    1.209e-10      35.82       1.1e-17       3.377e-12     0.99998
      0.8    9.84e-10     1.0000000000000    5.411e-13      34.55       5.3e-22       1.566e-14     0.99999
      0.84   1.61e-11     1.0000000000000    5.304e-15      33.70       8.5e-26       1.574e-16     0.99999

So (P2) holds to 14 digits, (P3) is nearly an equality (the slack 1 - delta/(lambda ||phi||^2) is 1e-4 to 6e-6), and the
lower bound h_a delta_a of (P1) is weak by a factor 8 at a = 0.6 and 1e9 at a = 0.84, because h_a itself collapses
(Section 4).

## 3. Does the margin control the boundary coefficient? (the next calculation)

(P5) Hadamard formula of the resolvent (PROVED for the K-mode form by implicit differentiation; for the exact form wherever
lambda is differentiable and m is differentiable in a; verified to 10 digits). Differentiating m(lambda(a), a) = 1,
        d_a m(lambda, a)|_{lambda = lambda(a)} = - lambda'(a) ||phi||^2,      d_a m = - <phi, H_a' phi> + 2 <phi, v_a'>,
and with the derivative identity lambda' = -2|A_0|^2 = -2|alpha|^2/||phi||^2:
        2|alpha(a)|^2 = d_a m(lambda, a)|_{lambda = lambda(a)} .
The squared boundary coefficient of the unit-source solution is the a-derivative of the secular function at fixed spectral
parameter. Both sides of the target are therefore derivatives of the one function m(lambda, a):
        TARGET  <=>  d_a m(lambda*, a) <= K(a) [ m(lambda*, a) - m(0, a) ] = K(a) int_0^{lambda*} d_lambda m(lambda, a) dlambda ,
and since d_lambda m = ||phi(lambda)||^2 is nearly constant on [0, lambda*] (table: 41.486 against 41.497 at a = 0.6),
this is the slope condition of the level curve m = 1: d log lambda*/da >= -K(a) (1 + o(1)).

Measured (same runs):

      a      2|alpha|^2 = d_a m(lambda*)   2|alpha|^2/delta_a   Phi' = -lambda'/lambda   2|alpha_0|^2 = d_a m(0)   2|alpha_0|^2/delta_a
      0.6    1.7245e-3                     69.666               69.657                   1.7550e-3                 70.901
      0.7    8.6751e-7                     90.880               90.877                   8.7545e-7                 91.712
      0.75   1.1826e-8                     97.777               97.775                   1.1932e-8                 98.656
      0.8    5.5972e-11                    103.446              103.445                  5.6297e-11                104.046
      0.84   6.4357e-13                    121.336              121.335                  6.4660e-13                121.907

Answer to the question. The margin controls the boundary coefficient with the constant K(a) = Phi'(a) = pi kappa(a) T*(a)
to one part in 10^4, and with nothing to spare: the chain (P4) loses only the factor delta/(lambda ||phi||^2) = 0.9999.
The target is the prime-side comparison itself, written on the resolvent at the polar source; it is not weaker. What the
separation buys is a change of object, not a change of difficulty.

(P6) The lambda = 0 reduction (HYPOTHESIS for the Hadamard step). The unit-source solution at lambda = 0, phi_a(0) =
H_a^{-1} v_a, needs no eigenvalue and no secular equation. Its boundary coefficient alpha_0 differs from alpha by
0.5 to 1.8 per cent at the five supports (ratio 1.018, 1.009, 1.009, 1.006, 1.005, decreasing), since
||phi(lambda*) - phi(0)|| <= lambda* ||phi||/h_a. If the Hadamard formula holds for phi_a(0),
        delta_a'(a) = - d_a m(0, a) = - 2|alpha_0(a)|^2        (HYPOTHESIS: Hadamard at lambda = 0),
then the entire deduction runs on the margin alone: delta_a > 0 <=> lambda(a) > 0 under (H+), and
        2|alpha_0(a)|^2 <= K(a) delta_a   gives   delta_a(b) >= delta_a(a_0) exp( - int_{a_0}^{b} K ) > 0 .
The margin's measured logarithmic rate is 70.9 and 104.4 at a = 0.6 and 0.8 (central differences over +-0.002), equal to
2|alpha_0|^2/delta_a above: consistent with (P6).

(P7) What the finite-tree representation gives (DEFINITIONS of the pieces; the inequality remains the HYPOTHESIS).
Below a_inf the effective symbol sigma_eff = sigma~_inf - R is explicit (Proposition prop:tree(iii), a finite continued
fraction along each chain), with Wiener-Hopf factors sigma_+-. The unit-source solution phi has the representation of
Corollary cor:edgeform with the polar source in place of the minimizer's self-consistent source: alpha is the constant
term of the amplitude A(t) = t P_-[ (H^_- - e^{-2ita} Phi_-)/sigma_+ ](t) - t P_+[ e^{-2ita} Phi_-/sigma_+ ](t), a LINEAR
functional of the source transform through 1/sigma_+; and m(0, a) = <v_a, phi_a(0)> = (1/2 pi) int conj(V_a) Phi^_a dt is
the PAIRING of the source with its solution, a quadratic functional of the source through the resolvent. So the target
compares the square of a residue-at-infinity functional of (sigma_+, V_a) with one minus a pairing functional of
(sigma_eff, V_a). Writing both out explicitly for one entry (a_2 < a < a_3, sigma_eff = sigma~_inf - c_2^2/(sigma~_inf))
and for the tree with the entries 2, 3, 4 is the next calculation; no inequality between the two functionals is known.

## 4. A finding on the prime-side operator (COMPUTATION)

The bottom h_a of H_a itself, with no polar term, collapses doubly exponentially: 3.1e-3, 6.5e-6, 9.2e-8, 9.8e-10, 1.6e-11
at a = 0.6, 0.7, 0.75, 0.8, 0.84 (second eigenvalue 0.066 and 4.7e-5 at 0.6 and 0.8), with lambda*/h_a = 1.9e-4,
3.9e-5, 3.7e-5, 1.6e-5, 9.8e-6. The ground state e_0 of H_a is nearly orthogonal to the polar source: <e_0, v>^2/h_a is
0.085 at a = 0.6 and 2.2e-8 at a = 0.8, and m_a = 1 - delta_a is built almost entirely from the higher modes. The
cancellation is two-stage: the archimedean part against the entries first (h_a tiny), then the polar source tuned to the
near-null space of H_a (delta_a tinier by 1e-4 to 1e-5 relative).

Remark (QUESTION, not a claim). (H+), the positivity of H_a, is the polar-free half of Weil positivity: by the explicit
formula <f, H_a f> = sum_rho g^(rho) + |<v_a, f>|^2 (the polar term moved to the right), so under RH it is a sum of
nonnegative terms plus the polar square and (H+) holds with h_a >= lambda(a). Whether (H+) can be proved without RH is
open and is not known to be equivalent to RH: a single off-line pair at displacement m contributes a term of size
m^2 e^{2am} in the window of support a, while the polar square can be as large as c e^{a} ||f||^2, and every nontrivial
zero has m < 1/2. If (H+) were a theorem, the Hypothesis beyond a_Z would be exactly the scalar statement
m_a = <v_a, H_a^{-1} v_a> <= 1 for every a, and the comparison would be the statement (P6) about the one function
delta_a(a) = 1 - m_a(a): its logarithmic derivative is bounded by an integrable K.

## 5. Status

PROVED: P1-P5 (P5 for the K-mode form; for the exact form under the derivative identity and differentiability of m).
HYPOTHESIS: (H+) where used; the Hadamard formula at lambda = 0 (P6); the TARGET 2|alpha|^2 <= K delta_a, equivalently
the prime-side comparison. NEXT: the explicit one-entry and three-entry forms of alpha and m(0,a) through sigma_eff and its
factors (P7), to see what structure, if any, bounds the residue functional by the pairing deficit.

## 6. The source transfer (the user's calculation of 2026-10-03, second part; verified)

DEFINITION. u_z = (H_a - zI)^{-1} v_a for z < h_a; u_0 = H_a^{-1} v_a; u_lambda at the root lambda = lambda(a).

(P8) Exact identities (PROVED, resolvent identity and (P3)):
        u_lambda = u_0 + lambda (H_a - lambda)^{-1} u_0,         delta_a = lambda <u_0, u_lambda> .
     The second is (P3) written on u_0: delta = lambda <v, H^{-1}(H-lambda)^{-1} v> = lambda <H^{-1} v, (H-lambda)^{-1} v>.

(P9) Transfer of the boundary coefficient (PROVED where the boundary representation applies to u_0, u_lambda and
     (H_a - lambda)^{-1} u_0 with a common linear boundary-coefficient functional):
        alpha_lambda = alpha_0 + lambda beta_lambda,        beta_lambda = boundary coefficient of (H_a - lambda)^{-1} u_0 .

(P10) From source bounds to the comparison (PROVED, given (P8)-(P9)). Suppose
        (S1)  |alpha_0(a)|^2 <= M(a) delta_a,          (S2)  |beta_lambda(a)| <= D(a),
     with M, D >= 0 locally bounded. Then, since ||u_lambda|| >= ||u_0|| (spectral calculus) and
     ||u_0|| >= <v_a, u_0>/||v_a|| = m_a/||v_a||,
        L(a) = 2|alpha_lambda|^2/||u_lambda||^2 <= 4(|alpha_0|^2 + lambda^2 |beta_lambda|^2)/||u_lambda||^2
             <= 4 M delta_a/||u_lambda||^2 + 4 lambda^2 D^2 ||v_a||^2/m_a^2 ,
     and delta_a = lambda <u_0, u_lambda> <= lambda ||u_0|| ||u_lambda|| <= lambda ||u_lambda||^2, so
        L(a) <= K(a) lambda(a),        K(a) = 4 M(a) + 4 lambda(a) D(a)^2 ||v_a||^2 / m_a^2 ,
     which is locally bounded (lambda <= lambda(a_Z) on the living range, m_a = 1 - delta_a >= 1/2 once delta_a <= 1/2,
     ||v_a||^2 = sinh(a) - a... explicit). This is the user's explicit coefficient up to the constants chosen here.

What (S1) and (S2) are.
  (S2) is an ABSOLUTE bound on the boundary coefficient of a resolvent applied to the smooth source u_0: the kind of
       bound the source norms of Corollary cor:edgeform are expected to give (paper, Lemma lem:A0bound, Proposition 9.15),
       locally uniform in a away from the entries. It carries no cancellation and is the accessible half.
  (S1) is RELATIVE: the squared boundary coefficient of H_a^{-1} v_a against the margin. By the Hadamard formula at
       lambda = 0 (P6, HYPOTHESIS), 2|alpha_0|^2 = d_a m(0,a) = -delta_a'(a), so (S1) is  -delta_a' <= 2M delta_a : the
       Gronwall statement for the margin itself. It carries the whole cancellation; it is the comparison in its lambda = 0
       form, and nothing in (P8)-(P10) weakens it. The measured values: 2|alpha_0|^2/delta_a = 70.9, 91.7, 98.7, 104.0,
       121.9 at a = 0.6, 0.7, 0.75, 0.8, 0.84 (Section 3), so M(a) must be at least Phi'(a)/2 (1 + 1%).

(P11) An integrated consequence (PROVED under Q_a >= 0 for all a, i.e. under RH, with (P6)). Since m(0,a) is
     non-decreasing in a and m(0,a) <= 1,
        int_a^infinity 2|alpha_0(a')|^2 da' = lim m - m(0,a) <= delta_a .
     The far-edge cancellation that (S1) asks for pointwise is forced, integrated over all larger supports, by positivity
     alone. A proof of (S1) must supply something that holds at every single support; positivity supplies it only on
     average, and the average is useless because delta decays doubly exponentially.

## 7. The structure of the cancellation (the explicit near-edge functional)

In the dictionary convention F(t) = int f(x) e^{-itx} dx, the polar source v_a = sqrt2 sinh(x/2) on [-a,a] has
        V_a(t) = sqrt2 [ e^{-ita} Sigma_R(t) + e^{+ita} Sigma_L(t) ],
        Sigma_R(t) = ( e^{a/2}/(1/2 - it) + e^{-a/2}/(1/2 + it) )/2,     Sigma_L(t) = -( e^{-a/2}/(1/2 - it) + e^{a/2}/(1/2 + it) )/2,
with Sigma_L(t) = -Sigma_R(-t) (oddness). In the right-edge frame (multiply by e^{ita}) the source is Sigma_R + e^{2ita} Sigma_L:
a near-edge part of size e^{a/2} and a far-edge part of the same size carrying the phase e^{2ita}. The exact
representation (paper, Proposition prop:whexact for the archimedean operator; Proposition prop:tree with sigma_eff below
a_inf) writes the solution's transform in that frame as P_-[H^_-/sigma_+]/sigma_-, where H^_- is the source plus the
force outside the FAR edge, with the real-zero corrections at +-t_lambda; and the boundary coefficient is the residue at
infinity, alpha = lim t P_-[ (H^_- - far)/sigma_+ ](t) - lim t P_+[ far/sigma_+ ](t) (Corollary cor:edgeform), that is,
up to the slowly varying factor 1/sigma_+ ~ (log t)^{-1/2}, a Cauchy integral of the source against 1/sigma_+.

So alpha_0 = ell_near + ell_far with
        ell_near = -(1/2 pi i) int Sigma_R(u)/sigma_+(u) du  (regularized at +-t_lambda)  =  O(e^{a/2}) ,
and ell_far the contribution of e^{2ita} Sigma_L and of the outside force beyond the far edge. For the odd problem the
force beyond the far edge is the mirror image of the force beyond the near edge, g_out(-a-v) = -g_out(a+v), and the near
force is itself determined by the solution: the far-edge term is not an independent datum but a self-consistency term,
a scalar equation for alpha_0 of the form
        alpha_0 = ell_near(a, lambda) + rho(a, lambda) conj(alpha_0) + (smaller),
with rho the return of the near edge through the far edge (the mirror chain of the paper, Section 9: the walk across the
window and back, weight O(1) through the real zeros +-t_lambda, not exponentially small, because the symbol has real
zeros). (S1) says |alpha_0|^2 <= M delta_a with delta_a doubly exponentially small while ell_near is e^{a/2}: the
cancellation between ell_near and rho conj(alpha_0) + (smaller) must be exact to a relative e^{-Phi/2 - a/2}. That is the
statement "the source representation retains opposite-edge feedback; its cancellation has not been established", made
explicit.

NEXT (analytical, no computation): write the scalar equation for alpha_0 exactly for the archimedean operator (prime-free,
Proposition prop:whexact, where sigma_+- are explicit through Lemma lem:factor) with the polar source, isolating rho and
the remainder; then with one entry (sigma_eff = sigma~_inf - c_2^2/sigma~_inf). The question is whether the equation,
combined with m(0,a) = <v_a, u_0> written in the same data, forces |alpha_0|^2 <= M (1 - m(0,a)) by an identity rather
than by a cancellation one has to prove separately. Nothing below a_inf involves the zeros.

## 8. Audit of the updated addendum (received/Prime_Side_Completion_Addendum_v2.pdf, Sections 8-10)

| statement | check | status |
|---|---|---|
| Proposition 4 (polar decomposition; Q >= 0 iff m <= 1; Q >= delta <f,Hf> >= h delta ||f||^2) | identical to (P1); necessity of m <= 1 by f = H^{-1}v with Q = m(1-m) | PROVED |
| (8.1): u_lambda = u_0 + lambda R_lambda u_0, delta = lambda <u_0, u_lambda>, the inner product real positive | resolvent identity; positivity by spectral calculus | PROVED |
| Proposition 5: delta <= lambda N <= (h/(h - lambda)) delta, delta >= lambda ||u_0||^2 | spectral measure of H at v; numerically delta/(lambda N) = 0.99987, 0.99998, 0.99998, 0.999990, 0.999994 against 1/(1 + lambda/h) = 0.99981, 0.99996, 0.99996, 0.999984, 0.999990 at a = 0.6..0.84 (within the bound), and delta = lambda ||u_0||^2 to 5 digits | PROVED; sharp here because lambda << h |
| (9.1) alpha_lambda = alpha_0 + lambda beta_lambda | linearity of the boundary map on its justified domain | PROVED under the domain hypothesis |
| Proposition 6: 2|alpha_lambda|^2 <= (4M + 4 D^2 lambda/n) delta, hence L <= K lambda | uses lambda^2 <= lambda delta/n from Proposition 5 and delta <= lambda N; cleaner than the K of (P10) | PROVED under (9.2) |
| Section 10: (10.1) |B H^{-1} v|^2 <= M (1 - <v, H^{-1} v>) is the target; the far-edge feedback depends on the solution; ker H = ker Q cap v^perp when Q >= 0 | correct; the kernel statement is immediate from <f,Hf> = Q(f) + |<v,f>|^2 | correct |

Nothing to correct. The two source bounds (9.2) are exactly (S1), (S2) of Section 6; the assessment there stands: (S2) is
absolute and accessible, (S1) = (10.1) carries the cancellation.

## 9. Where the cancellation lives: not in the archimedean operator

[computation already in the repository, data/dilation_a*_primefree.log, K = 48] The prime-free floor, the bottom of
K_infinity - |v><v| (archimedean part minus the polar term, no entries), is
        lambda_pf(0.3) = +0.2226,   lambda_pf(0.4) = -0.0780,   lambda_pf(0.5) = -0.3222,   lambda_pf(0.6) = -0.5359 ,
negative from about a = 0.35 on, while the full floor is +0.0147, +1.9e-4, +6.0e-7 at a = 0.4, 0.5, 0.6. So:
  (i) positivity of Q_a, and of H_a itself beyond a = 0.35, is a PRIME effect: at the minimizer the autocorrelation at the
      entries is negative and the term -2 sum c_n g(log n) is positive, lifting an indefinite archimedean form to a
      doubly-exponentially small positive floor;
 (ii) the archimedean-only scalar equation for alpha_0 (the NEXT of Section 7) is not the right first step: in the
      prime-free problem the margin 1 - <v, K_infinity^{-1} v> is negative beyond 0.35 and there is nothing to cancel;
(iii) the first informative case is ONE ENTRY, a_2 = 0.347 < a < a_3 = 0.549, where H_a = K_infinity - c_2 (tau_{log 2} +
      tau_{-log 2}), c_2 = (log 2)/sqrt2, the effective symbol is sigma_eff = sigma~_inf - c_2^2/sigma~_inf (Corollary
      cor:full), and the full floor already falls from 0.0147 to 1.9e-4 across the interval.

The forced problem in this case. (H_a - lambda) u = v with v = sqrt2 sinh(x/2): the representation of Corollary cor:full
applies verbatim with the minimizer's self-consistent polar coefficient 2P replaced by the fixed coefficient sqrt2
(Proposition prop:whfull is proved for sources in L^1 cap L^2 with y h in L^2, which v is). In its proof the near source
of H^'_- is the step -2iP sinh(a/2) - i c f(a - d) at the edge plus a smooth part, and "the far feedback is again
-i r(0^+) + O(1/t)". For the forced problem the step is -i sqrt2 sinh(a/2) - i c_2 u(a - log 2). So

        alpha = [ Wiener-Hopf functional of the step  -i sqrt2 sinh(a/2) - i c_2 u(a - log 2) ]  -  i r(0^+)  +  (rational corrections at +-t_lambda),

where r is the regular part of the force at the edge, determined by the solution u itself (the paper's Lemma lem:cont:
the regular part of K u is continuous across the edge, so for the minimizer g_reg(a) = 2P sinh(a/2), and for the forced
problem g_reg(a) = (H u)(a^-) + lambda u(a) = v(a) = sqrt2 sinh(a/2) since u(a) = 0). This is the exact form of the
"opposite-edge feedback": the regular part of the force at the near edge equals the source value there, by continuity,
and the boundary coefficient is what is left of the near-source step after that return is subtracted, together with the
echo term c_2 u(a - log 2) and the rational corrections. The cancellation that (10.1) requires is therefore between the
step of the source at the edge, the value of the solution at the echo point a - log 2 (one step inside), and the return
through the real zeros of sigma_eff. All three are explicit in the one-entry case.

NEXT (replacing the NEXT of Section 7): (a) write, for a_2 < a < a_3, the Wiener-Hopf solution of (H_a - lambda) u = v in
the right-edge frame with sigma_eff and its factors sigma_eff,+-, the near-source step -i sqrt2 sinh(a/2) - i c_2
u(a - log 2), the far feedback -i r(0^+) with r(0^+) = sqrt2 sinh(a/2) by continuity, and the rational corrections at
+-t_lambda; (b) read off alpha as the residue at infinity and m(0,a) = <v, u_0> as the pairing; (c) express u(a - log 2),
the only non-explicit datum, through the echo relation (eq:localecho with one echo: u near a - log 2 equals c_2 (k * u(a
+ .)) plus smooth); (d) check whether the resulting scalar relation for alpha_0 yields |alpha_0|^2 <= M (1 - m) with M
locally bounded on (a_2, a_3). Everything in (a)-(c) is in the paper's Section 8 for the minimizer; the forced problem
changes only the coefficient of the step.

## 10. Audit of the third version (received/Prime_Side_Completion_Addendum_v3.pdf, Sections 11-12)

| statement | check | status |
|---|---|---|
| (11.1): Q u_0 = delta v and Q(u_0) = m delta, with 0 < m < 1 at a positive secular root | H u_0 = v gives (H - v v*) u_0 = (1 - m) v; pair with u_0; m(0) < m(lambda) = 1 since m increases in lambda and lambda > 0 | PROVED, exact; no boundary input |
| (11.2): E_>(u;T) + E_<=(u;T) = m delta with the polar subtraction -m^2 booked in E_<= | (1/2pi) int Psi_a |U|^2 = <u, H u> = m; subtract m^2 | PROVED (bookkeeping) |
| Proposition 7: (11.3) E_> >= |alpha_0|^2/(pi T) and (11.4) E_<= >= -C m delta give |alpha_0|^2 <= pi T (1 + C) m delta | immediate | PROVED under (11.3), (11.4) |
| Proposition 8: |beta_lambda|^2 <= pi T V^2 h^{-2} ( (h - lambda)^{-1} + J (h - lambda)^{-2} ) under the forced tail bound (12.1) | pairing (H - lambda) w = u with w; low part >= -J ||w||^2 by Parseval; ||u|| <= V/h, ||w|| <= V/(h(h - lambda)) | PROVED under (12.1) and (H+) |
| (12.3): K = 4 pi T (1 + C) m + 4 lambda D^2/||H^{-1} v||^2; local boundedness under h_0, m >= 1/2 | from Propositions 6-8 | PROVED under the hypotheses |

Nothing to correct. Two remarks on content.

(R1) The identity Q u_0 = delta v is the sharpest exact statement so far: the zero-parameter forced solution is a near-null
vector of the window form with residual delta v, and its Weil energy is m delta. Everything about the comparison is a
statement about how a near-null vector of Q with residual delta v distributes its energy in frequency.

(R2) What (11.3) and (11.4) are, against the paper. The tail identity of Proposition prop:tailid and the sum rule (Lemmas
lem:overlap, lem:cont, Proposition prop:sumrule) are proved for a critical point K f = lambda f + 2 P s + g_out on the
line, and never use the self-consistency P = <f, s>; the Wiener-Hopf representation (Propositions prop:whexact,
prop:whfull, prop:tree) is proved for sources in L^1 cap L^2 with y h in L^2, which the polar source is. So, below
a_inf, the forced solution u_0 has the edge law with two derivatives and the asymptotic tail law by the same proofs with
the step coefficient 2P replaced by sqrt2 (status: expected from the same proofs; to be re-verified line by line before it
is called a theorem):
        E_>(u_0; T) = 2 |alpha_0|^2 / (pi T) (1 + o(1)),   T -> infinity,
a factor 2 stronger than (11.3). But then (11.4) at a cutoff T, together with this asymptotic, says exactly
        2 |alpha_0|^2 / (pi T) (1 + o(1)) <= (1 + C) m delta        at T = T_0(a),
which is the comparison with K(a) = pi T_0(a) (1 + C): the obligation (11.4) is the statement that the tail law has set
in by the height T_0(a) with its error controlled, i.e. the Fourier-tail onset for the forced solution, locally uniform
in a. For T beyond 2 |alpha_0|^2/(pi m delta) the inequality (11.4) with C = 0 is automatic and Proposition 7 returns
|alpha_0|^2 <= 2 |alpha_0|^2. So Sections 11-12 fold the three gates into one: the onset T_0(a) of the tail law for the
forced solution u_0 = H_a^{-1} v_a at a locally bounded multiple of the horizon, with E_<= not more negative than a
bounded multiple of the whole. That is item N1d.1 of the tree (the local law of the band) for the forced solution.

(R3) Where the sign of E_<= comes from. E_<= = (1/2pi) int_{|t| <= T} Psi_a |U|^2 - m^2, and Psi_a is negative on a set
reaching from t = 0 (Psi_infinity(0) = -5.4) up to the last sign change (33, 259, 397 horizons at a = 1, 1.25, 1.5 in the
paper's Computation comp:band; below the horizon for a < a_inf). (11.4) therefore asks that the forced solution's
transform carry little mass on the negative set of the symbol relative to m delta: the forced-solution form of the
observation that the minimizer's transform hides from the negative set of the symbol as it hides from the zeros below
the horizon (paper, discussion after Proposition prop:window; the fiction's "darkness"). In the finite-tree regime the
negative set lies below the horizon and is explicit: {t : sigma~_inf(t) < X_max}, with X_max = 2.25 at a = 0.6.

NEXT, analytical: in the one-entry regime, (a) carry the tail identity through for u_0 (step coefficient sqrt2) and
isolate the onset, i.e. the first T at which the two order-C L^{-1/2}/T pieces (the polar piece and the smooth overlap)
cancel to within a bounded multiple of 2|alpha_0|^2/(pi T), using the continuity lemma g_reg(a) = sqrt2 sinh(a/2) for the
forced problem; (b) bound E_<= from below on the explicit negative set of sigma_eff by the representation of U below the
horizon (Proposition prop:tree gives U = P_-[H^''_-/sigma_eff,+]/sigma_eff,- with all sources explicit). If (a) and (b)
give T_0(a) <= c T*(a) and C(a) bounded on (a_2, a_3), the comparison is a theorem there, and the first null is excluded
on (a_2, a_3) by Proposition 7 and the closing deduction, without zeros.

## 11. Audit of the fourth version (received/Prime_Side_Completion_Addendum_v4.pdf, Sections 13-15)

| statement | check | status |
|---|---|---|
| Lemma 9: for f supported in [-A, A], (1/2pi) int_{|t|>T} |f^|^2 >= eta(A,T) ||f||^2 with eta = 1 - ||K_{A,T}|| > 0 | the sinc operator on [-A,A] is compact, positive, a contraction; norm one would give a band-limited f of compact support, hence f = 0 (Paley-Wiener); this is the top Slepian eigenvalue lambda_0(c) < 1, c = AT | PROVED (classical) |
| (13.2): Psi_a >= Psi_inf - B(A), and beyond T_A with Psi_inf >= B(A) + 1 the high energy obeys E_> >= eta ||f||^2 | immediate | PROVED |
| Theorem 10: (14.1) E_<= >= -C_A Q_a(f) for all admissible odd f, a <= A, gives lambda(a) >= eta/(1 + C_A) | algebra | PROVED under (14.1) |
| Corollary 11: E_<=(a,u) >= -C m delta at a cutoff with (13.2) gives delta >= eta/((1 + C) V^2 + eta) | (1+C) m delta >= E_>(u) >= eta ||u||^2 >= eta m^2/V^2 | PROVED under its hypothesis |
| Section 15: a large support-dependent cutoff is not a locally controlled estimate; the crude absolute bound gives a coefficient with no proof that it is below one | correct | correct |

Nothing to correct. Three remarks on content, the last of which connects the route to the paper.

(R4) Scope of (14.1). With Lemma 9 and the crude bound E_<= >= -(J_A + V_A^2) ||f||^2, the hypothesis (14.1) for all f
and a <= A is EQUIVALENT to uniform coercivity lambda(a) >= c_A > 0 on (0, A]: one direction is Theorem 10; conversely
E_> = Q - E_<= <= Q + (J_A + V_A^2) ||f||^2 <= (1 + (J_A + V_A^2)/c_A) Q. So (14.1) is positivity with a margin on the
range, restated; it bypasses the derivative and the regularity but not the content.

(R5) The fixed cutoff is incompatible with the data. Corollary 11 gives a lower bound for the margin, delta >=
eta/((1+C)V^2 + eta), uniform in a <= A, whereas delta is doubly exponentially small (2.5e-5, 5.4e-13, 5.3e-15 at a =
0.6, 0.8, 0.84). Hence at a fixed cutoff T_A the constant must satisfy C >= eta (1 - delta)/(V^2 delta) - 1, which is
5e3, 1e11, 9e12 at those supports for eta = 0.01: the signed comparison with bounded C is false at any fixed cutoff. The
forced solution's energy above T_A is of order one while its total Weil energy is m delta; the low part is of order
minus one. The only cutoff at which C can be bounded is support-dependent, T_0(a) ~ 2|alpha_0|^2/(pi m delta) = kappa T*(a)
(Section 10), the horizon scale.

(R6) At the horizon cutoff the route reproduces the cooling form, and it is the paper's prolate count. With T = kappa
T*(a) and A = a, Slepian's asymptotic for the top eigenvalue, 1 - lambda_0(c) ~ 4 sqrt(pi c) e^{-2c}, c = a T, gives
-log eta(a, kappa T*) ~ 2 a kappa T*(a) = 4 pi kappa a e^{2a}, against the measured Phi = -log lambda:

      a      2 a kappa T*    Phi     ratio   (4a/pi)
      0.6        26.6        14.3    1.86     0.76
      0.8        52.7        31.8    1.66     1.02
      1.0       108.1        59.5    1.82     1.27
      1.25      229.8       116.2    1.98     1.59
      1.5       469.1       209.6    2.24     1.91
      1.75      926.3       372.3    2.49     2.23

So Theorem 10 at the horizon cutoff with bounded C would give lambda(a) >= exp(-4 pi kappa a e^{2a} (1 + o(1))), the
Cooling Theorem's doubly exponential form, weaker than the truth by the factor exp(2c - Phi) and consistent with C = 0;
this is the "prolate-type count" of the paper's Section 7, which "bounds Phi below by about half of what is observed"
(the ratio column). The obligation (14.1) at that cutoff, for the minimizer, is theta(kappa T*) <= 1 + C, true with C = 0
under the Born rule; for all f it is coercivity (R4). The route therefore sits exactly where the others sit: the local
law of the band at the horizon, now in the form "the window form dominates its own energy above kappa T*, for every odd
f of support a", with a coefficient 1 + C locally bounded. Nothing is lost and nothing is gained relative to Section 10,
except the form: no derivative identity, no absolute continuity, and no boundary asymptotic are needed on this route,
which is a real simplification of the deduction if the all-function estimate can be approached directly.

NEXT, unchanged: the one-entry regime, where the negative set of the symbol is explicit and below the horizon, and
where the all-function estimate at T = c T* can be attacked on the explicit representation of U below the horizon.

## 12. Audit of the fifth version (received/Prime_Side_Completion_Addendum_v5.pdf, Sections 16-17)

| statement | check | status |
|---|---|---|
| Proposition 12: the fixed-cutoff criterion (14.1) for all odd f, a <= A, is equivalent to lambda(A) > 0 | the converse uses nesting of the windows and the crude bound with C_A = (J_A + V_A^2)/gamma_A; this is (R4) of Section 11 | PROVED |
| (17.1): 0 <= v - s_{p,r} = N - u <= v^2/u for both path-norm cases | rationalize sqrt(u^2+v^2) - u and (sqrt(u^2+4v^2) - u)/2; denominators >= 2u | PROVED (the path-norm formulas themselves come from [T, CP20], not in this repository) |
| Proposition 13: S_sq(a) = a + O(1) | sum_p sum_r v^2/u = sum_p (log p) p^{-3/2}/(1 - p^{-3/2}) < infinity; r >= 2 terms bounded; sum_{p<e^a} log p/p = a + O(1) (Mertens; the PNT with its classical error gives the O(1) by partial summation) | PROVED; checked numerically: sum - a = 1.33, 1.33, 1.33 (constant) at a = 4, 6, 8 |
| (17.3): B_ref(a) = sum_{p^k < e^{2a}} log p / p^{k/2} ~ 2 e^a | partial summation from psi(X) ~ X; checked: ratio to 2e^a is 1.09, 1.03, 1.01 at a = 4, 6, 8 | PROVED |
| consequence: the scalar saving is a/(2e^a) of the reflection charge; the sign of the full form is decided by the retained positive path matrices, the continuous reserve and their interactions | correct as stated | correct |

Remark. B_ref(a) is the Schur bound B(a) of the addendum's Lemma 1 up to the factor 2 (the band top 2 sum c_d of the
paper); so Section 17 says that the prime-square reflection pairing of [T] saves O(a) out of a charge of order e^a:
nothing at the exponential scale. That is consistent with everything above: positivity is a cancellation to doubly
exponential precision between the archimedean part and the entries at the minimizer, not a bookkeeping of the charges.
[T] (Michaels_Consolidated_Theta_Local_Amplitude_Proof 2.pdf, CP20) is not in this repository; its formulas are taken as
quoted.

Status after the five versions. Proved and recorded: the polar decomposition, the source transfer, the quantitative
equivalence of the margins, the small-energy identity Q u_0 = delta v, the fixed-cutoff equivalence, the Mertens-scale
size of the reflection saving. Open and the same throughout: the all-support signed comparison, in any of its equivalent
forms (the prime-side comparison L <= K lambda; (10.1); (11.4) at a horizon cutoff; (14.1) on each window), which is
Weil positivity with a margin on the range in question. The structural facts established on the way (Sections 4, 9, 11:
the collapse of h_a, the prime-free form negative beyond 0.35, the natural cutoff at the horizon with the prolate scale)
locate where the content sits: the one-entry regime, below the horizon, on the explicit negative set of sigma_eff.
