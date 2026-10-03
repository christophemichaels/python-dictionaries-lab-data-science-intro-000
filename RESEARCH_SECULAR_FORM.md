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
