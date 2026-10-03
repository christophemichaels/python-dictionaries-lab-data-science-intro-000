# The prime-side form at a finite-tree support: the secular function of the polar source

Research note, 2026-10-03. Real mathematics; status of every statement marked. Companion to RH_IF_THEN_TREE.md, item
N1d.2, and to received/Prime_Side_Completion_Addendum.pdf (whose Section 6 states the target this note rewrites).

## 1. Audit of the addendum against the paper

| addendum | paper | status |
|---|---|---|
| Lemma 1, the Schur bound ||C_a|| <= B(a) = 2 sum_{n<e^{2a}} Lambda(n)/sqrt(n), resolvent cutoff S(a) | the band top 2 sum c_d of Computation comp:band (data/band_a*.log: "band top sigma_inf = 2 sum c") | known; elementary; correct |
| Lemma 2, bounded variation of f gives a Lipschitz autocorrelation, hence (H_fin) | the mechanism of Proposition 9.15 under (H_BV) | known; correct |
| the Cantor-function caution: an a.e. derivative identity does not integrate without absolute continuity | why Proposition 9.17 needs (H_fin) on the exceptional null set | correct; this is Step 1's remainder |
| Theorem 3, the first-zero argument with K integrable through every prospective zero | Proposition prop:AimpliesRH; Proposition 9.30 | known; correct |
| the measure inequality nu(E) <= int_E K lambda as a combined target | new formulation of (regularity + comparison) | useful; equivalent to AC plus (d) with integrable C |
| Section 6, the finite-tree target 2|A_0^(1)|^2 <= K lambda ||phi_1||^2 at the root of the secular equation | the paper's "prime-side form" paragraph after Theorem thm:reduction | the same statement; rewritten below as a level-set inequality |
| "an absolute bound |A_0|^2 <= M(a) does not establish it" | correct: K = 2M/lambda can diverge at a first zero | correct |

Nothing in the addendum is false; nothing in it goes beyond the paper except the two reformulations, which are adopted.

## 2. The secular function

Setting (odd sector, window [-a,a], the paper's normalization). Write the window form as

    Q_a(f) = <f, K_a f> - 2 <f, s_a>^2,        s_a(x) = sinh(x/2) on [-a,a],

where K_a is the prime-side operator (archimedean part minus the entries below e^{2a}; no polar term) and the polar term is
the rank-one -2 s s^T. Let mu_0(a) <= mu_1(a) <= ... be the eigenvalues of K_a on the window (discrete; Lemma lem:exist).
Define the SECULAR FUNCTION

    sigma(lambda, a) = < (K_a - lambda)^{-1} s_a , s_a >,        lambda < mu_0(a),

increasing in lambda from 0 (at -infinity) to +infinity (at mu_0^-), with d sigma / d lambda = ||phi||^2, phi = (K_a - lambda)^{-1} s_a.

(T1) [theorem, rank-one perturbation theory; verified to 20 digits at a = 0.6, 0.7, 0.75, 0.8, 0.84 on the K-mode form]
The floor lambda*(a) is the lowest root of sigma(lambda, a) = 1/2, and the minimizer is f_a = phi/||phi|| at lambda = lambda*.
(Euler-Lagrange: K f - 2 <f,s> s = lambda f, so f = 2<f,s> (K - lambda)^{-1} s and <f,s> = 2 <f,s> sigma, i.e. sigma = 1/2.)

(T2) [theorem, implicit differentiation; verified to 10 digits] At the root,

    d sigma / d a (lambda*, a) = - lambda'(a) ||phi||^2,         d sigma / d a = - <phi, K_a' phi> + 2 <phi, s_a'>.

(T3) [theorem where the boundary law is: finite-tree non-entry supports, Corollary cor:tree(v); a.e. in general]
With the boundary law lambda' = -2|A_0|^2 and f = phi/||phi||, A_0 = A_0^(1)/||phi|| (A_0^(1) the amplitude of the unit-source
solution phi), (T2) reads

    d sigma / d a (lambda*, a) = 2 |A_0^(1)(lambda*, a)|^2 :

the a-derivative of the secular function at fixed spectral parameter is twice the squared edge amplitude of the unit-source
solution. This is the Hadamard formula of the resolvent.

(T4) [theorem, from (T1)-(T3)] The inequality (d), 2|A_0|^2 <= c T* lambda, is equivalent at such supports to

    d sigma / d a  <=  c T*(a) lambda  d sigma / d lambda        on the level set sigma(lambda, a) = 1/2,

that is, the level curve lambda*(a) of the secular function has logarithmic slope at least -c T*: d log lambda*/da >= -c T*.
With an integrable C(a) in place of c T*, the same with C.

(T5) [theorem, from monotonicity of sigma in lambda] Where mu_0(a) > 0,

    lambda*(a) >= 0   <=>   sigma(0, a) = < K_a^{-1} s_a, s_a > <= 1/2 ,

and the DEFICIT D(a) = 1/2 - sigma(0, a) satisfies D(a) = lambda*(a) ||phi(0, a)||^2 (1 + O(lambda* ||phi||^2 ...)), so
D decays like the floor. So, beyond a_Z and given mu_0 > 0: RH <=> the quadratic form of K_a^{-1} at the polar source
never exceeds 1/2. Measured: sigma(0, a) = 0.49998762, 0.4999999952, 0.49999999994, 0.4999999999997, 0.499999999999997 at
a = 0.6, 0.7, 0.75, 0.8, 0.84; the deficit's logarithmic rate -D'/D is 70.9 and 104.4 at a = 0.6 and 0.8 against
Phi' = 69.7 and 103.4 (ratios 1.018, 1.009).

## 3. A finding on the prime-side operator

[computation, K-mode form, 600 bits, K = 48..72, 2026-10-03] The bottom of K_a itself, the prime-side operator WITHOUT the
polar term, collapses doubly exponentially:

      a       mu_0(K_a)      mu_1(K_a)     lambda*        lambda*/mu_0     <e_0,s>^2/mu_0    sigma(0)
      0.6     3.07e-3        0.066         5.97e-7        1.9e-4           0.0427            0.49998762
      0.7     6.54e-6                      2.56e-10       3.9e-5
      0.75    9.18e-8                      3.38e-12       3.7e-5
      0.8     9.84e-10       4.68e-5       1.57e-14       1.6e-5           1.08e-8           0.4999999999997
      0.84    1.61e-11                     1.57e-16       9.8e-6

So the cancellation of the archimedean part against the entries is already in K_a: its bottom is 1e-10 at a = 0.8 (and its
second eigenvalue 5e-5), and the polar term lowers the bottom by a further factor 1e-4 to 1e-5 to the floor. The ground
state e_0 of K_a is nearly orthogonal to the polar source (<e_0,s>^2 = 1e-17 at a = 0.8) and the secular value 1/2 is
built almost entirely from the higher modes of K_a (the "rest" column: 0.457 at 0.6, 0.49999999 at 0.8). Under RH the
K_a-form is sum_rho |F(gamma)|^2 + 2<f,s>^2 >= 0 (the explicit formula with the polar term moved to the right), so mu_0 >= 0
is the weaker, polar-free half of Weil positivity; what is new is that it is small to the same double-exponential order.

Reading. The floor is not "three O(1) quantities cancelling once" but a two-stage cancellation: K_a nearly has a null
vector, and the polar source is tuned to that near-null space to a further 1e-5. The inequality (d), in the form (T4), is a
statement about how the level curve of sigma moves as the window grows; both derivatives in it are explicit resolvent
quantities of K_a at s_a.

## 4. What this gives the attack (N1d.2), and what it does not

Gives: (i) an exact, polar-free object, K_a, whose resolvent at s_a carries the whole problem; (ii) the target as a
level-set slope, with both sides computable from (K_a, s_a) and their a-derivatives (the engine supplies K_a, K_a', s_a,
s_a' in one pass); (iii) the deficit D(a) as the quantity whose sign is the Hypothesis beyond a_Z and whose logarithmic
derivative is the relay.

Does not give: a bound. The equivalence (T4) is exact and therefore as hard as (d). The sub-questions it isolates:
  (Q1) why mu_0(K_a) collapses (the polar-free half of the cancellation), and at what rate: is -log mu_0 ~ pi kappa_K T*
       with its own kappa_K? (Data: -log mu_0 rises 61, 86, 90, 105 per unit a across the five supports, against
       Phi' = 70, 91, 98, 103, 121.)
  (Q2) the orthogonality <e_0, s>^2 << mu_0: a structural reason (the ground state of K_a hides from sinh(x/2) as it hides
       from the zeros) would be the polar analogue of darkness.
  (Q3) the slope inequality (T4) in the finite-tree regime a < a_inf = 0.843, where sigma_eff and its factors are finite
       continued fractions and A_0^(1) is explicit (Corollary cor:edgeform): the first place to try to prove it.

## 5. Not to drift into
More supports beyond what (Q1)-(Q3) need; refitting; the fiction. The computations above are the five finite-tree supports
and nothing else; the a = 2.0 tail-law run (user's earlier request) continues in the background.
