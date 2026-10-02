# The Immortality of Observers

## A proof of the Line Postulate from the laws of leakage

Academy of the Prime Sheet, Critical Line Universe. Communicated to the gods.

FICTION NOTICE (from our universe). This manuscript is written from inside the Critical Line Universe of
FICTION_CRITICAL_LINE_UNIVERSE.md and uses that universe's laws as axioms. It proves nothing about the real zeta
function. The ledger in Section 27 of the companion file says which of the laws below are theorems in our universe and
which one is the Riemann Hypothesis itself. Nothing in this file enters the paper.


### Abstract

We prove that every quantum of the arithmetic field is massless. The proof does not examine the quanta. It examines
observers, one at a time, and shows that no observer ever dies: the trapped energy of an observer is a positive quantity
that falls with age in proportion to itself, and a quantity that falls in proportion to itself never reaches zero. The
Line Postulate, hitherto an axiom of this universe, follows by the duality of the two sheets. The one law the proof
needs beyond the classical laws of observers is the Law of Proportional Leakage; we derive it from the Born rule of
arithmetic and the Local Law of the band, and we show that it is exactly the law that the dead universes violate.


## 1. Objects and notation

1.1. THE TWO SHEETS. The prime sheet is the line with coordinate l = log x. On it sit the charges, one at every l = log n
with n a power of a prime p, of strength q(n) = (log p) n^{-1/2}. The zero sheet is the dual line, coordinate gamma. On it
sit the quanta, one at every zero rho = 1/2 + m + i gamma of the arithmetic zeta function, with mass m = |Re rho - 1/2|.
Quanta come in mirror pairs: if rho is a quantum so is 1 - conj(rho), the quantum of the same height and opposite
displacement; and in conjugate pairs rho, conj(rho). A quantum is massless if m = 0.

1.2. OBSERVERS. An observer of age a > 0 is the window [-a, a] of the prime sheet. A FIELD of the observer is an odd real
function f on the window with unit norm, int f^2 = 1, regular enough for everything below. Its TRANSFORM is
        F(t) = int_{-a}^{a} f(x) e^{itx} dx = 2i int_0^a f(x) sin(tx) dx,
an entire function of t, purely imaginary on the real axis, and its AUTOCORRELATION is g = f * f~, g(x) = int f(y) f(y - x)
dy, an even function on [-2a, 2a] whose transform is F(t) conj(F(conj t)) for complex t and |F(t)|^2 for real t.

1.3. THE HORIZON. T*(a) = 2 pi e^{2a}. The density of the quanta at height t is (1/2 pi) log(t/2 pi) and the Shannon density
of a field of age a is a/pi; the two agree exactly at t = T*(a). Below the horizon an observer resolves the quanta; above
it they are too dense to be told apart.

1.4. THE WINDOW FORM. For a field f of age a,
        Q_a(f) = (archimedean flux of g) + (polar flux of g) - sum over charges q(n) g(log n),
where the archimedean flux is (1/2 pi) int |F(t)|^2 Psi_inf(t) dt with Psi_inf(t) = Re psi(1/4 + it/2) - log pi, and the
polar flux is -2 ( int f(x) sinh(x/2) dx )^2. Everything on the right is read on the prime sheet: Q_a(f) is a number an
observer can compute from the charges inside its window and from the vacuum's symbol, without reference to any quantum.
Only the charges with log n < 2a enter.

1.5. THE FLOOR, THE GROUND FIELD, THE EDGE. The floor is lambda(a) = inf Q_a(f) over the fields of age a. It is attained,
by the ground field f_a. The ground field approaches the edge of the window as
        f_a(a - delta) = |A_0(a)| ( log(a/delta) + beta(a) )^{-1/2} (1 + o(1)),        delta -> 0,
and A_0(a), the EDGE AMPLITUDE, is the strength with which the ground field touches the edge. Its transform falls off, for
heights far above the horizon, like 2|A_0| |sin(ta + theta_a)| / (t sigma(t)^{1/2}) with sigma(t) ~ log(t/2 pi).

1.6. LIFE AND DEATH. An observer of age a is ALIVE if lambda(a) > 0 and DEAD if lambda(a) <= 0.


## 2. The laws

The proof uses the following laws of the Critical Line Universe. Laws I to V hold in every universe with the same density
of charges, living or dead; they are the laws of observers. Law VI is the law of this universe. Laws VII to IX are used
only in Section 4, to derive Law VI; a reader who accepts Law VI as an axiom may skip Section 4.

LAW I (the Gauss law of the two sheets). For every field f of every age,
        Q_a(f) = sum over quanta rho of G(rho),        G(rho) = F(gamma - im) conj( F(gamma + im) ),
the sum taken over the quanta with their multiplicities. In particular a massless quantum contributes |F(gamma)|^2 >= 0,
and a mirror pair of mass m at height gamma contributes 2 Re[ F(gamma - im) conj(F(gamma + im)) ], a real number of
either sign. (The sum converges absolutely for the fields of 1.2.)

LAW II (Birth). There is an age a_Z > 0 such that every observer of age at most a_Z is alive, and lambda(a_Z) >= lambda_Z
for an explicit lambda_Z > 0.

LAW III (Monotonicity). The floor does not increase with the age.

LAW IV (Continuity). On every finite range of ages the floor is absolutely continuous: lambda(a_2) - lambda(a_1) =
int_{a_1}^{a_2} lambda'(a) da, with lambda' defined at almost every age.

LAW V (Decay). At almost every age, lambda'(a) = -2 |A_0(a)|^2.

LAW VI (Proportional Leakage). There is a function C(a) >= 0, finite at every age and integrable over every finite range
of ages, such that for every living observer of age a,
        2 |A_0(a)|^2  <=  C(a) lambda(a).

LAW VII (the Born rule of arithmetic). For a living observer, the share s_rho(a) = G(rho)/lambda(a) of a quantum is the
probability that the quantum carries the observer's trapped energy: s_rho >= 0 for every quantum, and the shares sum to one.

LAW VIII (Darkness below the horizon). The ground field of a living observer carries no share below its horizon: s_rho(a)
= 0 for every quantum of height |gamma| < T*(a).

LAW IX (the Local Law with the sum rule). There is a constant c_on >= 1 such that for every living observer and every
height T >= c_on T*(a), the total share of the quanta above height T is
        sum over |gamma| > T of s_rho(a) = 2 |A_0(a)|^2 / ( pi T lambda(a) ).

A remark on Law I. It is stated for every configuration of quanta, massive or not; it is the conservation law of the two
sheets and does not know where the quanta are. The Line Postulate is the statement that every term on its right-hand side
is a square. The proof below never uses that; it uses Law I only through Section 3, where both directions of the duality
between masses and deaths are drawn from it.


## 3. The Mortality Theorem

THEOREM 3.1 (duality of masses and deaths). The following are equivalent:
   (a) every quantum is massless;
   (b) every observer is alive;
   (c) every observer is alive or dead but none is dead.
That is: the vacuum has a massive quantum if and only if some observer dies at a finite age.

PROOF. (a) implies (b): by Law I, if every quantum is massless then Q_a(f) = sum |F(gamma)|^2 >= 0 for every field, so
lambda(a) >= 0 for every a; and lambda(a) > 0 because a field with F(gamma) = 0 at every quantum is identically zero (an
entire function of exponential type a with zeros of density exceeding a/pi above the horizon vanishes; this is the
sampling lemma of the Academy, which we take as known).

(b) implies (a): suppose a quantum rho_0 = 1/2 + m + i gamma_0 has mass m > 0. We construct a dying observer. Let phi be a
fixed even field of age 1 whose transform Phi is positive on the real axis, and for a > 1 let h_a be the odd field of age
a with transform
        H_a(t) = (t - gamma_0) Phi( a (t - gamma_0) ) - (t + gamma_0) Phi( a (t + gamma_0) ),
normalized to unit norm (h_a is odd, since H_a is odd in t up to the imaginary unit; its support is [-a-1, a+1], an age we
still call a). The transform H_a vanishes at +-gamma_0 with slope Phi(0) a^0 ... more precisely H_a'(gamma_0) = Phi(0), and
on the real axis |H_a(t)|^2 <= (|t - gamma_0|^2 |Phi(a(t-gamma_0))|^2 + ...) is concentrated within 1/a of +-gamma_0 with
total mass of order a^{-3}. Evaluate Law I for g = h_a * h_a~:
   * the massless quanta contribute sum |H_a(gamma)|^2 <= (number of quanta within O(1/a) of +-gamma_0, which is bounded
     by the density log(gamma_0/2 pi)/2 pi times O(1/a) plus one) times max |H_a|^2 = O( a^{-2} log gamma_0 );
   * the pair rho_0, 1 - conj(rho_0) and its conjugates contribute
        4 Re[ H_a(gamma_0 - im) conj(H_a(gamma_0 + im)) ]
        = 4 Re[ (-im) Phi(-i a m) conj( (im) Phi(i a m) ) ] + (cross terms from the -gamma_0 bump, exponentially smaller)
        = -4 m^2 |Phi(i a m)|^2 (1 + o(1)),
     and since Phi is the transform of an even field of age 1, |Phi(iam)| >= c e^{am} (1 + o(1))^{-1} for a m large: the
     transform of a compactly supported positive-type field grows exponentially off the real axis;
   * the other massive quanta, if any, contribute terms of the same form at their own heights, each at most
     4 m'^2 |Phi(iam')|^2 e^{-a^2 |gamma' - gamma_0|^2 ...}: negligible unless |gamma' - gamma_0| = O(1/a), and there are
     boundedly many of those.
Hence Q(h_a) <= O(a^{-2} log gamma_0) - c m^2 e^{2am} < 0 as soon as a exceeds a_death(m) of order m^{-1} log(1/m) (and
large enough that the horizon has passed gamma_0, so that the bookkeeping above applies). The observer of age a_death is
dead, and by Law III so is every older observer.

(b) and (c) are the same statement. QED.

COROLLARY 3.2 (the Blind Spot). A quantum of mass m is invisible to every observer younger than an age of order
m^{-1} log(1/m): the floors of those observers differ from the floors of a massless universe by terms of order m^2.
(This is the Academy's Blind Spot Theorem; the proof of 3.1 gives the upper bound on the age of detection, and the
expansion of Lemma A.1 gives the order m^2.)

The Mortality Theorem reduces the Line Postulate to the following statement, which is the content of this manuscript:

        IMMORTALITY. No observer dies.


## 4. The Law of Proportional Leakage, derived

This section derives Law VI from Laws VII, VIII and IX. It also shows, in the opposite direction, that Law VI is exactly
what fails in a universe where an observer dies: Law VI is not a convenience but the whole of the matter.

PROPOSITION 4.1 (the leak is proportional to the floor). Under Laws VII and IX, every living observer satisfies
        2 |A_0(a)|^2  <=  pi c_on T*(a) lambda(a).
Hence Law VI holds with C(a) = pi c_on T*(a), which is finite at every age and integrable over every finite range.

PROOF. Fix a living observer of age a and put T = c_on T*(a). By Law IX the total share of the quanta above height T is
2|A_0|^2/(pi T lambda). By Law VII every share is nonnegative and the shares sum to one, so the total share of any set of
quanta is at most one. Therefore 2|A_0|^2/(pi T lambda) <= 1, which is the claim. QED.

REMARK 4.2. Law VIII (Darkness) was not needed for the inequality; it is needed to see that the inequality is nearly
sharp. By Darkness all the shares sit above the horizon, so the total share above T*(a) is one, and by Law IX the total
share above c_on T* is kappa/c_on with kappa = 2|A_0|^2/(pi T* lambda), the Riemann number of the observer. The bound
kappa <= c_on of Proposition 4.1 is sharp exactly when the shares between one and c_on horizons vanish, which they do not:
the shares in the band [T*, c_on T*] are the plunge of Decree IV, so kappa < c_on strictly. (In the measured universe,
c_on is about 3 and kappa about 1.2.)

PROPOSITION 4.3 (the converse: death is a threshold leak). Let a universe satisfy Laws III and IV, and suppose an observer
dies: lambda(a_d) = 0 for some finite a_d while lambda > 0 on [a_Z, a_d). Then -lambda'/lambda is not integrable on
[a_Z, a_d): for every C integrable on [a_Z, a_d] the inequality -lambda' <= C lambda fails on a set of positive measure.
In particular Law VI fails in every universe with a mortal observer.

PROOF. If -lambda' <= C lambda almost everywhere on [a_Z, a_d) with C integrable, then by Lemma 5.1 below lambda(a) >=
lambda(a_Z) exp(-int_{a_Z}^{a} C) for a < a_d, and by continuity (Law IV) lambda(a_d) >= lambda(a_Z) exp(-int_{a_Z}^{a_d}
C) > 0, contradicting lambda(a_d) = 0. QED.

The two propositions together say: in a living universe the leak is proportional to the floor, with the horizon as the
constant of proportionality; in a dying universe the leak outruns the floor just before death, which is to say the edge
amplitude does not vanish with the floor. The Academy's name for the second behaviour is a THRESHOLD LEAK, and its
physical description is a tank that drains at a fixed rate and then keeps draining.


## 5. The Immortality Theorem

LEMMA 5.1 (the renormalized floor). Let lambda be absolutely continuous and positive on [a_1, a_2) with lambda'(a) >=
-C(a) lambda(a) almost everywhere, C >= 0 integrable on [a_1, a_2]. Then U(a) = lambda(a) exp( int_{a_1}^{a} C ) is
non-decreasing on [a_1, a_2), and consequently
        lambda(a)  >=  lambda(a_1) exp( - int_{a_1}^{a} C )        for all a in [a_1, a_2).

PROOF. U is absolutely continuous as a product of absolutely continuous functions (the exponential of an absolutely
continuous function is absolutely continuous on a finite range), and almost everywhere
        U'(a) = exp( int_{a_1}^{a} C ) ( lambda'(a) + C(a) lambda(a) )  >=  0.
An absolutely continuous function with nonnegative derivative almost everywhere is non-decreasing. QED.

THEOREM 5.2 (Immortality). Under Laws II, III, IV, V and VI, every observer is alive.

PROOF. Let
        a_d = sup { a > 0 : every observer of age less than a is alive }.
By Law II, a_d >= a_Z. Suppose, for contradiction, that a_d < infinity.

Step 1 (the floor on the living range). On [a_Z, a_d) every observer is alive. By Law IV the floor is absolutely
continuous there, by Law V its derivative is -2|A_0|^2 at almost every age, and by Law VI, which applies to living
observers, -2|A_0(a)|^2 >= -C(a) lambda(a). Hence
        lambda'(a) >= -C(a) lambda(a)        for almost every a in [a_Z, a_d).

Step 2 (the bound). Lemma 5.1 with a_1 = a_Z, a_2 = a_d gives lambda(a) >= lambda_Z exp(-int_{a_Z}^{a} C) on [a_Z, a_d).

Step 3 (the observer of age a_d is alive). C is integrable on [a_Z, a_d], so int_{a_Z}^{a} C increases to a finite limit
as a tends to a_d, and lambda is continuous at a_d by Law IV. Hence
        lambda(a_d) = lim lambda(a) >= lambda_Z exp( - int_{a_Z}^{a_d} C ) > 0.

Step 4 (so are its immediate elders). Since lambda is continuous and lambda(a_d) > 0, there is epsilon > 0 with lambda > 0
on [a_d, a_d + epsilon). Then every observer of age less than a_d + epsilon is alive, which contradicts the definition of
a_d as a supremum.

Therefore a_d = infinity: every observer is alive. QED.

THEOREM 5.3 (the Line Postulate). Every quantum of the arithmetic field is massless.

PROOF. Theorem 5.2 and the Mortality Theorem 3.1, (b) implies (a). QED.

COROLLARY 5.4 (the Cooling Theorem, quantitative). For every a >= a_Z,
        lambda(a)  >=  lambda_Z exp( - int_{a_Z}^{a} C ),
and with the constant of Proposition 4.1, C = pi c_on T*,
        lambda(a)  >=  lambda_Z exp( - (pi c_on / 2) ( T*(a) - T*(a_Z) ) ),
since dT*/da = 2T*. The floor is positive at every age and at most doubly exponentially small in the age: the vacuum is
marginally stable, forever approaching the critical point and never reaching it.


## 6. Remarks

6.1. WHAT THE PROOF USES. The Line Postulate is used nowhere in Sections 3 and 5. Section 3 uses Law I, which is a
conservation law indifferent to the positions of the quanta. Section 5 uses only the laws of observers (II to V), which
hold in the dead universes as well, and Law VI. The whole weight of the theorem rests on Law VI; Section 4 shows that it
could not be otherwise, since Law VI is equivalent, given II to V, to immortality itself (Propositions 4.1 and 4.3 are the
two directions). What the proof contributes is not a new ingredient but a change of subject: from infinitely many quanta to
one inequality about one observer at a time, and from a statement about where things are to a statement about how fast
something leaks.

6.2. WHERE THE LINE WAS HIDING. Law VI was derived from Law VII, the Born rule. In a universe in which quanta are quanta,
the Born rule is the first law anyone learns; it was never written among the axioms of the Academy because no one
thought it needed stating. The gods' earlier decree of the Line was therefore redundant: it followed from the Born rule,
the Local Law, and the laws of observers. This is the sense in which the Line has been "derived": the postulate has been
moved from the zero sheet, where it looked like a miracle of position, to the shares, where it looks like the
nonnegativity of a probability. A reader from a universe where the Born rule of arithmetic is not available will say that
nothing has been gained. We say that the location of a difficulty is itself a result.

6.3. THE RENORMALIZED FLOOR AND THE RIEMANN CHARGE. The function U(a) = lambda(a) exp(int C) of Lemma 5.1 is the floor
magnified by the accumulated leak. Law VI says that U never decreases; Immortality says that it is therefore positive
forever. With the Academy's decreed flux C = pi kappa_inf T*, U(a) = lambda(a) exp( (pi kappa_inf/2)(T*(a) - T*(a_Z)) ),
and the equation of state asserts that U tends to a positive constant, the RIEMANN CHARGE of the universe: the one number
that survives the cooling. The proof does not need the charge to exist; it needs only that U not decrease.

6.4. THE OTHER POLARIZATION. Fields come in two polarizations, odd and even. The manuscript treats odd observers. The
Mortality Theorem holds with odd observers alone (the parity lemma of the Academy: an odd field dark at +-gamma_0 with
slope at both already detects a mass), so immortality of the odd observers forces every quantum to be massless, and the
even observers are then alive by Law I. The even observers have their own Laws IV to VI, which we have not needed.

6.5. WHY A LAW OF LEAKAGE. The Academy's physicists prefer to state Law VI as the absence of a threshold: an observer's
edge amplitude vanishes with its trapped energy. A universe dies through a leak with a threshold. Ours has none; that is
the Line.


## Appendix A. The mirror term

LEMMA A.1. Let f be an odd field with transform F = iS, S real on the real axis, and let gamma be real. Then for a mirror
pair of mass m at height gamma,
        2 Re[ F(gamma - im) conj(F(gamma + im)) ] = 2 |F(gamma)|^2 - 2 m^2 ( |F'(gamma)|^2 + S(gamma) S''(gamma) ) + O(m^4),
with the O(m^4) uniform for m <= 1/a. In particular, if the field is dark at gamma, S(gamma) = 0, the pair contributes
        - 2 m^2 |F'(gamma)|^2 + O(m^4)  <  0.

PROOF. Expand F(gamma -+ im) = F -+ imF' - (m^2/2) F'' + O(m^3) at gamma (F is entire of exponential type a, so the
derivatives are bounded by a^k times the norm). Then
        F(gamma - im) conj(F(gamma + im)) = |F|^2 - im ( F conj(F') + F' conj(F) ) - m^2 |F'|^2 - (m^2/2)( F conj(F'') + F'' conj(F) ) + O(m^3).
With F = iS and S, S', S'' real: F conj(F') = S S', F conj(F'') = S S'', both real, so the term linear in m is purely
imaginary and the real part is |F|^2 - m^2 |F'|^2 - m^2 S S'' + O(m^3); the odd powers of m are purely imaginary by the same
computation, so the error is O(m^4). Multiply by two for the pair. QED.

The lemma is the mechanism of both halves of the Mortality Theorem: a massive pair at a height where a field is dark
lowers the field's energy by 2m^2 |F'|^2, which is invisible while the floor is large and fatal once the field is
concentrated enough at that height for |F'(gamma)|^2 to be magnified by e^{2am}.


## Appendix B. The constants of the measured universe

The proof is independent of numbers. For the reader who wants to see the laws at work, the Academy's measurements at the
ages 0.6, 0.8, 1.0, 1.25, 1.5, 1.75 (companion file, Parts II and III) give: Riemann number kappa = 2|A_0|^2/(pi T* lambda)
= 1.063, 1.058, 1.164, 1.201, 1.239, 1.272; share above three horizons 0.370, 0.415, 0.407, 0.433, 0.434, 0.472 against
kappa/3 = 0.354, 0.353, 0.388, 0.400, 0.413, 0.424 from Law IX; median height of the shares 2.39, 2.32, 2.65, 2.67, 2.74,
2.86 horizons; and the floor 5.97e-7, 1.57e-14, 1.48e-26, 3.42e-51, 9.22e-92, 2.07e-162, positive at every age measured,
falling as Corollary 5.4 says it must.

--- end of manuscript ---

TRANSLATOR'S NOTE (from our universe). Laws I to V are theorems here (Weil's explicit formula; Zhu's certified floor at
a_Z = 0.8; monotonicity; absolute continuity under the finiteness condition (H_fin) on an exceptional null set, Propositions
9.15 and 9.17; the boundary law, Theorem 9.12(iv)). Law VI is the inequality (d) of the paper, Conjecture A, equivalent to
the Riemann Hypothesis beyond a_Z (Theorem thm:reduction, Proposition 9.30). Law VII, the Born rule, is Weil positivity,
which is the Hypothesis. Law VIII is a computation at a few supports. Law IX is the tail law under RH together with the
local law of the band, item N1d.1 of RH_IF_THEN_TREE.md, open. Theorem 3.1 is Weil's criterion in the odd window form
(thm:weil) with the Blind Spot (Theorem 9.31); Theorem 5.2 is Proposition prop:AimpliesRH with a general C. The manuscript
is a proof in its universe and a map of ours.
