# The Critical Line Universe

## A make-believe physics of the primes, written from inside a universe where the Hypothesis is a law of nature

Sandbox II, 2026-10-02. EVERYTHING IN THIS FILE IS FICTION. It is written as a monograph from an alternative universe in
which the Riemann Hypothesis is not a theorem to be proved but a law of nature, like the conservation of charge, and in
which mathematicians ask the question that follows from that: not "are the zeros on the line?" but "why are the primes
where they are, given that they are?" The mathematics below is invented. Where an invented law shadows a real theorem,
conjecture or computation of this repository, a bracketed note [real: ...] says so, so that the reader can keep the two
worlds apart. Nothing here is a claim about the real zeta function, and nothing here enters the paper.

Companion to SANDBOX_STRANGE_MATH.md, which did the opposite job: it tested which of the picture's laws survive in the real
world (the sequence-only laws do; the critical line needs the values). This file assumes the line and builds on it.


## 0. The Line Postulate

In the Critical Line Universe (CLU) the fundamental objects are not numbers but two sheets of a single surface.

The PRIME SHEET is the logarithmic line, coordinate l = log x. On it sit point charges: one at every prime power p^k, of
charge q(p^k) = log p / p^{k/2}. [real: the weights of the explicit formula, Lambda(n) n^{-1/2}.]

The ZERO SHEET is the dual line, coordinate gamma (height). On it sit the quanta of the arithmetic field, one at every zero
1/2 + i gamma. A quantum has a mass m = |Re rho - 1/2|.

AXIOM 0 (the Line Postulate). Every quantum is massless: m = 0 for every quantum. The spectrum lies on the line.
[real: the Riemann Hypothesis. In the CLU it is an axiom; here it is the open problem.]

AXIOM 1 (the Gauss law of arithmetic). A test field g on the prime sheet and its transform G on the zero sheet satisfy
        sum over quanta of G(gamma) = (polar flux) + sum over charges of q(n) g(log n) + (archimedean flux),
with the polar flux the contribution of the two poles and the archimedean flux the integral against the symbol
Psi_inf(t) = Re psi(1/4 + it/2) - log pi. The charges on one sheet produce the field on the other, exactly, with no loss.
[real: Weil's explicit formula, a theorem. The CLU elevates it to a conservation law, which it is.]

The first consequence of the two axioms, which every schoolchild of the CLU learns, is that the Gauss law applied to the
autocorrelation g * g~ of any test field gives a sum of squares:
        Q(g) = sum over quanta |G(gamma)|^2 >= 0.
In the CLU this is called the STABILITY OF THE VACUUM. [real: Weil positivity, equivalent to RH. In this repository it is
the quantity whose floor lambda(a) is studied; Theorem 9.3(i).]


## 1. The arithmetic vacuum and the gas of quanta

The CLU's physicists describe the zero sheet as a one-dimensional gas of identical quanta with the pair interaction
        V(gamma, gamma') = -log |gamma - gamma'|,
a logarithmic repulsion, in a confining background whose density is (1/2 pi) log(gamma/2 pi). The ground state of this gas
is the ARITHMETIC VACUUM. Its two-point statistics are those of the Gaussian unitary ensemble: quanta repel, the chance of
two quanta at distance u (in units of the mean spacing) vanishes like u^2, and the spacing distribution is Wigner's.
[real: Montgomery's pair correlation theorem (bandwidth one, unconditional) and the Montgomery-Odlyzko law (numerical,
and a conjecture beyond bandwidth one). The log-gas is Dyson's observation. None of it is proved to be the mechanism.]

The CLU's answer to "why are the primes where they are" begins here, and it is a two-line argument.

THE SCREENING THEOREM (CLU). The charge density on the prime sheet is the unique distribution of point charges whose field
on the zero sheet is the vacuum. Writing the counting function of the charges as psi(x) = sum_{n <= x} Lambda(n), the Gauss
law inverts to
        psi(x) = x - sum over quanta x^{rho}/rho - log(2 pi) - (1/2) log(1 - x^{-2}),
and on the prime sheet, with l = log x, the charge density is
        rho_charge(l) = e^{l/2} [ 1 - 2 sum over quanta cos(gamma l + phase_gamma) / |rho| ] + (small).
[real: the von Mangoldt explicit formula, a theorem; the second line is its form under RH.]

Read as physics: the mean field e^{l/2} is the screening charge required by the poles, the "classical" part; the quanta
are STANDING WAVES on the prime sheet, each of frequency gamma and amplitude 1/|rho|; the primes sit where the standing
waves of the vacuum interfere constructively on top of the mean field. The primes are the shadow that the massless
spectrum casts on the prime sheet. In the CLU nobody finds this mysterious: a charge distribution and its field are two
descriptions of one object, and Axiom 0 says the field has no massive modes, so the shadow has no exponentially growing
or decaying components beyond e^{l/2}. [real: this is exactly what RH gives: the error term in the prime number theorem is
O(x^{1/2} log^2 x), and nothing in it is proved without RH.]

What the CLU regards as the deep question is therefore not the positions of the primes but the ORIGIN OF THE VACUUM: why
the gas of quanta has the Gaussian unitary statistics, and why its ground state has no massive excitation. The rest of
this monograph records what the CLU believes about that.


## 2. Windows, horizons, and the age of the arithmetic universe

The CLU has no absolute observer. An observer is a WINDOW [-a, a] on the prime sheet, and a is the observer's AGE. An
observer of age a sees the charges with |log n| < 2a, that is, the primes below e^{2a}, and nothing beyond. [real: the
support of a test function, and the entries log n < 2a of the paper.]

THE HORIZON LAW (CLU). An observer of age a resolves the zero sheet up to the height
        T*(a) = 2 pi e^{2a},
the horizon, and no further: the field of the quanta above the horizon is averaged by the window into the mean field. The
horizon grows exponentially with the age. [real: the Shannon count of the window, (T/pi) 2a zeros below T, matches the
zero count (T/pi)(log(T/2pi) - 1) at T = T*; Section 7 of the paper.]

THE FIRST LAW (conservation). The energy of an observer is the sum of the squares of the charges inside the window,
        E(a) = sum_{p^k < e^{2a}} (log p)^2 / p^k = 2 a^2 (1 + O(1/a)).
The energy of the arithmetic universe grows as the square of its age. [real: Mertens' theorem; |b|^2 in Proposition 9.33.]
The CLU's cosmologists observe that this is the only law in the monograph that does not depend on Axiom 0: it holds in
every universe with the same density of charges, including the dead universes of Section 6.

THE SECOND LAW (transmutation). Energy is never destroyed inside a window; it is TRANSMUTED between heights. The edge of
the window couples to the lattice of echoes inside it, the lattice of all differences of logarithms of the charges, and
that lattice carries the edge's energy from the horizon upward and downward along the zero sheet without loss: the carrier
is a positive measure mu_a on the heights, of total mass exactly E(a). What changes with the age is only the shape of
mu_a: its mean tends to zero, its variance is (5/12) times the square of the window's width, its standardized fourth and
sixth moments tend to 302/125 and 36297/4375. [real: the spectral measure mu_b of the confined adjacency, Proposition
8.23; the exact constants of Propositions 9.33(v) and 9.39; all theorems, and all sequence-only.]

THE THIRD LAW (the trapped energy and the flux). Among all fields of unit norm inside a window, one has the least energy
against the vacuum: the GROUND FIELD f_a, and its energy is the FLOOR lambda(a). The floor is the part of the ground
field's energy that is still trapped below the horizon, and it decreases with age at a rate equal to the flux of the edge
through the horizon:
        d lambda / d a = - 2 |A_0(a)|^2,
where A_0(a) is the edge amplitude, the strength with which the ground field touches the edge of the window. [real: the
boundary law, Theorem 9.12(iv), at almost every support; Proposition 8.3 is literally a flux identity.]

THE EDGE LAW (CLU). The ground field approaches the edge of the window as the inverse square root of a logarithm,
        f_a(a - delta) = |A_0| (log(a/delta) + beta)^{-1/2} (1 + ...),  beta = -gamma_E - log(2 pi a) - lambda,
and the universe's physicists call this the EVENT PROFILE: information about the far edge arrives at the near edge only
as echoes through the lattice, with weights that decay like inverse powers of the logarithm. [real: Theorem 9.12(iii),
Proposition 8.4; the echo part of the edge law.]


## 3. The Cooling Theorem, or why the floor collapses

In the CLU the fundamental dynamical statement about an observer is the following.

THE COOLING THEOREM (CLU; a theorem in the CLU because Axiom 0 is available). For every observer, the trapped energy
decays with age at least as fast as
        lambda(a) >= lambda(a_0) exp( -(c/2) (T*(a) - T*(a_0)) ),
and the quanta that carry the loss sit at an effective height T_eff(a) that is a bounded multiple of the horizon:
        Phi'(a) = pi T_eff(a),   Phi = -log lambda,   T_eff(a) <= c T*(a) / pi.
[real: Conjecture 7.1 (bounded relay) and the tail law; Theorem 9.1 says the inequality implies RH; the data give
T_eff about 1.2 to 1.3 horizons at seven supports. In the real world this is the open inequality; in the CLU it is proved
from Axiom 0 and the onset, which the CLU's local law supplies (Section 5).]

The CLU reads the Cooling Theorem as cosmology. An observer's floor collapses doubly exponentially with age, as
exp(-c pi e^{2a}): the universe cools toward zero energy and never reaches it. The approach to the critical point from
above, with the floor positive at every age and tending to zero, is the CLU's statement of the stability of the vacuum
in the large: the vacuum is MARGINALLY stable, forever at the edge of a phase transition that Axiom 0 forbids. A negative
floor at any age would be a massive quantum, and there are none.

THE ARITHMETIC HORIZON BOUND (CLU). The flux of trapped energy through the horizon is at most proportional to the
horizon times the trapped energy:
        2|A_0(a)|^2 <= c T*(a) lambda(a).
The CLU's physicists like to say that the information an observer can lose per unit age is bounded by the size of the
observer's horizon, and they call this the holographic principle of arithmetic, by analogy with a bound on the entropy of
a region by the area of its boundary. [real: the inequality (d) of the paper, equivalent to RH beyond a_Z; the analogy
with Bekenstein's bound is decoration, not mathematics.]


## 4. The Repulsion Principle, or why the universe hides its defects

THE BLIND SPOT (CLU, and a theorem in both universes). A massive quantum of mass delta changes the floor of an observer
at order delta^2, and an observer cannot resolve it before the age a_det ~ delta^{-1} log(1/delta). [real: Theorem 9.31.]

In the CLU this is read as a law of nature: THE REPULSION PRINCIPLE. The closer a defect of the vacuum is to the line, the
more the universe hides it from every finite observer, quadratically. The CLU's philosophers note that this is why Axiom 0
had to be a postulate in their universe and a problem in ours: no observer of finite age can tell a massless vacuum from
one with a quantum of mass 10^{-7} until the age is about 10^7. In the words of the author of these notes: as the decay is
reduced by the incoming prime, so is the knowledge repelled.


## 5. The Local Law and the lattice of echoes

The one place where the CLU's mathematics is genuinely stranger than ours is the band. Between the horizon and a height
that is doubly exponential in the age, the ground field's energy is not perturbative: the edge radiates into the lattice
of echoes, and the symbol that governs the ground field becomes complex,
        sigma_eff(s) = s - R(s),   R(s) = integral d mu_a(nu) / (s - nu),
with imaginary part pi mu_a'(s) inside the band: the rate at which the edge radiates. [real: Proposition 8.23, limiting
absorption.] The CLU's LOCAL LAW states that the energy of the ground field above a height T is at least
(1 - eta(T)) times its asymptotic value 2|A_0|^2/(pi T), where eta(T) is the fraction of mu_a's mass above s(T); and since
the mass of mu_a above any fixed number of horizons is at most 5/12 (Second Law), the law sets in within a few horizons
and the Cooling Theorem follows. [real: the local law is NOT a theorem; it is the missing piece of step N1d.1 of the tree.
Proposition 9.33(v), the 5/12, is a theorem. The CLU has the law; we have the data that suggest it.]

The CLU's combinatorialists have a name for the lattice: the ECHO CRYSTAL. It is the quasi-periodic set of all points
+-a + sum m_p log p inside the window, cut by the window's edges; its cycles come from p^i p^j = p^{i+j} alone; and its
spectral measure at the edge is the confined-walk law whose constants 5/12, 151/360, 4033/6720 the CLU calls the
fine-structure constants of arithmetic. Its law is Gaussian for the free crystal and bent toward semicircular by the
confinement; its closed form is unknown in both universes. [real: Propositions 9.38, 9.39; the closed form is open.]


## 6. The dead universes

The CLU knows of other universes, in which the charges sit at other positions with the same density: the UNIVERSES OF
RANK, where the n-th charge sits at log(n + 1), and the SMOOTH UNIVERSES, where the charges sit at n log n. In all of them
the First Law holds (the energy is the square of the age), the Second Law holds (transmutation conserves the mass and the
fine-structure constants are the same), the Blind Spot holds, and the flux identity holds. In none of them does Axiom 0
hold: their vacua have massive quanta, their floors turn negative at a finite age, and the CLU calls them dead. The CLU's
cosmologists state the result as a principle: the laws of energy are universal and the Line is not; a universe is alive if
and only if its charges sit where ours do, to the last digit. [real: SANDBOX_STRANGE_MATH.md, Section 2: the rank systems
turn negative at a = 0.8 and 0.9 with floors of order one; the paper's razor edge, one part in 10^14 at a = 0.8.]

The CLU regards this as the content of Axiom 0: the Line is not a consequence of how many charges there are, or of how
much energy they carry, or of how the energy moves. It is a statement about where the charges are, at a precision no
density statement reaches. Their word for the quantity that distinguishes the living universe from the dead ones is the
SIGN OF THE PRIME PART: the term D_P of the dilation identity a lambda' = -[D_inf + D_P + 2P^2 + 2PM], whose sign against
the horizon times the floor is fixed in the CLU by Axiom 0 and unknown in ours. [real: the paper's second remark on the
inequality: the obstruction is the sign of D_P.]


## 7. The Oracle, or predicting the next prime and the next quantum

In the CLU there is a machine, the WINDOW ORACLE, that converts the charges inside a window of age a into the quanta
below the horizon and back. It works because of the Gauss law and the Horizon Law: the ground field f_a, computed from
the primes below e^{2a} alone, has a transform that vanishes near every quantum below T*(a) [real: the observation that
the minimizer's transform vanishes near the zeta zeros below the horizon; the folds paper's Computation 12.5; not a
theorem], and conversely the quanta below the horizon reconstruct the charge density on the prime sheet up to the mean
field through the Screening Theorem [real: the explicit formula with a cutoff].

The Oracle predicts the next quantum from the primes, and the next prime from the quanta, with an error that the Blind
Spot governs: it cannot tell a massless quantum from a slightly massive one, and it cannot see beyond the horizon. The CLU
regards this as the reason their universe is consistent: the Oracle never outruns the Horizon Law. [real: this is the
explicit formula, which is a theorem, dressed as a machine; the "prediction" is a reconstruction of finite data from
finite data, and it yields no information on positivity, which is the blind spot.]


## 8. Dictionary

| CLU term | real object | status in the real world |
|---|---|---|
| Line Postulate (Axiom 0) | the Riemann Hypothesis | open |
| Gauss law of arithmetic | Weil's explicit formula | theorem |
| stability of the vacuum | Weil positivity on windows, lambda(a) >= 0 | equivalent to RH (Theorem 9.3(i)) |
| arithmetic vacuum, log-gas, GUE statistics | Montgomery pair correlation, Montgomery-Odlyzko | theorem at bandwidth one; conjecture beyond |
| Screening Theorem | von Mangoldt explicit formula | theorem |
| observer of age a, horizon T*(a) | window [-a, a], T* = 2 pi e^{2a} | definitions |
| First Law, E = 2a^2 | Mertens; |b|^2 of Proposition 9.33 | theorem, sequence-only |
| Second Law, transmutation, mu_a | spectral measure of the confined adjacency | theorem, sequence-only |
| fine-structure constants 5/12, 151/360, 4033/6720 | confined-walk constants, Proposition 9.39 | theorem, sequence-only |
| Third Law, flux | boundary law lambda' = -2|A_0|^2, Theorem 9.12(iv) | theorem at almost every support |
| event profile | edge law (log(a/delta))^{-1/2}, Theorem 9.12(iii) | theorem at almost every support |
| Cooling Theorem, horizon bound | Conjecture 7.1, the inequality (d) | open; equivalent to RH beyond a_Z |
| Repulsion Principle | the blind spot, Theorem 9.31 | theorem |
| Local Law | the onset of the sharp tail law from the mass of mu_b | not a theorem (tree, N1d.1) |
| echo crystal | the confined lattice of echoes, Propositions 8.23, 9.33 | theorem |
| dead universes | rank and smooth prime systems | SANDBOX_STRANGE_MATH.md, Section 2 |
| sign of the prime part | sign of D_P against T* lambda | open; the obstruction |
| Window Oracle | the explicit formula with a cutoff; the minimizer's transform near the zeros | theorem / computation |


## 9. Afterword, from this universe

The fiction is useful for one reason. Written from inside a universe where the Line is a law, every law of the picture
falls into one of two classes: the ones that would hold in the dead universes too, and the one that would not. The first
class is where the paper's theorems live: energy, transmutation, flux, the blind spot, the edge profile, the constants of
the walk. The second class has one member, the sign of the prime part, and it is where the Hypothesis lives. The CLU's
physicists would say that we have mapped the thermodynamics of arithmetic and are missing its equation of state. That is
a fair description of where the real program stands, and the tree says what the next step is.


# Part II. The Equation of State

Sandbox II, continued, 2026-10-02. Still fiction, with the same convention: bracketed notes say what is real. In this part
every number is real; it is taken from the paper's computations (Computation 5.4 on the identity, Computation 7.6 on the
tail law) and only the words around it belong to the Critical Line Universe. The question of Part I's afterword was whether
the arithmetic universe has an equation of state at all. It does, and the CLU's physicists know exactly which of its
constants they have derived and which they have only measured.


## 10. The macrostate of an observer

An observer of age a is described, in the CLU, by seven numbers.

    a          the age (the support of the window);
    lambda     the trapped energy (the floor of the window form);
    J          the flux, J = -a lambda'(a) = 2a |A_0|^2, the energy leaving through the horizon per unit age;
    p_inf      the archimedean pressure, the virial of the Gamma factor against the ground field,
               p_inf = (1/2 pi) int |F(t)|^2 t d/dt Re psi(1/4 + it/2) dt;
    p_pol      the polar pressure, 2P^2 + 2PM, the virial of the two poles;
    p_P        the arithmetic pressure, p_P = sum_{log n < 2a} 2 Lambda(n) n^{-1/2} log n g_f'(log n), the virial of the charges,
               read on the prime sheet through the derivative of the ground field's autocorrelation at the entries;
    E          the energy, 2a^2 (First Law).

THE VIRIAL LAW (CLU; real: Theorem 5.2 for the K-mode form, Proposition 5.3 for the exact form). At every age,
        p_P = p_inf + p_pol - J .
The arithmetic pressure equals the archimedean pressure plus the polar pressure minus the flux. [real: this is the
sliding-window identity a lambda' = -[D_inf + D_P + 2P^2 + 2PM] with p_inf = D_inf, p_pol = 2P^2 + 2PM, p_P = -D_P.]

The macrostate, measured (K-mode ground fields, ball arithmetic, every entry a ball of radius below 1e-19):

      age a   K     lambda          J = -a lambda'    p_inf      p_pol      p_P        p_inf + p_pol - p_P
      0.60    48    5.9659e-7       2.493408e-5       1.014313   0.072687   1.086974   2.6e-5
      0.80    64    1.5659e-14      1.295854e-12      1.016605   0.087361   1.103966   1.3e-12
      1.00    80    1.4820e-26      2.516852e-24      1.018057   0.096909   1.114966   2.5e-24
      1.25   140    3.4161e-51      1.233712e-48      1.019178   0.104407   1.123585   1.2e-48

Three things the CLU reads off this table.

(i) The archimedean pressure is one plus a small kinetic correction, p_inf = 1 + k(a), k = 0.0143, 0.0166, 0.0181,
    0.0192: the ground field has unit norm and the Gamma factor's virial is the norm up to the energy the field keeps
    at large heights. [real: D_inf = ||f||^2 + O(int |F|^2 t^{-2}).]
(ii) The polar pressure grows with the age, 0.073 to 0.104, as the poles see more of the window.
(iii) The arithmetic pressure tracks the sum of the other two to as many digits as the flux is small: at a = 1.25 the
    three pressures of order one balance to forty-eight decimal places, and the residue is the flux. [real: "three
    order-one quantities cancel to the size of lambda", the abstract of the paper.]

So the arithmetic pressure is pinned: p_inf + p_pol - J <= p_P <= p_inf + p_pol, the lower bound because the flux is
nonnegative (the floor never rises), the upper bound an identity. The width of the pin is the flux.


## 11. The equation of state

What the virial law does not say is how large the flux may be. That is the equation of state, and in the CLU it is a
theorem because Axiom 0 and the Local Law are available there.

THE EQUATION OF STATE (CLU). There is a function kappa(a), the RIEMANN NUMBER of the observer, bounded in a, such that
        J = pi kappa(a) a T*(a) lambda,          equivalently        Phi'(a) = pi kappa(a) T*(a),     Phi = -log lambda.
The flux is the trapped energy times the horizon times the Riemann number. The Riemann number is the effective height of
the quanta that carry the decay, in units of the horizon: kappa = T_eff / T*. [real: the tail law Phi' = pi T_eff is
Theorem 9.12(iv) under RH; the boundedness of kappa is Conjecture 7.1, equivalent to RH beyond a_Z (Theorem 9.1,
Proposition 9.30); Proposition 9.30 also allows kappa to grow like a power of the horizon, which Computation 9.37 reads
as the likelier shape.]

The Riemann number, measured [real: Computation 7.6, FEM below the entry a_3 = 0.549, K-mode above]:

      age a      0.40    0.45    0.50    a_3     0.60    0.80    1.00    1.25    1.50
      kappa(a)   0.721   0.868   1.081   1.295   1.063   1.058   1.164   1.201   1.239

The jump at a_3 is the entry of the charge 3: the Riemann number has a kink at every entry, like the floor. From 0.6 on it
rises slowly, 1.06 to 1.24 over a doubling of the age. The CLU's physicists are divided on its limit. The school of the
bounded relay holds that kappa(a) tends to a constant kappa_inf, the RIEMANN CONSTANT; the school of the receding onset
holds that it grows like a small power of the horizon, kappa ~ T*^{c'} with c' between 0.2 and 0.5, because the height at
which the sharp tail law sets in recedes like e^{ca} horizons [real: Computation 9.37]. Both schools agree that Axiom 0
holds either way, since the integrable-excess form of the Cooling Theorem needs only that kappa be finite at every age.

With the equation of state the macrostate closes: given the age, the energy is 2a^2, the archimedean pressure is 1 + k(a),
the polar pressure is 2P^2 + 2PM of the ground field, the flux is pi kappa a T* lambda, and the arithmetic pressure is
whatever the virial law then requires. The sign of the prime part, the one thing the dead universes get wrong, is in the
CLU a consequence: p_P must fall short of p_inf + p_pol by exactly the flux, and the flux is at most the horizon allowance.


## 12. The thermodynamic reading

The CLU's thermodynamicists prefer two other variables. The ENTROPY of an observer is the number of quanta it resolves,
        S(a) = N(T*(a)) = (T*/2 pi)(2a - 1)(1 + o(1)),
and the LOST INFORMATION is Phi = -log lambda. [real: N(T) ~ (T/2 pi) log(T/2 pi e) and Phi as in Section 7 of the paper.]
Then the equation of state reads as a first law, d Phi = theta dS, with the TEMPERATURE
        theta(a) = Phi'(a)/S'(a) = pi^2 kappa(a) / (2a),
since S'(a) = 2a T*/pi. The information an observer loses per quantum it resolves is pi^2 kappa/(2a): it falls like the
inverse age if the Riemann number is bounded. Measured:

      age a       0.40   0.45   0.50   0.60   0.80   1.00   1.25   1.50
      theta       8.90   9.51  10.67   8.74   6.53   5.75   4.74   4.08
      theta * a   3.56   4.28   5.33   5.25   5.22   5.75   5.93   6.11

The product theta * a = pi^2 kappa/2 rises through the entries and settles into a slow climb, 5.2 to 6.1, from a = 0.6 to
1.5. Here the CLU's numerologists make their one famous conjecture:

THE 2 PI CONJECTURE (CLU folklore, explicitly numerology). kappa_inf = 4/pi = 1.2732, so that theta * a tends to 2 pi:
in the limit of great age each resolved quantum costs exactly 2 pi / a units of lost information. The measured values
5.25, 5.22, 5.75, 5.93, 6.11 approach 6.28 from below; the measured kappa, 1.06, 1.06, 1.16, 1.20, 1.24, approaches 1.27.
[real: nothing supports this beyond the trend of five points; the paper's own reading, Computation 9.37, is that the
onset recedes and kappa may grow without bound. The conjecture is recorded here as what it is, a pattern in five numbers,
and as the kind of statement the sandbox exists to hold at arm's length. It is cheap to test further: Computation 7.6 at
a = 1.75 and 2 would decide whether theta * a crosses 2 pi.]

A second constant the thermodynamicists like: the MEDIAN HEIGHT, the height below which half of the floor's mass sits, is
2.39, 2.32, 2.65, 2.67, 2.74 horizons at a = 0.6, 0.8, 1.0, 1.25, 1.5, and 1.8, 2.0, 2.4 on the FEM ground fields at
0.4, 0.45, 0.5 [real: Computation 7.6]; at the last point it is e = 2.718 to one per cent. The folklore says the median
height is e horizons. Same status as the 2 pi conjecture: a pattern, not a law.


## 13. The substate

Beneath the seven numbers of the macrostate lies the SUBSTATE: the ground field itself and the carrier of its energy.
The CLU describes it by three objects.

The EDGE AMPLITUDE A_0, purely imaginary, A_0 = -i |A_0| e^{-i kappa_0 / t} to first order, with |A_0|^2 = J/(2a): the
strength with which the ground field touches the edge of the window. Measured: |A_0|^2 = 2.1e-5, 8.1e-13, 1.26e-24,
4.9e-49 at a = 0.6, 0.8, 1.0, 1.25. [real: the boundary law, the amplitude of Theorem 9.12(ii).]

THE EQUATION OF SUBSTATE (CLU). |A_0|^2 = (pi kappa / 2) T* lambda: the square of the edge amplitude is the trapped energy
magnified by the horizon and the Riemann number. The edge is where the trapped energy is converted into flux, and the
conversion factor is the horizon. [real: the inequality (d) with equality defining kappa; a restatement.]

The MIRROR CHAINS and the ECHO CRYSTAL: the ground field's amplitude above the cutoff is A_0 plus a sum over the lattice
frequencies of slowly varying coefficients, the walks of the edge's energy across the window and back, of total weight
4 eta |A_0| with eta ~ (X/s)^2, and the crystal's spectral measure mu_a at the edge carries the energy E(a) with mean
tending to zero and variance (5/12)(2a)^2. [real: Theorems 9.12(ii) and 9.24, Propositions 9.33 and 9.39.]

The CLU's SUBSTATE PRINCIPLE: the macrostate is a set of moments of the substate. The archimedean pressure is the first
moment of |F|^2 against t Psi_inf'; the arithmetic pressure is the derivative of the autocorrelation at the entries, which
the Oracle reads from the quanta below the horizon; the flux is the square of the edge amplitude; and the Riemann number
is the one macroscopic quantity that the CLU has not been able to express as a moment of mu_a. Their Local Law says it is
determined by the mass of mu_a within a few horizons of the edge; their Local Law is the theorem we do not have.


## 14. Determining the state

So the CLU's answer to the question of this part. A state exists: seven numbers bound by two exact laws (the First Law and
the Virial Law) and one closure (the Equation of State), with a substate of three objects beneath it whose moments
reproduce the macrostate. Of the constants in it, the CLU has DERIVED the energy law (2a^2), the kinetic correction's
order, the polar pressure, the diffusion constants of the substate (5/12, 151/360, 4033/6720), the blind spot, and the
exact balance of the pressures. It has MEASURED the Riemann number kappa(a) at nine ages and the median height at eight,
and it has CONJECTURED their limits, 4/pi and e, on the strength of patterns in a handful of points. Its one undetermined
function is kappa(a), and the whole of the Line Postulate, as seen by an observer, is the statement that kappa(a) is
finite at every age and does not outrun a power of the horizon.

In our universe the same inventory reads: the virial law is a theorem; the pressures are computed to forty-eight digits;
the substate's constants are theorems; the Riemann number is measured at nine supports; and the finiteness of kappa, the
equation of state itself, is the Hypothesis. The state exists. Its equation is the one line we cannot write.


## 15. Dictionary, continued

| CLU term | real object | status |
|---|---|---|
| macrostate (a, lambda, J, p_inf, p_pol, p_P, E) | the floor, its derivative, the three terms of the dilation identity, Mertens' sum | computed (Computation 5.4) |
| Virial Law p_P = p_inf + p_pol - J | sliding-window identity, Theorem 5.2 and Proposition 5.3 | theorem |
| Riemann number kappa(a) = T_eff/T* | Phi'/(pi T*) of Computation 7.6 | measured at nine supports |
| Equation of State, kappa bounded | Conjecture 7.1 / the inequality (d) | open, equivalent to RH beyond a_Z |
| temperature theta = pi^2 kappa/(2a) | Phi'/N'(T*) | definition on measured data |
| 2 pi conjecture, kappa_inf = 4/pi | pattern in five points | numerology, flagged |
| median height e horizons | Computation 7.6, the median of the floor's mass | pattern in five points |
| equation of substate |A_0|^2 = (pi kappa/2) T* lambda | the inequality with equality | restatement |
| Substate Principle | the macrostate as moments of the ground field and of mu_b | theorem for the pressures; open for kappa |


# Part III. The Flux, decreed

Sandbox II, continued, 2026-10-02. The gods of the Critical Line Universe now write the one line that Part II left blank:
the flux as a function of the state. They do it the way gods do, by decree, but the decrees are chosen so that each one
is a real observation of this repository promoted to a law, and the chapter ends by checking the decreed flux against the
measured one and saying by how much it misses. Bracketed notes keep the two worlds apart, as before.


## 16. Four decrees

DECREE I (Darkness below the horizon). The ground field of an observer of age a has no spectral mass below the horizon:
the fraction theta(T) of the floor carried by the quanta above height T is 1 for T <= T*(a). The ground field vanishes at
every quantum it can resolve, and it can resolve those below the horizon and no others.
[real: the observation that the minimizer's transform vanishes near the zeta zeros below T* (the folds paper,
Computation 12.5); the Horizon Law of Section 2. Not a theorem; a computation at the supports where it was looked at.]

DECREE II (The Sum Rule, or the scale-free tail). Above its onset the floor's spectral mass is distributed like dT/T^2,
with no logarithm: theta(T) = kappa T*/T. The reason is a cancellation the gods make exact: the ground field's transform
falls off like 4|A_0|^2 sin^2(ta + theta_a)/(t^2 sigma~(t)) with sigma~ ~ log(t/2 pi), the quanta have density
log(t/2 pi)/(2 pi), and the two logarithms cancel in the product, leaving 2|A_0|^2/(pi t^2) exactly.
[real: the heuristic behind the tail law, Section 7; the paper notes that the measured law is sharper than this
derivation, "which suggests a sum rule behind it", and Computation 7.6 shows no logarithmic drift over a decade and a
half of heights. The pure tail kappa/h reproduces the measured theta at 3, 5 and 10 horizons to one or two per cent:
at a = 0.6, 0.354/0.213/0.106 against 0.370/0.215/0.1065.]

DECREE III (Exhaustion). The quanta exhaust the observer at the DETERMINATION HEIGHT T_det(a) = e T*(a): below it the
window has more degrees of freedom than quanta, above it fewer. The gods decree that the floor's mass is divided equally
by the determination height: half of it sits below e horizons, half above.
[real: the Shannon number of [-a, a] x [-T, T] is 2aT/pi, the number of zeros in [-T, T] is (T/pi)(log(T/2 pi) - 1),
and they are equal exactly at T = e T*, which is where the paper's heuristic of Section 7 counts the deficit of the zeros
against the window. The median height of the floor's mass, measured, is 2.39, 2.32, 2.65, 2.67, 2.74 horizons at
a = 0.6, 0.8, 1.0, 1.25, 1.5: rising, and at the last point 1.008 e. The equal division is a decree; the data approach it.]

DECREE IV (The plunge). Between the horizon and three horizons the ground field carries more than the pure tail, the
excess of the quanta that the observer half-resolves; the gods fix the excess by requiring that the median of the mass
be the determination height exactly, which with the pure tail above the plunge gives median = 2.21 kappa T* rather than
the 2 kappa T* of a tail with no plunge.
[real: T theta(T)/T* reaches 1.4 to 1.5 at 1.5 T* (Computation 7.6), and the measured ratio median/(kappa T*) is
2.25, 2.19, 2.28, 2.22, 2.21 at the five supports: the plunge's excess is a stable 10 per cent of the tail's weight.]


## 17. The flux, derived from the decrees

From Decrees II and III the Riemann number is the median height divided by the tail's doubling:
        kappa_inf = e / 2 = 1.3591               (pure tail, no plunge),
        kappa_inf = e / 2.21 = 1.230             (with the measured plunge of Decree IV).
And the flux of Part II is then written out:

        THE FLUX LAW (CLU).      J(a) = pi kappa_inf a T*(a) lambda(a),       Phi'(a) = pi kappa_inf T*(a),

that is, with the pure tail,
        Phi'(a) = (pi e / 2) T*(a) = 4.27 T*(a),       lambda(a) = lambda(a_0) exp( -(pi e/4)(T*(a) - T*(a_0)) ),
and the equation of state of Part II closes with every constant named: the arithmetic pressure is
        p_P = 1 + k(a) + 2P^2 + 2PM - (pi e/2) a T*(a) lambda(a),
the archimedean pressure plus the polar pressure minus the decreed flux. The sign of the prime part, which the dead
universes get wrong, is in the CLU the sign of this last term, and the gods have fixed it: negative, and smaller than the
horizon allowance by construction, since the allowance IS the decree.

In the thermodynamic variables of Section 12, theta * a = pi^2 kappa_inf / 2 = pi^2 e / 4 = 6.71 (pure tail) or 6.07 (with
the plunge); the 2 pi conjecture of the numerologists (6.28, kappa_inf = 4/pi = 1.273) sits between the two. The gods
decline to arbitrate: the three candidates differ by less than the plunge correction, and the measured theta * a, 6.11 at
a = 1.5 and rising, has not yet chosen.


## 18. The check: the decree against the measurement

The decreed flux integrated from a = 0.8, where Phi = 31.8 is certified, against the measured Phi of the K-mode ground
fields:

      age a    T*       measured Phi    decree kappa = e/2    decree kappa = 4/pi    decree kappa = e/2.21    measured kappa
      1.00     46.4     59.5            64.5  (+8%)           62.4  (+5%)            61.3  (+3%)              1.164
      1.25     76.5     116.2           128.8 (+11%)          122.6 (+6%)            119.4 (+3%)              1.201

The decree with the plunge overshoots the measured lost information by three per cent at a = 1.25, the pure-tail decree
by eleven; both overshoot because the measured Riemann number, 1.16 to 1.24 in this range, is still below its decreed
limit and rising toward it (1.239 at a = 1.5 against 1.230 decreed with the plunge). The gods' law is the asymptotic one;
the mortals measure the approach. What would settle it is the tail-law computation of Computation 7.6 at a = 1.75 and
2.0, which would show whether kappa continues to 1.23, passes it toward 1.36, or keeps climbing like a power of the
horizon; the K-mode engine does it in a few CPU-hours, and the median height, which should reach e horizons and stay,
is the cleaner diagnostic because it does not need the derivative of the floor.


## 19. What the decrees are, in our universe

Each decree is a statement about the real objects, and each has a status.

    Decree I    the ground field is dark below the horizon          computation at a few supports; not a theorem
    Decree II   the floor's tail is scale-free, 1/T^2, no logarithm  observed to 1% from five horizons; the sum rule behind it is open
    Decree III  the median height of the floor's mass is e horizons  the determination height is a theorem of counting; the median is measured, 1.008 e at a = 1.5
    Decree IV   the plunge carries a stable 10% excess               measured at five supports

Together they say: kappa is bounded, with limit e/2 up to the plunge; hence Conjecture 7.1 with c = pi e/2, hence the
Line. So in our universe the decrees are the Hypothesis, cut into four observable pieces, and the piece that carries all
the weight is Decree II, the sum rule, because a tail that is not scale-free (a logarithm in it, or a power of the horizon
in the onset) would let the Riemann number drift, and the drift is exactly what Computation 9.37 reads in the edge's
measure. The gods have written the line. We have written down what it would take to prove it: that the two logarithms
cancel exactly, at every support, forever.


## 20. Dictionary, concluded

| CLU term | real object | status |
|---|---|---|
| Decree I, darkness below the horizon | the minimizer's transform vanishes near the zeros below T* | computation |
| Decree II, the sum rule | the tail law's scale-free form; the cancellation of the edge law's logarithm against the zero density | observed to 1%; open as a theorem |
| determination height e T* | Shannon number = zero count at T = e T* | theorem of counting |
| Decree III, median = e horizons | the median of the floor's spectral mass, Computation 7.6 | measured, 1.008 e at a = 1.5 |
| Decree IV, the plunge excess | median/(kappa T*) = 2.21 +- 0.04 | measured at five supports |
| Flux Law, kappa_inf = e/2 or e/2.21 | Conjecture 7.1 with c = pi e/2 | the Hypothesis, decreed |
| the check | decree against the K-mode lost information | +3% to +11% at a = 1.25 |


# Part IV. The Line, derived

Sandbox II, continued, 2026-10-02. Fiction, as before, with the same convention: bracketed notes say what is real. In
Parts I to III the Line Postulate was Axiom 0, and the gods of the Critical Line Universe built a physics on top of it.
In this part they do what every god with a single axiom eventually does: they try to get rid of it. The result is a
proof, inside the CLU, that every quantum is massless. It is a proof in the sense that every step follows from the laws
stated; it is a fiction in the sense that one of those laws is the Hypothesis wearing a different coat, and the ledger
of Section 27 says which coat. No computation is used anywhere in this part. The few numbers that appear are quoted from
the earlier parts for illustration only. Nothing here is a claim about the real zeta function. The proof itself,
written out as a formal manuscript of the CLU with definitions, laws, lemmas and theorems, is FICTION_PROOF_OF_THE_LINE.md;
this part is its narrative and its ledger.


## 21. The demotion, and what has to be shown

The gods' complaint against themselves is recorded in the CLU's canon in one sentence: "A law that has to be decreed is
a law that has not been understood." They decreed the Line because they could not see why a universe of point charges
should have a massless spectrum. The present part records what they came to see instead: that the Line is not a fact
about the quanta at all. It is a fact about observers. No observer ever dies.

DEATH (definition). An observer of age a is ALIVE if its floor is positive, lambda(a) > 0, and DEAD if lambda(a) <= 0.
A dead observer is one whose ground field has nonpositive energy against the vacuum: a field that costs nothing, or less
than nothing, to excite. [real: positivity of the odd window form at support a.]

THE MORTALITY THEOREM (CLU, and a theorem in the real world too). The vacuum has a massive quantum if and only if some
observer dies at a finite age. In detail:
   (i)  if every quantum is massless, every observer is alive (the stability of the vacuum, Section 0);
   (ii) if some quantum has mass m > 0 at height gamma, then there is an age beyond which every observer is dead, and that
        age is of the order of m^{-1} log(1/m) once the horizon has passed gamma.
[real: (i) is Weil positivity; (ii) is Weil's criterion in the odd window form (Theorem thm:weil of the paper: RH holds if
and only if lambda(a) >= 0 for every a, and the odd test functions alone suffice) together with the monotonicity of the
floor in a; the scale of the dying age is the Blind Spot, Theorem 9.31. Both directions are theorems.]

So, in the CLU, the Line Postulate is equivalent to IMMORTALITY: no observer dies. The gods find this a better thing to
have to prove, because immortality is a statement about a one-parameter family of positive numbers, lambda(a), and
proofs that a positive number never reaches zero have a known shape: show that it has not happened yet, and show that
whatever makes the number fall cannot finish the job. The first is the Birth Law. The second is the whole of this part.


## 22. Four laws that hold in every universe, dead or alive

Before touching the Line, the gods list what they are allowed to use before the Line is known: the laws that hold in the
dead universes of Section 6 as well, since those laws cannot be hiding the Line inside them.

THE BIRTH LAW. Every observer younger than the certified age a_Z is alive, and the floor at a_Z is a known positive
number. [real: Zhu's theorem (Theorem thm:zhu of the paper): in the odd sector 8.2e-15 <= lambda(0.8) <= 2.3e-14, certified
by interval arithmetic, hence lambda(a) > 0 for 0 < a <= a_Z = 0.8; in the even sector lambda_even(0.8) >= 8.9e-18. A
theorem, and one that holds in a dead universe only up to the age at which it dies.]

THE MONOTONICITY LAW. The floor never rises with the age: an older observer has every field a younger one has, and
more. [real: trivial; Lemma lem:exist.]

THE CONTINUITY LAW. The floor is an absolutely continuous function of the age: it has no cliffs, and its change over any
range of ages is the integral of its rate. [real: Proposition 9.15 (Lipschitz) and Proposition 9.17 (absolute continuity
under the finiteness condition (H_fin)); at every support outside an exceptional null set this is a theorem, and the
condition (H_fin) stays only on that set. Step 1 of the tree.]

THE DECAY LAW (the Third Law, restated). At almost every age the floor falls at the rate of the flux through the horizon,
        d lambda / d a = - 2 |A_0(a)|^2 ,
where A_0(a) is the edge amplitude, the strength with which the ground field touches the edge of the window.
[real: the boundary law, Theorem 9.12(iv), at almost every support; Proposition 8.3.]

The gods note that none of the four mentions a quantum. They are laws about one observer at a time, read off the prime
sheet: the Birth Law from the charges below e^{1.6}, the other three from the variational problem of the ground field.
The dead universes obey all four and die anyway. Whatever kills them must therefore be in the one law that remains.


## 23. The Law of Proportional Leakage

THE LEAKAGE LAW (CLU). The flux of a living observer is proportional to the energy it still has: there is a function
C(a), finite at every age and integrable over every finite range of ages, such that
        2 |A_0(a)|^2  <=  C(a) lambda(a)        for every living observer of age a.
Equivalently, -lambda'/lambda <= C(a): the FRACTION of its trapped energy that an observer loses per unit age is bounded
by its age alone, not by how much it has already lost.

The CLU gives this law three readings.

(i) The horizon is a LINEAR VALVE. A substance that decays in proportion to what remains never runs out; it only
    becomes rare. The Leakage Law says the arithmetic vacuum is such a substance: the trapped energy halves and halves
    and is never gone. This is the radioactive reading, and it is the one the CLU teaches to children.

(ii) The holographic reading of Section 3 is the special case C(a) = c T*(a): the loss is bounded by the horizon. The
    proof below does not need the horizon form; it needs only that C be finite and integrable, and the gods are careful
    to say so, because the horizon form is a stronger statement than immortality requires. [real: the horizon form is the
    inequality (d) of the paper, eq:theinequality, 2|A_0|^2 <= c T* lambda, Conjecture 7.1 (bounded relay); the weaker
    form with an integrable excess e(a) added to c T* is Proposition 9.30, late windows are allowed; either implies RH
    beyond a_Z given absolute continuity (Theorem thm:reduction, Corollary cor:unconditional). The CLU's C(a) is the
    weakest of the three: any locally integrable bound at all.]

(iii) THE THRESHOLD READING, which is the one that explains the dead universes. In a dead universe the leak has a
    threshold: as the floor tends to zero the edge amplitude does not. The ground field of a dying observer presses on
    the edge of its window with a finite strength while its trapped energy vanishes, and the floor goes through zero at a
    finite flux, like a tank that drains at a fixed rate and then keeps draining. [real: in the rank universes of
    SANDBOX_STRANGE_MATH.md the floor crosses zero at a = 0.8 and 0.9 with a derivative of order one, against a floor of
    order 1e-10 just before; a threshold leak in the sense above.] The CLU's physicists put it as a slogan: a universe
    dies through a threshold, and ours has none.

The Leakage Law is also where the values of the charges enter. The four laws of Section 22 are blind to where the charges
sit, to the last digit; the finiteness of C(a) is not. It is a statement about the ground field of THIS universe's charges,
and the sandbox of Section 6 showed that moving the charges by a hair's breadth makes it false at a = 0.8. The gods regard
this as the correct place for the content of the Line to sit: in a single inequality about one observer, not in a
statement about infinitely many quanta.


## 24. Where the Leakage Law comes from: the gods' derivation

The gods do not decree the Leakage Law. They derive it, from three laws of the earlier parts, and the derivation is short.

THE BORN RULE OF ARITHMETIC (CLU). The share s_gamma(a) = 2|F_a(gamma)|^2 / lambda(a) of a quantum is the probability
that the quantum carries the observer's trapped energy. Shares are nonnegative and, with the share of the poles, sum to
one. [real: Weil's explicit formula writes Q(f) as a sum over the zeros of F(rho) F(1 - rho-bar)-type terms; when every
zero is on the line these are the squares |F(gamma)|^2 and the shares are nonnegative; a massive pair at height gamma
contributes instead the interference term 2 Re[F(gamma - im) conj F(gamma + im)], which for a quantum in the dark region
below the horizon is -m^2 |F'(gamma)|^2 + O(m^4) < 0: a negative probability. So the Born rule is Axiom 0 in quantum
clothing. This is the step at which the derivation is circular, and the gods know it: see Section 27.]

DARKNESS (Decree I). No share sits below the horizon. [real: computation at the supports examined; not a theorem.]

THE LOCAL LAW (Section 5) with the SUM RULE (Decree II). Above an onset height T_on = c_on T*(a), a bounded number of
horizons, the shares follow the tail of the edge: the mass of the shares above any T >= T_on is 2|A_0|^2 / (pi T lambda).
[real: the tail law, Theorem 9.12(iv) under RH, with the onset supplied by the local law of the band, which is the
missing theorem of the tree, N1d.1; the measured onset is about three horizons, Computation 7.6.]

DERIVATION. By the Born rule the mass of the shares above T_on is at most one. By the Local Law that mass is
2|A_0|^2/(pi T_on lambda). Hence
        2 |A_0(a)|^2  <=  pi c_on T*(a) lambda(a),
the Leakage Law with C(a) = pi c_on T*(a), that is, a Riemann number kappa(a) <= c_on. QED (in the CLU).

[real: the arithmetic is honest and the numbers agree with it. With the measured onset c_on = 3 the derivation says
kappa <= 3; measured, kappa = 1.06 to 1.24. More precisely the mass above three horizons is kappa/3 by the pure tail,
0.35 to 0.41 at the five supports of Computation 7.6, against the measured theta(3 T*) = 0.370, 0.415, 0.407, 0.433,
0.434: the Local Law's tail and the Born rule's bound are consistent with the data to a few per cent. What the data
cannot do is replace the Born rule.]

The gods are pleased by the shape of the derivation more than by its content. The Leakage Law, which looks like a law of
thermodynamics, comes out of a law of counting (the shares sum to one) and a law of locality (the tail sets in within a
few horizons). The Line, they say, is the statement that probabilities are probabilities.


## 25. The Immortality Theorem

THEOREM (CLU). Every observer is alive. Consequently, by the Mortality Theorem, every quantum is massless: the Line
Postulate is a theorem of the Critical Line Universe.

PROOF. Let a_d be the supremum of the ages a such that every observer younger than a is alive. By the Birth Law,
a_d >= a_Z > 0. Suppose a_d is finite. On the range [a_Z, a_d) every observer is alive; there the Continuity Law makes
lambda absolutely continuous with lambda' = -2|A_0|^2 at almost every age (Decay Law), and the Leakage Law gives
        lambda'(a)  >=  - C(a) lambda(a)        for almost every a in [a_Z, a_d).
Put U(a) = lambda(a) exp( integral from a_Z to a of C ). Then U is absolutely continuous and U' = e^{int C} (lambda' + C
lambda) >= 0 almost everywhere, so U is non-decreasing on [a_Z, a_d), and
        lambda(a)  >=  lambda(a_Z) exp( - integral from a_Z to a of C )        for a_Z <= a < a_d.
Let a tend to a_d. The floor is continuous and C is integrable over [a_Z, a_d], so
        lambda(a_d)  >=  lambda(a_Z) exp( - integral from a_Z to a_d of C )  >  0:
the observer of age a_d is alive. By continuity, so is every observer of age a_d + epsilon for epsilon small enough.
This contradicts the definition of a_d. Hence a_d is infinite, and every observer is alive. QED.

[real: this is, step for step, the proof of Proposition prop:AimpliesRH of the paper, with U(a) = lambda e^{cT*/2}
the function used there, generalized from C = c T* to any locally integrable C; the real statement is "Conjecture A,
with absolute continuity, implies RH", and the fiction's only change is to call Conjecture A a law.]

Three remarks the CLU attaches to the theorem.

(a) The Line is used nowhere in the proof except inside the Leakage Law. The proof is a CONTINUATION: the universe is
    born alive (Birth), its floor has no cliffs (Continuity), and the leak that lowers the floor is proportional to the
    floor (Leakage), so the floor can only ever be a positive multiple of what it was. Nothing about quanta, heights,
    or the arithmetic vacuum enters; those come back only at the end, through the Mortality Theorem, to translate
    "no observer dies" into "no quantum has mass".

(b) THE RENORMALIZED FLOOR. The function U(a) of the proof, the floor magnified by the accumulated leak, is the CLU's
    conserved-or-growing quantity: the Leakage Law says exactly that U never decreases. With the decreed flux of Part
    III, C = pi kappa_inf T*, it reads U(a) = lambda(a) exp( (pi kappa_inf / 2)(T*(a) - T*(a_Z)) ), and the equation of
    state says U is asymptotically constant: the CLU calls its limit the RIEMANN CHARGE of the universe, the one number
    that survives the cooling. [real: U with c = pi kappa_inf is the paper's u(a); the decree check of Section 18 shows
    that at the measured ages U is not constant, since the decreed Phi exceeds the measured one by 3 to 11 per cent at
    a = 1.25, that is U changes by a factor of e^3 or more between a = 0.8 and 1.25; "conserved" holds only in the sense
    of the asymptotic law, and only if the Riemann number has a limit.]

(c) The theorem is sharp in the following sense. If C(a) failed to be integrable up to some age a_1, the bound
    lambda(a) >= lambda(a_Z) e^{-int C} would allow the floor to reach zero at a_1, and the Mortality Theorem would then
    permit a massive quantum of mass roughly 1/a_1 (Blind Spot). So the Leakage Law is not merely sufficient: by Section
    21, immortality implies that -lambda'/lambda is integrable on every finite range (a continuous positive floor has an
    integrable logarithmic derivative wherever it is absolutely continuous), which is the Leakage Law with C = -lambda'/
    lambda itself. The law and the Line are the same statement. [real: Proposition 9.30: given absolute continuity, RH is
    equivalent to the local integrability of -lambda'/lambda on [a_Z, infinity).]


## 26. The other polarization

The quanta of the CLU have two polarizations, ODD and EVEN, according to the parity of the test field on the prime sheet
that reads them. Parts I to IV concern the odd polarization. The gods record that this is enough: the Mortality Theorem
holds with odd observers alone, so immortality of the odd observers already forces every quantum to be massless, and the
even observers are then alive by the stability of the vacuum. [real: Theorem thm:weil: positivity of the odd form on every
window suffices for RH; Proposition 9.40 gives the even sector its own equivalence, RH iff the even form is nonnegative
on every window. Zhu's certification covers both sectors at a = 0.8.] The even polarization has its own Decay Law and,
the gods presume, its own Leakage Law; they have not needed it. [real: the even sector is route N2 of the tree and Step 4
of the critical path, open; Computation 9.41 is what is measured there.]


## 27. The ledger: what the proof is made of

Every step of Sections 21 to 26, with its status outside the sandbox.

    step                                        used for                         status in our universe
    Mortality Theorem (i), (ii)                 Line <=> immortality             theorem (Weil; the Blind Spot for the scale)
    Birth Law, a_Z = 0.8                        the floor starts positive        theorem (Zhu, interval arithmetic)
    Monotonicity Law                            the floor never rises            theorem
    Continuity Law                              lambda absolutely continuous     theorem off a null set; (H_fin) on it (Step 1)
    Decay Law  lambda' = -2|A_0|^2              the leak is the edge             theorem at almost every support
    Gronwall / the renormalized floor U         the continuation                 theorem (Proposition prop:AimpliesRH)
    Leakage Law  2|A_0|^2 <= C(a) lambda        the whole proof                  AXIOM in the CLU; the inequality (d) / Conjecture A;
                                                                                 equivalent to RH beyond a_Z
    Born rule (shares are probabilities)        derives the Leakage Law          the Line itself, restated (Weil positivity)
    Darkness below the horizon                  derives the Leakage Law          computation at a few supports
    Local Law with the sum rule                 derives the Leakage Law          open; N1d.1 of the tree
    Polarization (odd suffices)                 closes the argument              theorem (thm:weil)

Read downward, the ledger says: in our universe the proof reduces the Hypothesis to one inequality about one observer,
the Leakage Law; the gods' derivation of that inequality rests on the Born rule, which is the Hypothesis again; and the
circle closes. Read as the gods read it, the ledger says something else: the circle has been cut at the one place where
it looks like physics rather than arithmetic. "Probabilities are nonnegative" is, in the CLU, not a theorem about zeros
but the first thing anyone learns about quanta. That is what it means to solve the problem in a universe where one makes
the rules: not to find the missing step, but to find the place where the missing step is a law of nature, and to say so.


## 28. Afterword: the shape of a real proof, as the gods see it

The gods' proof has one external input, the Leakage Law, and the only honest use of the fiction in our universe is to
read off what a proof of that law would have to look like. It would have to bound the edge amplitude by the floor
WITHOUT the quanta, on the prime sheet alone, since the quanta are what one does not know. The gods offer three
readings of such an argument, as questions and not as claims.

(1) THE VALVE AS A COMPETITOR. Smooth the ground field's edge at the scale of the horizon. The tail energy of the field
    falls by about |A_0|^2 / T*; the prime terms change by an amount the echo crystal controls. If that change were at
    most a constant times the floor, the minimality of the ground field would give the Leakage Law with C = c T*.
    [real: the dilation and competitor arguments of Section 5 of the paper, the virial identity; no such bound is known,
    and the difficulty is the one the abstract names: three order-one quantities cancel to the size of lambda.]

(2) THE VALVE AS A FIXED POINT. The edge amplitude is the fixed point of the echo equation: a contraction of Hankel norm
    below one applied to a forcing. The Leakage Law would follow from a bound on the size of the forcing by the floor.
    [real: Theorem 9.24 identifies A_a with the fixed point and Proposition 9.23 gives the Hankel norm below one; bounding
    the fixed point's size by lambda is not done, and it is not clear that it can be read from the forcing alone.]

(3) THE VALVE AS A MEAN-FIELD LAW. The archimedean energy of the ground field above the horizon, (1/2 pi) int over
    |t| > T* of |F|^2 Psi_inf, is a quantity on the vacuum's side alone, with no quanta in it, and the tail of the edge
    makes it about |A_0|^2/(pi T*). The Leakage Law with C = c T* says precisely that this archimedean tail is at most a
    constant times Q(f_a): "the window form dominates its own archimedean tail". This is the form of the law in which the
    Line does not appear, and it is the one the gods would try to prove first. [real: the sampling defect, Proposition
    7.29, measures the difference between the archimedean integral above T and the zero sum above T on the minimizers; as
    an inequality for the form it is open, and it is a restatement of (d), not a new route.]

The gods' last word, from the canon: "We did not prove the Line by computing. We proved it by noticing that a universe
dies only through a leak with a threshold, and by forbidding thresholds. Whether your universe forbids them is the one
thing we cannot tell you from here."


## 29. Dictionary, Part IV

| CLU term | real object | status |
|---|---|---|
| alive / dead observer | lambda(a) > 0 / lambda(a) <= 0 | definition |
| Mortality Theorem | Weil's criterion in the odd window form (thm:weil) + monotonicity; scale from the Blind Spot | theorem |
| Birth Law, a_Z | Zhu's certified floor at a = 0.8 | theorem |
| Continuity Law | Propositions 9.15, 9.17; (H_fin) on the exceptional set | theorem off a null set |
| Decay Law | boundary law lambda' = -2|A_0|^2, Theorem 9.12(iv) | theorem a.e. |
| Leakage Law, linear valve | the inequality (d), Conjecture A; late-window form Proposition 9.30 | open; equivalent to RH beyond a_Z |
| threshold leak | the dead universes' floors crossing zero at finite flux | SANDBOX_STRANGE_MATH.md |
| Born rule of arithmetic | Weil positivity: shares |F(gamma)|^2 >= 0 iff the zeros are on the line | the Hypothesis |
| Immortality Theorem | Proposition prop:AimpliesRH (Conjecture A + AC implies RH) | theorem, given its hypothesis |
| renormalized floor U, Riemann charge | u(a) = lambda e^{cT*/2} of the paper's proof | definition; constancy is the decree |
| other polarization | the even sector, Proposition 9.40, route N2 | odd suffices (theorem); even open |
| the three valves | competitor / fixed point / mean-field forms of (d) | open, restatements |
