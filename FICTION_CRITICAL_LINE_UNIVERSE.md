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
