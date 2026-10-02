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
