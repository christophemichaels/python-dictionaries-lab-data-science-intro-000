# Strange math: a sandbox

A side project, 2026-10-02. Nothing in this file is a claim about the zeta function, and nothing here enters the paper. The
rule of the sandbox: we may invent any law we like, but every invented law is then held against the real object, and the
file records which laws survive the contact and which do not. The surviving laws turn out to be the ones that depend only
on the sequence of the primes (their ranks and their density). The one that dies is the critical line.

Vocabulary from the paper (weil_window.pdf, 2026-10-02): the window [-a, a]; the horizon T* = 2 pi e^{2a}; the floor
lambda(a) and its logarithm Phi = -log lambda; the entries log n < 2a, the prime powers inside the window; the edge
amplitude A_0; the edge's spectral measure mu_b of mass |b|^2 on the height variable s (s = 2a is one horizon, s = 2a + log k
is k horizons); the constants kappa_2 = 5/12, kappa_4 = 151/360, kappa_6 = 4033/6720 of the confined walk.


## 1. The inversion: two is one, three is two, five is three

Write p_n for the n-th prime and r(p_n) = n for its rank. The inversion replaces the value by the rank. Three things are
true and provable about the window picture, and depend on the primes only through their ranks and density.

LAW 1 (the rank clock). A prime enters the window at the window-time a_n = (1/2) log p_n, and these are the only events:
    the floor is smooth between consecutive entries and has a kink at each (the entering term, Section 5). By the prime
    number theorem a_n = (1/2) log(n log n) (1 + o(1)): the clock runs logarithmically in the rank, and the inter-arrival
    time of the primes in window-time is a_{n+1} - a_n = (p_{n+1} - p_n)/(2 p_n)(1 + o(1)), on average 1/(2n).

LAW 2 (the energy law). Call E(a) = sum_{p^k < e^{2a}} (log p)^2 / p^k the prime energy of the window: the sum of the squares
    of the weights with which the primes enter the form. Mertens' theorem gives E(a) = 2a^2 + O(a). This is |b|^2 in the
    paper, the total mass of the edge's spectral measure (Proposition 9.33(iii)). The energy density in window-time is
    dE/da = 4a: linear. It depends on the density of the primes and on nothing else.

LAW 3 (conservation and transmutation). The edge couples to the lattice of echoes through the vector b, and the lattice
    turns b into the positive measure mu_b on the heights, of total mass exactly E(a) (Proposition 8.23). The mass is
    conserved; what the lattice does is move it among the heights. Its mean tends to zero (Proposition 9.33(iii)); its
    variance is (5/12)(2a)^2 (1 + o(1)) (Proposition 9.33(v)); its fourth and sixth cumulant-type constants are exact
    rationals (Proposition 9.39). The proofs use Mertens' theorem and the combinatorics of walks on the lattice of
    logarithms, nothing about which numbers are prime: the same constants hold for any Beurling system of generalized
    primes with the same counting function. The first alien constants, 5/12, 151/360, 4033/6720, are sequence-only.

LAW 4 (the flux law). The boundary law lambda' = -2|A_0|^2 is a flux identity (Proposition 8.3): the decrease of the floor
    is the flux of the edge's energy through the horizon. In the tail law the same quantity is pi T_eff times the energy of
    the minimizer above the height T_eff, which is a bounded multiple of the horizon at every support computed. Read as a
    conservation law: the floor is the part of the minimizer's energy that has not yet left through the horizon, and
    nothing is destroyed, only radiated outward and recorded in the zeros. The law itself holds at almost every support
    (Theorem 9.12); its proof uses the arithmetic of the primes only through the Diophantine condition on 2a.

LAW 5 (the repulsion). An off-line zero at distance delta from the critical line changes the window's floor at order
    delta^2 (Theorem 9.31). "As the decay is reduced by the incoming prime, so the knowledge is repelled": the lattice hides
    its defects at second order, and the detection scale a_det ~ delta^{-1} log(1/delta) recedes as the defect shrinks.

So four of the five laws the picture suggests are theorems, and three of them hold for any sequence with the density of the
primes. What depends on the values is what the sandbox now tests.


## 2. The test: feed the window the ranks

The K-mode odd Weil form is computed (rh_sandbox_rank_primes.py, log data/sandbox_rank_primes.log) with four prime systems
in the place of the primes, each entering the form as -2 c g_f(d) at positions d with weights c:

    T   the true primes and prime powers: d = k log p, c = log p / p^{k/2};
    R1  "every integer is a prime": d = log n, c = log n / sqrt(n) for every n >= 2 (the rank system: 2 is 1, 3 is 2, ...);
    R2  true positions, rank energies: d = k log p_n, c = log(n+1)/(n+1)^{k/2};
    R3  rank positions, true energies: d = k log(n+1), c = log p_n / p_n^{k/2}.

Below a = (log 4)/2 = 0.693 the systems T and R1 coincide (the integers below 4 are the primes 2, 3), and for a < 0.805 they
differ only in the weight at 4: log 4 as a prime against log 2 as the square of 2. The results:

    K = 30 modes, 40 digits (rh_sandbox_rank_primes.py 30 40; the K-mode floor is an upper bound for the true floor; Zhu's
    certified value at a = 0.8 is 8.2e-15 <= lambda <= 2.3e-14, and the row below sits inside it).

      a     2a^2  |  T  terms  E(a)    lambda          |  R1 terms  E(a)   lambda      |  R2 terms  E(a)   lambda      |  R3 terms  E(a)   lambda
     0.50   0.50  |  1  0.240  +1.94e-4   |  1  0.240  +1.94e-4    |  1  0.240  +1.94e-4    |  1  0.240  +1.94e-4
     0.60   0.72  |  2  0.643  +5.99e-7   |  2  0.643  +5.99e-7    |  2  0.643  +5.99e-7    |  2  0.643  +5.99e-7
     0.70   0.98  |  3  0.763  +2.58e-10  |  3  1.123  +3.79e-10   |  3  0.763  +2.58e-10   |  4  1.281  +4.20e-10
     0.80   1.28  |  3  0.763  +1.63e-14  |  3  1.123  -0.129  (2) |  3  0.763  +1.63e-14   |  4  1.281  -0.378  (2)
     0.90   1.62  |  4  1.281  +8.52e-20  |  5  2.176  -0.216  (2) |  4  1.243  -0.0104 (3) |  6  2.344  -0.501  (3)
     1.00   2.00  |  5  1.822  +5.12e-25  |  6  2.717  -0.509  (4) |  5  1.761  -0.0231 (4) |  7  2.850  -0.738  (4)

    (the number in parentheses is the count of negative eigenvalues of the K-mode form; "terms" is the number of positions
    d < 2a; E(a) = sum c^2 is the prime energy; log data/sandbox_rank_primes.log).

What the table says. The true primes give a floor that is positive at every support and collapses as the paper says, from 2e-4 at a = 0.5
to 5e-25 at a = 1.0 (and the K = 30 value at 0.8 sits inside Zhu's certified interval). Each rank system agrees with the
truth exactly as long as it coincides with it (R1 and R3 to a = 0.6, R2 to a = 0.8, because the ranks 1, 2 of the primes 2, 3
are the primes themselves shifted by one), and the first time it differs the floor collapses to a negative number OF ORDER
ONE, not to a slightly negative number: R1 at a = 0.8, where the only change is that 4 counts as a prime of weight log 4
instead of the square of 2 with weight log 2, has floor -0.13 with two negative eigenvalues, against +1.6e-14 for the
truth; R3, which puts the prime 5 at the position of 4, has -0.38. R2 is the sharpest test: it keeps every position and
changes only the energies of 5 and 7 by about four per cent (log 4/2 = 0.693 for log 5/sqrt 5 = 0.720; log 5/sqrt 5 for
log 7/sqrt 7 = 0.735), and that is enough to take the floor from +8.5e-20 to -0.010 at a = 0.9. One more detail: just
after a modification, before the sign change, the modified floor is LARGER than the true one (a = 0.7: 3.8e-10 and
4.2e-10 against 2.6e-10): extra energy at an entry lifts the floor for a moment and then breaks it.

So the sequence of the primes, their ranks, their density and their energy law, is not what holds the floor above zero.
The positions and the energies of the individual primes are, at the precision of their actual values; four per cent in
one weight is already fatal. This is the razor edge of the paper (positivity fails under a shift of one part in 10^14 at
a = 0.8) seen from the other side: the window does not read "a prime-like sequence", it reads the primes.


## 3. The alien axioms, and what the real object says to each

The sandbox proper. Each axiom is stated as if it were a law of an unknown physics; after it, in brackets, what is actually
true of the real object.

AXIOM A (the prime is a quantum of window-time). A prime is not a number but an event: the moment a_n at which a new
    reflector appears inside the window. Its value is the time of its arrival; its rank is its index in the sequence of
    events. The clock is logarithmic in the index.
    [True as a description of the floor's kinks (Law 1). False as a replacement: the test of Section 2 shows that moving
    the events to rank positions destroys positivity.]

AXIOM B (energy is conserved and transmuted). Each event deposits energy (log p)^2/p; the window's energy is 2a^2; the
    lattice of reflections redistributes it over heights without loss; the floor is the trapped part and its decrease is
    the flux through the horizon.
    [Laws 2, 3, 4: theorems. The one quantitative thing the sandbox cannot supply is the rate: that the flux is at most
    a bounded multiple of the horizon times the trapped energy, which is Conjecture 7.1, and which fails exactly where
    the floor is negative.]

AXIOM C (the diffusion constant of the primes is 5/12). The energy spreads over heights like a walk confined to the window,
    with variance (5/12)(2a)^2, kurtosis 302/125, sixth standardized moment 36297/4375; the walk's law is not Gaussian and
    not semicircular, and its closed form is unknown.
    [Theorem (Proposition 9.39), and sequence-only: any Beurling system with the primes' density has the same constants.
    The alien constants are real; they are the constants of a combinatorial object, the confined walk, not of the primes.]

AXIOM D (the inversion). The values of the primes are not important; the sequence is. Replace p_n by n + 1.
    [Refuted by the test. The floor at a = 0.8 is 8 x 10^{-15} for the true primes (Zhu, certified), and the paper records
    (memo RH_TOP3_PROOF_ARCHITECTURE.md, Section 2.4; and the deleted form of Section 5 of the paper, rh_deleted_form.py: removing the term of the prime 2 makes the floor negative 0.025 beyond its entry) that moving one prime by one part in 10^{14} already makes it negative: the
    critical line is a statement about the values at that precision, and the rank system is nowhere near it. The sequence
    carries the energy law and the diffusion constants; it does not carry positivity.]

AXIOM E (the repulsion of knowledge). The closer a defect is to the line, the less the window sees it, quadratically.
    [Theorem 9.31. The sandbox adds nothing; the real object is already this strange.]

AXIOM F (the strangeness is confined to the line by one sign). All of the above is universal, true for the primes and for
    any sequence like them. The single statement that is not universal is the sign of the prime part D_P of the dilation
    identity against T* lambda: a lattice of reflectors with the primes' density but the wrong positions has the same
    energy, the same diffusion constants, the same flux law, and a negative floor. Whatever confines the zeros to the
    line is in the positions of the reflectors at a precision of fourteen digits and beyond, and in nothing coarser.
    [This is the paper's second remark on the inequality, restated. It is the content of the sandbox.]


## 4. What the sandbox gives back to the program

Three things, none of them a proof.

(i) A sharpened picture of what is universal and what is arithmetic. The confined-walk constants, the energy law, the
    flux law and the blind spot are sequence-level facts; they will hold for every Beurling system with the primes'
    density, and one could prove them in that generality to make the point formal. The critical line is the one non-
    universal fact, and any proof of it must use the positions of the primes at a precision no density statement reaches:
    the memo's "structure theorem for the minimizer stated in terms of the primes alone" is a theorem about the values.

(ii) A test that any proposed "law of the primes" must pass before it is taken seriously here: does it change when the
    primes are replaced by the rank system R1, or by the smooth density model? If it does not change, it cannot by itself
    prove positivity. The energy law and the alien constants fail this test honourably; the inequality passes it, because
    it fails for the rank system.

(iii) A side quest worth one scoped computation later: how far below the horizon the zeros of the minimizer's transform
    track the true zeros of zeta, as a function of the support, since the minimizer at support a is a predictor of the
    zeros below the horizon from the primes inside the window. The blind spot says this cannot decide positivity; it can
    still measure how much of the explicit formula a finite window resolves.

Not to be done: to treat any of the alien axioms as evidence for RH, or to let this file pull on the critical path of
RH_IF_THEN_TREE.md, where the next step is the prime-side form of the inequality at a support with one prime inside.
