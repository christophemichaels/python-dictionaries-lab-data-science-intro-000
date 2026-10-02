# The if-then tree: from the present state of the Weil-window program to the Riemann Hypothesis

Written 2026-10-02, after Theorem 9.24 closed the identification point. Statement numbers refer to `weil_window.pdf`
of that date (93 pages). This file is a PLAN. It records the proven status of every node as it stands in the paper,
and it is updated after each round, never before: the order of work is mathematics and computation first, verification
second, the paper third, PROGRAM.md and this tree fourth, build, commit and PDF last. Nothing enters the paper that has
not been established; the only forward-looking texts in the repository are the NEXT lines of PROGRAM.md and this tree.

Status tags.  [THM] proved in the paper, or classical and cited.  [THM a.e.] proved at every support that is not an entry
and satisfies the Diophantine condition (D_a) of Lemma 9.6(iii): the complement of a null set, which contains every
resonant support.  [COND] proved under a stated hypothesis.  [OPEN] a precise statement, neither proved nor disproved.
[REFORM] equivalent to RH: proving it proves RH and nothing less.  [FALSE] shown false as stated.  [NUM] numerical
evidence only.  [PLAN] a step not begun.  [DEAD END] a branch closed with the reason.  [DO NOT] a branch that cannot
help a proof.

Leaves.  P = RH PROVED.  F = RH FALSE.  R = REFORMULATION (an equivalence, publishable, not a proof).  C = CONDITIONAL
RESULT (RH under hypotheses).  Every node reads: statement, status, IF ... THEN -> node, IF NOT ... THEN -> node.


## 0. The root: Weil's criterion on windows

N0. [THM: Theorem 9.3(i), Proposition 9.40(i)]
    RH  <=>  lambda(a) >= 0 for every a > 0        (odd sector: the floor of Q(f) = (1/2pi) int |F|^2 Psi - 2P^2 on [-a,a])
    RH  <=>  lambda^even(a) >= 0 for every a > 0   (even sector: Q^even = <f, K_a f> + 2<f, c>^2, c = cosh(y/2))
    Either sector alone decides RH (the odd test functions suffice because a zero and its mirror 1 - conj(rho) share
    |gamma|; the even ones by the same symmetry).
    Known: lambda(a) > 0 for 0 < a <= a_Z = 0.8 and lambda^even(0.8) > 0, certified by interval arithmetic
    (Theorem 1.3, Zhu). Beyond a_Z the floor collapses doubly exponentially (under Conjecture 7.1,
    lambda(a) >= lambda(a_Z) exp(-(c/2)(T*(a) - T*(a_Z))), T* = 2 pi e^{2a}), and the window reads an off-line zero at
    distance delta from the line only at second order delta^2 (Theorem 9.34, the blind spot), with detection scale
    a_det ~ delta^{-1} log(1/delta). No finite computation decides N0 in either direction.

    IF a certified computation finds lambda(a_0) < 0 or lambda^even(a_0) < 0 at some a_0
       THEN RH is false  -> leaf F.  Theorem 9.3(iii) then reads the distance of the farthest zero off the line from the
            rate, lambda(a) <= -e^{2 delta_max a}/(2 delta_max)(1+o(1)); Proposition 9.31 transfers the negative floor
            to a certified Suzuki defect. Status: every computed floor is positive. Not a route to a proof (blind spot).
    IF NOT  THEN a proof that lambda >= 0 for ALL a is needed. Four routes:
            N1 the paper's reduction (odd sector, Section 9.1);  N2 the even-sector criterion (Proposition 9.40);
            N3 a direct positivity structure (Section 9.4, the Atlas bridge);  N4 programs outside the window.


## 1. Route N1: the reduction (Theorem 9.1, Corollary 9.13, Proposition 9.30)

N1. [COND] Four inputs:
      (a) lambda(a_Z) > 0                                                   [THM, Theorem 1.3]
      (b) lambda locally absolutely continuous on [a_Z, inf)                [OPEN on the null set E only: (H_fin), Proposition 9.17]
      (c) lambda'(a) = -2|A_0(a)|^2 at almost every a                       [THM a.e.: Theorems 9.12 and 9.24]
      (d) 2|A_0(a)|^2 <= (c T*(a) + e(a)) lambda(a) at almost every a >= a_Z, e >= 0 locally integrable   [REFORM]
    THEN (Groenwall) lambda(a) >= lambda(a_Z) exp(-int_{a_Z}^a (cT* + e)) > 0 on [a_Z, inf), hence RH  -> leaf P1.
    Sharp form (Proposition 9.30): given (a), (b), (c),
                RH  <=>  lambda > 0 on [a_Z, inf)  <=>  -lambda'/lambda = 2|A_0|^2/lambda is locally integrable on [a_Z, inf).
    So the weakest sufficient form of (d) is "2|A_0|^2/lambda locally integrable"; Conjecture 7.1 (c T*, e = 0) is the
    quantitative version the data suggest (T_eff a bounded multiple of the horizon at seven supports to a part in a
    thousand, Computation 7.6). (a) and (c) are done; (b) is removable analysis; (d) is the Hypothesis.
    (H_inf), the boundedness of the minimizer, is assumed in the statements and is a corollary at every (D_a) support
    (paragraph after Proposition 8.6, Theorem 9.24): it is not an input any more, only a name.

### 1.1  Node (c): the boundary law.  DONE at almost every support.

N1c. [THM a.e.] lambda'(a) = -2|A_0(a)|^2 at every support that is not an entry and satisfies (D_a).
    Chain: Proposition 8.18 (the window's symbol, the Wiener-Hopf representation of the transform at every support)
    -> Theorem 9.12 (the factors, the amplitude, the edge law with its echo part, the two laws)
    -> Proposition 9.21 (the amplitude as a fixed point of the mirror equation, a contraction above T_2(a))
    -> Proposition 9.22 (uniqueness on the whole line, by energy) and Proposition 9.23 (the Hankel operator of the
       window's symbol has norm below one) -> Theorem 9.24 (the amplitude is the fixed point)
    -> Proposition 8.3 (edge law => boundary law, through the flux pairing) and Lemma 9.26 (the echo part contributes
       nothing to the flux).
    Below a_inf = 0.843 at every non-entry support without (D_a) (Corollary 8.15(v)); verified at every support
    computed (Computation 8.11, a part in a thousand).
    IF an error is found in 9.12-9.24 THEN repair it (exposed parts: the constants of Lemmas 9.9-9.11; the contraction
       bracket of 9.21(ii) with the new frequency-zero term O(X^2/s_*); the kernel conventions and the odd-part
       hypothesis of Lemma 9.20; the limit theta_a -> 0 in 9.23(i), which rests on Step 1 of 9.12 and on
       theta_inf = -pi/(4 l) + ...).
    IF NOT THEN nothing more is needed here for a.e.; the exceptional null set E does not matter for (c).
    Written into the paper as Corollary 9.25 (every support), 2026-10-02:
    (D_a) enters Theorem 9.12 only through the far feedback of Step 2 (the sums sum_u |gamma_u|/|u - 2a - nu|).
    Step 1 (the factors), Proposition 9.23 and Theorem 9.24(i)-(ii) hold at EVERY support that is not an entry: the
    whole-line amplitude equation is uniquely solvable in H = {B : B/t in L^2} with ||A_a/t||_2 <= (1 + sqrt2/(1 - ||H_Z||))
    ||N_ex/t||_2 everywhere. This is the handle for (b).

### 1.2  Node (b): absolute continuity of the floor.  The removable analytic gap.

N1b. [THM] lambda is continuous and non-increasing (Lemma 1.4). The Groenwall step needs lambda(a_2) - lambda(a_1) =
    int lambda'; a non-increasing function may carry a singular decrease (Cantor), and then the lower bound fails.
    Four ways in, from the weakest to the strongest conclusion.

N1b.1 Keep (H_BV) [HYP: Proposition 5.6]: near-minimizers of bounded variation with ||f||_inf Var f <= M locally give
    a Lipschitz floor. Every computed minimizer has Var f = 4 sup|f| exactly (unimodal on each half window).
    IF kept THEN the final statement is CONDITIONAL -> leaf C1 (RH under (H_BV) and (d)). This is Theorem 9.1 today.

N1b.2 The boundary law off a countable set, with the height T_0(a) locally bounded, gives Lipschitz (Proposition 9.15,
    through Lemma 9.14, |A_0|^2 <= pi T_0 (2P^2 + max|sigma_a|), and Saks' theorem on Dini derivates).
    [DEAD END as stated]: Theorem 9.12 gives the law off E, which is uncountable and dense, and T_0, T_2 blow up near E
    (c_D -> 0). Replacing (D_a) by a condition that fails only on a countable set is impossible for this proof: the
    Hilbert errors need sum_w (X/s_*)^{|w|}/|u_w - 2a| < inf, whose divergence set is of Liouville type, uncountable
    for every exponential rate.

N1b.3 Local boundedness of the dilation form of the minimizer at EVERY support.  [OPEN]
    Proposition 5.3 gives a D^- lambda(a) >= -D(f_a) at every non-entry support at which the autocorrelation g_f is
    differentiable at the entries, D(f) = D_inf(f) + D_P(f) + 2P^2 + 2PM. D_inf, P, M are bounded by the L^2 norm and
    the energy, locally bounded in a. D_P(f) = -sum_n 2 Lambda(n) n^{-1/2} log n g_f'(log n) is bounded by
    ||f||_inf Var f: this is (H_BV) for the minimizer alone. At (D_a) supports D(f_a) = 2a|A_0(a)|^2 (Proposition 8.3)
    and f_a is continuous with f' in L^1 (edge profile plus echo cusps with summable weights, proof of Proposition
    9.15); the summability of the cusp weights is exactly what (D_a) gives, so near E boundedness is not available.
    IF a ↦ D(f_a) is shown locally bounded over all supports THEN lambda is Lipschitz (Saks) -> (b) holds.
    IF NOT THEN N1b.4.

N1b.4 Finiteness instead of boundedness.  [DONE 2026-10-02 as analysis: Lemma 9.16, Proposition 9.17; what remains
    is the condition (H_fin) at the exceptional set E]
    Lemma 9.16 (proved in the paper, three steps, no citation needed): a continuous F whose upper left derivate is
    > -inf at every point off a countable set and >= -m almost everywhere with m integrable satisfies
    F(y) - F(x) >= -int_x^y m; a non-increasing such F is absolutely continuous. Here m = |lambda'| = 2|A_0|^2 a.e.,
    integrable because lambda is monotone. Proposition 9.17: (i) at EVERY non-entry support, with no hypothesis,
    a D_- lambda(a) >= -D-bar(f_a), the dilation form with g_f' replaced by the lower right Dini derivate D_+ g_f at
    the entries (the dilate f_b, b < 1, is admissible at ab; only the prime terms need the derivate); (ii) under
        (H_fin): D_+ g_{f_a}(log n) > -inf at every entry log n < 2a, at every support off a countable set,
    lambda is absolutely continuous, and Theorem 9.1, Corollary 9.13 and Proposition 9.30 hold with (H_fin) in place
    of (H_BV); (iii) (H_fin) holds at every (D_a) support (f_a continuous, f_a' in L^1, so g_f Lipschitz), hence it
    is a condition on the null set E alone, implied there by f_a in L^inf of bounded variation.
    (H_BV) is therefore replaced by a statement about the minimizer alone, at a null set of supports, about finitely
    many one-sided derivates. What is NOT done: (H_fin) on E itself. At a in E what is known: f_a lies in the form domain
    int |F|^2 sigma~_a < inf, a logarithmic Sobolev space below H^{1/2}, which gives nothing pointwise; the whole-line
    amplitude equation is uniquely solvable in H at every support (N1c, observation), with the explicit bound; and
        g_f'(log n) = -(1/2pi) int |F(t)|^2 t sin(t log n) dt,
    with |F|^2 = |A_a|^2/(t^2 sigma~_a) + ..., converges iff the component of |A_a(t)|^2 at the frequency log n has a
    coefficient integrable against dt/(t sigma~_a). At (D_a) supports that component comes from mirror chains of
    weight O(s^{-2}) and converges; the frequency log n is a lattice point, so the component exists at every support.
    IF the frequency-(log n) component of |A_a|^2 is controlled at every support (through the H-bound, or through a
       formula for it that does not involve 1/|u_w - 2a - nu|) THEN (H_fin) holds on E, (b) holds, and
       Theorem 9.1 reads "the inequality a.e. implies RH" with no analytic hypothesis -> the critical path is (d) alone.
       The energy alone does not give it: a component decaying like 1/(t log^3 t) is compatible with the energy and
       makes D_+ g_f = -inf; what excludes it at (D_a) supports is the boundedness of A_a (the component then decays
       like 1/(t^2 log t) and the cusp of g_f is |h|/log(1/|h|), with derivate zero). So (H_fin) on E is a statement
       about the boundedness of the amplitude, or the decay of the log n component of |A_a|^2, at the exceptional
       supports (paragraph after Proposition 9.17).
    IF NOT THEN N1b.3 is harder still; fall back to N1b.1 (leaf C1) and pursue N2, whose regularity problem is its own.
    Three cheap tests, decidable before any writing:
      (i)  is c_Z(a) = ||T_{e^{2i theta_a}}^{-1}||^{-1} continuous in a between consecutive entries? (there theta_a
           changes only through lambda(a) and the real zeros t_j; compute theta_a from the FEM symbol);
      (ii) does the frequency-(log n) coefficient of |A_a|^2 have a closed expression in the lattice data alone?
           (at a finite-tree support everything is explicit: start at a_2 < a < a_3);
      (iii) does Proposition 8.3's limit lim B_T exist on a model lattice with a Liouville 2a (the Diophantine sum made
           to diverge by hand), or does B_T oscillate without a limit? If it oscillates, the Dini derivatives of lambda
           at such a support are the limsup/liminf of B_T, and finiteness is still enough for N1b.4.

N1b.5 Integrated inequality instead of absolute continuity.  [OPEN, lower priority]
    The Groenwall step uses only lambda(a_2) >= lambda(a_1) - int_{a_1}^{a_2} (cT* + e) lambda. IF (d) is proved in this
    integrated form directly (for instance from an integrated dilation identity) THEN (b) is not needed. No candidate
    argument is known: the integrated dilation identity needs the same regularity.

### 1.3  Node (d): the inequality.  The Hypothesis itself.

N1d.0 Status and shape.  [REFORM]
    2|A_0|^2 >= 0, so (d) fails wherever lambda < 0: (d) implies RH beyond a_Z, and (d) is implied by RH plus the onset
    (N1d.1). Any proof must use the minimality of f in a way that excludes a negative floor; no argument treating
    lambda as a given number can succeed (the paper's second remark on the inequality). Equivalent forms:
      - Conjecture 7.1 (bounded relay): T_eff <= c T*, the zeros carrying the decay of the floor sit within a bounded
        multiple of the horizon; with e(a): a power of the horizon above one is allowed (Proposition 9.30).
      - Pohozaev form: the dilation identity a lambda' = -[D_inf + D_P + 2P^2 + 2PM] (Theorem 5.2 for the K-mode form,
        Proposition 5.3 for the exact form) is the Pohozaev identity of the Weil form; (d) asks for a sign of its
        remainder against T* lambda; the obstruction is D_P, whose sign is not fixed. For the Dirichlet Laplacian the
        analogue u'(a)^2 <= C lambda is unconditional because the remainder is a sum of squares.
      - Prime-side form: at a finite-tree support, with phi_1 the solution of the eigen-equation with unit polar source
        and lambda as a parameter, f = P phi_1, P^2 = 1/||phi_1||^2, lambda(a) the lowest root of the secular equation
        <phi_1, s> = 1, and A_0 = P A_0^{(1)}(lambda, a): (d) reads 2|A_0^{(1)}|^2 <= c T* lambda ||phi_1||^2 at the
        root. Only the archimedean symbol and the entries below e^{2a} appear, through sigma_eff and the spectral
        measure mu_b of the echo lattice at the edge (Proposition 8.23, Proposition 9.33).
    Two things that do not work: a bootstrap in height (the verified zeros give positivity to a ~ 9; the gain per
    step is O(1/T*), which does not sum); numerical extension (N1d.3).

N1d.1 The reachable target: RH => (d).  [OPEN, provable with effort; second priority]
    Under RH, 2|A_0|^2 = lim pi T sum_{|gamma|>T} |F(gamma)|^2 (Theorem 9.12(iv), the analytic tail law, no zeros
    involved in its proof) and the zero sum above kappa T* is at most lambda minus the sum below, every term a square
    (Propositions 7.27, 7.29, 7.32): (d) with e = 0 says that the minimizer's tail energy above kappa horizons is at
    most c lambda/(pi kappa), the ONSET. Known: the edge's measure mu_b keeps at least 7/12 of its mass below any fixed
    number of horizons uniformly in a (kappa_2 = 5/12, Proposition 9.33(v)), the moments cap the provable mass bound
    near one third (kappa_4 = 151/360, kappa_6 = 4033/6720 exact, Proposition 9.39), the true mass above five horizons
    is 1e-6 to 1e-4 at a <= 1.25 (Computation 9.36); the uniform onset is FALSE as stated (the mass above one horizon
    tends to about 6%, the 1% height recedes like e^{ca}, Computation 9.37), the natural rate being a power of the
    horizon above one, which Proposition 9.30 absorbs into e(a). The free lattice law is Gaussian (Proposition 9.38);
    the confined law's closed form is open.
    IF the shape of mu_b (not its moments) gives mu_b([s, inf)) <= C s^{-1-epsilon} uniformly in a, by a resolvent bound
       or the Perron vector at the top of the component's band THEN RH => (d), and the paper states
       RH <=> Conjecture 7.1 with integrable excess  -> leaf R1 (an exact reformulation; publishable; not a proof).
    IF NOT THEN the equivalence stays one-directional, and N1d.2 proceeds without the exact shape of the target.

N1d.2 Unconditional attack on (d).  [OPEN; as hard as RH]
    The missing piece (memo RH_TOP3_PROOF_ARCHITECTURE.md, Section 3.2): a structure theorem for the minimizer stated
    in terms of the primes alone, sensitive at relative resolution e^{-cT*} (positivity fails under a perturbation of
    size 1e-14 of the prime positions at a = 0.8). Three vectors, each with its likely fatal flaw (memo Section 4):
    V1 Fourier interpolation (Bondarenko-Radchenko-Seip: bases with the zeros as nodes on one side and +-log n on the
       other). Flaw: run from the log n side, the dual nodes are the zeros only if RH holds.
       IF the duality can be run from the prime side with Diophantine input only THEN a prime-side description of
          f_a, which is the structure theorem -> N1d.2 closes.
       IF NOT THEN a reformulation ("RH <=> the dual basis is real-noded") -> leaf R. Survives either way: rigorous
          interpolation of the zeros below (1 - epsilon)T* assuming RH below T*.
    V2 Deformation invariance (move the primes continuously; look for an invariant that cannot change).
       Flaw: there is no open neighbourhood of the primes in which positivity holds, so the invariant cannot be
       topological in the prime positions; it must be algebraic, exactly preserved by "the primes are the primes of Z"
       (Fourier quasicrystal rigidity, Lev-Olevskii, Kurasov-Sarnak).
       IF such an invariant exists and forces the sign of D_P against T* lambda THEN (d).
       IF NOT THEN [DEAD END for topological arguments]; survives: the experiment a_fail(epsilon).
    V3 Traces to the spectral edge (Alpoege-Furman): RH <=> Hankel positivity of the moment sequence m_k of the
       compressed Weil matrix, m_k a k-point correlation of Lambda at bandwidth <= 1; m_2 is known (Montgomery), m_3
       needs triple correlations, and finitely many moments never pin the edge.
       IF the moment problem can be closed by a structural constraint (the confined-walk law gives the moments of
          mu_b exactly, Proposition 9.39; its limit law in closed form would be a candidate for "the full sequence")
          THEN (d) in prime-side form.
       IF NOT THEN a reformulation -> leaf R.
    The clean target, to be written out first [PLAN, third priority]: at a finite-tree support (a_2 < a < a_3, one
    prime, then a < a_inf) write 2|A_0^{(1)}(lambda, a)|^2 - c T* lambda ||phi_1(lambda)||^2 as an explicit function of
    lambda through the Wiener-Hopf data, locate the lowest root of the secular equation, and test numerically whether
    |A_0^{(1)}|^2/||phi_1||^2 is monotone in lambda with an explicit crossing.
    IF it is THEN (d) reduces to locating the root: an arithmetic statement about the entries below e^{2a}.
    IF it is not THEN the inequality depends on the root and the function jointly, and the structural attack must
       produce both at once.

N1d.3 Numerical extension: more supports, finer grids, larger K, refitted constants.  [DO NOT]
    The blind spot (Theorem 9.34) and the doubly exponential collapse make it uninformative for a proof. Computations
    stay as checks of theorems (every theorem of Section 9 has one).

### 1.4  Reading N1 as a whole
    P1 needs (b) and (d). (b) is analysis and is decidable; N1b.4 is the step to take. (d) is the Hypothesis; N1d.1
    turns the reduction into an equivalence and tells exactly what must be proved; N1d.2 is the attack.
    IF (b) is done and (d) is not THEN the result is "the inequality a.e. => RH", unconditional -> a sharper leaf R.
    IF (d) is done and (b) is not THEN leaf C1 with (H_BV) alone, and N1b.5 becomes worth a try.


## 2. Route N2: the even sector (Proposition 9.40, Computation 9.41)

N2. [THM for the equivalence] RH <=> Q^even >= 0 on every window. Folded to L^2(0,a), K_a^even has a kernel negative
    off the diagonal, so its semigroup is positivity improving and the lowest eigenvalue is simple with a positive
    eigenfunction of positive overlap with c = cosh(y/2) (Proposition 9.40(ii)); with exactly one negative eigenvalue,
    Q^even >= 0 <=> 1 + 2<c, K_a^{-1} c> <= 0 (interlacing and a determinant, Proposition 9.40(iii)). The "exactly one
    negative eigenvalue" is [NUM] at the computed supports (Computation 9.41), where the secular number is just below
    zero and collapses doubly exponentially like the odd floor. No sum of squares appears.
    IF "exactly one negative eigenvalue" is proved for all a THEN RH <=> a single scalar inequality in a -> leaf R2.
    IF a sign structure for the secular number is found (monotone in a with a nonpositive limit; a sum of squares; a
       Pohozaev identity for the even form with a signed remainder) THEN RH -> leaf P2.
    IF NOT THEN build the even analogue of Section 9.2 (edge law, boundary law; the polar term is +2<f,c>^2 instead
       of -2P^2, so the edge mechanism differs) and look for a cross-sector identity, for instance that the two floors
       collapse at the same rate -> a structural statement with no odd-sector counterpart [PLAN, fourth priority].


## 3. Route N3: direct positivity structures

N3.1 The Atlas bridge (Propositions 9.28-9.32): the window form inside the Atlas's Q_a^-, two readings of delta_max
    (Theorem 9.29), the witness in the Suzuki metric (Proposition 9.31), the paired-prime matrices and the echo tree
    (Proposition 9.32). These transfer a NEGATIVE floor into other invariants: detection, the F branch.
    IF a positivity structure of the Atlas (the ground-state transform, a Krein extension) applies to the window form
       THEN N3 is a proof route -> leaf P3.
    IF NOT THEN N3 serves leaf F only.
N3.2 Rellich-Pohozaev with a signed remainder: a sum-of-squares decomposition of D_P plus lower-order terms against
    T* lambda.  [OPEN; none found in either sector.]


## 4. Route N4: outside the window

N4. The Moebius Green energy and the Nyman-Beurling-Baez-Duarte criterion (mobius_modifier_v2), the Atlas (Theta
    kernels), the hundred routes of RH_ROUTES.md, each with its named obstruction. Independent of N1; a proof there
    would make N1 a corollary. Their status is in their own documents; none is a proof, and none is on this path.


## 5. The critical path, in order

  1. N1b.4  [DONE as analysis, 2026-10-02] Lemma 9.16 and Proposition 9.17 reduce (H_BV) to (H_fin) on the
            exceptional null set E. Remaining: (H_fin) on E, through the three cheap tests, then the boundedness of
            the amplitude at exceptional supports.
  2. N1d.1  RH => (d) through the shape of mu_b: makes the reduction an equivalence and fixes the exact target.
  3. N1d.2  the prime-side form at a finite-tree support: the structural attack with the target in explicit
            Wiener-Hopf terms.
  4. N2     the even sector as the second reading, in parallel and at low cost.
  Not on the path: more grids, more supports, refitting the decay constant; the bootstrap in height; topological
  deformation arguments; any pointwise claim in the band without an interior-resolving solver.
  What would stop the program: a certified negative floor (leaf F); a proof that the constant c of Conjecture 7.1 must
  grow faster than any power of the horizon while RH holds (then only the integrable form of (d), Proposition 9.30,
  remains, and that form is the one to prove).


## 6. The leaves, collected

  P1  RH proved through N1: needs (b) and (d).            P2  RH proved through N2 (the even secular number).
  P3  RH proved through a positivity structure (N3).      F   RH false: a certified negative floor (not a proof route).
  R1  RH <=> Conjecture 7.1 with integrable excess (from N1d.1).     R2  RH <=> one scalar inequality in a (from N2).
  C1  RH under (H_fin) on E (formerly (H_BV)) and (d): Theorem 9.1 as it stands.
  Dead ends: N1b.2 as stated; topological invariants in V2; finitely many moments in V3; the height bootstrap.


## 7. Dependency table

| statement | assumes | gives | status |
|---|---|---|---|
| Theorem 1.3 (Zhu) | interval arithmetic | lambda(0.8) > 0, lambda^even(0.8) > 0 | THM |
| Lemma 1.4 | none | lambda continuous, non-increasing, differentiable a.e. | THM |
| Theorem 5.2, Proposition 5.3 | g_f differentiable at the entries | a lambda' = -[D_inf + D_P + 2P^2 + 2PM] (K-mode exact; exact form two-sided Dini) | THM |
| Proposition 5.6 | (H_BV) | lambda locally Lipschitz | COND |
| Lemma 9.16 | none | absolute continuity from an integrable derivate finite off a countable set | THM |
| Proposition 9.17 | (H_fin) on the exceptional set E | lambda absolutely continuous; Theorem 9.1 with (H_fin) for (H_BV) | THM / COND on (H_fin) |
| Corollary 9.25 | none | Proposition 9.23 and Theorem 9.24(i)-(ii) at every non-entry support | THM |
| Proposition 8.3 | edge law with two derivatives | D(f) = 2aC^2, lambda' = -2C^2 | THM |
| Proposition 8.6 | none | prime-free minimizer bounded: (H_inf) prime-free | THM |
| Proposition 8.18 | f in L^2 | Wiener-Hopf representation with the window's symbol, at every support | THM |
| Theorem 9.1, Corollary 9.13 | (H_inf), (H_BV) or (H_fin) for AC, (d) | Conjecture 7.1, RH | COND |
| Theorem 9.3 | (i) none; (ii)-(iii) (H_inf) and TV at most exponential | RH <=> lambda >= 0; the rate of negativity = delta_max | THM / COND |
| Theorem 9.12 | (H_inf) [a corollary], (D_a) | factors, amplitude, edge law with echo part, boundary law, tail law under RH | THM a.e. |
| Lemma 9.14, Proposition 9.15 | law at every support of an interval, T_0 bounded | lambda Lipschitz there | COND (hypothesis not available) |
| Lemma 9.20 | odd part of the kernel integrable at 0 | Hilbert transform of a symbol, one log per derivative | THM |
| Proposition 9.21 | setting of 9.12, (D_a) | fixed point on Y_+, contraction 1/2, atom functional c_1 | THM a.e. |
| Proposition 9.22 | A in the class A | uniqueness of the whole-line equation | THM (every support) |
| Proposition 9.23 | theta_a -> 0 (Step 1 of 9.12, no (D_a)) | ||H_Z|| < 1, ||H_{chi Z}|| < 1 | THM (every support) |
| Theorem 9.24 | 9.21(i)-(ii), 9.23 | I - T invertible on H and on Y_+; A_a = fixed point | THM a.e. ((i)-(ii) every support) |
| Proposition 9.30 | (a), (b), (c) | RH <=> -lambda'/lambda locally integrable | COND on (b) |
| Proposition 9.33 | none | band and mass of mu_b; m_2 <= (5/12)(2a)^2 (1 + O(1/a)) | THM |
| Theorem 9.34 | none | the blind spot: second order in the distance from the line | THM |
| Propositions 9.38, 9.39 | none | free law Gaussian; kappa_2, kappa_4, kappa_6 exact | THM |
| Proposition 9.40 | none; (iii) one negative eigenvalue | even-sector criterion | THM / NUM for the eigenvalue count |
| Conjecture 7.1 | | T_eff <= c T* | REFORM (with Proposition 9.30: exactly RH) |
| Conjecture 7.7 | RH | tail law | THM under RH at (D_a) supports (Theorem 9.12(iv)) |
