# Dictionary of the program

Fixed 2026-10-03. Each term has one mathematical meaning throughout the proof. Three status labels are used everywhere:
DEFINITION (fixes a meaning), PROVED (a relationship between meanings, established, with its source), HYPOTHESIS (a
relationship still to be established). A dictionary entry fixes meaning; a theorem establishes a relationship.

## 1. The window, the form, the floor

| Term | Fixed mathematical meaning |
|---|---|
| Support parameter a | Half-width of the interval [-a, a] supporting the test functions. |
| Admissible space | Real odd functions in the form domain of the window form, supported in [-a, a], extended by zero to the real line. |
| Window form Q_a | The number-theoretic quadratic form on the admissible space for that interval: Q_a(f) = <f, H_a f> - |<v_a, f>|^2 (Section 2). |
| Floor lambda(a) | The variational infimum inf_{||f||_2 = 1} Q_a(f) over the admissible space. |
| Minimizer f_a | A normalized admissible function attaining lambda(a), wherever attainment is established. |
| Edge amplitude A_0(a) | The boundary asymptotic coefficient of f_a, with the normalization of the representation theorem (paper, Theorem 9.12(ii), Corollary cor:edgeform). Not an ordinary endpoint value. |
| Absolute leak L(a) | The nonnegative quantity 2|A_0(a)|^2. |
| Derivative identity | The assertion lambda'(a) = -L(a), at supports where it has been proved (paper, Corollary cor:tree(v) below a_inf = 0.843 off the entries; Theorem thm:everysupport almost everywhere). |
| Relative leak r(a) | L(a)/lambda(a), defined where lambda(a) > 0. Where the derivative identity holds, r = -(log lambda)'. |
| Prime-side comparison | A proved estimate L(a) <= K(a) lambda(a) with a nonnegative coefficient K integrable on every finite support interval, including across any proposed first null. |
| Integration regularity | Local absolute continuity of lambda, sufficient to recover its increments by integrating lambda'. |
| First null a_* | The first finite support after a positive starting point at which lambda(a_*) = 0, if such a support exists. |
| Exceptional supports | Parameters excluded by the hypotheses of a particular representation or regularity theorem. |
| Horizon T*(a) | The reference height 2 pi e^{2a}. |

## 2. The prime-side operator and the polar source

| Term | Fixed mathematical meaning |
|---|---|
| Prime-side operator H_a | The self-adjoint operator of the archimedean contribution and every signed prime-power term on the admissible space, so that <f, H_a f> = (1/2 pi) int |F_a|^2 [Re psi_Gamma(1/4 + it/2) - log pi] dt - 2 sum_{n >= 2} Lambda(n) n^{-1/2} g_a(log n). It is Q_a without the polar term. (Earlier notes wrote K_a; H_a is the fixed name.) |
| Polar source v_a | v_a(x) = sqrt(2) sinh(x/2) on [-a, a], so that the polar term is -|<v_a, f>|^2. |
| h_a | The bottom of the spectrum of H_a, inf sigma(H_a). Strict positivity H_a >= h_a I > 0 is a HYPOTHESIS where it is used. |
| Secular function m(lambda, a) | <v_a, (H_a - lambda)^{-1} v_a>, for lambda < h_a; m_a = m(0, a) = <v_a, H_a^{-1} v_a>. |
| Margin delta_a | 1 - m_a. |
| Unit-source solution phi | phi_a(lambda) = (H_a - lambda)^{-1} v_a; phi_a = phi_a(lambda(a)) at the root, phi_a(0) = H_a^{-1} v_a at zero. |
| Boundary coefficient alpha(a) | The boundary asymptotic coefficient of the unnormalized phi_a(lambda(a)), in the same normalization as A_0, wherever the boundary representation is established; alpha_0(a) the same for phi_a(0). |
| Polar branch | A minimizer with <v_a, f_a> = 0; its equation is H_a f_a = lambda(a) f_a and the unit-source representation does not apply. |

## 3. Fourier side

| Term | Fixed mathematical meaning |
|---|---|
| Fourier transform F_a(t) | F_a(t) = int_R f_a(x) e^{-itx} dx, with zero extension outside the window. |
| Fourier tail | The portion with |t| > T of a specified Fourier integral; weight and normalization written explicitly each time. |
| Tail onset T_0(a) | A threshold beyond which the stated tail estimate holds with its stated error bound. Local uniformity in a requires proof. |
| Autocorrelation g_a(d) | int_R f_a(x + d) conj(f_a(x)) dx. |
| Prime adjacency C_a | The weighted operator on the confined echo graph, with allowed steps +- log n and weights Lambda(n)/sqrt(n) (paper, Proposition 8.21). |
| Spectral band | The spectrum of C_a in its spectral parameter. |
| Resolvent cutoff | A threshold beyond which (sI - C_a)^{-1} exists and satisfies the proved quantitative bounds (addendum, Lemma 1: s >= max{1, 2B(a), 2 sqrt(beta(a))}). |

## 4. The closing deduction (PROVED, given its three inputs)

If lambda(a_0) > 0, integration regularity holds on [a_0, b], the derivative identity holds almost everywhere there, and
the prime-side comparison L <= K lambda holds with K integrable on [a_0, b], then
        lambda(b) >= lambda(a_0) exp( - int_{a_0}^{b} K(u) du ) > 0
for every finite b >= a_0 (paper, Proposition prop:AimpliesRH; addendum, Theorem 3). With a_0 = a_Z = 0.8 (Zhu) the
Riemann Hypothesis follows (paper, Theorem thm:weil).

The remaining gates are the prime-side comparison, the Fourier-tail onset, and regularity. RESEARCH_SECULAR_FORM.md
carries the prime-side comparison in this notation.
