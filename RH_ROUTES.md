# 100 Routes Toward the Riemann Hypothesis, Ranked by Credibility

*Prepared 2026-09-29 as a companion to C. Michaels, "Primes, Folds, and the One Dot" (draft of 2026-09-28).*

## Read this first

- **Nobody knows a route that proves RH.** Every entry below is a research program with a named obstruction. None of them is a proof strategy anyone can vouch for.
- **"Credible" here means three things:** serious mathematicians work on it; the obstruction is understood; and there is an intermediate result that would be publishable on its own.
- **Tiers:**
  - **[A]** A live program with theorems from 2019–2026. A partial result would be publishable now.
  - **[B]** Serious and well studied, but stalled at a known obstruction.
  - **[C]** An equivalent reformulation or a speculative framework. Useful for partial cases and intuition; low odds as a complete route.
- **Two results don't count as progress toward "all zeros":** a proportion (even 100%) and a finite verification (any height). Both are valuable. Neither implies RH.

## Seven filters every credible route must pass

1. **Euler-product filter.** The Davenport–Heilbronn function and Epstein zeta functions of class number > 1 have a functional equation of the same shape as ζ, and they have zeros off the line. An argument that uses only the functional equation, symmetry, growth, or the *positions* of prime powers proves something false. The argument must use multiplicativity somewhere essential.
2. **Beurling filter.** Diamond–Montgomery–Vorhauer (2006) built generalized prime systems whose integers are extremely well distributed, but whose zeta functions have zeros close to Re s = 1. Counting axioms alone can't give RH.
3. **Function-field filter.** Does the route specialize to a proof (or at least a sensible statement) for curves over F_q, where RH is Weil's theorem? If it can't even be stated there, ask why it would work over ℤ.
4. **Uniformity filter.** Finite checks don't reach infinity. For the Weil form, resolution costs double-exponentially (see the calibration below and Zhu's abstract).
5. **Proportion filter.** "100% of zeros" is compatible with infinitely many off-line zeros.
6. **Equivalence filter.** Restating RH is progress only if the new form has a proved partial case that the old form lacked.
7. **Sign filter (spectral routes).** In the explicit formula the prime terms enter with the sign opposite to periodic orbits in Gutzwiller/Selberg trace formulas. This is why Connes' realization is an *absorption* spectrum. A spectral route must say how it handles this sign.

---

## Calibration: what your paper's computations can and can't reach

### 1. Weil positivity at support *a* is worth RH only up to height ≈ T\*(a)

Your Heuristic 12.4 and Remark 12.7 point this way, and Zhu's abstract says so explicitly: positivity alone "cannot resolve RH ... due to exponential frequency resolution thresholds".

| probe half-width *a* | T\*(a) = 2πe^{2a} | λ(a) |
|---|---|---|
| 1.1 (your range) | 57 | 7.7 × 10⁻³⁵ |
| 1.19 (Zhu, exploratory, 950 modes) | 68 | between 10⁻⁴⁸ and 10⁻⁴⁶ |
| 13.4 | 3 × 10¹² (Platt–Trudgian) | ≈ 10^(−2.3 × 10¹²) (extrapolated) |

Matching the existing verification height would need an eigenvalue with about 2.3 trillion digits. So the computational side of this route can't compete with direct zero verification. Its value is in suggesting a mechanism that is uniform in *a*, and in testing controls (route 4).

### 2. Question 12.8: the data can't yet support a constant C

I refit your Table 8 values plus a value at a = 1.2 that an earlier reading attributed to Zhu (`rh_decay_fit.py`, numpy least squares on a ≥ 0.6). Correction: Zhu's paper reports an exploratory, non-certified floor between 10⁻⁴⁸ and 10⁻⁴⁶ at half-width 1.19 (950 modes, 70 digits) and certified two-sided bounds at 0.8 only (odd sector 8.2×10⁻¹⁵ ≤ λ ≤ 2.35×10⁻¹⁴); the fit below should be read with the last point replaced by that interval, which does not change its conclusion that the data cannot separate T\* from T\* log T\*. Zhu's own proposed law is −log λ ≃ 2π²·N(T\*)/log N(T\*), fitted on his upper bounds for 0.5 ≤ L ≤ 2.0.

| a | T\* | −log λ | ratio −log λ / T\* | local slope d(−log λ)/dT\* |
|---|---|---|---|---|
| 0.6 | 20.9 | 14.3 | 0.687 | 1.53 (0.5→0.6) |
| 0.8 | 31.1 | 31.8 | 1.022 | 1.70 |
| 1.0 | 46.4 | 59.5 | 1.281 | 1.81 |
| 1.1 | 56.7 | 78.6 | 1.385 | 1.86 |
| 1.2 | 69.3 | 102.1 | 1.474 | 1.87 |

**Model comparison** (RMS residual in nats):

- −log λ ≈ 1.816·T\* − 24.2 fits with RMS 0.53.
- −log λ ≈ 0.381·T\* log T\* − 9.2 fits with RMS 0.57.
- −log λ ≈ c·a·T\* fits worst (RMS 1.24).

Fitting on a ≤ 1.1 and predicting a = 1.2 gives 1.8 × 10⁻⁴⁴ under the linear law and 7.7 × 10⁻⁴⁶ under the T\* log T\* law. These bracket Zhu's value of 4.76 × 10⁻⁴⁵, so the data **can't distinguish the two laws**.

**Consequences:**

- **If growth is ≍ T\* log T\*,** Question 12.8 is **false even assuming RH**. No constant C works.
- **"Consistent with C = 1.5" should come out.** Under your own linear fit, the ratio crosses 1.5 near a ≈ 1.25. The local slope (1.87 and still rising) says the crossing may come sooner.
- **Safer phrasings:**
  - *"Does lim sup −log λ(a)/T\*(a) exist and is it finite?"*
  - *"Prove, assuming RH, two-sided bounds of the form exp(−C₁·Φ(a)) ≤ λ(a) ≤ exp(−C₂·Φ(a)), and identify Φ."*

  The second is a well-posed, publishable problem (route 3).

### 3. Your Open Question 5 has a ceiling

The Alpöge–Furman paper (arXiv:2608.13637) states:

- Its argument "was discovered and written by Claude", and it comes with a Lean 4 verification.
- With the Montgomery–Taylor window it reaches ≈ 0.6725.
- The *bandwidth-one ceiling* of the method is ≈ 0.682.
- "Improving on 2/3 by this route would require pair-correlation information beyond Fourier support 1."

**Follow-up work:**

- Lamzouri (arXiv:2609.02882) reproved 67.25% (simple and on the line) and 83.62% (distinct) with a single Hilbert-space inequality.
- Wang (arXiv:2609.24167) improved both by about 10⁻⁷.

So a window search inside the odd Weil matrices can gain at most about 0.01. The real lever is pair correlation past support 1 (routes 9 and 38).

### 4. Connes–van Suijlekom requires an even eigenfunction

Their real-zeros theorem (arXiv:2511.23257, CMP 2025) assumes the lowest eigenvalue is simple and isolated *with an even eigenfunction*. You work in the odd sector. An **odd-sector analogue** would be a clean, self-contained lemma, and it would make your Computation 12.5 rigorous in part. See route 6.

---

## Top five routes for your toolkit

You have high-precision Weil-form engines, certified linear algebra, and AI-assisted drafting.

1. **Route 4.** The Davenport–Heilbronn control of the horizon heuristic. It is falsifiable and computable at your precision: if zeta's decay law carries over with the conductor-5 horizon (2π/5)e^{2a}, λ is about 10⁻⁵⁶ near a ≈ 2.1. The prime side has to be rederived, since DH has no Euler product and has zeros in Re s > 1.
2. **Route 3.** Settle the decay law, and restate Question 12.8 in a form that could be true.
3. **Route 2.** Turn the horizon into a theorem in both directions. It's unconditional, new, and doesn't claim RH.
4. **Route 6.** An odd-sector version of Connes–van Suijlekom, plus the Hurwitz convergence question.
5. **Route 5.** Compute λ(a) for curves over F_q, where RH is known, and see which mechanism sustains the relay there.

---

## The 100 routes

Tier counts: 19 [A], 42 [B], 39 [C].

Format: **Route** [tier]: idea. *Blocker:* the known obstruction. *Next step:* a credible intermediate result.

### I. Weil positivity and the explicit formula (routes 1–13)

1. **Uniform relay** [B]: Prove Q_a ≥ 0 by continuation in *a*, turning Computation 2.1 into a lemma: each prime power entering at ½ log n restores the margin that the archimedean term erodes. *Blocker:* positivity for all *a* is equivalent to RH. λ(a) falls double-exponentially. By the Euler-product filter, the argument must use multiplicativity of Λ, not just the positions and sizes of its terms. *Next step:* a computer-free proof of the first hand-off (at 2), stating exactly which property of Λ is used.
2. **Localized Weil criterion** [A]: Prove Heuristic 12.4 in both directions. First, an off-line zero at height H forces λ(a) < 0 once a ≥ ½ log(H/2π) + E(H). Second, RH up to height H plus a zero-density bound forces λ(a) > 0 for a ≤ ½ log(H/2π) − E′(H). *Blocker:* sampling and interpolation estimates for Paley–Wiener functions (Cartwright–Levinson, Plancherel–Pólya) with explicit constants. *Next step:* either direction with any explicit E. Check Zhu (arXiv:2608.24827) for overlap first.
3. **The decay law of λ(a)** [A]: Determine Φ(a) = −log λ(a) up to constants, assuming RH. Is Φ ≍ T\* or ≍ T\* log T\*? This is a time–frequency concentration problem for functions of exponential type *a* sampled on the zeta zeros (Landau–Widom). *Next step:* bounds in both directions assuming RH. Then Question 12.8 becomes a statement that could be true.
4. **Davenport–Heilbronn and Epstein controls** [A]: Your Open Question 3. The horizon picture predicts that the DH floor turns negative near a ≈ 2.1. *Care:*
   - DH has zeros with Re s > 1 (Saias–Weingartner), so −f′/f has non-multiplicative coefficients and the "prime side" must be derived from scratch.
   - DH is entire, so there are no polar terms.
   - Epstein zetas of class number > 1 forms give a second control.
5. **Function-field Weil forms** [A]: Compute the analog of λ(a) for L-functions of curves over F_q, where RH is a theorem and the prime side is periodic in log q. See which mechanism sustains positivity as the support grows, then ask which ingredient is missing over ℚ. It's one of the few ways to test a proposed mechanism on a case where the answer is known.
6. **Ground states and Hurwitz** [A]: Connes–van Suijlekom (arXiv:2511.23257): if the lowest eigenvalue of such a form is simple and isolated with an *even* eigenfunction ξ, then ξ̂ has only real zeros. Groskin (arXiv:2605.20224) and your §12 show those zeros match zeta zeros to hundreds of digits. If normalized ground-state transforms converged locally uniformly on ℂ to Ξ, Hurwitz's theorem would give RH. *Blocker:* convergence off the real axis needs control of a spectral gap as small as λ(a). *Next step:* the odd-sector analogue of Connes–van Suijlekom.
7. **Zeta spectral triples** [A]: Connes–Consani–Moscovici (arXiv:2511.22755) build self-adjoint rank-one perturbations of scaling operators from Euler products. Their spectra approach the zeros numerically, and the authors note that proving convergence would prove RH. *Next step:* convergence of the regularized determinants to Ξ on compact subsets of a fixed strip.
8. **Rank–trace compressions (Alpöge–Furman)** [A]: Your Open Question 5. *Ceiling:* ≈ 0.682 at bandwidth one (see calibration §3). *Next step:* make sharp the conditional statement. The paper says a Hardy–Littlewood-type higher-moment hypothesis would give proportion 1.
9. **Lamzouri's inequality with wider pair correlation** [A]: Montgomery's pair-correlation theorem is now unconditional (Baluyot–Goldston–Suriajaya–Turnage-Butterbaugh) and feeds a single Hilbert-space inequality. Every further gain comes from Fourier support beyond 1. *Blocker:* support beyond 1 is equivalent to variance estimates for primes in short intervals (Goldston–Montgomery). Even full pair correlation gives 100%, not all.
10. **Computer-free small support** [B]: Yoshida proved positivity for small support analytically. Connes–Consani (Selecta 2021) analyze the range where only the archimedean place contributes, support [2^(−1/2), 2^(1/2)], i.e. your a ≤ ½ log 2, through the Sonin trace and prolate functions. Zhu certified positivity by computer to a = 0.8. *Next step:* push the analytic method past the first prime.
11. **Characterize the minimizer** [B]: Computation 12.5 shows the minimizer interpolates zeta zeros below T\*(a) and is exponentially small above. Prove this structurally, for example as a perturbed prolate eigenfunction. A uniform-in-*a* argument almost certainly needs a description of the minimizer, not just its eigenvalue. Related: Bondarenko–Radchenko–Seip, *Fourier interpolation with zeros of zeta and L-functions* (Constr. Approx. 2023) — the rigorous framework for "interpolating on zeros and log n".
12. **Weil forms over families** [B]: Alpöge–Furman extends to primitive Dirichlet L-functions. Averaging over a family gives "positivity on average", in the way Bombieri–Vinogradov is GRH on average. *Blocker:* averages never isolate one L-function.
13. **Li's criterion** [C]: λ_n = Σ_ρ [1 − (1 − 1/ρ)^n] ≥ 0 for all n is equivalent to RH (Li 1997; Bombieri–Lagarias 1999). *Next step:* explicit theorems of the form "RH up to height T ⇒ λ_n > 0 for n ≤ f(T)".

### II. Spectral realizations: Hilbert–Pólya (routes 14–24)

14. **Connes' adelic trace formula** [B]: RH is equivalent to positivity of the Weil distribution. The semi-local trace formula is proved, and the zeros appear as an absorption spectrum (Connes 1999; Meyer 2005). *Blocker:* the missing positivity is RH.
15. **Prolate wave operator** [B]: Connes–Moscovici (PNAS 2022): the UV spectrum of the prolate operator matches the zeros' counting asymptotics. *Blocker:* matches density, not individual zeros. *Next step:* connect to routes 6 and 7.
16. **Berry–Keating xp and Sierra's models** [B]: The semiclassical xp Hamiltonian reproduces the smooth zero count. Sierra–Townsend (Landau levels) and Sierra (2014, Rindler Dirac fermion) insert log p through boundary conditions. *Blocker:* the fluctuating part has to be derived, not inserted, and the sign filter applies.
17. **Selberg analogy** [B]: For compact hyperbolic surfaces, the Selberg zeta satisfies "RH" apart from finitely many small eigenvalues, because the Laplacian is self-adjoint. So: look for a space whose primitive orbits have lengths log p. *Blocker:* the sign filter; no such space is known.
18. **Scattering on the modular surface** [B]: Faddeev–Pavlov and Lax–Phillips: the zeros of ζ are resonances at ρ/2 of the scattering matrix ξ\*(2s−1)/ξ\*(2s). RH says they all lie on Re s = 1/4. Burnol recast this as a causality property. *Blocker:* self-adjointness doesn't constrain resonances.
19. **Pseudo-Laplacians** [B]: Colin de Verdière (1982) truncated Eisenstein series to get self-adjoint operators with zeta-related eigenvalues. Bombieri–Garrett ("Designed pseudo-Laplacians") study what such operators can and can't capture. Zero-spacing statistics limit them. *Next step:* find the obstruction's exact form.
20. **Mayer's transfer operator** [B]: The Selberg zeta of PSL₂(ℤ) is a Fredholm determinant of Mayer's Gauss-map operator, and it vanishes at s = ρ/2. So RH becomes a statement about when an explicit operator has eigenvalue ±1 on Re s = 1/4 (Lewis–Zagier period functions give the spectral side). *Blocker:* the operator isn't self-adjoint, and the ρ/2 zeros come from the continuous spectrum.
21. **de Branges spaces** [B]: de Branges' positivity conditions fail for ζ (Conrey–Li 2000). *Next step:* find the right Hilbert space of entire functions; Lagarias and Burnol's Sonine spaces are the starting points. *Blocker:* no candidate space is known to have the needed positivity.
22. **Nyman–Beurling–Báez-Duarte** [B]: RH is equivalent to an explicit L² distance d_N → 0 (dilations of the fractional-part function). Burnol showed d_N² ≥ (C + o(1))/log N, so any convergence is slow. *Next step:* a partial mechanism for restricted approximants.
23. **Bost–Connes and the primon gas** [C]: ζ(β) is the partition function of a quantum statistical system with a phase transition at β = 1. *Blocker:* the system sees ζ at real β, not the zeros.
24. **PT-symmetric Hamiltonians** [C]: Bender–Brody–Müller (2017). *Blocker:* no rigorous self-adjoint (or similar-to-self-adjoint) realization, and the eigenvalue condition presupposes the zeros (Bellissard's critique).

### III. Entire functions, Laguerre–Pólya and heat flow (routes 25–35)

25. **de Bruijn–Newman upper bounds** [A]: Known results:
    - Λ ≥ 0 (Rodgers–Tao 2020).
    - Λ ≤ 0.22 (Polymath 15) and Λ ≤ 0.2 (Platt–Trudgian 2021).
    - For every t > 0, H_t has only finitely many non-real zeros (Ki–Kim–Lee 2009).

    RH ⇔ Λ = 0. *Blocker:* each upper bound is bought with verified height and zero-free regions, and the method can't reach 0. *Next step:* lower the upper bound.
26. **Jensen polynomials** [B]: Griffin–Ono–Rolen–Zagier (2019): for each degree d, J^{d,n} is hyperbolic for n ≥ N(d). RH ⇔ hyperbolic for all d and n. *Next step:* effective N(d), uniform in d.
27. **Turán and Laguerre inequalities** [C]: Csordas–Norfolk–Varga (1986) proved the Turán inequalities for Ξ's Taylor coefficients; Dimitrov–Lucas (2011) proved higher-order ones. RH ⇔ an infinite family of such inequalities (the Laguerre–Pólya class). This connects to your §9.1.
28. **Multiplier sequences (Pólya–Schur)** [C]: Find a real-zero-preserving operator that carries a kernel with known real zeros (e.g. Pólya's ∫ e^{−λ cosh u} cos(zu) du) to Riemann's Φ. *Blocker:* Φ lies in no known closed class.
29. **Lee–Yang and Newman** [B]: Newman (1976) linked Fourier transforms with only real zeros to Lee–Yang measures of ferromagnets. Realizing Ξ as a limit of ferromagnetic partition functions would give RH by the Lee–Yang circle theorem. *Blocker:* no such model is known, and Newman's conjecture, now a theorem, says RH is at best "barely true".
30. **Real-zero approximants plus Hurwitz** [B]: Find explicit entire functions with only real zeros that converge to Ξ locally uniformly on ℂ. Approximations like Σ_{n≤X} n^{−s} + χ(s)Σ_{n≤X} n^{s−1} (studied by Gonek–Montgomery and others) have most zeros on the line. *Blocker:* convergence inside the strip is the whole problem.
31. **Lagarias positivity and horizontal monotonicity** [C]: RH is equivalent to each of:
    - Re ξ′/ξ(s) > 0 for Re s > 1/2 (Lagarias 1999);
    - |ξ(σ + it)| strictly increasing in σ > 1/2 (Sondow–Dumitrescu; Matiyasevich–Saidak–Zvengrowski).

    These are exactly your Figure 20 "valleys". *Blocker:* a pointwise equivalent.
32. **Speiser: zeros of ζ′** [B]: RH ⇔ ζ′ has no zeros in 0 < Re s < 1/2 (Speiser 1934). Levinson–Montgomery (1974) quantified this, and that work led to Levinson's method. *Next step:* sharper counts of zeros of ζ′ left of the line.
33. **High derivatives of Ξ** [B]: Conrey (1983): the proportion of zeros of Ξ^(k) on the line tends to 1 as k → ∞. *Route:* descend to k = 0. *Blocker:* Rolle's theorem only runs one way, so integrating can create non-real zeros.
34. **Balazard–Saias–Yor** [C]: RH ⇔ ∫ log|ζ(1/2 + it)| dt/(1/4 + t²) = 0, a Jensen-formula identity. *Blocker:* equivalent to RH.
35. **Heat-flow dynamics of zeros** [B]: Zeros of H_t move by a Calogero-type ODE, and off-line zeros are pulled onto the line forward in time. *Route:* control the backward flow. *Blocker:* the backward flow is ill-posed, and Lehmer pairs show zeros can nearly collide.

### IV. Proportions, correlations and zero statistics (routes 36–46)

36. **Levinson–Conrey mollifiers** [B]: 1/3 (Levinson 1974), 2/5 (Conrey 1989), 5/12 (Pratt–Robles–Zaharescu–Zeindler 2020). Farmer's θ = ∞ conjecture would give 100%. *Blocker:* the proportion filter. This line is now behind the Weil-form methods.
37. **Longer mollifiers** [B]: Mollifier length beyond Conrey's 4/7 needs off-diagonal bilinear Kloosterman bounds (Deshouillers–Iwaniec; Bettin–Chandee–Radziwiłł 2017). *Blocker:* the proportion filter.
38. **Primes in short intervals ⇒ pair correlation** [A]: Goldston–Montgomery: Montgomery's F(α) for 1 ≤ α ≤ A is equivalent to the variance of ψ in short intervals. After routes 8 and 9, this is the lever for every proportion result. *Next step:* any unconditional information on F(α) for α slightly above 1.
39. **Higher correlations** [B]: Rudnick–Sarnak n-level correlations with restricted support. *Next step:* feed them into a Lamzouri-type inequality.
40. **Ruling out the Alternative Hypothesis** [B]: The Alternative Hypothesis says zeros are asymptotically spaced at half-integer multiples of the mean spacing. Lagarias–Rodgers (2020) showed higher correlations constrain it. Ruling it out is a natural checkpoint beyond pair correlation.
41. **Simple and distinct zeros** [B]: 67.25% simple and on the line, 83.62% distinct (Lamzouri; Wang 2026). *Next step:* bounded multiplicity unconditionally. *Blocker:* says nothing about location.
42. **Selberg's sign-change method optimized** [B]: Pearce-Crump (arXiv:2609.15329) gets ≥ 7% via positive-semidefinite mollifier families. Of structural interest, not a record.
43. **Low-lying zeros in families** [B]: Katz–Sarnak and Iwaniec–Luo–Sarnak one-level density with larger Fourier support. *Blocker:* family averages.
44. **Number variance and Berry saturation** [C]: Your Computation 11.1. Berry's saturation follows from the explicit formula plus prime-pair correlations. *Blocker:* needs Hardy–Littlewood; heuristic value only.
45. **Gaps between zeros** [C]: Small and large gaps (Conrey–Ghosh–Gonek; Bui–Milinovich and others) connect to Lehmer pairs and to Λ. *Blocker:* statistical.
46. **Moments** [C]: The Keating–Snaith and CFKRS conjectures; sharp upper bounds under RH (Soundararajan, Harper) and unconditional ones (Heap–Radziwiłł–Soundararajan). *Blocker:* moments don't locate zeros.

### V. Zero-free regions and zero density (routes 47–54)

47. **The Vinogradov–Korobov exponent** [B]: Zero-free region σ > 1 − c/((log t)^{2/3}(log log t)^{1/3}). The exponent has not changed since 1958, even after Bourgain–Demeter–Guth proved Vinogradov's mean value theorem. Beating 2/3 would be major.
48. **Density hypothesis via large values** [A]: Guth–Maynard (2024): N(σ, T) ≪ T^{30(1−σ)/13+ε}. *Next step:* the density hypothesis, N ≪ T^{2(1−σ)+ε}. *Blocker:* the large-values problem; and density bounds never exclude isolated zeros.
49. **Positivity kernels in the explicit formula** [B]: Every unconditional zero-free region comes from a nonnegative kernel (3 + 4cos θ + cos 2θ ≥ 0 and its descendants; Mossinghoff–Trudgian–Yang), with kernels optimized by Beurling–Selberg methods (Carneiro–Milinovich–Soundararajan). *Blocker:* these only give regions of width ≈ 1/log t.
50. **Turán's power-sum method** [B]: An independent route to zero density near σ = 1, with explicit versions. *Blocker:* same regime as route 47.
51. **Quasi-RH** [B, as a goal]: No zeros with Re s > 1 − δ for a fixed δ > 0. This would already be a landmark, giving ψ(x) = x + O(x^{1−δ+ε}). No route is known. It's listed because it's the natural stepping stone.
52. **Lindelöf and subconvexity** [C]: Bourgain's 13/84. Lindelöf doesn't imply RH, but it's equivalent to having few zeros in boxes σ > 1/2 + ε near each height (Backlund).
53. **Pretentious methods** [B]: Granville–Soundararajan and Koukoulopoulos: zero-free statements through multiplicative functions "pretending" to be n^{it}. *Blocker:* so far the same barrier as route 47.
54. **Landau–Siegel zeros** [B]: Eliminate exceptional real zeros of L(s, χ), the real-zero special case of GRH. Y. Zhang's 2022 claim has not, to my knowledge, been accepted.

### VI. Geometry: transplanting the function-field proof (routes 55–64)

55. **Weil's Hodge-index proof** [C]: For curves, RH follows from positivity of the intersection form on C × C. Over ℤ that needs "Spec ℤ ×_{F₁} Spec ℤ". This is the central dream; no such object has been constructed.
56. **Connes–Consani–Marcolli dictionary** [B]: "Weil's proof and the geometry of the adeles class space" (2007) translates each step of Weil's proof. *Blocker:* the positivity step.
57. **Arithmetic site, scaling site, Riemann–Roch for Spec ℤ** [B]: Connes–Consani's topos-theoretic "curve", and their Riemann–Roch for the compactified Spec ℤ. *Next step:* a Hodge-index analogue.
58. **Deninger's foliated dynamical systems** [C]: A conjectural cohomology carrying a flow with periodic orbits of length log p, where RH would follow from Kähler-type positivity. *Blocker:* the space hasn't been constructed.
59. **F₁ geometry** [C]: Soulé, Borger (λ-rings), Lorscheid (blueprints), Durov, Kapranov–Smirnov, Manin. *Blocker:* none yet gives a cohomology of the right size.
60. **Deligne's Weil I strategy** [C]: Monodromy plus the tensor-power (Rankin–Selberg squaring) trick. *Blocker:* ζ isn't a member of a family with large monodromy.
61. **Bombieri–Stepanov polynomial method** [C]: An elementary proof of RH for curves via auxiliary polynomials. *Blocker:* it uses Frobenius, and there is no analogue of x ↦ x^q over ℤ.
62. **Grothendieck's standard conjectures** [C]: The Hodge standard conjecture gives RH for varieties over finite fields. *Blocker:* the number-field analogue needs a Weil cohomology for Spec ℤ.
63. **Arakelov intersection theory** [C]: The Faltings–Hriljac Hodge index theorem and Gillet–Soulé arithmetic Riemann–Roch give positivity on arithmetic surfaces. *Blocker:* Spec ℤ ×_ℤ Spec ℤ is just Spec ℤ.
64. **Prismatic and condensed settings** [C]: Bhatt–Scholze prismatic cohomology and Clausen–Scholze analytic geometry have been discussed as possible homes for a cohomology over ℤ. Speculative: I know of no concrete RH mechanism proposed there.

### VII. Automorphic forms, and toy RHs that are theorems (routes 65–72)

65. **Weng zeta functions** [B]: Lagarias–Suzuki (2006) proved RH for Weng's rank-2 zeta function of ℚ, which is built from integrals of Eisenstein series; Ki and Suzuki–Weng extended this. Some functions built from ζ provably satisfy RH. *Route:* find which feature transfers.
66. **RH for period polynomials** [B]: The zeros of period polynomials of Hecke eigenforms lie on the unit circle (Conrey–Farmer–Imamoğlu 2013; El-Guindy–Raji 2014; Jin–Ma–Ono–Soundararajan 2016). These are proven toy RHs from modular forms.
67. **Nonnegativity of central values** [B]: GRH implies L(1/2, π) ≥ 0 in self-dual families. Period formulas prove this without GRH (Waldspurger, Kohnen–Zagier, Katok–Sarnak): geometry forces an RH-shaped positivity. *Route:* find period formulas for Weil-form-type quantities.
68. **Zagier's horocycles** [B]: RH is equivalent to long closed horocycles on SL₂(ℤ)\H equidistributing at rate O(y^{3/4−ε}) (Zagier 1981; Sarnak). *Blocker:* spectral methods give y^{1/2}; the missing 1/4 is RH.
69. **Epstein families as a deformation laboratory** [B]: Epstein zeta functions vary continuously with the lattice. At Euler-product points (class number 1) GRH is expected; elsewhere there are off-line zeros (Davenport–Heilbronn). *Route:* track how zeros leave the line as the lattice deforms. This is concrete, computational, and a quantitative form of the Euler-product filter.
70. **Selberg class** [C]: Kaczorowski–Perelli: there are no functions of degree in (0, 1) or (1, 2). *Route:* axioms that imply RH. *Blocker:* the axioms must make the Euler product do real work.
71. **Rankin–Selberg positivity** [B]: Hoffstein–Lockhart and Goldfeld–Hoffstein–Lieman used positive combinations of L-functions to kill Siegel zeros for GL(2). *Blocker:* this only gives regions near σ = 1.
72. **Functoriality and symmetric powers** [C]: Gives continuation, functional equations and Ramanujan-type bounds, but nothing directly about where zeros lie. Listed because any analogue of Deligne's trick (route 60) would need it.

### VIII. Dynamics, probability and the Beurling boundary (routes 73–79)

73. **Beurling generalized primes** [B]: Diamond–Montgomery–Vorhauer (2006) is the second filter. *Route:* identify the minimal extra axiom (the additive structure of ℤ?) that forces RH. That would be a valuable meta-theorem.
74. **Bagchi's strong recurrence** [C]: RH is equivalent to ζ approximating itself in the strip in Voronin's universality sense (Bagchi 1981).
75. **Möbius randomness** [C]: RH ⇔ M(x) = O(x^{1/2+ε}). Matomäki–Radziwiłł, Tao's logarithmic Chowla and Sarnak's disjointness give cancellation, but always with savings far from √x.
76. **Random multiplicative functions and chaos** [C]: Harper (better than square-root cancellation for random multiplicative functions) and Saksman–Webb (ζ on the line converges to Gaussian multiplicative chaos). *Blocker:* these model statistics, not locations.
77. **Farey fractions** [C]: Franel–Landau: RH is equivalent to an L² discrepancy bound for Farey fractions. Boca–Cobeli–Zaharescu Farey statistics are the partial results.
78. **Redheffer matrix** [C]: det R_n = M(n), so RH is equivalent to a determinant bound. There is spectral work (Barrett–Jarvis, Vaughan). *Blocker:* the determinant is not a spectral quantity you control.
79. **Cramér's mean square** [C]: Your Computation 13.3. The log-mean of ((ψ − x)/√x)² is bounded if and only if RH holds. *Blocker:* equivalent to RH; your 10¹⁰ data is a check, not evidence of mechanism.

### IX. Elementary and arithmetic equivalents (routes 80–86)

80. **Robin's inequality** [C]: σ(n) < e^γ n log log n for n > 5040 is equivalent to RH (Robin 1984). It is proved for odd n, squarefree n and other classes (Choie–Lichiardopol–Moree–Solé 2007), and for t-free n with small t (Morrill–Platt). It is verified to astronomically large n. *Next step:* new infinite classes.
81. **Lagarias's harmonic-number inequality** [C]: σ(n) ≤ H_n + e^{H_n} log H_n (Lagarias 2002).
82. **Nicolas's primorial criterion** [C]: N_k/φ(N_k) > e^γ log log N_k for all primorials N_k (Nicolas 1983).
83. **Solé–Planat criterion** [C]: The Dedekind ψ-function on primorials (Solé–Planat 2011).
84. **Schoenfeld's explicit π(x) bound** [C]: |π(x) − li(x)| < √x log x/(8π) for x ≥ 2657.
85. **Riesz's series** [C]: Σ (−1)^{k+1} x^k/((k−1)! ζ(2k)) = O(x^{1/4+ε}) (Riesz 1916). Hardy–Littlewood gave a variant.
86. **The Π⁰₁ form** [C]: RH is equivalent to an explicit Diophantine equation having no solutions (Davis–Matiyasevich–Robinson). So if RH is independent of ZFC, it is true. *Blocker:* there are no tools for proving such independence. Listed for completeness.

### X. Physics-inspired routes (routes 87–90)

87. **Modular bootstrap** [C]: Benjamin–Chang (JHEP 2022) relate crossing equations for Narain CFTs to ζ zeros. *Blocker:* bootstrap bounds are inequalities, not locations.
88. **Adelic string amplitudes** [C]: The Freund–Witten product formula over all places. Speculative.
89. **Semiclassical quantum chaos** [C]: Bogomolny–Keating and Berry derive the zero statistics from prime-pair correlations. *Blocker:* heuristic, and relies on Hardy–Littlewood.
90. **Fourier quasicrystals and crystalline measures** [B]: The explicit formula says that Σ δ_γ has a Fourier transform supported on ±log n plus a smooth part, and RH makes that measure real-supported. Kurasov–Sarnak built crystalline measures from stable (real-rooted) polynomials, and Olevskii–Ulanovskii and Lev–Olevskii characterize Fourier quasicrystals. *Route:* characterize the measures with this Fourier support and positivity, and show the zeta zeros must be one of them.

### XI. Computation, certification and formalization (routes 91–100)

None of these prove RH. All are credible, citable contributions.

91. **Extend rigorous verification past 3 × 10¹²** [A]: Platt–Trudgian. This feeds route 25 and explicit prime bounds.
92. **High-height statistics** [A]: Odlyzko-style computations near 10²³ and beyond: Lehmer pairs, GUE fit, Gram-law failures.
93. **Formalize the Weil criterion and explicit formula in Lean** [A]: Alpöge–Furman already ships Lean 4 proofs, and the Kontorovich–Tao PNT+ project is building the analytic background.
94. **Certify your relay** [A]: Recompute Table 1 in interval arithmetic (Arb/FLINT) so that Computation 2.1 becomes a certified statement, as Zhu did for positivity.
95. **Certified λ(a) beyond a = 1.2** [A]: Worth doing, but calibrate claims: a = 1.2 corresponds to height ≈ 70.
96. **Machine search with formal verification** [A]: The Alpöge–Furman argument was found by Claude and checked in Lean. Searching window, mollifier and test-function families with machine help plus formal checking is now a demonstrated route to proportion results, not to RH.
97. **Certified Li coefficients** [C]: λ_n to large n, with rigorous error bounds tied to verified height (pairs with route 13).
98. **Robin and Lagarias verification ranges** [C]: Push the certified ranges for routes 80 and 81. These make useful constants, not evidence for RH.
99. **The DH/Epstein Weil-form atlas** [A]: The computational half of routes 4 and 69: a public, certified table of Weil floors for functions *with* off-line zeros. It would be the first benchmark that any claimed positivity mechanism has to fail on.
100. **Explicit horizon constants** [A]: The computational half of route 2: measure E(H) empirically for synthetic "zeta plus one off-line quartet" models, as in your Computation 13.3(4). This gives the constants a theorem should aim for.

---

## Traps: routes that look promising and are known dead ends

- **The Mertens conjecture** |M(x)| < √x is false (Odlyzko–te Riele 1985). RH needs only x^{1/2+ε}.
- **The Pólya conjecture** is false (Haselgrove 1958).
- **de Branges' positivity conditions** fail for ζ (Conrey–Li 2000).
- **Symmetry-only arguments** fail the Euler-product filter (Davenport–Heilbronn; see your §8).
- **"Proof by 100%"** fails the proportion filter.
- **"Proof by verification"** fails the uniformity filter, whatever the height.
- **Operators that encode the zeros by fiat**, where the boundary condition or domain presupposes ζ(ρ) = 0.
- **Probabilistic "RH holds with probability 1" arguments** (Denjoy's heuristic) are heuristics, not proofs.

## Sources checked for this note (2026-09-29)

These 2026 items postdate my training, so I checked each one on arXiv:

- Alpöge–Furman, [arXiv:2608.13637](https://arxiv.org/abs/2608.13637). Confirmed: 2/3 simple and on the line; Lean 4; the argument was "discovered and written by Claude"; ceiling ≈ 0.682 at bandwidth one.
- Zhu, [arXiv:2608.24827](https://arxiv.org/abs/2608.24827). Confirmed: certified two-sided bounds; the "frequency resolution" limitation.
- Connes–Consani–Moscovici, [arXiv:2511.22755](https://arxiv.org/abs/2511.22755). Confirmed: convergence would prove RH.
- Connes–van Suijlekom, [arXiv:2511.23257](https://arxiv.org/abs/2511.23257). Confirmed: real zeros for simple, isolated, *even* ground states.
- Groskin, [arXiv:2605.20224](https://arxiv.org/abs/2605.20224). Confirmed.
- Lamzouri, [arXiv:2609.02882](https://arxiv.org/abs/2609.02882). Confirmed: 67.25% and 83.62%.
- Wang, [arXiv:2609.07918](https://arxiv.org/abs/2609.07918) and [arXiv:2609.24167](https://arxiv.org/abs/2609.24167).
- Pearce-Crump, [arXiv:2609.15329](https://arxiv.org/abs/2609.15329). Confirmed: 7% via Selberg's method.
- Connes–Consani, [arXiv:2006.13771](https://arxiv.org/abs/2006.13771) (archimedean place).
- Bondarenko–Radchenko–Seip, [arXiv:2005.02996](https://arxiv.org/abs/2005.02996).

The other attributions come from the standard literature. Check any one before citing it in the paper. Where I wasn't certain of a detail, I described the result rather than quoting it.
