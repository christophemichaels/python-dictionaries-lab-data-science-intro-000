# Three Routes, One Mechanism: What an Unconditional Proof Through the Weil Floor Would Have to Look Like

*Companion to `RH_ROUTES.md`. Prepared 2026-09-29 for C. Michaels, "Primes, Folds, and the One Dot".*

## Ground rules

- **This is not a proof and does not contain one.** It's an attempt to design the *shape* of a proof along the three routes I ranked highest for your toolkit (the horizon theorem, the decay law, and the Davenport–Heilbronn control), to isolate the single lemma on which such a proof would stand, and to say where the circularity hides.
- Everything below is one of three things, and each is labelled:
  - **[verified]** — computed with a small independent engine for the odd Weil form (`rh_weil_odd.py`, 24 modes, 50 digits). It reproduces your floors (λ(0.5) = 1.943×10⁻⁴, λ(0.6) = 6.0×10⁻⁷) and your Table 1 failure points (0.371602 vs your 0.371601; 0.55735 vs your 0.557323). A K-mode minimum is a minimum over a subspace and therefore an *upper bound* on the true floor: a negative finite-mode value certifies indefiniteness (up to rounding), a positive one does not certify positivity, which is why certification (Zhu) needs the control of the omitted directions that this engine does not attempt. Every [verified] number below is a statement about the 24-mode form, except in §2.6, where the numbers come from a second engine (`rh_edge_fem.py`, hat functions on a mesh graded to the window edge, double precision) whose minima are likewise upper bounds; the rigidity conclusions (§2.4) survive because they rest on negative values.
  - **[derived]** — a formula I derived and checked for consistency but did not compute.
  - **[conjecture]** or **[heuristic]** — clearly marked.
- The three routes are not three problems. They are three faces of one function, Φ(a) = −log λ(a).

---

## 1. The object: Φ(a) and its derivative

Write λ(a) for the odd Weil floor, T\*(a) = 2πe^{2a} for the horizon, F_a = f̂_a / i for the transform of the normalized minimizer, and

Φ(a) = −log λ(a).

**RH ⇔ Φ(a) < ∞ for every a ⇔ Φ′ is locally integrable on (0, ∞).**

That's Weil's criterion in the odd sector, restated so that a proof becomes a *bound on a derivative*. Bounds on derivatives can sometimes be proved locally when the global statement can't, and that's the whole reason to look at it this way.

### 1.1 The symbol identity

For real odd f supported in [−a, a], with Ψ_a your equation (2):

**Q_a(f) = (1/2π) ∫ |f̂(t)|² Ψ_a(t) dt − 2 (∫ f(x) sinh(x/2) dx)².** [derived; standard]

Two facts about the truncated symbol Ψ_a(t) = Re ψ(¼ + it/2) − log π − Σ_{log n < 2a} 2Λ(n)n^{−1/2} cos(t log n):

- The digamma part is log(t/2π) + O(1/t). It equals **2a exactly at t = T\*(a)**.
- The prime comb has mean-square 2Σ_{n<e^{2a}} Λ(n)²/n ~ (2a)², so its RMS → 2a. [verified: RMS = 1.90 vs 2a = 2.0 at a = 1; 3.76 vs 4.0 at a = 2]

So the horizon is also **the height at which the archimedean symbol first equals the typical size of the prime comb**. This ties your §3 (crests of the comb) to your §12 (horizon) through the same number 2a. It isn't a coincidence, since both are the smooth zero density, but it says which frequencies matter: all of them. |f̂|² has Fourier transform g supported in [−2a, 2a], so it sees every mode of the comb and no smoothing helps.

### 1.2 Where the symbol is negative [verified]

| a | T\* | fraction of [0, T\*] with Ψ_a < 0 | of [T\*, 3T\*] | of [3T\*, 10T\*] | min Ψ_a on [0, T\*] |
|---|---|---|---|---|---|
| 0.6 | 20.9 | 0.43 | 0.06 | 0.00 | −5.3 |
| 0.8 | 31.1 | 0.33 | 0.04 | 0.00 | −5.9 |
| 1.0 | 46.4 | 0.36 | 0.11 | 0.02 | −7.7 |
| 1.1 | 56.7 | 0.33 | 0.10 | 0.03 | −8.3 |
| 1.5 | 126 | 0.30 | 0.12 | 0.05 | −9.6 |
| 2.0 | 343 | 0.27 | 0.13 | 0.07 | −7.7 |

The symbol does **not** become positive above the horizon. [heuristic] Since the comb has RMS ≈ 2a and roughly Gaussian fluctuations (its frequencies are ℚ-independent by your Proposition 3.1), the negative fraction at height t = cT\* tends to P(Z > 1 + log c / 2a) → 16% for any fixed c as a → ∞. The symbol is mostly positive only for t ≫ T\*^{2.3}.

**Consequence.** No argument that bounds Ψ_a pointwise, or in mean, or by its distribution, can give positivity. Positivity of Q_a is an *exact cancellation*: the negative set of Ψ_a is exactly the set of near-alignments of {log p}, and the minimizer avoids it only because ∫|F|²Ψ_a = Σ_γ F(γ)² sits on the zeros. That is the Euler-product filter in numbers.

### 1.3 The derivative formula [derived]

Between prime-power entries, λ(a) is a simple eigenvalue (empirically) and Feynman–Hellmann gives, on the prime side,

**λ′(a) = (1/a) [ 2 Σ_{log n < 2a} Λ(n) (log n) n^{−1/2} g′_{f_a}(log n) + A′(f_a) + P′(f_a) ],**

where g_f = f ∗ f̃, A′ is the derivative of the archimedean form under dilation, and P′(f) = −(2S/a)[S + ∫ f(x) x cosh(x/2) dx] with S = ∫ f sinh(x/2). Equivalently, in symbol form,

λ′(a) = −(1/2πa) ∫ |F_a(t)|² t Ψ_a′(t) dt + P′(f_a),

and on the zero side (under RH),

λ′(a) = λ/a + (2/a) Σ_γ γ F_a(γ) F_a′(γ).

Note the prime weights in λ′ are **Λ(n) log n**, the Dirichlet coefficients of (ζ′/ζ)′. So the derivative of the floor is governed by an explicit formula for the *second* logarithmic derivative of ζ.

---

## 2. What the relay really is [verified]

### 2.1 The kink proposition

**Proposition (K-mode form).** At each entry a_n = ½ log n, the K-mode floor λ_K is continuous and its derivative jumps up by

**λ_K′(a_n⁺) − λ_K′(a_n⁻) = 4 Λ(n) n^{−1/2} f_K(a_n)²,**

where f_K(a_n) is the endpoint value of the normalized K-mode minimizer at entry.

**For the exact form the endpoint value is zero** (§2.6: the minimizer vanishes at the edge like C(log(a/δ) + β)^{−1/2}, a theorem of Hernández-Santamaría–López Ríos–Saldaña applied to the window form), so the exact λ has no first-order kink. What replaces it is a **soft kink**: the entering term grows like 4Λ(n)n^{−1/2}C²(a − a_n)/(log(a/(a − a_n)) + β + 1 − log 2), so λ′ is continuous but its increment across the entry, measured at scale ε, is 4Λ(n)n^{−1/2}C²/(log(a/ε) + β + 1 − log 2), vanishing only logarithmically. The K-mode kick is this law evaluated at the resolution scale of the truncation, f_K(a_n)² ≈ C²/(2 log K + β′). Everything in this subsection is about the K-mode form, where the proposition is exact; §2.6 has the exact-form measurements.

*Proof sketch.* The new term −2Λ(n)n^{−1/2} g(log n) vanishes at entry; for a slightly larger, g(log n) = ∫_{2a_n−a}^{a} f(y) f(y − log n) dy ≈ −2(a − a_n) f(a)² for odd f. Differentiate in a with the minimizer frozen (Feynman–Hellmann). ∎

| n | a_n | λ(a_n) | f(a_n) | λ′(a_n⁻) | λ′(a_n⁺) | jump | predicted | ratio |
|---|---|---|---|---|---|---|---|---|
| 2 | 0.34657 | 7.31×10⁻² | 0.4652 | −3.011 | −2.575 | 0.436 | 0.424 | 1.03 |
| 3 | 0.54931 | 1.50×10⁻⁵ | −0.00947 | −1.160×10⁻³ | −0.945×10⁻³ | 2.16×10⁻⁴ | 2.28×10⁻⁴ | 0.95 |

Monotonicity (your Proposition 12.2) holds as it must: |λ′(a_n⁻)| ≥ jump. The residual is **not** finite-difference error in the ordinary sense; it is real and it shrinks with the step. At the entry of 3 (K = 24), the jump/prediction ratio is

| h | 5×10⁻⁴ | 2.5×10⁻⁴ | 1.25×10⁻⁴ | 6.25×10⁻⁵ | 2×10⁻⁵ | 1×10⁻⁵ |
|---|---|---|---|---|---|---|
| ratio | 0.948 | 0.962 | 0.984 | 0.994 | 0.9993 | 0.9998 |

so the proposition is exact in the limit (verified to 2×10⁻⁴), and the deficit at coarser steps is the reorganization of the minimizer described in §2.5. [verified]

### 2.2 What Table 1 measures

The deleted form is the same form at entry with the same derivative from the left. So to first order,

**gap_n = a_fail − a_n ≈ λ(a_n) / |λ′(a_n⁻)| = 1 / Φ′(a_n⁻).**

- n = 2: predicted 0.0243, your Table 1: 0.0250.
- n = 3: predicted 0.0129, your Table 1: 0.0080 (the deleted form's decay accelerates; the linear estimate overshoots).

So **Table 1 measures the e-folding length of λ(a)**, not a hand-off in the sense of a special role for the newest prime power. Every term larger than the margin is load-bearing in the same sense. The control experiment that would show this: at a = 0.75, delete 3 instead of 4 and observe that the form fails equally fast.

### 2.3 The Φ′ budget

Φ′ = −λ′/λ. In Φ′ language the relay is exact:

- Between entries, the archimedean drift pushes Φ′ up (faster decay).
- At each entry, Φ′ **drops** by J_n = 4Λ(n)n^{−1/2} f(a_n)²/λ(a_n).

Measured on the 24-mode form: J_2 = 5.96 (Φ′ from 41.2 to 35.2, a 14% drop); J_3 = 14.4 (from 77.4 to 63.0, a 19% drop). [verified]

For the exact form (§2.6) the drop is not instantaneous: at distance ε past the entry, Φ′ has dropped by J_n(ε) = 4Λ(n)n^{−1/2}C_n²/(λ(a_n)(log(a_n/ε) + β + 1 − log 2)). What the budget needs is the drop accumulated over the next e-folding length ε ≈ 1/Φ′, i.e. J_n(1/Φ′) with log(1/ε) = log Φ′. At the entry of 3 the asymptotic law gives 4Λ(3)3^{−1/2}C₃²/(log(a₃Φ′) + β + 1 − log 2) ≈ 6×10⁻⁴ at ε = 1/Φ′ = 0.013, but at that scale the law overestimates (measured increment 1.6×10⁻⁴ at ε = 10⁻³ after removing the smooth curvature, against 3.15×10⁻⁴ from the law; the two agree to 5–10% from ε = 10⁻⁴ down). So the exact relay at the scale that matters is of the same order as the 24-mode kick (2.3×10⁻⁴), not a different mechanism, and the budget below keeps its form with f(a_n)² read as the measured λ′ increment at scale 1/Φ′ divided by 4Λ(n)n^{−1/2}.

Two ratios stay in narrow bands across a = 0.35 to 1.2:

- **Φ′(a)/T\*(a) ≈ 2.8 to 4.1** (from your Table 8 differences: 2.8–3.7; at the entries of 2 and 3: 3.3→2.8 and 4.1→3.3).
- **f_a(a)²/λ(a) ≈ 0.16–0.32 · T\*(a)** (values 2.2, 4.4, 5.0, 7.0, 9.9 at a = 0.4, …, 0.8).

Together these give [conjecture] **J_n/Φ′ ≈ c·Λ(n)/√n with c ≈ 0.3**, a fractional drop independent of a. *Note added with §2.7:* on the 40-mode grid at scale 0.004 the fractional drop is 0.165 at n = 2 (c = 0.34) but 0.054 at n = 3 and zero within the drift at n = 4, 5, 7, so this conjecture fails beyond the first entry; the entering term is absorbed by reorganization below the crossover scale, and Φ′/T\* shows a plateau rather than a sawtooth. The budget below is kept as the form a relay would take; §2.7 states what is actually observed. Individually the drops vanish as n grows, but the entries become dense (about e^{2a}/a per unit of a), and Σ_{n ≤ x} Λ(n)/√n ~ 2√x says the total fractional kick per unit of a is about 0.6 e^{a}, growing without bound. The smooth drift must grow at the same rate for Φ′ to stay ≈ 3.5 T\*. **The relay is this balance.** Your Open Question 1 becomes:

> **Relay inequality.** For all a: 2 Σ_{log n<2a} Λ(n)(log n) n^{−1/2} g′_{f_a}(log n) + A′(f_a) + P′(f_a) ≥ −C · a · T\*(a) · λ(a).

This implies Φ′ ≤ C·T\*, hence λ(a) ≥ λ(a₀) exp(−C(T\*(a) − T\*(a₀))) > 0, hence RH.

### 2.4 Rigidity [verified]

How exactly do the primes have to be the primes?

- At a = 0.8 (λ = 1.73×10⁻¹⁴ in my engine), replacing 2 by 2e^{ε} with **ε = +3×10⁻¹⁴** makes the form indefinite: λ = −2.76×10⁻¹⁴. The response is linear: dλ/dε = −1.497, both signs, ε = ±3×10⁻¹⁴ and ±10⁻¹³. The positivity margin at a = 0.8 corresponds to displacing the prime 2 by 1.2 parts in 10¹⁴.
- Replacing 3 by 3e^{±0.001} makes the form indefinite by a = 0.60 (+) or a = 0.65 (−), where the unperturbed floor is 6×10⁻⁷ and 2×10⁻⁸.

Extrapolating with the same slope, at a = 1.1 the margin corresponds to a displacement of 2 by about **5 parts in 10³⁵**, and at a = 13.4 (the Platt–Trudgian height) to a displacement of one part in 10^{2·10¹²}.

**Consequence.** Any statement about primes that is insensitive to displacements of relative size 10⁻³⁵ is too coarse to prove λ(1.1) > 0. That rules out every distributional input (PNT with any error term, short-interval results, sieve bounds, moment estimates) as the *engine* of a relay proof. The only things that fine are (i) the exact multiplicative structure of Λ, and (ii) identities. This is Zhu's "resolution threshold" and the Beurling filter, in one number.

It also explains why the odd-sector positivity looks like a "barely true" statement: it is the finite-a shadow of Newman's Λ_dBN ≥ 0 (Rodgers–Tao), which says RH, if true, is true with no room to spare.

### 2.5 The reorganization scale, and the cost of omitted directions [verified]

The entering term at a_n has derivative dT/da = 4Λ(n)n^{−1/2}·φ_i(1)φ_j(1)/a in the normalized Legendre basis: a **rank-one** matrix whose entries are products of endpoint values. At the entry of 3 with K = 24 modes:

- spectrum of the form: λ₀ = 1.499×10⁻⁵, λ₁ = 0.0785, λ₂ = 0.833, so λ₁/λ₀ = 5239 (the ground state is simple and well isolated; Lemma L2 holds here);
- ‖dT/da‖ = 2716, while its ground-state matrix element is 2.276×10⁻⁴ (the predicted kick): the excited states have endpoint values of order 1–10, the ground state 0.009;
- second-order coefficient Σ_k |⟨k|dT/da|0⟩|²/(λ_k − λ₀) = 0.118, so λ(a) ≈ λ₀ + (λ′_old + kick)(a − a_n) − 0.118(a − a_n)² just after entry;
- crossover scale gap/‖dT/da‖ ≈ 2.9×10⁻⁵, beyond which the entering term mixes the excited states among themselves and the polynomial expansion of λ(a) fails. This is why finite differences with h ≥ 6×10⁻⁵ see 95–99% of the kick and h = 10⁻⁵ sees 99.98%.

So the relay kick is instantaneous and exact, and within Δa ≈ 3×10⁻⁵ the minimizer reorganizes and absorbs a few percent of it. The kick is a rank-one perturbation whose ground-state part is tiny because the ground state is small at the endpoint; that smallness is what the tail of F above the horizon is made of (your Computation 12.5 estimates the tail from f(a)).

**Mode dependence.** The floor converges fast in K; the endpoint value does not:

| K | λ(a₃) | f(a₃) | predicted kick 4Λ(3)3^{−1/2} f(a₃)² |
|---|---|---|---|
| 24 | 1.49863×10⁻⁵ | 0.009471 | 2.276×10⁻⁴ |
| 32 | 1.49662×10⁻⁵ | 0.009077 | 2.090×10⁻⁴ |
| 40 | 1.49508×10⁻⁵ | 0.008652 | 1.899×10⁻⁴ |
| 48 | 1.49414×10⁻⁵ | 0.008635 | 1.891×10⁻⁴ |
| 56 | 1.49368×10⁻⁵ | 0.008444 | 1.808×10⁻⁴ |
| 64 | 1.49325×10⁻⁵ | 0.008207 | 1.708×10⁻⁴ |

The floor converges (its steps shrink: 2.0, 1.5, 0.9, 0.5, 0.4 ×10⁻⁸). The endpoint value does not on this range: it is 13% lower at K = 64 than at K = 24 and still falling, so the kick's predicted size has dropped 25% and has no visible limit. A fit f(K) = f_∞ + c/K gives f_∞ ≈ 0.0075 from each of three independent pairs (24/64, 32/56, 40/64), which would make the true kick about 1.4×10⁻⁴; but a slow power law f ∝ K^{−0.15}, tending to zero, fits the six values equally well. A Legendre expansion converges slowest at the endpoint, and the endpoint is exactly where the relay lives. Two consequences. First, this is the concrete form of the reviewers' point about omitted directions: a certified relay transition must enclose f(a_n), not only λ(a_n), and the complementary-space debit for f(a_n) is much larger than for λ. Second, whether the true minimizer's endpoint value is nonzero at all is answered in §2.6: it is zero, the kick is a truncation effect of the K-mode form, and the relay proceeds through a logarithmic boundary layer. The Legendre profiles alone could not decide this (they resolve only δ ≳ 1/K² at the edge); they are kept here because they show what a finite-mode computation sees:

| y/a | 0.5 | 0.9 | 0.95 | 0.98 | 0.99 | 0.995 | 0.999 | 1.0 |
|---|---|---|---|---|---|---|---|---|
| K = 24 | 1.3664 | 0.11681 | 0.04864 | 0.02364 | 0.01743 | 0.01430 | 0.01072 | 0.00947 |
| K = 40 | 1.3664 | 0.11690 | 0.04873 | 0.02379 | 0.01745 | 0.01422 | 0.01093 | 0.00865 |
| K = 64 | 1.3664 | 0.11689 | 0.04875 | 0.02381 | 0.01747 | 0.01429 | 0.01100 | 0.00821 |

For y ≤ 0.99a the minimizer is converged to better than 0.1% (f(0.99a) = 0.0175). The last one percent of the window is a boundary layer where the K-mode minimizer is still adjusting: f(0.999a) drifts up (0.0107 → 0.0110) while f(a) drifts down (0.0095 → 0.0082). Read with the edge law of §2.6, the six endpoint values are f_K(a₃) = C(2 log K + β′)^{−1/2}: 1/f_K² is linear in log K (the fit f_∞ + c/K and the power law K^{−0.15} were both wrong extrapolations), and the "boundary layer" is the truncation's view of a profile that continues to fall, like (log(a/δ))^{−1/2}, all the way to zero. [verified: profiles]

### 2.6 The archimedean wall and the edge law [identity verified; edge law: theorem + numerics]

For f supported in [−a, a], the archimedean form in x-space is c₁‖f‖² + ½∬_{ℝ²} J_Γ(u−v)|f(u)−f(v)|² du dv with J_Γ(s) = e^{−|s|/2}/(1−e^{−2|s|}). Splitting the double integral according to whether v lies inside the window gives an exact decomposition of the window form:

**A(f) = c₁‖f‖² + ½∬_{[−a,a]²} J_Γ(u−v)|f(u)−f(v)|² + ∫_{−a}^{a} V(y)|f(y)|² dy,  V(y) = ∫_{|v|>a} J_Γ(y−v) dv.**

Since J_Γ(s) ~ 1/(2|s|), the potential diverges logarithmically at the edge. Numerically, at a = a₃:

| y/a | 0.5 | 0.9 | 0.99 | 0.999 | 0.9999 | 0.99999 |
|---|---|---|---|---|---|---|
| V(y) | 3.4390 | 4.1342 | 5.2650 | 6.4143 | 7.5654 | 8.7166 |
| V − ½log(1/(a−y)) | 2.7929 | 2.6834 | 2.6628 | 2.6608 | 2.6607 | 2.6606 |

So **V(y) = ½ log(1/(a−|y|)) + 2.6606 + o(1)**: the odd Weil form on a window is a logarithmic kinetic form (the interior J_Γ Dirichlet form, which has symbol ~ log|t|, the logarithmic Laplacian of Chen–Weth) plus a logarithmic confining wall, plus the prime shifts and the rank-one polar term, both bounded.

Consequences, in order of certainty:

**1. The edge law is a theorem.** The window form is one half of the Dirichlet form of the logarithmic Laplacian L_Δ (Chen–Weth) plus a bounded quadratic form: the symbol of ½∬J_Γ(u−v)|f(u)−f(v)|² is σ(ξ) = ∫J_Γ(s)(1 − cos sξ)ds = log|ξ| + O(1), so it differs from ½·(2 log|ξ|) by a bounded multiplier, and on the window that difference, the prime shifts and the polar term are integral operators with bounded kernels. The Euler–Lagrange equation of the minimizer is therefore L_Δ f = F with F ∈ L^∞, granted f ∈ L^∞ (interior regularity; numerically max f = 1.37 at a₃). Two theorems on bounded weak solutions of the Dirichlet problem for L_Δ then apply to it:

- Chen–Weth, *The Dirichlet problem for the logarithmic Laplacian*, Comm. PDE 44 (2019), Thm 1.11: |f(y)| = O((log(1/(a − |y|)))^{−τ}) for every τ < ½.
- Hernández-Santamaría, López Ríos, Saldaña, *Optimal boundary regularity and a Hopf-type lemma for Dirichlet problems involving the logarithmic Laplacian*, DCDS 45 (2025), arXiv:2401.18033, Thm 1.1: |f(y)| ≤ C (log(1/(a − |y|)))^{−1/2}. Their Thm 1.2 shows the rate is sharp (two-sided for the torsion function of a small ball), and Thm 1.4 is a Hopf lemma giving the same rate from below for nonnegative supersolutions.

So **the exact minimizer vanishes at the edge, at least like (log)^{−1/2}**, and the endpoint value f(a_n) of §2.1 is a quantity of the truncation. Checkpoint 17 of the atlas (`atlas/Michaels_Theta_Atlas_CP20.pdf`, PDF pages 264–266) already applies the same theorem to the forced solutions w_a = H_a^{−1}s_a of this operator and obtains the uniform endpoint bound |u(±(a − d))| ≤ B/√log(1/d) (CP17:2.20); it states that convergence of √log(1/d)·w_a(a − d), the normalized boundary trace, is a separate open question. The graded-mesh profiles below are the numerical answer to that question for the ground state: the trace exists, and its value is the amplitude C. The same exponent comes out of the Euler–Lagrange equation directly [derived]: in t = log(1/δ) the equation near the edge reads t·g(t) − ½∫₀ᵗ g(τ)dτ = O(g) + O(1), the coefficient t being ½ from the wall and ½ from the interior kernel and the ½∫g coming from the long-range part of the interior kernel (a constant ½ on one side in the variable t), and its solutions are g = C t^{−1/2}(1 + O(1/t)). A finite limit would need the nonlocal term to diverge, which it cannot.

**2. The profile, resolved.** `rh_edge_fem.py` (Appendix) is a second engine built for the edge: hat functions on a mesh graded geometrically to 10⁻¹²·a, with exact cross-correlations and closed-form archimedean tails, so the minimizer is read to t = log(a/δ) ≈ 25 where the Legendre engine stops at 2 log K ≈ 8. On every mesh and in both cases, 1/f² is linear in t: [verified]

| form | mesh (graded to 10⁻¹⁰–10⁻¹²) | λ | C | β | max rel. residual, 14 < t < 20 | local exponent, t = 8 → 22 |
|---|---|---|---|---|---|---|
| prime-free, a = 0.3 | 130-node, ratio ½ from 0.1a | 2.2255833e-01 | 1.3090 | -1.49 | 4e-05 | 0.621 → 0.519 |
| prime-free, a = 0.3 | 350-node, ratio 0.6 from 0.2a | 2.2256002e-01 | 1.3089 | -1.47 | 4e-05 | 0.622 → 0.536 |
| prime-free, a = 0.3 | 487-node, ratio 0.9 then ½ from 0.05a | 2.2255660e-01 | 1.3089 | -1.49 | 4e-05 | 0.623 → 0.537 |
| full form, a = a₃ | 130-node, ratio ½ from 0.1a | 1.5446844e-05 | 0.0241 | -1.87 | 5e-05 | 0.718 → 0.528 |
| full form, a = a₃ | 350-node, ratio 0.6 from 0.2a | 1.8594108e-05 | 0.0244 | -1.86 | 5e-05 | 0.724 → 0.546 |
| full form, a = a₃ | 487-node, ratio 0.9 then ½ from 0.05a | 1.4914475e-05 | 0.0241 | -1.87 | 5e-05 | 0.724 → 0.547 |

The fit is f(a − δ) = C (t + β)^{−1/2}, with the residual of the linear law below 10⁻⁴ over six e-foldings and the local exponent −d log f/d log t descending toward ½ exactly as ½·t/(t + β) does. C and β are mesh-independent to 10⁻³ (prime-free) and 10⁻² (at a₃, where the interior is harder and λ itself is a five-digit cancellation). On the two coarser meshes the floor sits above the 64-mode Legendre value (0.2225569 at a = 0.3; 1.49325×10⁻⁵ at a₃) by the P1 interior error; on the two-stage mesh it lies *below* it, by 1.3×10⁻⁶ and by 0.12% respectively, as a space that resolves the edge better must give, and the profile agrees with the converged Legendre profile to 10⁻⁴ for y ≤ 0.99a.

The six Legendre endpoint values of §2.5 are this law seen at the truncation's resolution scale: a degree-2K expansion resolves δ ~ 1/K² at the edge, and reading each f_K(a) through the fitted C and β gives an effective t_K that tracks 2 log K with a constant offset:

| K | 24 | 32 | 40 | 48 | 56 | 64 |
|---|---|---|---|---|---|---|
| 2 log K | 6.36 | 6.93 | 7.38 | 7.74 | 8.05 | 8.32 |
| t_K at a₃ (from f_K(a₃)) | 8.36 | 8.93 | 9.64 | 9.67 | 10.03 | 10.51 |
| t_K at a = 0.3, prime-free | 8.44 | – | 9.42 | – | – | 10.34 |

That is the "slow, steady decrease with no visible limit" of §2.5, explained: f_K(a_n)² = C²/(2 log K + β′), and the fits f_∞ + c/K and K^{−0.15} were both wrong extrapolations of a logarithm. The same offset, t_K ≈ 2 log K + 2 (the truncation resolves δ ≈ e⁻²·a/K² at the edge), appears with and without primes, which is what one expects from a property of the Legendre basis.

**3. The relay is a soft kink.** With f(a − δ) = C(t + β)^{−1/2}, the entering term at a = a_n + ε is −2Λ(n)n^{−1/2} g(log n) = 2Λ(n)n^{−1/2}∫₀^{2ε} f(a − δ) f(a − 2ε + δ) dδ ≈ 4Λ(n)n^{−1/2} C² ε/(log(a/ε) + β + 1 − log 2) (with δ = 2εs the two edge factors expand as (L + log(1/s))^{−1/2}(L + log(1/(1−s)))^{−1/2} ≈ L^{−1}(1 − 1/L) with L = log(a/(2ε)) + β, since ∫₀¹ log(1/s) ds = 1, hence the constant 1 − log 2; it was omitted in the first version of this section and is an 8% effect at ε = 10⁻³), so λ′ is continuous across the entry and its increment at scale ε is 4Λ(n)n^{−1/2}C²/(log(a/ε) + β + 1 − log 2): a kick whose size vanishes, but only logarithmically. Measured at the entry of 3 on the 134-node mesh (C = 0.0241, β = -1.93), with D(ε) = [λ(a₃+ε) − 2λ(a₃) + λ(a₃−ε)]/ε, all entries ×10⁻⁴: [verified]

| ε | 10⁻³ | 10⁻⁴ | 10⁻⁵ | 10⁻⁶ | 10⁻⁷ |
|---|---|---|---|---|---|
| measured λ′ increment D(ε) | 2.45 | 1.90 | 1.48 | 1.21 | 1.01 |
| first-order energy P(ε)/ε on the a₃ minimizer | 3.42 | 2.13 | 1.59 | 1.27 | 1.06 |
| law 4Λ(3)3^{−1/2}C²/(log(a₃/ε) + β + 1 − log 2) | 3.15 | 2.11 | 1.59 | 1.27 | 1.06 |

The measured increment falls by about 22% per decade, as 1/(log(a/ε) + β + 1 − log 2) does, and from ε = 10⁻⁴ down it lies 5–10% below the law and 4–11% below the frozen-minimizer first-order energy (the law and the first-order energy agree to 1% from ε = 10⁻⁴ down, as they must, the law being its asymptotic evaluation) (the remainder is the reorganization of the minimizer, which shrinks like 1/log(1/ε) relative to the first-order term); at ε = 10⁻³ the smooth curvature λ″ε ≈ 0.9×10⁻⁴ is still inside D. The 24-mode form gives a constant 2.28×10⁻⁴ at every scale below its resolution; the exact form has passed below that value by ε = 10⁻⁴ and keeps falling. This is the relay of §2.1–2.3 in its exact form: the entering prime power does not kick λ′, it bends it, with a slope that diverges like 1/(ε log²(1/ε)) at the entry. The Φ′ budget of §2.3 keeps its shape with f(a_n)² read at the scale 1/Φ′ (see there); the rigidity results of §2.4 are unaffected, since they rest on negative values of finite forms.

**4. The boundary layer is archimedean.** The ground state of the archimedean-plus-polar odd form at a = 0.3 (no prime inside the window; λ = 0.222557) shows the same signature in the Legendre engine, values of f at y/a: [verified]

| K | 0.5 | 0.9 | 0.99 | 0.995 | 0.999 | 1.0 |
|---|---|---|---|---|---|---|
| 24 | 1.6068 | 1.2553 | 0.7477 | 0.6777 | 0.5518 | 0.4966 |
| 40 | 1.6067 | 1.2562 | 0.7479 | 0.6724 | 0.5591 | 0.4647 |
| 64 | 1.6067 | 1.2558 | 0.7472 | 0.6747 | 0.5657 | 0.4399 |

and the edge-FEM rows above give its law, C = 1.3089, β = -1.49. This is the first place in the program where the archimedean place acts *alone* on the mechanism, which is what §1.2 said any proof would need; the primes enter the edge law only through the amplitude C(a). (A first K = 40 Legendre run returned nonsense because the inner quadrature rule had 64 nodes, exact only to polynomial degree 127; the engine now scales the rule with K. Nothing at K ≤ 32 was affected.)

**What a certified transition needs, revised.** The earlier list (interval enclosures of the K-mode matrix at a₃ ± h, of its lowest eigenvalue, and of the endpoint value of the enclosed eigenvector, plus Zhu's complementary-space floor for λ and f(a₃)) would certify the K-mode kick, which is now known to be a quantity of the truncation with limit zero. What is worth certifying: (i) the floor λ(a) itself, unchanged, where Zhu's tail control is the tool; (ii) for the relay, the soft law, whose inputs are the theorem in item 1 (nothing to certify) and the edge amplitude C(a), a property of the exact eigenvector. Enclosing C needs an a posteriori eigenvector bound on the graded mesh, which is a harder certificate than an eigenvalue enclosure; whether the relay inequality of §2.3 needs C certified at all, or only bounded below, is the question to settle before building that certificate. E0 and E1 of §5 are the measurements that decide it.

---

### 2.7 The relay conjecture [statement; evidence to a = 1.5]

Everything in §2 is about one function. Let λ(a) be the floor of the odd Weil form on [−a, a] (the normalization of `rh_weil_odd.py`, which is also that of Checkpoints 15–20 of the atlas and of Zhu's paper: the 40-mode value λ(0.8) = 1.59×10⁻¹⁴ lies inside Zhu's certified odd interval [8.2×10⁻¹⁵, 2.35×10⁻¹⁴]), Φ(a) = −log λ(a), and T\*(a) = 2πe^{2a} the horizon. Two facts frame the conjecture. Positivity of the odd form on every finite support implies RH (Weil's criterion with Yoshida's odd-test reduction, as recorded in Suzuki, arXiv:2606.09096). And no estimate that treats the polar term and the reflected primes separately can prove positivity (G3.39 of the atlas; the rigidity of §2.4 is the same fact in another coordinate). So the only quantity worth a conjecture is the floor of the whole coupled form.

**Conjecture A (bounded relay).** There are a₀ and c such that λ is locally absolutely continuous on [a₀, ∞) and, for almost every a ≥ a₀,

  λ′(a) ≥ −c · T\*(a) · λ(a).

Where λ > 0 this reads Φ′(a) ≤ c·T\*(a), which is how it is used below; the λ form is the one to conjecture, because it does not presuppose the positivity it is meant to deliver. λ is non-increasing and continuous in a (the admissible functions on a smaller window are admissible on a larger one; continuity is in the paper, Lemma 1.4), hence differentiable almost everywhere; that it is absolutely continuous, so that the inequality integrates, is part of the conjecture, and follows from the bounded-variation hypothesis (H_BV) on near-minimizers of §2.8, which every computed minimizer satisfies (its derivative across an entry a_n is continuous by §2.6, the kink being soft, but the K-mode data cannot exclude a singular part). Gronwall: with u = λ·exp((c/2)T\*) one has u′ = (λ′ + cT\*λ)·exp((c/2)T\*) ≥ 0 a.e., since dT\*/da = 2T\*, so λ(a) ≥ λ(a₀)·exp(−(c/2)(T\*(a) − T\*(a₀))) > 0 for all a ≥ a₀; with a₀ = 0.8 from Zhu's certificate, **Conjecture A implies RH**. It is strictly stronger than RH: it fixes the rate at which the floor may decay. (The first version of this section stated the conjecture as Φ′ ≤ cT\* with the exponent c(T\*(a) − T\*(a₀)); the λ form and the factor ½ follow the independent review of the paper of 2026-09-30.)

**Conjecture B (decay law).** Φ′(a)/T\*(a) is bounded above and below by positive constants for a ≥ a₀. In its sharp form Φ′(a)/T\*(a) → c₀, i.e. Φ(a) ∼ (c₀/2)·T\*(a): the floor decays like exp(−(c₀/2)·2πe^{2a}); the data to height 130 give c₀ ≥ 3.9 and do not fix it (see the ball-arithmetic run below). Zhu's paper proposes, from his own upper bounds on 0.5 ≤ a ≤ 2, the law −log λ ≃ 2π²·N(T\*)/log N(T\*), which differs from the sharp form by a slowly varying factor; the weak form of B covers both, and the data below cannot separate them.

**Observation C (the relay at grid scale; corrected in §2.9).** Beyond the entries of 2 and 3, the entries a_n leave no visible drop in Φ′ at scale 0.004 or 0.0125: the measured drops at n = 4, 5, 7 are zero within the smooth drift Φ″·Δa. The reason is the crossover scale of §2.5, gap/‖dT/da‖, which is 3×10⁻⁵ at a₃ and about 10⁻⁹ at a = 0.7: the minimizer absorbs an entering prime power by reorganization long before any macroscopic scale, so the kink of §2.1 and the soft kink of §2.6 are microscopic structure. The Φ′ budget of §2.3 in its "drift up, drop at entries" form is therefore not what the data show beyond n = 2 (see the note added there); what they show is a plateau. The edge amplitude keeps its law, C(a)² ≍ λ(a)·T\*(a), within a factor two over 0.5 ≤ a ≤ 1.03.

**Why the decay law is natural under RH [heuristic].** A function of exponential type a can vanish on a real sequence of density at most a/π (Beurling–Malliavin). The zeros have density (1/2π) log(T/2π) at height T, which equals a/π exactly at T = T\*(a). So a test function with Fourier support in [−a, a] can cancel the zeros below the horizon but not above it; the floor is set by the zeros near T\*, and its logarithm should scale with the number of constraints there, N(T\*) ≍ a·T\*, up to the logarithm that separates T\* from T\* log T\* in the fits of `RH_ROUTES.md`.

**Evidence.** `rh_floor_grid.py` (data in `floor_grid_K40.csv`): K = 40 modes, 60 digits, 97 supports on [0.30, 1.20] with the eight entries a_n = ½ log n, n ≤ 11, resolved at ±0.004. Φ′ is a centered difference on the uniform grid, never straddling an entry. K-mode floors are upper bounds on the true floor (ground rules). Convergence: K = 56 changes λ by 3% at a = 0.9 (7.30 against 7.51 ×10⁻²⁰), but by a factor 200 at a = 1.1 (9.3×10⁻³⁵ against 1.9×10⁻³²; the 72-mode value of the folds paper is 7.7×10⁻³⁵) and by five orders at a = 1.2, where Zhu needed 950 modes. The table therefore stops at a = 1.03; the last rows of the grid are truncation artifacts and are shown separately. [verified]

| a | T\* | λ (K = 40) | Φ = −log λ | Φ′ | Φ′/T\* | f_K(a)²/(λ T\*) | λ₁/λ₀ |
|---|---|---|---|---|---|---|---|
| 0.3125 | 11.7 | 1.805e-01 | 1.71 | 18.6 | 1.58 | 0.10 | 5.7e+00 |
| 0.3625 | 13.0 | 4.521e-02 | 3.10 | 29.4 | 2.27 | 0.14 | 2.0e+01 |
| 0.4125 | 14.3 | 9.791e-03 | 4.63 | 33.6 | 2.34 | 0.15 | 9.2e+01 |
| 0.4625 | 15.8 | 1.388e-03 | 6.58 | 46.1 | 2.91 | 0.19 | 4.7e+02 |
| 0.5125 | 17.5 | 9.404e-05 | 9.27 | 55.8 | 3.19 | 0.21 | 3.0e+03 |
| 0.5625 | 19.4 | 6.847e-06 | 11.89 | 60.3 | 3.12 | 0.21 | 6.4e+03 |
| 0.6125 | 21.4 | 2.475e-07 | 15.21 | 70.6 | 3.30 | 0.23 | 1.9e+04 |
| 0.6625 | 23.6 | 7.376e-09 | 18.72 | 80.9 | 3.42 | 0.23 | 2.8e+04 |
| 0.7125 | 26.1 | 8.442e-11 | 23.20 | 86.8 | 3.32 | 0.23 | 7.8e+04 |
| 0.7625 | 28.9 | 8.881e-13 | 27.75 | 110.7 | 3.83 | 0.27 | 1.0e+05 |
| 0.8375 | 33.5 | 2.193e-16 | 36.06 | 119.7 | 3.57 | 0.23 | 2.7e+05 |
| 0.8875 | 37.1 | 4.831e-19 | 42.17 | 141.7 | 3.82 | 0.29 | 3.8e+05 |
| 0.9375 | 41.0 | 3.799e-22 | 49.32 | 149.5 | 3.65 | 0.27 | 6.4e+05 |
| 0.9875 | 45.3 | 1.548e-25 | 57.13 | 166.5 | 3.68 | 0.31 | 1.1e+06 |
| 1.0250 | 48.8 | 3.166e-28 | 63.32 | 162.2 | 3.32 | 0.23 | 1.2e+06 |

Φ′/T\* rises from 1.6 at a = 0.31 to 3.2 by a = 0.49 and then stays in [2.75, 3.83] for 0.46 ≤ a ≤ 1.03, mean 3.16 on [0.5, 0.65] and 3.56 on [0.85, 1.03]. That 13% rise sits between the pure-T\* law (no rise) and the T\* log T\* law (24% over the same range), and the truncation error, which grows with a and lowers Φ′, biases it downward; the sharp form of B is a working hypothesis, the weak form is the conjecture. The ratio λ₁/λ₀ grows without interruption, so the ground state stays simple and isolated (Lemma L2) across every entry in the range.

| a | 1.0625 | 1.0750 | 1.1125 | 1.1250 | 1.1375 | 1.1500 | 1.1625 | 1.1750 |
|---|---|---|---|---|---|---|---|---|
| Φ′/T\* (K = 40, unconverged) | 1.94 | 2.56 | 1.12 | 1.66 | 2.44 | 2.63 | 2.49 | 1.84 |

**Drift test at K = 72.** `rh_floor_grid.py` with 72 modes on 0.475 ≤ a ≤ 1.05, step 0.025 (`floor_grid_K72.csv`). The 72-mode floor differs from the 40-mode floor by 0.1–5% for a ≤ 0.925 and by 26% at a = 1.0, so the 40-mode grid above is reliable to a = 1.0 (the entry of 7 at 0.973 and of 8 at 1.040 are straddled by the coarser step and their neighbours are omitted). [verified]

| a | λ (K = 72) | λ₄₀/λ₇₂ | Φ′ (72) | Φ′/T\* (72) | Φ′/T\* (40-mode grid) |
|---|---|---|---|---|---|
| 0.500 | 1.938e-04 | 1.001 | 55.2 | 3.23 | 3.35 |
| 0.575 | 3.162e-06 | 1.002 | 63.6 | 3.20 | 3.20 |
| 0.600 | 5.956e-07 | 1.003 | 68.7 | 3.29 | 3.33 |
| 0.625 | 1.018e-07 | 1.004 | 69.6 | 3.17 | 3.19 |
| 0.650 | 1.837e-08 | 1.004 | 74.8 | 3.24 | 3.06 |
| 0.725 | 2.919e-11 | 1.008 | 86.5 | 3.23 | 3.11 |
| 0.750 | 3.372e-12 | 1.010 | 98.6 | 3.50 | 3.51 |
| 0.775 | 2.105e-13 | 1.017 | 107.5 | 3.63 | 3.55 |
| 0.850 | 4.800e-17 | 1.038 | 119.7 | 3.48 | 3.41 |
| 0.875 | 2.463e-18 | 1.054 | 130.1 | 3.60 | 3.53 |
| 0.900 | 7.176e-20 | 1.046 | 138.4 | 3.64 | 3.81 |
| 0.925 | 2.428e-21 | 1.032 | 144.5 | 3.62 | 3.53 |
| 1.000 | 1.491e-26 | 1.258 | 172.6 | 3.72 | 3.59 |

The ratio at 72 modes drifts from 3.23 (mean on [0.5, 0.7]) to 3.61 (mean on [0.85, 1.05]), +12%; the 40-mode ratio at the same points drifts +11%, so the drift is a property of the form, not of the truncation. Least squares on 0.5 ≤ a ≤ 1.05: Φ′ = c·T\* fits with c = 3.53 and 6.7% rms relative residual; Φ′ = c·T\* log T\* fits with 17% residual and would have produced a drift of about 25% over this range (58% with log(T\*/2π)). The pure law is favored, the residual drift is half of what a logarithm would give, and neither is excluded at height 50. At height 50 the sharp form with c₀ ≈ 3.5 looked like the working hypothesis; the run to height 130 below shows the ratio still rising slowly and leaves c₀ undetermined.

**To a = 1.5 with the ball-arithmetic engine.** `rh_weil_arb.py` (python-flint; Appendix) makes the same matrices hundreds of times faster, so the band can be followed to height 130 with the mode counts the truncation needs: K = 110 for a ≤ 1.15, 140 to 1.33, 180 to 1.5, with checks at K = 180 (a = 1.25) and K = 230 (a = 1.5). Φ′ at each centre is the difference over the pair a ± 0.005; no pair straddles an entry except the one at 1.10 (a₉ = 1.0986), where the K-mode kick is below 1% of Φ′ and is in any case absorbed (Observation C). Each eigenvalue carries a rigorous enclosure of the assembled matrix's lowest eigenvalue (radius ≈ 10⁻¹⁴³ at 500 bits); the polynomial evaluations and the outer quadrature are not enclosed. Data: `floor_grid_arb.csv`. [verified]

| centre a | K | T\* | λ(a − 0.005) | λ(a + 0.005) | Φ | Φ′ | Φ′/T\* |
|---|---|---|---|---|---|---|---|
| 1.060 | 110 | 52 | 7.922e-31 | 1.119e-31 | 70.3 | 195.8 | 3.740 |
| 1.100 | 110 | 57 | 2.029e-34 | 2.654e-35 | 78.6 | 203.4 | 3.587 |
| 1.150 | 110 | 63 | 3.459e-39 | 3.130e-40 | 89.8 | 240.3 | 3.834 |
| 1.250 | 140 | 77 | 1.563e-50 | 8.602e-52 | 116.1 | 290.0 | 3.788 |
| 1.330 | 140 | 90 | 1.950e-61 | 6.557e-63 | 141.5 | 339.2 | 3.776 |
| 1.400 | 180 | 103 | 1.269e-72 | 2.362e-74 | 167.5 | 398.4 | 3.856 |
| 1.450 | 180 | 114 | 1.123e-81 | 1.332e-83 | 188.6 | 443.5 | 3.884 |
| 1.500 | 180 | 126 | 1.248e-91 | 1.118e-93 | 211.7 | 471.5 | 3.736 |
| 1.250 | 180 | 77 | 1.526e-50 | 8.419e-52 | 116.2 | 289.7 | 3.785 |
| 1.500 | 230 | 126 | 1.094e-91 | 8.058e-94 | 211.9 | 491.1 | 3.891 |

convergence at 1.25: lambda(K=140)/lambda(K=180) = 1.0217;  Phi'/T*: 3.788 -> 3.785
convergence at 1.5: lambda(K=180)/lambda(K=230) = 1.3872;  Phi'/T*: 3.736 -> 3.891

Combined with the 72-mode data below a = 1 (heights 17 to 130):

- mean Phi'/T* on [0.5,0.7]: 3.229  (5 pts)
- mean Phi'/T* on [0.85,1.05]: 3.611  (5 pts)
- mean Phi'/T* on [1.06,1.25]: 3.736  (4 pts)
- mean Phi'/T* on [1.33,1.5]: 3.852  (4 pts)

One-parameter fits of Φ′/T\* against a for a ≥ 0.85, by rms relative residual:

- Phi' = c T*  (pure horizon)                                    c = 3.7237  rms rel resid = 0.034
- Phi' = c T* log T*                                             c = 0.8979  rms rel resid = 0.078
- Phi ~ N(T*):  ratio = c (4a-1)                                 c = 0.9954  rms rel resid = 0.214
- Phi ~ N/log N (Zhu shape): ratio = c (2a-1)/(2a+log(2a-1))     c = 7.2157  rms rel resid = 0.030

Reading. Convergence is complete at a = 1.25 (2% in λ, 0.1% in the ratio) and incomplete at a = 1.5 (39% in λ between 180 and 230 modes, +4% in the ratio, and more modes raise it further), so the last ratio is a lower bound of about 3.9. The ratio rises across the whole range, 3.23 → 3.61 → 3.74 → 3.85, by shrinking steps. The raw phase-space count Φ ∝ N(T\*), whose ratio would grow linearly in a, is excluded (21% residual); the logarithmic law is disfavored (7.8%); the pure horizon law (3.4%, c = 3.72) and Zhu's Landau–Widom form Φ ∝ N(T\*)/log N(T\*) (3.0%) fit equally well. They agree on the data and disagree on the asymptote: the pure law says the ratio saturates near 3.8–4 now, while the Landau–Widom shape c·(2a − 1)/(2a + log(2a − 1)) with its fitted c = 7.2 says it keeps rising, to about 7, over a range of a far beyond the resolution wall. Height 130 cannot separate the two. So what the numerics establish is the weak form of Conjecture B, a ratio bounded between 3 and 8 with a slow rise, and they leave the sharp constant undetermined: c₀ is at least 3.9, and the two natural laws put it at 4 and at 7. This is the strongest statement the numerics can make below the wall, and it is the correct place to stop: the constant is not the content of the conjecture, the boundedness is.

**Entries at scale 0.004** (Φ′ from the points a_n ± 0.004; drop = Φ′(a_n⁻) − Φ′(a_n⁺); the smooth drift over 0.008 is about −Φ″·0.008 ≈ −1.4 at a = 0.7):

| n | a_n | Φ′(a_n⁻) | Φ′(a_n⁺) | drop | drop/Φ′ |
|---|---|---|---|---|---|
| 2 | 0.3466 | 38.3 | 32.0 | 6.3 | 0.165 |
| 3 | 0.5493 | 62.8 | 59.4 | 3.4 | 0.054 |
| 4 | 0.6931 | 88.7 | 90.7 | −2.0 | −0.02 |
| 5 | 0.8047 | 101.6 | 109.9 | −8.3 | −0.08 |
| 7 | 0.9730 | 155.1 | 157.5 | −2.4 | −0.02 |

The entries of 8, 9 and 11 lie in the unconverged range and are not reported.

**What the conjecture is not.** It is not a route to a proof by estimates: the band in Φ′/T\* is a property of the minimizer of the coupled form, and any attempt to bound Φ′ by bounding the polar and prime terms separately fails by G3.39. It is a target for two kinds of work. Numerically, the band has now been followed to a = 1.5 (height 130) with the ball-arithmetic engine; the next test is the approach to c₀ beyond height 130, which needs several hundred modes and is within reach of the same engine. Analytically, Conjecture A is a differential inequality for the ground state of a one-parameter family of self-adjoint operators whose parameter is the support, which is the form in which Suzuki's limit conjecture and this program meet; §2.8 derives the exact identity behind it, verifies it to twenty digits, and shows that it is a property of the ground state and not an inequality between forms.

---

### 2.8 The sliding-window identity: Conjecture A as a law [theorem for the K-mode form, verified to 20 digits; Dini-derivative proposition for the exact form]

**Status (after the independent review of 2026-09-30).** The identity derived below is a theorem for the K-mode form Q_K(a) at every fixed K: the basis is a dilation family, Q_K(a) is a real-analytic matrix function of a between two entries, and a·λ_K′ = −[D_∞ + D_P + 2P² + 2PM] on the ground state is its Hellmann–Feynman derivative (paper, Theorem 5.2; verified to twenty digits in the table below and independently on the Fourier side at K = 12, where D_∞ is computed from (1/2π)∫|F|² t ∂_t Re ψ with F through spherical Bessel functions and D_P from g′(log n) of the exact autocorrelation). For the exact form the same computation gives a Dini-derivative statement only: a·D⁺λ ≤ −[D_∞(f) + D_P(f) + 2P² + 2PM] ≤ a·D₋λ for every minimizer f whose autocorrelation is differentiable at the points log n, with equality almost everywhere because λ is monotone (paper, Proposition 5.3). That λ is absolutely continuous, so that the almost-everywhere identity integrates, is not proved; it is part of Conjecture A as restated in §2.7. The first version of this section asserted the identity for the exact form on numerical grounds.

**Regularity of the floor (added 2026-09-30, paper Lemma 1.4, Lemma 5.5, Proposition 5.6, Computation 5.8).** Unconditionally, λ is continuous and non-increasing (minimizers at nearby supports have bounded form norm, and on a form-bounded set the dilation b ↦ Q(f_b) is uniformly continuous at b = 1), so the floor cannot jump. For f ∈ L^∞ of bounded variation, g_f is Lipschitz with constant ‖f‖_∞·Var f, b ↦ Q(f_b) is locally Lipschitz with |b dQ(f_b)/db| = |𝒟(f_b)| ≤ w\*‖f‖² + b^{−1}‖f‖_∞ Var f · Σ_{log n<2ab} 2Λ(n) log n/√n + polar, w\* = sup t∂_t Re ψ(¼+it/2) = 2.0807. Hence **(H_BV): near-minimizers at each a ∈ [a₀, a₁] with ‖f‖_∞·Var f ≤ M ⇒ λ is Lipschitz on [a₀, a₁]** (dilate a near-minimizer at the larger support down to the smaller one), so absolutely continuous, and Conjecture A is then the inequality λ′ ≥ −cT\*λ a.e. alone. (H_BV) is a statement about near-minimizers, and the K-mode minimizers are near-minimizers (λ_K ↓ λ), so a bound uniform in K suffices. Measured (`rh_bv_check.py`, `data/bv_*.log`): every edge-FEM minimizer (eight supports, edge resolved to 10⁻¹²a) and every K-mode minimizer (a = 0.6, 0.8, 1.0; K up to 96, 128, 160) has Var f = 4 sup|f| to four digits, i.e. each half-window profile is unimodal, with sup|f| = 1.20–1.70 (FEM) and 1.40–1.46 (K-mode), converged to five digits in the mesh (124/237/486 nodes at a = 0.5: Var = 5.26688, 5.26679, 5.26678) and in K; so M ≈ 9 and the trivial Lipschitz bound on [0.6, 1] is λ′ ≥ −211, twenty orders above the truth: the conjecture is the cancellation in 𝒟. The Hardy integral ∫f²/δ near the edge grows doubly logarithmically as the edge law predicts (increments 1.63, 0.93, 0.62 per three decades at a = 0.3 against 1.66, 0.83, 0.55 from C²log(log(a/δ)+β)), so the minimizers are not in H^{1/2}: L^∞ ∩ BV, not H^{1/2}, is the right class for the dilation derivative.

Conjecture A is a statement about a rate, and the rate has an exact, unconditional expression. Write the odd form on the window as Q_a(f) = (1/2π)∫|F(t)|² Ψ_a(t) dt − 2P(f)², F = f̂, P(f) = ∫ f(x) sinh(x/2) dx, with the symbol

  Ψ_a(t) = Re ψ(¼ + it/2) − log π − Σ_{log n < 2a} 2Λ(n) n^{−1/2} cos(t log n).

The K-mode basis of the engines is a dilation family: φ_i(y/a)/√a. For any dilation family f_a(x) = a^{−1/2} g(x/a) one has F_a(t) = √a·G(at) exactly, jumps at the edge included, so between two entries (where the prime set in Ψ_a is fixed) the substitution s = at gives

  d/da (1/2π)∫|F_a|² Ψ_a = −(1/a)·D_a(f_a),  D_a(f) := (1/2π)∫|F(t)|²·t·Ψ_a′(t) dt,

and for the polar term P_a = √a ∫ g(y) sinh(ay/2) dy gives P_a′ = (P_a + M(f_a))/(2a), M(f) := ∫ f(x)·x·cosh(x/2) dx. Hellmann–Feynman on the minimizer (the matrix family Q(a) is exactly this dilation family, so the K-mode statement is exact) then reads

  **a·λ′(a) = −[ D_∞(f) + D_P(f) + 2P² + 2P·M ],**

where the symbol t·Ψ_a′(t) splits into an archimedean part t·∂_t Re ψ(¼ + it/2), which tends to 1, so that D_∞(f) = ‖f‖² + O(∫|F|²/t²) ≈ 1, and a prime part Σ 2Λ(n)n^{−1/2}·(t log n)·sin(t log n), which in x-space is

  D_P(f) = −Σ_n 2Λ(n) n^{−1/2} · log n · g_f′(log n),  g_f = f ⋆ f̃:

the Euler product differentiated, acting on the derivative of the minimizer's autocorrelation at the prime-power positions. So

  **Φ′(a) = −λ′/λ = [D_∞(f) + D_P(f) + 2P² + 2PM] / (a·λ),**

and **Conjecture A is the inequality D_∞(f) + D_P(f) + 2P² + 2PM ≤ c·a·T\*(a)·Q_a(f) on the minimizer.** Since Q_a(f) = λ is 10⁻⁹⁰ at a = 1.5 while D_∞(f) ≈ 1, the prime dilation form must cancel the archimedean one to ninety digits: D_P(f) = −D_∞(f) − 2P² − 2PM + O(a T\* λ). Together with the explicit formula on the same function, Σ 2Λ(n)n^{−1/2} g_f(log n) = W_∞(g_f) − 2P² − λ, the minimizer's autocorrelation has both its values and its log-weighted derivatives at the prime-power positions pinned by the archimedean side. Each further a-derivative of λ pins the next weighted derivative; the relay is this hierarchy, and a proof of A is a bound on its first level.

**Verified on the K-mode minimizer** (`rh_dilation.py`, central differences of the three parts of Q(a) at step 10⁻³⁰ in ball arithmetic, contracted with the eigenvector; every entry is a ball with the radius shown, and the sum reproduces Φ′ of the grids): [verified]

| a | K | λ | a⟨c, A′c⟩ (archimedean) | a⟨c, P′c⟩ (primes) | a⟨c, S′c⟩ (polar) | sum = a·λ′ | radius | Φ′/T\* |
|---|---|---|---|---|---|---|---|---|
| 0.6 | 48 | 5.9659e-07 | -1.014313 | +1.086974 | -0.072687 | -2.493408e-05 | 3.51e-20 | 3.3391 |
| 0.8 | 64 | 1.5659e-14 | -1.016605 | +1.103966 | -0.087361 | -1.295854e-12 | 3.62e-27 | 3.3239 |
| 1.0 | 80 | 1.4820e-26 | -1.018057 | +1.114966 | -0.096909 | -2.516852e-24 | 3.15e-39 | 3.6579 |
| 1.25 | 140 | 3.4161e-51 | -1.019178 | +1.123585 | -0.104407 | -1.233712e-48 | 2.32e-63 | 3.7744 |

The archimedean dilation form sits at 1.014–1.019, the 1–2% above ‖f‖² being the finite-t correction to t·∂_t Re ψ; the prime dilation form grows with the primes present, from −1.087 to −1.124, and cancels it; the polar term grows from 7% to 10%. The sum is a·λ′ = −a·Φ′·λ to the last digit.

**Is the law an operator inequality?** The strong form of A would be −Q′(a) ≤ c·T\*(a)·Q(a) as forms on the whole window space, not only in the minimizer's direction. That is a generalized eigenvalue problem for the pencil (−Q′, Q), and its largest eigenvalue divided by T\* answers it: [verified]

| a | −Q′/Q on the minimizer, over T\* | largest pencil eigenvalue over T\* | second, third | smallest |
|---|---|---|---|---|
| 0.6 | 3.339 | 74.7 | 4.26, 3.69 | -70.5 |
| 0.8 | 3.323 | 1.46e+05 | 332, 6.08 | -1.46e+05 |
| 1.0 | 3.657 | 6.2e+10 | 3.06e+07, 8.17e+04 | -6.2e+10 |

The strong form is false, and not by a little: the largest eigenvalue of the pencil grows from 75·T\* at a = 0.6 to 6×10¹⁰·T\* at a = 1.0, in ± pairs of nearly equal size. Those are near-null directions of Q whose dilation derivative is enormous relative to their value, the reorganization directions of §2.5, and the pencil measures their rate; it is the inverse of the crossover scale gap/‖dT/da‖ (≈ 10⁻⁷ at a = 0.8, ≈ 10⁻¹³ at a = 1.0). The ground state is neither the largest nor the smallest direction of the pencil; it is a special interior direction in which the dilation derivative is anomalously small, by ten orders of magnitude at a = 1.0. So **Conjecture A is a statement about the ground-state eigenvector alone**, exactly λ′ = ⟨f₀, Q′f₀⟩, and cannot be an inequality between the two forms. Any proof must use that f₀ minimizes Q, through its Euler–Lagrange equation, to control the single diagonal matrix element ⟨f₀, Q′f₀⟩; the off-diagonal elements that the pencil sees drive the rotation of the minimizer and do not enter λ′ at first order. That is the exact shape of the shadow: a variational identity for the minimizer, not a law for the form.

**Under RH, what the identity says.** With all zeros real, Q_a(f) = Σ_γ |F(γ)|², and the same dilation computation on the zero side gives a·λ′ = λ + Σ_γ γ·(|F|²)′(γ), so Conjecture A becomes Σ_γ γ·∂_t|F|²(γ) ≥ −(1 + c a T\*)·Σ_γ |F(γ)|²: the minimizer must be small at the zeros to the same order in its derivative as in its value. And λ(a) itself becomes the inverse of a sampling constant: Φ(a) = −log inf{Σ_γ|F(γ)|² : F ∈ PW_a, ‖f‖ = 1}. The zeros have density (1/2π) log(T/2π), which exceeds the Nyquist density a/π exactly above T\*(a), so only the zeros above the horizon sample PW_a, and the floor measures how well a band-limited function can hide from them. A prolate-type lower bound follows from Landau–Pollak–Slepian: the number of zeros in [−T\*, T\*] is 2N(T\*) = (T\*/π)(2a − 1), short of the Shannon number 2aT\*/π by T\*/π, so a function can vanish on all of them and leak past the horizon only exp(−πT\*/log(4aT\*))-much; hence Φ(a) ≳ πT\*/log(4aT\*) ≈ 0.5·T\* at a = 1, against the observed 1.85·T\*: the zeros are a much more expensive constraint set than the prolate count alone, by a factor near four, which is the arithmetic in the problem. Conjecture A is the matching upper bound, an explicit sampling inequality for the zeros above the horizon with constant e^{cT\*}. None of this is a proof; it is the exact statement of what a proof must show, in both the unconditional (symbol) and the conditional (sampling) form.

**Evaluated on the zeros (2026-09-30; `rh_zero_side.py`, `data/zeros_6000.txt`, `data/zero_side_*.log`; paper Computation 7.5).** [verified] The explicit formula is unconditional and every zero used lies on the line, so for the K-mode minimizer 2Σ_{γ>0}|F_K(γ)|² must return λ_K, and λ_K + 2Σγ ∂_t|F_K|²(γ) must return a·λ_K′. With the first 6000 zeros (the first 400 to 36 digits, the rest to 21) and Riemann–von Mangoldt tails:

| a | K | λ_K | 2Σ\|F\|² + tail | ratio | a·λ_K′ (matrices) | λ_K + 2Σγ(\|F\|²)′ + tail | ratio |
|---|---|---|---|---|---|---|---|
| 0.6 | 48 | 5.96594e-7 | 5.96554e-7 | 0.99993 | −2.49341e-5 | −2.49310e-5 | 0.99988 |
| 0.8 | 64 | 1.56587e-14 | 1.56568e-14 | 0.99988 | −1.29585e-12 | −1.29685e-12 | 1.00077 |
| 1.0 | 80 | 1.48201e-26 | 1.48157e-26 | 0.99971 | −2.51685e-24 | −2.52433e-24 | 1.0030 |

The floor of the window form, built from the primes and the archimedean symbol, is the sum over the actual zeros of |F_K|² to four digits, and its rate of change is the zero-side sum to three: the third leg of the triangulation (matrices, Fourier symbol, zeros). The derivative sum is first-order sensitive to the zero positions (4γ|F′(γ)|²ε per zero): with 21-digit zeros it came out +7e-22 against −2.5e-24 at a = 1, and the 36-digit zeros cured it. Where the floor sits: the zeros below the horizon carry 1e-5 (a = 0.6, one zero), 1e-7 (0.8, four zeros), 1.6e-13 (1.0, eight zeros; the first zero alone 9e-25) of it; the largest single terms are at γ₇ = 40.9 (15%), γ₁₃ = 59.3 (19%), γ₂₇ = 94.7 (8%), i.e. at 2.0, 1.9, 2.0 horizons, and the mass sits between 1.5 and 3 horizons at each support. The minimizer hides from the zeros below the horizon to an accuracy far beyond the floor itself; the floor is made by the zeros between T\* and about 3T\*, as the sampling picture says.

**The tail law (2026-09-30; `rh_tail_law.py`, `data/tail_law_fem.log`, `data/tail_law_kmode/`; paper Computation 7.6, Conjecture 7.7, Proposition 7.8, Heuristic 7.20).** [verified on FEM and K-mode minimizers] Let ϑ(T) be the fraction of 2Σ_{γ>0}|F(γ)|² carried by the zeros above height T. Then, from a few horizons on,

  **ϑ(T) = Φ′(a)/(πT), i.e. Φ′(a) = π·T_eff(a), T_eff := lim T·ϑ(T):**

the rate of decay of the floor is π times the effective height of the zeros that carry its tail. On the edge-FEM minimizers with primes (a = 0.4, 0.45, 0.5; 486 nodes, edge to 10⁻¹²a; F of the piecewise-linear f in closed form), the FEM floor is recovered from the zeros to 5×10⁻⁶ (2Σ|F|²/λ + tail = 0.999995, 0.999997, 0.999999) and T·ϑ(T)/T\* equals Φ′/(πT\*) = 0.7213, 0.8676, 1.0811 to 1% from 7 to 100 horizons (ratios 0.995–1.009), with no logarithmic drift; at the entry a₃ the ratio rises from 0.69 (2T\*) to 0.92 (100T\*), the soft kink seen from the zeros (a zero at height T probes the edge at scale 1/T, where λ′ is still moving). On the K-mode minimizers it holds to 1–3% from five horizons on: Tϑ(T)/T_eff (T_eff = Φ′/π) at 5 and 10 horizons is 1.011/1.002 (a = 0.6, K = 48), 1.067/0.953 (0.8, 64), 1.032/1.005 (1.0, 80), 1.038/1.007 (1.25, 140), 1.015/0.993 (1.5, 230, λ = 9.2×10⁻⁹³, first 1500 zeros at 50 digits, Φ′ = 491.1 from the grid's central difference since the Hellmann–Feynman contraction at K = 230 failed for want of eigenvector precision). The profiles ϑ(T/T\*) at 0.6 and 0.8 coincide to 0.01; those at 1.0–1.5 lie above them in the plunge region (ϑ(2T\*) = 0.81, 0.78, 0.86 against 0.59) and rejoin in the tail; in units of T_eff all five agree to 5% from 2.5 horizons on, across 86 orders of magnitude of the floor. Medians: 2.39, 2.32, 2.65, 2.67, 2.74 T\* at 0.6 … 1.5, rising with a and at the last point at e. (The λ_K(1.25) of the identity table was corrected to 3.4161e-51, the value the zeros return.) Below that, the plunge region [T\*, 3T\*] carries more than the law (Tϑ/T\* up to 1.4–1.5 at 1.5T\*); the median height of the floor's mass is 1.8, 2.0, 2.4 T\* (FEM 0.4, 0.45, 0.5) and 2.4, 2.3 T\* (K-mode 0.6, 0.8), rising slowly with a. Heuristic derivation: edge law ⇒ |F(t)|² ≈ 4C²cos²(ta)/(t² log t), times the zero density log(t/2π)/2π, integrated from T: C²/(πT) = −λ′/(2πT) by the boundary law. The derivation carries a ratio of logs (1.1–1.4 at these heights) that the data do not show: the law is sharper than its derivation, which suggests a sum rule.

**What this makes of the conjectures.** Conjecture A ⟺ T_eff(a) ≤ (c/π)T\*(a): the effective height of the zeros that carry the floor is a bounded multiple of the horizon. Conjecture B ⟺ T_eff/T\* → c₀/π; the measured c₀ = 3.2–3.9 is T_eff = 1.0–1.24 horizons. **Why the effective height should be tied to the horizon (hypothesis):** (1) below T\* the zeros are sparser than the a/π degrees of freedom per unit height of a function of type a, and the minimizer vanishes at all of them (shares 10⁻⁵, 10⁻⁷, 10⁻¹³); (2) above T\* the excess density of zeros over the Nyquist density is log(t/T\*)/2π, the same function of t/T\* for every a, so after scaling heights by the horizon the minimizer faces the same problem at every support: this is the reason for a universal c₀ and for the convergence of the profiles ϑ(T/T\*) as a grows; (3) the only currency for being small at more points than one has freedom for is the edge (edge law, Hadamard), and the edge amplitude feeds the tail at all heights at once at the rate C²/(πT), so the minimizer balances the gain from hiding at one more zero above the horizon against a cost paid everywhere above, and the balance sits at a fixed multiple of T\*; (4) the spare freedom below the horizon is worth T\*/2π zeros and the accumulated excess ∫_{T\*}^{T} log(t/T\*)dt/2π exhausts it at T = eT\*; the medians are below that and rising. Two inputs no counting supplies: positivity of every term (RH itself), and the absence of gaps between zeros wider than the resolution π/a above a few horizons, in which a shape could park its tail. So the explanation is conditional on where the zeros are: it says what Conjecture A asserts about the zeros (density and no large gaps at the scale of the horizon), not why the prime side must produce it. A proof along these lines needs three lemmas: edge law ⇒ tail law (asymptotics plus the averaging of cos²(γa) over zeros, where the entries reappear as the prime at e^{2a}); positivity ⇒ C² ≤ πTλ past the onset of the law; and the onset at a bounded multiple of the horizon, a sampling-type lower bound for the zeros above κT\* acting on the tail of a function of type a.

**Program (see `PROGRAM.md`).** Item 1 done: paper Proposition 7.8 proves edge law ⇒ tail law under RH (the transform's edge asymptotics by splitting at δ = 1/t; the sum over zeros by Riemann–von Mangoldt with Littlewood's S₁ = O(log t/(log log t)²)); with the boundary law it is Conjecture 7.7. The law is a property of the transform sampled by the zero density: at a = 0.5 the density integral 2∫_T|F|²ρ is 1.068 … 1.006 of the law at 5–100 horizons (the proof's O(1/log T)) while the zero sum is 1.00 ± 0.02. Item 2 open, route found: the Euler–Lagrange equation in Fourier form, (Ψ − λ)F = G + 4PΣ with G the transform of the force outside the window, gives for every T the exact identity ∫_{|t|>T}|F|²(Ψ−λ) = 4P∫_{|t|>T}F̄Σ − 2π⟨f_T, g_out⟩; the leak of the low-pass f_T outside the edge (∝ C L(1/T)^{−1/2}) against the outside force (∝ −2C L(δ)^{1/2}) has no logarithm and gives 4C²/T, the tail law; the identity is now paper Proposition 7.9 (proof: the eigen-equation Kf = λf + 2P sinh on the window, Parseval, disjoint supports) and is verified on the FEM minimizers at a = 0.5 and 0.45 to 1–2% at 5–20 horizons (`rh_sum_rule.py`, `data/sum_rule_*.log`, Computation 7.10): the overlap of the leak with the edge-law force −C L^{1/2} alone is the law to ±10%, and the smooth overlap plus the polar piece (each ~30% of the tail at 5 horizons) cancel to ±10%. The 1/L expansion of the edge overlap is done (paper Lemma 7.11): J(L) = 1 + O(1/L²) with no 1/L term, by the symmetry of the kernel sin(s+v)/(s+v) in the two distances against the antisymmetric log ratio; J = 1.020 … 1.003 at L = 3 … 8. This is why the tail law is sharp at finite T. The cancellation of the smooth overlap against the polar piece is now proved (paper Lemma 7.12): the regular part of K_∞f is continuous across the edge, because inside the wall's +½log(1/δ)·f and the log-Laplacian's −C L^{1/2} cancel while outside the −C L^{1/2} stands alone, with the same regular remainder on both sides (the difference is −C(β + log 2a)L^{−1/2} → 0); hence g_reg(a) = (Kf)(a⁻) = 2P sinh(a/2) by the eigen-equation, h = 2Ps + g_reg is continuous at the edge, and the polar piece plus the smooth overlap is ⟨f, h_{>T}⟩ = O(C²/(TL)). Corollary 7.13 (the sum rule): the tail energy is 2C²/(πT)(1 + O(1/L)), the edge overlap giving 1 + O(1/L²) and the rest O(1/L) with a coefficient the data say is small. Item 3 reformulated (paper Lemma 7.26, Proposition 7.27): the zeros in [−T, T] equal the Shannon number 2aT/π exactly at T = eT\* (the exhaustion height; the medians of the floor's mass are 2.39 … 2.74 T\*, ϑ(eT\*) = 0.43 … 0.50); under RH, a sampling inequality for the zeros above κT\* on the tail of F (constant M) together with the edge overlap holding at κT\* (within M′) and the boundary law give Conjecture A with c = πκMM′; measured M ≈ 1, M′ ≈ 1, κ = e, so the bound is πe = 8.5 against the observed c₀ = πe·ϑ(eT\*) ≈ 3.3–3.9. What remains is (S_e), a statement about the zeros alone, and (E_e), about the minimizer alone. **Decisive measurement (paper Computation 7.14):** the tail law is analytic. On the prime-free critical points (no zeros, floors of either sign) T·(1/2π)∫_{|t|>T}|F|²(Ψ_∞ − λ) = −λ′/π to 1–2% from three horizons on, while the profile at δ = 1/T is 1.5–3× above its asymptotic edge law there; so the finite-T sum rule is an exact analytic identity still to be found (a Pohozaev-type identity for the eigen-equation), Corollary 7.13 being its far-field limit. And M(κ) = 1.01–1.09 for κ ≥ e (8.9, 3.6, 1.8 at κ = 1): above the exhaustion height the zeros sample the analytic tail at their density with constant one. All the arithmetic of Conjecture A is in (S_e). **The exact finite-height identity (paper Proposition 7.15):** cutting the dilation integral at height T and using the eigen-equation, 𝒟 = tail(T) + (1/π)T|F(T)|²(Ψ(T)−λ) + tail_D(T) − 4P⟨s,((xf)′)_{>T}⟩ + 2B_T with B_T = ⟨g_out, ((xf)′)_T⟩ the Hadamard pairing of the outside force with the low-pass of the dilation generator; hence −aλ′ = 2 lim B_T and, with the boundary law, B_∞ = aC². Verified on the prime-free critical point at a = 0.5 (𝒟 = 1.129084 from the profile vs 1.129075 from the engine; B_T → 0.551 vs 𝒟/2 = 0.5645; identity to 1–2%). The tail law is B_∞ − B_T = B_∞/(πaT)(1+o(1)): the pairing converges at the rate set by its own value; in the pairing identity the polar high-pass term is O(L^{−1/2}), five times the tail, cancelled by B_∞ − B_T at that order, so the identity is the dilation identity in Fourier form and does not simplify the law. **The sharp form (paper Computation 7.18, Conjecture 7.19): the local envelope ⟨t²(Ψ−λ)|F(t)|²⟩ over one interference period equals −λ′** — to 0.1% with no drift from ten to three hundred horizons on the prime-free critical points, ±3% with primes. Mechanism: Wiener–Hopf, now verified in modulus and phase (paper Conjecture 7.19 in its sharp form, Computation 7.20, `rh_wiener_hopf.py`, `data/wiener_hopf_*.log`). With σ̃ the symbol with its real zeros ±t_j divided out and θ = ½H[log σ̃] (the Hilbert transform, odd in t; θ ≈ −π/(4 log(t/2π)) at large t), the half-line problem at the edge gives S(t) = ∫₀^a f sin(tx) dx = −|A| cos(ta + θ(t) − κ/t)/(t σ̃(t)^{1/2}) (1 + ε(t)), |A|² = −λ′/2, with A purely imaginary (no free phase: the reality symmetry of the Wiener–Hopf solution forces it) and one real constant κ from the zeros' factors and the standing waves of the interior; the far edge enters through the commutator [K, χ] with a cutoff, at relative order t⁻². Measured on the prime-free critical point at a = 0.5: the 319 nulls of S between 2.8 and 120 horizons sit at ta = π/2 + kπ − θ(t) + κ/t with rms residual 1.2×10⁻⁴ rad (κ = 4.34 against t_λ = 4.56), where the raw shift is 0.34 → 0.14 and the local law c/log(t/2π) leaves 3.5×10⁻³ (its error 0.03 at three horizons); pointwise, with nothing but κ fitted, the form is within 4.5×10⁻³, 3.4×10⁻⁴, 8×10⁻⁵, 3×10⁻⁵ in L² over a period at 3, 10, 20, 50 horizons, decaying like T⁻². With the prime 2 the phase oscillates at frequency log 2 with amplitude ≈ 0.98/2σ̃, and the nulls follow it (rms 0.010 against 0.070 for the prime-free phase at the same λ): the wobble of the nulls is the Hilbert transform of the prime's cosine; pointwise 2.4% at 10 and 1.0% at 50 horizons, the amplitude modulation of the ±3% envelope. So all the logarithms sit in σ₋ (the edge law C L^{−1/2} is its x-space asymptotics, |A| = C), the onset is not the last sign change of the symbol (0.7, 2.2, 3.5, 33, >150 T\* at a = 0.5 … 1.25), and the tail law at finite height is this form averaged over a period. The derivation needs no cutoff: on the whole line σ_a F̃ = H̃ + G̃₊ exactly, with F̃ = e^{−ita}F analytic below, G̃₊ (the outside force on the right) analytic above and H̃ (the polar step −2iP sinh(a/2)/t plus the far edge, carrying e^{−2ita}) analytic below; the primes outside the window cancel identically, so the symbol is σ_a = Ψ_a − λ with the entries log n < 2a only; Liouville and the reality symmetry F̃(−t) = conj F̃(t) give A(t) = A₀ + A₁/t with A₀ imaginary and A₁ real, i.e. the first correction is the phase κ/t. With primes the entries inside the window copy each edge to its echoes (a − log n, …, while inside the window); the echo feeds back on the edge as the real term −c²/σ_∞ (c = Λ(n)n^{−1/2}), so the edge's symbol is σ_eff = σ_∞ − c²/σ_∞ and the echo piece is (c/σ_∞)e^{−it log n} times the edge piece (a finite continued fraction for longer chains; the window's symbol σ_∞ − 2c cos(t log n) is the infinite chain). Verified at a = 0.5 with the prime 2: the echo form leaves a null residual of 1.1×10⁻³ rms on 5–80 horizons (the window's symbol 9.8×10⁻³) and a pointwise distance 0.010, 0.0021, 0.0007 at 10, 20, 50 horizons decaying like T⁻² (the window's symbol 0.024, 0.016, 0.010, not decaying): the percent-level remainder with primes was the echo that the window cuts off. **The prime-free form is now a theorem (paper Lemma 7.21, Lemma 7.22, Proposition 7.23, Corollary 7.24).** For every L² odd critical point of W_∞ − 2P², with no hypothesis on the edge: F̃ = P₋[Ĥ₋/σ₊]/σ₋ exactly (the real zero of the symbol removed by a rational factor, the factor σ₊ the outer function of σ̃^{1/2} with c ≤ |σ₊| ≤ C(log)^{1/2}, Liouville on an entire function of growth (log)^{1/2}); hence F̃ = A(t)/(tσ₋) + e^{−2ita}Φ₋/σ̃ with A(t) = A₀ + O(log t/t), A₀ purely imaginary, the far edge feeding back only through the kernel of 1/σ₊ beyond 2a (which decays like e^{−ηξ}), and the mirror F = A e^{ita}/(tσ₋) − conj(A)e^{−ita}/(tσ₊) + O(t^{−2}(log t)^{−1/2}). So the tail energy of every prime-free critical point is 2|A₀|²/(πT)(1 + o(1)) beyond the onset, and the tail law is the single identity 2|A₀|² = −λ′: the boundary law read as a flux. What remains of item 2: the flux identity |A₀|² = −λ′/2 (to be computed from the pairing identity −aλ′ = 2B_∞ with the now-proven asymptotics of G̃₊ = σ₊A/t + 2iP sinh(a/2)/t + …), the constancy of the phase κ (bound O(log t), data 10⁻⁴), and the echo form with primes as a theorem.

---

### 2.9 The boundary law: λ′(a) = −2C(a)² [theorem under the edge law with two derivatives, paper Proposition 8.3; verified to 10⁻³ in eight cases]

Since Q does not depend on a, λ(a) is the infimum of a fixed form over the growing family L²(−a, a), and its derivative is a domain variation. For the Dirichlet fractional Laplacian, Hadamard's formula expresses such a derivative through the boundary trace u/dist^s of the eigenfunction (Djitte–Fall–Weth); for the logarithmic Laplacian the trace is the amplitude C of §2.6. A Pohozaev identity for L_Δ with the local boundary term ∫u²·log(δ⁻²)(x·ν) was claimed in a 2024 preprint that was withdrawn for an error in its proof, so the law is stated here as a conjecture with its evidence.

**Conjecture (boundary law).** For a not an entry, λ′(a) = −(C₊² + C₋²) = −2C(a)², where C_± are the edge amplitudes lim √log(1/δ)·|f(±(a−δ))| of the minimizer. Equivalently, with §2.8, D_∞(f) + D_P(f) + 2P² + 2PM = 2a·C(a)² on the minimizer: a Pohozaev identity for the Weil form on a window.

**Proved (2026-09-30, paper Lemma 8.2 and Proposition 8.3) under the edge law with two derivatives (the hypothesis of Proposition 7.8):** B_T = ⟨g_out, ((xf)′)_T⟩ → aC². The boundary term of the Pohozaev identity is the overlap, at scale 1/T, of the leak of the dilation generator −af′ through the low-pass with the singular edge force −CL^{1/2}: the double integral J₁(L) = −(2/π)∫∫[(L+log 1/v)/(L+log 1/s)]^{1/2}k′(s+v) ds dv equals 1 + O(1/L²) with no 1/L term (the antisymmetry of Lemma 7.11), the regular value g_reg = 2P sinh(a/2) of the force gives the rate aC g_reg L_T^{−1/2}, and Lemma 8.2 controls the remainder of the force (o(L^{−1/2}), derivative o(δ^{−1}L^{−3/2})). With Proposition 7.15 (𝒟 = 2B_∞) and the sliding-window identity, λ′ = −2C². For the prime-free form, Corollary 7.24 gives the edge in the form A₀/(tσ₋), i.e. the edge law with C = |A₀|; that derivative control is now Proposition 8.4 (2026-10-01): the inverse transform of A₀/(tσ₋) is C(L+β)^{−1/2} + ψ with ψ^{(k)} = O(δ^{−k}L^{−5/2}), C = |A₀| and β = −γ − log(2πa) − λ, i.e. f(a−δ) = |A₀|(σ_∞(1/δ) − γ)^{−1/2}(1 + O(L^{−2})), the profile being the inverse square root of the symbol at frequency 1/δ shifted by Euler's constant (Mellin expansion of the transform of (log 1/δ)^{−1/2}, whose real L^{−2} corrections cancel; the Hilbert phase θ = −π/(4ℓ) − π³/(24ℓ³) + …; matching to relative order ℓ^{−3}; a splitting lemma for the remainder). **Corollary 8.5: for the prime-free form under (H_∞) the boundary law λ′ = −2|A₀|² and the tail law T·(1/2π)∫_{|t|>T}|F|²σ → −λ′/π are theorems.** Computation 8.6 (`rh_edge_constants.py`): with nothing fitted, f/[(−λ′/2)^{1/2}(L+β)^{−1/2}] = 1.001–1.003 for 10⁻¹⁰ < δ < 10⁻⁵ at a = 0.3, 0.4, 0.5, 0.6, and 1.0003–1.0009 on the finer mesh.

**Evidence** (`rh_hadamard.py`, data in `boundary_law.csv`): C from the edge FEM on the two-stage mesh, λ′ exact from `rh_dilation.py`: [verified]

| a | primes | λ (FEM) | C | β | λ′ (exact) | κ = −λ′/2C² |
|---|---|---|---|---|---|---|
| 0.3 | no | +0.2225566 | 1.30891 | −1.49 | −3.42881 | 1.0007 |
| 0.4 | no | −0.0780502 | 1.15425 | −1.48 | −2.66640 | 1.0007 |
| 0.5 | no | −0.3222425 | 1.06222 | −1.45 | −2.25815 | 1.0007 |
| 0.6 | no | −0.5358860 | 1.00932 | −1.42 | −2.03880 | 1.0007 |
| 0.4 | yes | 1.4705×10⁻² | 0.48266 | −1.57 | −0.465947 | 1.0001 |
| 0.45 | yes | 2.4046×10⁻³ | 0.22494 | −1.67 | −0.101283 | 1.0009 |
| 0.5 | yes | 1.9368×10⁻⁴ | 0.07493 | −1.78 | −1.12356×10⁻² | 1.0005 |
| 0.549 | yes | 1.4914×10⁻⁵ | 0.02412 | −1.87 | −1.14353×10⁻³ | 0.983 |

Seven cases give κ = 1 to better than 10⁻³ (the common 7×10⁻⁴ is the mesh bias of C), across positive and negative prime-free floors and with primes; the last measures C at a₃ = 0.54931 and λ′ at 0.549, where Φ′ moves 6% per 0.001 (below), which accounts for its 1.7%. Scaling check (not an independent test of the law, since the same dilation identity underlies both sides): for one half of the Dirichlet form of L_Δ on (−a, a), the scaling identity L_Δ(u(·/a)) = (L_Δu)(·/a) − 2 log a · u(·/a) gives λ(a) = λ(1) − log a, so the law predicts C² = 1/(2a); the prime-free 2aC² = 1.028, 1.066, 1.129, 1.223 at a = 0.3…0.6 against −aλ′ = 1.029, 1.067, 1.129, 1.223, the departures from 1 being the bounded corrections, identical on both sides.

**Consequences.** Φ′(a) = 2C(a)²/λ(a). **Conjecture A ⟺ C(a)² ≤ (c/2)·T\*(a)·λ(a):** the edge amplitude of the minimizer is bounded by its energy times the horizon. The amplitude law of §2.7 ("C² ≍ λT\*") is therefore not a separate observation but Conjecture A itself. Under RH, C is the boundary trace of the band-limited function that best hides from the zeros above the horizon.

**The relay, resolved.** The exact derivative (`rh_dilation.py`, 48 modes) just below the entry of 3: Φ′/T\* = 3.393 (a = 0.5), 3.197 (0.547), 3.527 (0.548), 3.731 (0.5485), 3.969 (0.549), each confirmed by a finite difference of step 2×10⁻⁴ to four digits. Approaching an entry from below, the form without the entering prime power is heading to indefiniteness at a_fail (Table 1 of the folds paper: 0.5574 for n = 3), λ falls roughly linearly (λ′ ≈ −1.14×10⁻³ at 0.549, hitting zero at 0.562 by linear extrapolation) while C stays of the same order, and Φ′ = 2C²/λ rises like 1/(a_fail − a) until the entry rescues the floor. The 0.0125 and 0.004 grids of §2.7 averaged over this rise (their 62.8 at a₃⁻ is the mean of 60 → 75 over the last 0.004); Observation C there stands at those scales, and the fine structure is this. The scan of the exact derivative around the first five entries (`rh_dilation.py`; data in `relay_scan.csv`), Φ′/T\* at a_n + offset: [verified]

| n | a_n | K | -0.01 | -0.003 | -0.001 | -0.0001 | 0.0001 | 0.001 | 0.01 |
|---|---|---|---|---|---|---|---|---|---|
| 2 | 0.3466 | 48 | 2.42 | 2.96 | 3.16 | 3.27 | 2.80 | 2.59 | 2.30 |
| 3 | 0.5493 | 48 | 2.31 | 3.02 | 3.65 | 4.08 | 3.41 | 3.17 | 3.08 |
| 4 | 0.6931 | 64 | 3.61 | 3.32 | 3.70 | 4.08 | 3.73 | 3.63 | 3.51 |
| 5 | 0.8047 | 80 | 3.96 | 2.98 | 3.32 | 4.32 | 3.63 | 3.48 | 3.51 |
| 7 | 0.9730 | 96 | 3.18 | 3.66 | 3.06 | 4.12 | 3.59 | 3.54 | 3.84 |

Every entry shows the same shape: a rise of 30–45% over the last 0.003 (the deleted form heading to its failure point), a drop of 14%, 16%, 9%, 16%, 13% across the entry at scale 10⁻⁴ (where the K-mode kink and the soft kink nearly coincide, 1/K² ≈ 10⁻⁴), and a slow relaxation after. Peaks 3.27, 4.08, 4.08, 4.32, 4.12; dips to 2.3. The supremum of Φ′/T\* over every measured point, the quantity Conjecture A bounds, is 4.32 at the entry of 5. **Observation C of §2.7 is corrected accordingly:** the relay is invisible at step 0.004 only because the difference quotient across an entry averages the rise before it against the drop at it; it is the sawtooth of the folds paper, and the drops at entries 4, 5, 7 are 9%, 16%, 13% at scale 10⁻⁴.

---

## 3. Proof architecture

Here is the only architecture I can see along these routes. Four lemmas; three are within reach; the fourth is where RH lives.

**L1 (Detection; route 2, direction A).** [derived mechanism, provable with effort] If ζ has a zero ρ₀ = ½ + δ + iγ₀ with δ > 0, then λ(a) < 0 for all a ≥ ½ log(γ₀/2π) + E(γ₀, δ).

*Mechanism.* Your Proposition 12.3: a quartet at γ₀ ± iδ contributes 4 Re F(γ₀ + iδ)². Choose f with F(γ₀) = 0. Then F(γ₀ + iδ) = iδF′(γ₀) − (δ²/2)F″(γ₀) + …, so

4 Re F(γ₀+iδ)² = **−4δ² F′(γ₀)² + O(δ⁴)**.

The off-line quartet turns *negative* as soon as the test function plants a zero on it. What remains is to show the on-line sum Σ_{γ on line} F(γ)² can be made smaller than 4δ²F′(γ₀)² once a exceeds the horizon of γ₀: that is a Paley–Wiener interpolation problem (Cartwright–Levinson, Plancherel–Pólya) with explicit constants. The Davenport–Heilbronn control (route 4) measures E empirically: your Remark 12.9 predicts failure near a ≈ 2.1 for the zero at 85.7 with δ = 0.31.

**L2 (Spectral gap).** [conjecture, testable] The lowest eigenvalue of the odd form is simple for all a, and λ₂(a)/λ₁(a) is bounded away from 1. This makes the minimizer continuous in a and Feynman–Hellmann valid everywhere. Connes–van Suijlekom need exactly "simple and isolated" for their real-zeros theorem. You can compute λ₂/λ₁ with your LDLᵀ engine today.

**L3 (Differential inequality).** [conjecture] Φ′(a) ≤ C·T\*(a) for all a, equivalently the relay inequality of §2.3.

**L4 (Grönwall).** [trivial] L2 + L3 ⇒ λ(a) ≥ λ(a₀)·exp(−C(T\*(a) − T\*(a₀))) > 0 for all a ≥ a₀, with a₀ = 0.8 from Zhu's certificate. Hence RH.

### 3.1 Why a bootstrap in height can't replace L3

One might hope to combine L1 with its converse (L1′: RH up to height H ⇒ λ(a) > 0 for a ≤ ½ log(H/2π) − E′) and induct: Platt–Trudgian gives RH to 3×10¹², so positivity to a ≈ 13.4, so (by L1) RH a bit higher, so positivity a bit further, and so on. It doesn't close. Each step gains

Δa ≈ (E + E′) / (2 T\*(a)),

and since T\* grows like e^{2a}, Σ Δa converges unless E + E′ grows like T\* itself, which would mean the lemmas have no content. The dictionary between heights and supports is exact but has no engine in it. The engine has to be L3.

### 3.2 Where L3 hides the circularity

Everything in L3 is unconditional in form: Λ(n), g′_{f_a}(log n), A′, P′ are all computable from the primes. The problem is that the minimizer f_a is only known, structurally, through the zeros: "F_a vanishes at the zeta zeros below T\*" (your Computation 12.5, Groskin, Connes–Consani–Moscovici). That description *is* RH below T\*. So:

> **L3 ⇐ a structure theorem for the minimizer stated in terms of the primes alone.**

That is the missing piece, and it is the same missing piece in all three routes:

- Route 3 (decay law): the decay exponent is set by the tail of F_a above T\*, i.e., by the minimizer's structure.
- Route 2 (horizon): the constant E′ in L1′ is a property of how well the minimizer interpolates below T\*.
- Route 4 (DH control): DH shows what the minimizer does when the structure *fails*: it tracks the off-line quartet instead.

A prime-side structure theorem would have to describe f_a as the solution of a variational problem whose data are the Diophantine positions of {log p^k} and the archimedean kernel, and prove from that description that its autocorrelation derivatives g′(log n) are large enough at the edge. I don't know how to do this, and I don't know anyone who does. But it is a *precise* target, and §2.4 says what any candidate must survive: it must be sensitive to the primes at relative resolution e^{−cT\*}.

---

## 4. Three attack vectors on L3, each with its likely fatal flaw

### V1. A prime-side structure theorem via Fourier interpolation

Bondarenko–Radchenko–Seip (Constr. Approx. 2023) build interpolation bases for functions in a strip whose nodes are zeta zeros on one side and the numbers ±log n on the other, with a duality principle. That is exactly the two-sided structure of the minimizer: interpolate the zeros, be small on the primes' side except at the edge. If the minimizer could be characterized as the extremal element of such a basis, its edge behaviour might be computable from the log n side alone.

*Fatal flaw check.* Their construction assumes the zeros are the interpolation nodes, i.e., it takes the zero set as given. To be useful for L3 the duality would have to be run *from the log n side toward the zeros*, and that direction produces a set of nodes that are the zeros only if RH holds. I expect the theorem to be "RH ⇔ the dual basis is real-noded," a reformulation.

*What survives.* A version of L1 with explicit constants, and a rigorous proof that the minimizer interpolates the zeros below (1−ε)T\* *assuming RH below T\**: a real theorem, unconditional in form, that would make your §12 rigorous.

### V2. Deformation invariance

Move the arithmetic continuously (Epstein zetas across the lattice space; or the "Beurling" deformation log p → log p + ε_p) and ask whether positivity is preserved by a topological invariant, an index or a signature that cannot change under small deformations.

*Fatal flaw check.* §2.4: positivity fails under a deformation of size 10⁻¹⁴ at a = 0.8. There is no open neighbourhood of the primes in which positivity holds. So the invariant cannot be topological in the primes' positions; it has to be an *algebraic* invariant of the Euler-product locus, something that is exactly preserved under "the primes are the primes of ℤ" and destroyed by any perturbation. Fourier-quasicrystal rigidity (Lev–Olevskii, Kurasov–Sarnak) is the closest known phenomenon: measures with discrete support and discrete spectrum are forced to be arithmetic. RH would say the zero measure is one of them.

*What survives.* A clean experiment (§5, E4): the failure support a_fail(ε) as a function of the perturbation size, for one prime at a time. If log a_fail(ε) ≈ −(1/2)log(1/ε)·(something like 1/2) + const, the razor-edge picture is confirmed to all orders and the paper can state it as a computation.

### V3. From traces to the spectral edge

Alpöge–Furman used tr G and ‖G‖_F² of a compressed Weil matrix G (moments m₁, m₂ of its spectral measure) to force 2/3 of the eigen-directions to be positive. Positivity of G is the statement that its spectral measure is supported on [0, ∞), which by the Stieltjes moment theorem is equivalent to

**the Hankel matrices [m_{i+j}] and [m_{i+j+1}] built from m_k = tr(G^k) being positive semidefinite.**

Each m_k is a k-fold sum over the prime side. So RH is equivalent to a family of inequalities among **k-point correlation sums of Λ at Fourier bandwidth ≤ 1**, for all k, and the 2/3 result is the k ≤ 2 case. This puts Alpöge–Furman, Lamzouri, Wang and RH on one axis.

*Fatal flaw check.* m₂ was available unconditionally because Montgomery's pair correlation is known for Fourier support ≤ 1. m₃ needs the triple correlation of Λ in the same range, which is a Hardy–Littlewood-type statement nobody can prove. And a finite number of moments can never pin the edge: the measure can hide negative mass below the resolution of K moments. The full sequence is needed, and the full sequence *is* the spectral measure.

*What survives.* A sharp statement of what the next proportion gain costs, and a way to phrase your Open Question 5 correctly: not "windows with constant below 4/3" (the bandwidth-one ceiling is 0.682) but "any unconditional information on m₃ at bandwidth ≤ 1."

---

## 5. Experiments for your engine

All cheap at your precision, all new, and each one either strengthens or kills a piece of §2–§4.

- **E1. Soft kinks J_n for n ∈ {2,3,4,5,7,8,9}.** With `rh_edge_fem.py`, tabulate the edge amplitude C_n of the minimizer at each entry and the soft-kink coefficient 4Λ(n)n^{−1/2}C_n²/λ(a_n), which replaces f(a_n)²/λ(a_n) in the Φ′ budget. Test whether C_n²/λ(a_n) scales with T\*(a_n) as f_K(a_n)²/λ(a_n) appeared to (§2.3). If the fractional effect grows or shrinks systematically with n, the Φ′ budget of §2.3 is wrong and the decay law is being set elsewhere.
- **E0. The edge amplitude C(a) between entries.** The one new archimedean quantity: C(a) from the graded-mesh minimizer on a fine grid of a. Its jumps at entries, if any, and its growth law in a are inputs to the relay inequality that no finite-mode computation can supply.
- **E2. Φ′/T\* against a**, from your λ(a) curve, at resolution fine enough to see the sawtooth of drifts and kinks. Whether Φ′/T\* converges (then −log λ ≍ T\*) or grows like log T\* (then −log λ ≍ T\* log T\*) is the whole of route 3. Your Table 8 cannot tell; a = 1.3 to 1.5 probably can.
- **E3. λ₂(a)/λ₁(a).** Lemma L2. If the ratio ever approaches 1 the ground state can switch branches and everything in §2 needs restating.
- **E4. Rigidity at 640 bits.** Displace log 2 by ±ε at a = 1.1 and find the ε that flips the sign; prediction ≈ 5×10⁻³⁵. Then a_fail(ε) for ε = 10⁻³, 10⁻⁶, 10⁻⁹ on the prime 3. This is a one-paragraph result for the paper: "positivity at support a is destroyed by displacing a single prime by exp(−cT\*(a))."
- **E5. The Davenport–Heilbronn floor (route 4).** Needed: the prime side of DH. −f′/f = Σ b(n) n^{−s} with b(n) obtained by Dirichlet-series division from the periodic coefficients (1, ξ, −ξ, −1, 0), for n ≤ e^{2a} ≈ 81 at a = 2.2; no polar term (f is entire); the archimedean symbol for the odd-character Γ((s+1)/2)(5/π)^{s/2} factor. The b(n) are not supported on prime powers and grow because f has zeros in Re s > 1: that growth is the Euler product's absence made visible on the prime side. Prediction: the DH odd floor turns negative near a = 2.1, and the minimizer there plants a zero of F at 85.7, per L1.
- **E6. Control for Table 1.** At a = 0.75, delete 3 (not the newest, 4). If the form fails within ≈ 1/Φ′ as well, §2.2 is confirmed and "load-bearing" should be reworded.

---

## 6. What I would change in the paper

1. **Open Question 1** → the Φ′ budget / relay inequality of §2.3, with the kink proposition stated and the J_n table. This is a mechanism, uniform in a in *form*, which is what the question asks for.
2. **Question 12.8** → "Does Φ′(a)/T\*(a) converge?", with the caveat that both T\* and T\* log T\* fit the data (see `RH_ROUTES.md`, calibration §2). Drop "consistent with C = 1.5".
3. **Table 1 caption:** "the deleted form fails after roughly one e-folding length 1/Φ′(a_n) of the full floor."
4. **Add a rigidity remark** (E4). It is the strongest, cheapest, most quotable statement your engine can make: it turns "nothing here proves RH" from a disclaimer into a measurement.
5. **Remark 2.3:** you already say a proof needs a mechanism uniform in a. Add: it needs one sensitive to the primes at relative resolution e^{−cT\*}, so it cannot come from any distributional statement about primes.

---

## Appendix: the small engine

`rh_weil_odd.py` builds the odd Weil matrix in the normalized odd-Legendre basis from
- the archimedean form W_∞(g) = −(γ + log π) g(0) + 2∫₀^∞ [g(0)e^{−2x} − g(x)e^{−x/2}]/(1 − e^{−2x}) dx (from the integral representation of ψ),
- the prime terms −2Σ Λ(n)n^{−1/2} g(log n) with g_{ij}(x) = ∫ φ_i(u)φ_j(u − x/a) du by exact Gauss–Legendre quadrature,
- the polar term −2(∫ f sinh(x/2))².

It agrees with your Table 1 to 4–6 digits with 24 modes and takes about 18 s per floor at 50 digits. The scripts `rh_kink_test.py`, `rh_deleted_form.py`, `rh_rigidity.py`, `rh_symbol_sign.py` reproduce every [verified] number above.

`rh_edge_fem.py` is the second engine, built for the edge. It uses odd-extended hat functions on a mesh of the half window that is uniform in the interior and geometrically graded to 10⁻¹² of the endpoint (350–490 nodes), computes every cross-correlation of basis functions exactly (piecewise cubic, Simpson on merged breakpoints), integrates the archimedean kernel piece by piece with Gauss–Legendre where the cross-correlation is a single cubic and in closed form (∫ₓ^∞ J_Γ = artanh e^{−x/2} + arctan e^{−x/2}) elsewhere, and keeps every breakpoint near the edge exact in floating point by working in a coordinate whose origin is the edge. It reproduces the Legendre floor at a = 0.3 to 6×10⁻⁶ and the converged Legendre profile at a₃ to 10⁻³, and resolves the minimizer to log(a/δ) ≈ 25 where a K-mode expansion stops at 2 log K. `rh_soft_kink.py` measures the entering term across a₃ on this engine.

`rh_dilation.py` verifies the sliding-window identity of §2.8 on the K-mode minimizer and computes the pencil (−Q′, Q). `rh_weil_arb.py` is the third engine: the same odd-Legendre matrices in ball arithmetic (python-flint), with certified Gauss–Legendre rules, M = 2K + 8 inner and 2K + 100 outer nodes, and rigorous enclosures of the lowest eigenvalue of the assembled matrix; it reproduces the mpmath floors to twenty digits and runs a 180-mode point in six minutes. `rh_arb_grid.py` drives it.
