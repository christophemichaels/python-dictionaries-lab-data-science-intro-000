# Three Routes, One Mechanism: What an Unconditional Proof Through the Weil Floor Would Have to Look Like

*Companion to `RH_ROUTES.md`. Prepared 2026-09-29 for C. Michaels, "Primes, Folds, and the One Dot".*

## Ground rules

- **This is not a proof and does not contain one.** It's an attempt to design the *shape* of a proof along the three routes I ranked highest for your toolkit (the horizon theorem, the decay law, and the Davenport–Heilbronn control), to isolate the single lemma on which such a proof would stand, and to say where the circularity hides.
- Everything below is one of three things, and each is labelled:
  - **[verified]** — computed with a small independent engine for the odd Weil form (`rh_weil_odd.py`, 24 modes, 50 digits). It reproduces your floors (λ(0.5) = 1.943×10⁻⁴, λ(0.6) = 6.0×10⁻⁷) and your Table 1 failure points (0.371602 vs your 0.371601; 0.55735 vs your 0.557323). A K-mode minimum is a minimum over a subspace and therefore an *upper bound* on the true floor: a negative finite-mode value certifies indefiniteness (up to rounding), a positive one does not certify positivity, which is why certification (Zhu) needs the control of the omitted directions that this engine does not attempt. Every [verified] number below is a statement about the 24-mode form; the rigidity conclusions (§2.4) survive because they rest on negative values.
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

**Proposition.** At each entry a_n = ½ log n, λ is continuous and its derivative jumps up by

**λ′(a_n⁺) − λ′(a_n⁻) = 4 Λ(n) n^{−1/2} f_{a_n}(a_n)²,**

where f_{a_n}(a_n) is the endpoint value of the normalized minimizer at entry.

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

Measured: J_2 = 5.96 (Φ′ from 41.2 to 35.2, a 14% drop); J_3 = 14.4 (from 77.4 to 63.0, a 19% drop). [verified]

Two ratios stay in narrow bands across a = 0.35 to 1.2:

- **Φ′(a)/T\*(a) ≈ 2.8 to 4.1** (from your Table 8 differences: 2.8–3.7; at the entries of 2 and 3: 3.3→2.8 and 4.1→3.3).
- **f_a(a)²/λ(a) ≈ 0.16–0.32 · T\*(a)** (values 2.2, 4.4, 5.0, 7.0, 9.9 at a = 0.4, …, 0.8).

Together these give [conjecture] **J_n/Φ′ ≈ c·Λ(n)/√n with c ≈ 0.3**, a fractional drop independent of a. Individually the drops vanish as n grows, but the entries become dense (about e^{2a}/a per unit of a), and Σ_{n ≤ x} Λ(n)/√n ~ 2√x says the total fractional kick per unit of a is about 0.6 e^{a}, growing without bound. The smooth drift must grow at the same rate for Φ′ to stay ≈ 3.5 T\*. **The relay is this balance.** Your Open Question 1 becomes:

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

The floor converges (its steps shrink: 2.0, 1.5, 0.9, 0.5, 0.4 ×10⁻⁸). The endpoint value does not on this range: it is 13% lower at K = 64 than at K = 24 and still falling, so the kick's predicted size has dropped 25% and has no visible limit. A fit f(K) = f_∞ + c/K gives f_∞ ≈ 0.0075 from each of three independent pairs (24/64, 32/56, 40/64), which would make the true kick about 1.4×10⁻⁴; but a slow power law f ∝ K^{−0.15}, tending to zero, fits the six values equally well. A Legendre expansion converges slowest at the endpoint, and the endpoint is exactly where the relay lives. Two consequences. First, this is the concrete form of the reviewers' point about omitted directions: a certified relay transition must enclose f(a_n), not only λ(a_n), and the complementary-space debit for f(a_n) is much larger than for λ. Second, it is an open question whether the true minimizer's endpoint value is nonzero at all; if it tends to zero the kick is a truncation effect and the relay proceeds through the interior of the window instead. The question is analytical, not numerical: the Euler–Lagrange equation of the minimizer is an integral equation on [−a, a] whose kernel has the singularity J_Γ(s) ~ 1/(2|s|), and the endpoint behaviour of its solutions (a boundary layer or a finite limit) is a Wiener–Hopf question about that kernel. Numerically, the direct test is the profile of the minimizer near y = a at increasing K:

| y/a | 0.5 | 0.9 | 0.95 | 0.98 | 0.99 | 0.995 | 0.999 | 1.0 |
|---|---|---|---|---|---|---|---|---|
| K = 24 | 1.3664 | 0.11681 | 0.04864 | 0.02364 | 0.01743 | 0.01430 | 0.01072 | 0.00947 |
| K = 40 | 1.3664 | 0.11690 | 0.04873 | 0.02379 | 0.01745 | 0.01422 | 0.01093 | 0.00865 |
| K = 64 | 1.3664 | 0.11689 | 0.04875 | 0.02381 | 0.01747 | 0.01429 | 0.01100 | 0.00821 |

For y ≤ 0.99a the minimizer is converged to better than 0.1% and is nonzero (f(0.99a) = 0.0175). The last one percent of the window is a boundary layer where the K-mode minimizer is still adjusting: f(0.999a) drifts up (0.0107 → 0.0110) while f(a) drifts down (0.0095 → 0.0082), the two closing toward each other. A vanishing endpoint value would require the drop to occur inside the last 0.1% of the window, against that trend. Taken alone, the profiles suggest a finite limit f(a₃) ≈ 0.008–0.011. The structural fact below argues the other way, and it is the more reliable guide. [verified: profiles]

### 2.6 The archimedean wall [verified identity; reading conjectural]

For f supported in [−a, a], the archimedean form in x-space is c₁‖f‖² + ½∬_{ℝ²} J_Γ(u−v)|f(u)−f(v)|² du dv with J_Γ(s) = e^{−|s|/2}/(1−e^{−2|s|}). Splitting the double integral according to whether v lies inside the window gives an exact decomposition of the window form:

**A(f) = c₁‖f‖² + ½∬_{[−a,a]²} J_Γ(u−v)|f(u)−f(v)|² + ∫_{−a}^{a} V(y)|f(y)|² dy,  V(y) = ∫_{|v|>a} J_Γ(y−v) dv.**

Since J_Γ(s) ~ 1/(2|s|), the potential diverges logarithmically at the edge. Numerically, at a = a₃:

| y/a | 0.5 | 0.9 | 0.99 | 0.999 | 0.9999 | 0.99999 |
|---|---|---|---|---|---|---|
| V(y) | 3.4390 | 4.1342 | 5.2650 | 6.4143 | 7.5654 | 8.7166 |
| V − ½log(1/(a−y)) | 2.7929 | 2.6834 | 2.6628 | 2.6608 | 2.6607 | 2.6606 |

So **V(y) = ½ log(1/(a−|y|)) + 2.6606 + o(1)**: the odd Weil form on a window is a logarithmic kinetic form (the interior J_Γ Dirichlet form, which has symbol ~ log|t|, the logarithmic Laplacian of Chen–Weth) plus a logarithmic confining wall, plus the prime shifts and the rank-one polar term, both bounded.

Consequences:

- The minimizer is pushed off the edge by a divergent potential. The Euler–Lagrange equation [c₁ + V(y) − λ] f(y) + (kinetic term)(y) = (bounded) forces either f(y) → 0 as |y| → a or a compensating divergence of the nonlocal kinetic term. The simplest consistent behaviour is f(y) ≈ C/log(1/(a−|y|)), i.e. **f(a) = 0 with logarithmic approach**. The K = 64 profile is consistent with this on the converged range (f·log(1/δ) = 0.081, 0.076, 0.076 at δ = 0.01, 0.005, 0.001 in units of a), and a Legendre truncation resolving scales ~1/K² would then show f_K(a) ~ c/log K, which is exactly the slow, non-shrinking decrease observed from K = 24 to 64. A finite limit f(a₃) ≈ 0.008 is not excluded by the numerics, but it would require the kinetic term to diverge at the edge.
- If f(a) = 0 for the exact form, the kink of §2.1 is a **truncation feature of the K-mode form**: exact for every K, vanishing in the limit. The entering prime power then enters at *zero* first-order rate, with the term −2Λ(n)n^{−1/2} g(log n) growing like (a − a_n)·C²/log²(1/(a − a_n)): a soft kink, a derivative that is continuous but whose slope diverges. The relay is real but gentler than a kick, and the Φ′ budget of §2.3 must be re-derived with this growth law. The rigidity results (§2.4) are unaffected, since they rest on negative values of the finite form.
- The boundary layer seen in §2.5 is then **archimedean**, and the prime-free form confirms it. Ground state of the archimedean-plus-polar odd form at a = 0.3 (no prime inside the window; λ = 0.222557, converged to six digits), values of f at y/a: [verified]

| K | 0.5 | 0.9 | 0.99 | 0.995 | 0.999 | 1.0 |
|---|---|---|---|---|---|---|
| 24 | 1.6068 | 1.2553 | 0.7477 | 0.6777 | 0.5518 | 0.4966 |
| 40 | 1.6067 | 1.2562 | 0.7479 | 0.6724 | 0.5591 | 0.4647 |
| 64 | 1.6067 | 1.2558 | 0.7472 | 0.6747 | 0.5657 | 0.4399 |

  Same signature with no primes present: interior converged, endpoint value falling steadily with K (0.4966, 0.4647, 0.4399: −6.4% then −5.3%, not converging), f(0.999a) rising toward it (0.5518, 0.5591, 0.5657). The edge law differs between the two cases (f(0.99a)/f(0.999a) is 1.32 here at K = 64 against 1.59 at a₃), so the edge behaviour is not a universal power; that is what one expects from logarithmic corrections and is a further reason to settle it analytically rather than by fitting.

This is the point at which the numerics stop being informative and the analysis takes over: the edge asymptotics of the first eigenfunction of "log-Laplacian + ½log(1/dist) wall" on an interval is a well-posed problem, and settling it settles the kink. It is also the first place in this program where the archimedean place acts *alone* on the mechanism, which is what §1.2 said any proof would need. (A first K = 40 run returned nonsense because the inner quadrature rule had 64 nodes, exact only to polynomial degree 127; the engine now scales the rule with K. Nothing at K ≤ 32 was affected.)

**What a certified transition needs.** Interval enclosures of (i) the K-mode matrix at a₃ ± h for h ≤ 10⁻⁵ (the archimedean integrals reduce to incomplete-gamma series, so this is elementary-function interval arithmetic), (ii) the lowest eigenvalue by interval LDLᵀ from both sides, (iii) the endpoint value of the enclosed eigenvector, and (iv) Zhu's complementary-space floor for both λ and f(a₃). Items (i)–(iii) are within reach of `rh_weil_odd.py` rewritten over `mpmath.iv`; item (iv) is your team's certificate machinery.

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

- **E1. Kinks J_n for n ∈ {2,3,4,5,7,8,9}.** You have all the minimizers. Tabulate f(a_n)²/λ(a_n) and the fractional drop J_n/Φ′(a_n⁻). Test the conjecture J_n/Φ′ ≈ 0.3 Λ(n)/√n. If the fractional drop is instead growing or shrinking systematically with n, the Φ′ budget of §2.3 is wrong and the decay law is being set elsewhere.
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
