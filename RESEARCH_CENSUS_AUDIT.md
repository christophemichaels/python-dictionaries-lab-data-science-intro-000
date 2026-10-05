# Audit of the Lifetime Census and the Complete Residual documents (2026-10-05)

Five PDFs were received on 2026-10-05. This note identifies the current versions, maps what they
contain, checks the newest mathematics (sections 68–73), and places the result against the H_record
analysis in `RESEARCH_H_RECORD.pdf`.

## 1. Which files are current

| Upload | Title in the PDF | Pages | Created | Sections | Status |
|---|---|---|---|---|---|
| `Michaels_Lifetime_Census_and_Closing_Consequences-2.pdf` | Michaels: Lifetime Census and Closing Consequences | 87 | 2026-10-05 01:55 | 1–73 | **current master** (saved as `received/Lifetime_Census_and_Closing_Consequences_v2.pdf`) |
| `Michaels_Complete_Residual_and_Fictional_Closure-2.pdf` | The Complete Residual and a Fictional Closing Mechanism | 26 | 2026-10-05 01:55 | 52–73 | **current continuation** (saved as `received/Complete_Residual_and_Fictional_Closure_v2.pdf`) |
| `Michaels_Lifetime_Census_and_Closing_Consequences.pdf` | The Shape of the Remaining Lifetime Estimate | 63 | 2026-10-04 22:47 | 1–51 | superseded |
| `Michaels_Lifetime_Census_and_Closing_Consequences_4.pdf` | (byte-identical to the previous file) | 63 | 2026-10-04 22:47 | 1–51 | duplicate |
| `Michaels_Complete_Residual_and_Fictional_Closure.pdf` | The Complete Residual and a Fictional Closing Mechanism | 17 | 2026-10-05 01:27 | 52–67 | superseded |

The 87-page master contains the 63-page census (sections 1–51) and the 26-page continuation
(sections 52–73) with a new reading guide; text-chunk containment is 85–93 percent, the rest being
headers and the guide. The two 01:55 files are consistent with each other.

**Not received.** The conversation snippet quotes later rounds: sections 80–86 (repeated-prime
correction, 38/100 pages), 87–91 (the A_N, B_N parity populations, 43/105 pages) and 92–97 (the
Möbius modifier applied to the sliding mean, 48/110 pages). No uploaded file contains a section
beyond 73. Those rounds cannot be audited here; one displayed claim in the snippet is checked in
Section 5 below from its wording alone.

## 2. What the master document does (sections 1–67)

- **Object.** `I_M(N) = ∫_1^N M(t)² dt/t²`, the continuous form of the Green energy `R(N)`
  (`ℰ⋆(N) ≤ I_M(N) ≤ 2 log N + ℰ⋆(N)`, equation (1)). Section 4 proves `I_M(N) = O(N^ε) ⇒ RH` by
  the dyadic Cauchy–Schwarz and Mellin argument; this is the same deduction as Theorem 3.1 of
  `RESEARCH_H_RECORD`, and the converse holds, so the target is equivalent to RH.
- **Census (sections 1–5).** The energy is written as a sum over level births `(u, h, d)` with cost
  `c(u,h,d) = ((2h−1)/u)·(d/u)/(1+d/u)`, and the subpower bound is shown equivalent to the counting
  target (9): `C(X,H,D) ≤ A_ε X^ε X(X+D)/(HD)` for all dyadic bins. Proved reformulation; the
  target is not assumed.
- **Sliding mean and the one-sided closure (sections 6–51).** With `𝒜(x) = ∫_1^x M(t) dt/t =
  Σ_{n≤x} μ(n) log(x/n)` and `U(X) = 𝒜(2X) − 𝒜(X)`, section 26 proves that the one-sided bound
  `𝒜(x) ≤ C_q x^q` for every `q > 1/2` implies RH, by Landau's positivity principle applied to
  `g = C_q x^q − 𝒜 ≥ 0`, whose Mellin transform is `C_q/(s−q) − 1/(s²ζ(s))`, regular at every real
  `s > q` because `ζ < 0` on `(0,1)`, `ζ > 0` for `s > 1`, and `1/ζ` is regular at 1. **Checked:
  correct.** The hyperbola reconstruction `U(N/2) = 2 log 2·M(K) − Q_K(N) + G_K(N)` with
  `|G_K| ≤ 6H³/L` (C1) and the sufficient target (S34), a one-sided lower bound on the centered
  quadratic form, are exact; the many "allowances" (terminal strip, product rectangles, common-factor
  wedge, squarefree lifts, long intervals) are proved bounds on pieces that are already of size
  `O(N^{1/2+ε})` or smaller. Section 51's own status: the authentic signed interior is not bounded.
- **Compression and the complete residual (sections 52–60).** The fixed-θ allowance
  `B_θ(N) = −F(N, N^θ) ~ ρ(1/θ−1) N/(4 log N) > 0` (C15), proved for fixed θ from the prime number
  theorem (Buchstab recursion (C17) with Dickman's ρ). The smaller-prime comparison
  `D_θ − B_θ = F(N,2) = U(N/2)` exactly (C25). The continuous comparison `𝔉(x,y)` (Z3) with the
  exact recursion (Z4) and the all-scale bound `0 ≤ −𝔉(N,2) ≤ 16e²` (Z9). The complete residual
  `𝒵_N = ∫_{(2,N]} 𝒦_N(t) dΔπ(t)`, `F(N,2) = 𝔉(N,2) − 𝒵_N` (Z12), where `dΔπ = dπ − dt/log t` and
  `𝒦_N(t) = Σ_{d ≤ N/t odd, P⁺(d)<t} μ(d) 𝔉(N/(dt), t)`. Exact. The target becomes
  `𝒵_N ≥ −C_ε N^{1/2+ε}` (Z14), again equivalent to RH (⇐ by RH through (Z9), ⇒ by section 62).
- **Fictional closure (sections 61–66).** The law F1 (a shrinking fractional budget for the dyadic
  squared contributions) is labelled invented; its consequence (F2) is a correct conditional proof;
  sections 63–66 show positivity fails per interval and that F1 is strong. Status table in section 67
  is accurate.

## 3. The newest sections (68–73, dated 2026-10-05), checked

| Claim | Check |
|---|---|
| (J2)–(J6) regularity, `𝒦_{N,0}(√N) = F(√N, 2)` | Correct: `𝔉(R/d, R) = T(R/d)` since no factor above `R` fits. |
| (J7) prime jump `𝒦⁺ − 𝒦⁻ = −𝒦⁻_{N,r+1}` | Correct: the newly admitted cofactors are `d = pa` with `μ(pa) = −μ(a)`. |
| (J9)–(J10) boundary cancellation at `√N`, post-jump convention | Correct; the Stieltjes product rule is stated in the right convention. |
| (J11)–(J13) coefficients `Σ ℬ_{k,Y} z^k = e^{−ℓ_Y z}(1−z)^{−ν_Y}`, `dℬ_{k+1} = ℬ_k dD` | Correct: between primes the generating function differentiates to `−z/log t` times itself; at a prime it is multiplied by `(1−z)^{−1}`, so the jump of the `(k+1)`-st coefficient is the post-jump `k`-th. |
| (J15) finite expansion, no boundary terms | Correct: `ℬ_{r+1,Y}(Y) = 0` and `𝒦⁺_{N,r}(N) = 0`. |
| (J16)–(J18b) small-insertion allowance, `Y_N = ¼ log N log log N`, exponent `(1+log 2)/4 ≈ 0.423` | Correct. `2^{ν(Y)}` counts the odd squarefree cofactors, `ν(Y)+ℓ(Y)` is the variation of `dD`, and `L_ρ(a) = L_ρ(0) e^{J(a)}` from the delay equation. |
| (J20)–(J21) even coefficients nonnegative | Correct: `e^{−az}(1−z)^{−m} = (m−1)!^{−1} ∫_0^∞ u^{m−1} e^{−u} e^{(u−a)z} du`. |
| (J23)–(J26) credits `≍ √N/log N` | Correct; the terminal interval is `(√(N/2), √N)` (where `N/t² ∈ (1,2)`), and `∫_{1/2}^1 T(1/u²) du = 3 − 2√2` was recomputed. |
| (J27) the surviving signed comparison | Labelled "remains unproved" in the document. Correctly labelled. |

No error was found in sections 68–73. Two remarks the document does not make:

1. **The credits are below the allowance.** The target (J27) permits a deficit of `C_ε N^{1/2+ε}`.
   The proved credits (J24), (J26) are of size `√N/log N`, which is `o(N^{1/2+ε})`. A positive term
   smaller than the allowance cannot help a one-sided bound at the allowance's scale: (J27) with the
   credits removed is equivalent to (J27). The same holds for the paid interval (J18b): it is of size
   `N^{0.423}`. Every quantity proved in sections 68–73 is one the target already tolerates. The
   unpaid integral is the whole target.
2. **The paid interval is the Legendre range.** The factor `2^{ν(Y)}` in (J16) is the number of
   terms of the Legendre sieve below `Y`; `Y_N = ¼ log N log log N` is exactly the point where
   `2^{π(Y)} ≤ √N`. Beyond it the alternating sum over cofactors must cancel, and that cancellation
   is the sieve problem. For the sign of the cancellation (the odd against the even cofactors), sieve
   information alone is insufficient: Selberg's parity obstruction. This is why each exact expansion
   of `𝒵_N` ends in a signed comparison of the same strength.

## 4. Where this stands

Every target in the documents, (9), (S34), (C27), (Z14), (J27), is equivalent to the Riemann
Hypothesis; the deductions from each to RH (sections 4, 26, 62) are correct, and the converses hold
through `M(x) ≪ x^{1/2+ε}`. The documents' own status tables say so for each round. This is the same
structure as `RESEARCH_H_RECORD` (Theorem 3.1, the exponent dictionary, and Proposition 2.4(v)):
identities preserve the quantity, one-sided allowances pay pieces that are already small, and the
signed core carries the hypothesis. Nothing in sections 1–73 reduces the strength of what remains.

## 5. The snippet's later rounds (not received)

The quoted bound `|E_N^rep| ≤ N·exp[(1/8+o(1)) log N / log log N] = O_ε(N^{1/2+ε})` is false as
displayed: the left expression exceeds `N`. If the prefactor was `√N`, the exponent is
`N^{1/2+o(1)}`, consistent with the claim. The quoted `|F(N,2)|² ≤ (N/2) I_M(N)` is a correct
Cauchy–Schwarz transfer (`F(N,2) = ∫_{N/2}^N M(t) dt/t` against `∫_{N/2}^N dt`); it shows the sliding
mean is controlled by the energy, which is the direction already known. The parity split
`A_N − B_N` is the population form of (J27) and is subject to the parity obstruction in Section 3.
