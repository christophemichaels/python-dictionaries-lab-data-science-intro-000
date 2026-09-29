# Errata and edits for "Primes, Folds, and the One Dot" (draft of 2026-09-28)

Collected from the external review of 2026-09-29 and from this session's checks. Each item names the place, the problem, and the minimal edit.

## Must fix

1. **Theorem 13.2 (Cramér).** The statement assumes RH only. The limit Σ_ρ 1/|ρ|² = 2 + γ − log 4π requires RH *and simple zeros*; with a multiple zero the constant is Σ m_ρ²/|ρ|² over distinct zeros. Edit: "Assume RH and that all zeros are simple."

2. **Computation 2.1 / Table 1 / Figure 2, "the full form stays positive definite throughout 0.1 ≤ a ≤ 1.1".** A minimum over K modes is a minimum over a subspace, hence an upper bound on the true floor λ(a). A positive K-mode value does not certify positivity of the form; a negative one does certify indefiniteness (up to rounding). So: the "deleted form fails at" column is a certified statement about the full form; "the full form stays positive definite" is a statement about the 40- or 72-mode form unless Zhu's tail control is invoked. Edit: say "the K-mode form" where positivity is asserted, and cite Zhu for certified positivity on a ≤ 0.8 (and CP22 for a ≤ 3/5).

3. **Computation 12.5 and Table 8.** Same point: λ(a) and the zero sums are 72-mode quantities. The agreement "zero sum/λ = 0.99998" is an end-to-end check of the engine, not a certificate. Add one sentence.

4. **Heuristic 12.4 and Remark 12.7.** The horizon T*(a) is a heuristic scale. Nothing in the paper proves that an off-line zero is invisible below it or becomes detectable when its height is crossed; Proposition 12.3 gives the contribution of a quartet but not its size relative to the on-line sum. The current labels are right; make sure no sentence elsewhere (§14 (2), Remark 12.9) reads as if the threshold were proved. Suggested wording for Remark 12.9: "the heuristic predicts".

5. **Question 12.8.** Replace by: "Does Φ′(a)/T*(a) converge, where Φ = −log λ?" and drop "consistent with C = 1.5". Both −log λ ≍ T* and ≍ T* log T* fit the data to a = 1.2 (RMS 0.53 vs 0.57 in the fit of the companion memo). Under the second law Question 12.8 as stated is false even assuming RH.

## Should add

6. **Table 1 reading.** The gap a_fail − a_n equals λ(a_n)/|λ′(a_n⁻)| = 1/Φ′(a_n) to first order (0.0243 predicted vs 0.0250 observed for n = 2). Table 1 measures the e-folding length of λ, not a special role of the newest prime power. Every term larger than the margin is load-bearing in the same sense. Suggested caption: "the deleted form fails after roughly one e-folding length 1/Φ′(a_n) of the full floor."

7. **Kink proposition.** At each entry a_n = ½ log n, the derivative of the *K-mode* floor jumps up by 4Λ(n)n^{−1/2} f_K(a_n)², where f_K is the normalized K-mode minimizer (Feynman–Hellmann). Verified at the entry of 3 with 24 modes to 2×10⁻⁴ at step h = 10⁻⁵. The *exact* form has no first-order kink. On the window the archimedean form is one half of the Dirichlet form of the logarithmic Laplacian plus a bounded form (memo §2.6: interior log-kernel plus the wall V(y) = ½ log(1/(a − |y|)) + 2.66 + o(1)), and bounded weak solutions of such Dirichlet problems vanish at the boundary like (log(1/dist))^{−1/2} [Hernández-Santamaría, López Ríos, Saldaña, *Optimal boundary regularity and a Hopf-type lemma for Dirichlet problems involving the logarithmic Laplacian*, DCDS 45 (2025), Thm 1.1; the rate is sharp, Thm 1.2; Chen–Weth had (log)^{−τ} for τ < ½]. The minimizer computed on a mesh graded to 10⁻¹² of the edge (memo §2.6) has f(a − δ) = C (log(a/δ) + β)^{−1/2} (1 + o(1)) with 1/f² linear in log(a/δ) over fifteen e-foldings; C ≈ 0.024 at a₃ and the K-mode endpoint value is C(2 log K + β′)^{−1/2}, which is the slow decrease seen from K = 24 to 64. The entering term therefore grows like 4Λ(n)n^{−1/2} C² (a − a_n)/(log(1/(a − a_n)) + β): a **soft kink**, λ′ continuous with a jump at scale ε of 4Λ(n)n^{−1/2}C²/(log(1/ε) + β). Measured at the entry of 3: the finite-difference jump at scales ε = 10⁻⁴, 10⁻⁵, 10⁻⁶ is 1.90, 1.48, 1.21 ×10⁻⁴ (predicted 2.02, 1.54, 1.24), against a constant 2.28×10⁻⁴ for the 24-mode form. State the proposition for the K-mode form, where it is exact, and the soft law for the exact form; the relay is real but logarithmically gentler than a kick.

8. **Rigidity.** At a = 0.8, displacing log 2 by 3×10⁻¹⁴ makes the 24-mode form indefinite (certified direction: negative value). Displacing log 3 by 10⁻³ makes it indefinite before a = 0.65. Suggested remark: "positivity at support a is destroyed by displacing a single prime by roughly exp(−cT*(a))"; this is the quantitative content of Remark 2.3's "mechanism uniform in a" and rules out any distributional input as the engine of such a mechanism.

9. **Symbol picture.** The truncated symbol Ψ_a is negative on 27–43% of [0, T*] and on about 10% of [T*, 3T*]; the negative set's density equals the Haar measure of the chord region on the torus of prime phases (torus model vs line: 0.091 vs 0.098 at a = 1.1). Ties §3 to §12 and explains why no pointwise or statistical bound on Ψ_a can give positivity.

10. **§14 open questions.** Item (5): the Alpöge–Furman method has a bandwidth-one ceiling of ≈ 0.682, and the Montgomery–Taylor window already gives 0.6725; a window search gains at most ≈ 0.01. Rephrase as "any unconditional information on the third moment of the compressed Weil matrix at Fourier support ≤ 1."

## Minor

11. Table 8 caption: state that the last column is the distance from γ₁ to the nearest zero of F_f (it is in the text but not the caption).
12. Reference [46] (Zhu) now has a version of 2 September 2026; check the certified range quoted (a ≤ 0.8) against the current version.
13. Remark 12.6: Groskin's arXiv:2605.20224 is at v4 (14 August 2026); the "about 300 digits" figure should be checked against v4 (307–329 matching digits for the first ten zeros).
