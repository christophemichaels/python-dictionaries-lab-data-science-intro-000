# Riemann Hypothesis: Weil positivity, the Möbius Green energy, and the prime relay

Research repository of Christophe Michaels. Everything here surrounds the Riemann Hypothesis; nothing here proves it. Each document labels its statements as theorem (cited), proposition (proved here), computation (with precision and cross-checks), or heuristic.

## Papers and notes

| File | What it is |
|---|---|
| `mobius_modifier_v2.pdf` / `.tex` | *The Michaels Möbius Modifier and the Michaels Dynamic DNA Sieve*, version 2. The Green energy of the Möbius vector, its block decomposition, and (new in v2) its identification with the Nyman–Beurling–Báez-Duarte norm on the critical line; RH ⇔ R(N) = O(N^ε); closure statistics to 5×10⁷; the Möbius mean-square constant. |
| `atlas_potential_set.pdf` / `.tex` | Two remarks for the Mathematical Theory Atlas: the potential set 𝒫(a) of the light cone (Krein extensions of the truncated Weil distribution), and the ground-state transform of the Weil form with its spectral corollary (RH ⇔ E_Γ + E_P ≥ ½ Var_ν). |
| `weil_window.pdf` / `.tex` | *The Odd Weil Form on a Window: the sliding-window identity, the logarithmic edge law, and the bounded-relay conjecture.* The paper drawn from the memo: the wall and the (log)^{−1/2} edge law, the soft kink, the exact dilation identity, the boundary law λ′ = −2C², the floor to height 130, the conjectures and what a proof must be. |
| `RH_ROUTES.md` | One hundred research routes toward RH, tiered by credibility, each with its known obstruction and a publishable next step; seven filters any route must pass. |
| `RH_TOP3_PROOF_ARCHITECTURE.md` | The Weil-floor program as one object Φ(a) = −log λ(a): the derivative formula, the kink proposition for finite-mode forms and the soft kink of the exact form (the archimedean wall, the (log)^{−1/2} edge law, and its consequence for the relay), the Φ′ budget, rigidity of positivity under displacement of a prime, the relay conjecture (§2.7: Φ′ ≤ c·T\*, which implies RH, with its evidence to a = 1.2), and the four-lemma architecture with the lemma equivalent to RH isolated. |
| `FOLDS_ERRATA.md` | Corrections and additions for *Primes, Folds, and the One Dot* (draft of 2026-09-28). |
| `atlas/Michaels_Theta_Atlas_CP20.pdf` | *Theta Kernels, Weil Positivity, and the Mathematical Theory Atlas*, complete research compilation through Checkpoint 20 (23 September 2026), 679 pages, bookmarked. Foundations (Sections 1–20, Appendices A–G), Checkpoints 9–20, the 301 atlas entries, ten earlier notes, and the research record. |
| `ATLAS_INDEX.md` | Index of the atlas: contents with PDF page numbers, the 301 entries grouped as the atlas groups them, provenance of the recovered file, and where the atlas meets the work in this repository. |

## Code

All scripts are Python 3 and need only `numpy` and `mpmath` (plus `pymupdf` for PDF handling).

| Script | Purpose |
|---|---|
| `rh_weil_odd.py` | 24–64-mode engine for the odd-sector Weil form on [−a, a] at 50 digits. Reproduces the floors and the Table 1 failure points of the folds paper to 4–6 digits. Finite-mode minima are upper bounds on the true floor; negative values certify indefiniteness, positive values do not certify positivity. |
| `rh_kink_test.py`, `rh_kink_steps.py`, `rh_kink_modes.py`, `rh_kink_endpoint.py`, `rh_kink_profile.py` | The relay transition at the entry of 3: the derivative jump 4Λ(n)n^{−1/2} f(a_n)², its convergence in the step size and in the mode count, and the minimizer's profile near the endpoint. |
| `rh_edge_fem.py` | Edge-adapted finite-element solver for the same form (hat functions on a mesh graded geometrically to 10⁻¹² of the window edge; exact cross-correlations, closed-form archimedean tails). Resolves the minimizer to log(a/δ) ≈ 25 and shows the edge law f(a−δ) = C (log(a/δ) + β)^{−1/2}. Double precision, numpy only. |
| `rh_floor_grid.py` | The floor λ(a), the endpoint value and the two next eigenvalues on a grid of 97 supports on [0.30, 1.20] with the prime-power entries resolved at ±0.004 (K = 40, 60 digits); the data behind the relay conjecture of memo §2.7, saved in `floor_grid_K40.csv`, with the 72-mode drift test in `floor_grid_K72.csv`. Runs in parallel chunks. |
| `rh_weil_arb.py`, `rh_arb_grid.py` | Ball-arithmetic engine for the same form (python-flint), hundreds of times faster, with rigorous eigenvalue enclosures; the drift test to a = 1.5 (height 130) in `floor_grid_arb.csv`. |
| `rh_hadamard.py` | The boundary law λ′(a) = −2C(a)² (memo §2.9): the edge amplitude from the FEM against the exact derivative; data in `boundary_law.csv`. `relay_scan.csv`: the exact derivative at ten offsets around each of the first five entries (memo §2.9, paper Computation 4.4). |
| `rh_dilation.py` | The sliding-window identity (memo §2.8): the dilation derivative of the floor split into archimedean, prime and polar parts on the minimizer, verified to twenty digits, and the pencil (−Q′, Q) showing that the relay is a ground-state property, not an inequality between forms. |
| `data/` | The inputs of `rh_hadamard.py` and of the soft-kink table: the edge-FEM minimizers (`fem_a*.json`, prime-free at a = 0.3–0.6 and with primes at 0.4, 0.45, 0.5, a₃) with their edge fits, the ball-arithmetic dilation logs at the same supports (`dilation_a*.log`), and the 134-node soft-kink finite differences (`softkink_134node.json`). `python3 rh_hadamard.py data/fem_a0.3_primefree.json data/dilation_a0.3_primefree.log data/fem_a3_full.json data/dilation_a0.549_full.log` reproduces the boundary-law table. |
| `rh_fourier_check.py` | Independent Fourier-side check of the sliding-window identity at K = 12, a = 0.6 (paper Computation 5.9): D_∞ from (1/2π)∫|F|² t ∂_t Re ψ with F through spherical Bessel functions and the Parseval tail, D_P from g′(log n) of the exact autocorrelation, against the matrix-side a·λ_K′; twelve digits. |
| `rh_bv_check.py` | Bounded variation of the minimizers (paper Computation 5.8): sup, total variation and the Hardy integral of the edge-FEM minimizers (`data/fem_*.json`, including the 124- and 237-node refinements at a = 0.5) and of the K-mode minimizers; outputs in `data/bv_fem.log`, `data/bv_fem_refinement.log`, `data/bv_kmode.log`. Behind the hypothesis (H_BV) under which the floor is Lipschitz. |
| `rh_zero_side.py` | The floor and the sliding-window identity on the zero side (paper Computation 7.5): for the K-mode minimizer, 2Σ_{γ>0}|F(γ)|² against λ_K and λ_K + 2Σγ(|F|²)′(γ) against a·λ_K′, over the first 6000 zeros (`data/zeros_6000.txt`, from `mpmath.zetazero`, the first 400 to 36 digits), with Riemann–von Mangoldt tails, the share of the floor carried by the zeros below the horizon, and the zeros carrying most of it; logs in `data/zero_side_*.log`. |
| `rh_tail_law.py` | The tail law (paper Computation 7.6, Conjecture 7.7): the fraction of the floor carried by the zeros above height T is Φ′/(πT), i.e. Φ′ = π·T_eff. Modes `fem` (edge-FEM minimizers on the zeros, closed-form transform of the piecewise-linear function), `coefs`/`shares` (per-zero shares of the K-mode minimizers), `overlay` (ϑ(T/T*) across supports, medians). Data: `data/tail_law_fem.log`, `data/tail_law_kmode/`, `data/zeros_1500_50digits.txt` (first 1500 zeros to 50 digits, needed at a = 1.5). |
| `rh_alignment.py` | The zeros' rotation at the entries (paper Computation 7.10): the tail-weighted average of cos(2γa) over the first 6000 zeros anti-aligns to −0.5…−0.6 exactly at the entries of the primes and is a cusp there; output in `data/alignment.log`. |
| `rh_soft_kink.py` | The soft kink at the entry of 3: finite differences of the edge-FEM floor across a₃ at scales 10⁻³ … 10⁻⁷, against the first-order entering energy and the law 4Λ(3)3^{−1/2}C²/(log(1/ε) + β). |
| `rh_deleted_form.py` | Failure points of the form with the newest prime power deleted (Table 1 of the folds paper). |
| `rh_rigidity.py` | Sensitivity of positivity to displacing a single prime. |
| `rh_symbol_sign.py` | Where the truncated Weil symbol is negative, relative to the horizon. |
| `rh_decay_fit.py` | Fits of −log λ(a) against T\*(a) and T\* log T\*. |
| `rh_mobius_green.py`, `rh_mobius_zeros.py`, `rh_mobius_analyze.py` | The Möbius Green energy to 5×10⁷, the zero-side constants from the first 4000 zeros, and the comparison. |

## Status

Reviewed items are marked in the documents. External reviews of 2026-09-29 (independent referee on §12 of the Möbius paper; a four-branch review of the atlas and the folds paper) and of 2026-09-30 (independent referee on the Weil-window paper, whose report is applied in the revised `weil_window.pdf` and in memo §§2.6–2.9) have been applied; see the commit history.

## Building the papers

`tectonic mobius_modifier_v2.tex` or `pdflatex` (packages: amsmath, amssymb, amsthm, booktabs, hyperref, enumitem, graphicx).
