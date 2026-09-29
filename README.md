# Riemann Hypothesis: Weil positivity, the Möbius Green energy, and the prime relay

Research repository of Christophe Michaels. Everything here surrounds the Riemann Hypothesis; nothing here proves it. Each document labels its statements as theorem (cited), proposition (proved here), computation (with precision and cross-checks), or heuristic.

## Papers and notes

| File | What it is |
|---|---|
| `mobius_modifier_v2.pdf` / `.tex` | *The Michaels Möbius Modifier and the Michaels Dynamic DNA Sieve*, version 2. The Green energy of the Möbius vector, its block decomposition, and (new in v2) its identification with the Nyman–Beurling–Báez-Duarte norm on the critical line; RH ⇔ R(N) = O(N^ε); closure statistics to 5×10⁷; the Möbius mean-square constant. |
| `atlas_potential_set.pdf` / `.tex` | Two remarks for the Mathematical Theory Atlas: the potential set 𝒫(a) of the light cone (Krein extensions of the truncated Weil distribution), and the ground-state transform of the Weil form with its spectral corollary (RH ⇔ E_Γ + E_P ≥ ½ Var_ν). |
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
| `rh_soft_kink.py` | The soft kink at the entry of 3: finite differences of the edge-FEM floor across a₃ at scales 10⁻³ … 10⁻⁷, against the first-order entering energy and the law 4Λ(3)3^{−1/2}C²/(log(1/ε) + β). |
| `rh_deleted_form.py` | Failure points of the form with the newest prime power deleted (Table 1 of the folds paper). |
| `rh_rigidity.py` | Sensitivity of positivity to displacing a single prime. |
| `rh_symbol_sign.py` | Where the truncated Weil symbol is negative, relative to the horizon. |
| `rh_decay_fit.py` | Fits of −log λ(a) against T\*(a) and T\* log T\*. |
| `rh_mobius_green.py`, `rh_mobius_zeros.py`, `rh_mobius_analyze.py` | The Möbius Green energy to 5×10⁷, the zero-side constants from the first 4000 zeros, and the comparison. |

## Status

Reviewed items are marked in the documents. External reviews of 2026-09-29 (independent referee on §12 of the Möbius paper; a four-branch review of the atlas and the folds paper) have been applied; see the commit history.

## Building the papers

`tectonic mobius_modifier_v2.tex` or `pdflatex` (packages: amsmath, amssymb, amsthm, booktabs, hyperref, enumitem, graphicx).
