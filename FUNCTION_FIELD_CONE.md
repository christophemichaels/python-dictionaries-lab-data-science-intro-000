# The cone over a curve: the paperwork of the curves projected onto the magenta cone (2026-10-05)

The question was whether the information in the curves over finite fields, where the Riemann Hypothesis is
a theorem (Weil 1948), can be projected onto the light cone of the folds paper and the potential set of
`atlas_potential_set.pdf`. It can, exactly, and the projection says one thing: **over a curve the cone
closes at a finite half-width, `a = g log q`, where `g` is the genus; over the integers it never closes.**
At the closing support the floor of the Weil form is exactly zero, the potential set collapses to the
single true zero measure, the null vector of the form *is* the polynomial whose roots are the zeros, and
positivity of the form at that one support is already the Riemann Hypothesis for the curve. Everything
below is classical mathematics in the cone's language, plus a computed demonstration on five curves
(`rh_ff_cone.py`, `data/ff_cone_*.json`, figure `FUNCTION_FIELD_CONE.png`). Nothing here bears on RH for
`ζ` beyond the comparison in Section 6.

## 1. Dictionary

| Over the integers (the magenta cone) | Over a curve `X / 𝔽_q` of genus `g` |
|---|---|
| prime powers `n = p^k`, weight `Λ(n)` | closed points `P` of degree `d`, weight `d` (norm `N(P) = q^d`) |
| position `log n` on the prime line | position `log N(P) = d log q` |
| cone at half-width `a`: `log n < 2a` | cone at half-width `a = D log q`: `d ≤ 2D` |
| entry of `n` at `a = ½ log n` | entry of degree `d` at `a = ½ d log q` |
| zeros `½ + iγ`, infinitely many | Frobenius eigenvalues `α_i = √q e^{iθ_i}`, exactly `2g` of them |
| zero measure `ν_ζ = Σ δ_γ` on `ℝ` | zero measure `ν_X = Σ_{i≤2g} δ_{θ_i}` on the circle |
| explicit formula with `Σ Λ(n) n^{−1/2} g(log n)` | `ν̂_X(m) = Σ_i e^{imθ_i} = q^{m/2} + q^{−m/2} − N_m q^{−m/2}`, `N_m = #X(𝔽_{q^m}) = Σ_{d|m} d b_d` |
| Weil form at support `a` on `C_c^∞(−a, a)` | Toeplitz form `Q_D(h) = Σ_{|j|,|k|≤D} h(j) h(k) ν̂_X(j−k) = Σ_i |ĥ(θ_i)|²` on `h : [−D, D] → ℂ` |
| floor `λ(a)` | floor `λ(D)` = least eigenvalue of `T_{2D+1} = [ν̂_X(j−k)]` |
| potential set `𝒫(a)` (Krein) | `𝒫_X(D)` = positive symmetric measures on the circle with the moments `ν̂_X(m)`, `|m| ≤ 2D` (Carathéodory–Toeplitz) |
| horizon `T*(a) = 2π e^{2a}`, a height | horizon `a = g log q`, a degree: the record ends at degree `2g` |
| RH: `∩_a 𝒫(a) = {ν_ζ}` | RH for `X` (Weil): `𝒫_X(D) = {ν_X}` for every `D ≥ g` |

The identification of the two explicit formulas is exact: for a curve the "prime side" of the Weil form
at support `D` uses only the counts `N_1, …, N_{2D}`, i.e. the closed points of degree `≤ 2D`, which is the
cone `log N(P) < 2a` with `a = D log q`, with the same edge `log n = 2a` as the magenta cone.

## 2. What the cone does over a curve

Write `z_i = e^{iθ_i} = α_i / √q`, `i = 1, …, 2g`. The multiset `{α_i}` is closed under complex conjugation
(`P(T)` has integer coefficients) and under `α ↦ q/α` (the functional equation), so `{z_i}` is closed under
`z ↦ z̄` and `z ↦ 1/z̄`, and `ν̂_X(−m) = conj ν̂_X(m)`: `T_{2D+1}` is Hermitian for every `D`.

**Proposition.** Let `X / 𝔽_q` have genus `g` with Frobenius eigenvalues `α_1, …, α_{2g}` and `z_i` as above.

1. *(Open cone.)* If the `z_i` are distinct, `T_{2D+1}` is positive definite for `D ≤ g − 1`: `λ(D) > 0`
   and `𝒫_X(D)` is infinite. In general `T_{2D+1}` is positive definite exactly while `2D + 1 ≤ r`, `r` the
   number of distinct `z_i`.
2. *(The horizon.)* For every `D ≥ g`, `T_{2D+1}` is singular: the vector `p` of coefficients of
   `Π_i (w − z̄_i)` (degree `2g`, extended by zeros) satisfies `T p = 0`. Hence `λ(D) = 0` for `D ≥ g`,
   `𝒫_X(D) = {ν_X}` for `D ≥ g`, and the roots of the null polynomial are the zeros: the counts
   `N_1, …, N_{2g}` inside the cone reconstruct `θ_1, …, θ_{2g}` at the horizon.
3. *(Sectors.)* With distinct conjugate pairs and no real `α_i`: the even sector (`h(−n) = h(n)`,
   `ĥ` a cosine polynomial of degree `D`) has floor `> 0` for `D < g` and `= 0` for `D ≥ g`; the odd sector
   (`ĥ = sin θ ·` a polynomial of degree `D − 1` in `cos θ`) has floor `> 0` for `D ≤ g` and `= 0` for
   `D ≥ g + 1`. The even sector closes the cone; the odd sector follows one degree later.
4. *(Finite-support Weil criterion.)* Let `{z_i}_{i≤2g}` be any multiset in `ℂ^*` closed under `z ↦ z̄` and
   `z ↦ 1/z̄`, with `c_m = Σ_i z_i^m`. Then `T_{2g+1} = [c_{j−k}]` is singular, and it is positive
   semidefinite if and only if every `|z_i| = 1`. So for a curve, positivity of the Weil form at the single
   support `D = g` is equivalent to the Riemann Hypothesis for `X`, and an impostor configuration (a pair
   `r z, z/r` off the circle, `r ≠ 1`) is rejected by the cone at some `D ≤ g`.

*Proof.* (1) `T_{2D+1} = Σ_i v_i v_i^*` with `(v_i)_j = z_i^j`, `|j| ≤ D`; it is positive definite iff the
`v_i` span `ℂ^{2D+1}`, i.e. iff the Vandermonde matrix in the distinct nodes `z_i` has rank `2D + 1`, i.e.
iff `2D + 1 ≤ r`. A positive definite Toeplitz matrix has infinitely many positive extensions (the
nondegenerate Carathéodory–Toeplitz problem), so `𝒫_X(D)` is infinite.
(2) `(T p)_j = Σ_k c_{j−k} p_k = Σ_i z_i^j Σ_k z_i^{−k} p_k = Σ_i z_i^j p(1/z_i)` where `p(w) = Σ_k p_k w^k`;
since `{1/z_i} = {z̄_i}`, `p(1/z_i) = 0` for every `i`. For `D ≥ g` pad `p` with zeros. A singular positive
semidefinite Toeplitz matrix of rank `r` has a unique positive extension, a measure with `r` atoms
(Carathéodory–Fejér), so `𝒫_X(D) = {ν_X}`, and `ν_X` is that measure.
(3) The even and odd sectors are orthogonal and `T`-invariant because `ν_X` is symmetric; the vanishing
conditions are the stated interpolation conditions at the `g` nodes `cos θ_i`.
(4) Singularity is (2). Suppose `T_{2g+1} ⪰ 0`. By Carathéodory–Fejér, `c_m = Σ_{l≤r} ρ_l w_l^m` for
`|m| ≤ 2g`, with `ρ_l > 0`, `|w_l| = 1`, `r = rank T_{2g+1} ≤ 2g`, and the null space of `T_{2g+1}` consists
of the polynomials of degree `≤ 2g` vanishing at every `w_l`. The null vector `p` of (2) therefore vanishes
at every `w_l`, so each `w_l` is one of the `1/z_i = z̄_i`, which is one of the `z_i`. Now
`e_m = Σ_i z_i^m − Σ_l ρ_l w_l^m = 0` for `0 ≤ m ≤ 2g`; grouping equal nodes, `e_m` is an exponential sum
over at most `2g` distinct nodes (the distinct `z_i`), and a nonzero exponential sum over `s` distinct nodes
cannot vanish at `s` consecutive integers (Vandermonde). So every coefficient is zero: every `z_i` is some
`w_l`, of modulus `1`. The converse is `Q_D(h) = Σ |ĥ(θ_i)|² ≥ 0`. For the impostor, `T_{2g+1}` is still
singular by (2) and is not positive semidefinite by the "only if", so some principal `T_{2D+1}`, `D ≤ g`, has
a negative eigenvalue. ∎

Only (2) and (4) need the exact count `2g`. In the number-field cone the zeros are infinitely many, and
nothing in this proposition survives: see Section 6.

## 3. The computation

`rh_ff_cone.py q "c_0,…,c_{2g+1}"` takes a hyperelliptic curve `y² = f(x)`, `f` monic squarefree of degree
`2g + 1` over `𝔽_q` (`q` an odd prime), counts `N_m = #X(𝔽_{q^m})` for `m ≤ 2g` by brute force in
`𝔽_{q^m}` with the quadratic character (one point at infinity), forms `P(T)` from `exp(Σ N_m T^m / m)`,
checks the functional equation `P(T) = q^g T^{2g} P(1/(qT))` and `|α_i| = √q`, derives the closed-point
counts `b_d`, checks the counted `ν̂_X(m)` against the zeros and the explicit formula on a random `h`, and
computes the floors, the sectors, the ranks, the null-vector reconstruction at `D = g`, and the first `D`
at which an impostor pair `(√q r, √q / r)` replacing `α_1, α_2` is rejected.

| curve | `N_1, …, N_{2g}` | `P(T)` | checks |
|---|---|---|---|
| `g = 1`, `𝔽_3`, `f = 1 + x² + x³` | 6, 12 | `1 + 2T + 3T²` | FE ✓, `|α| = √3` ✓, `ν̂` 2e−16, EF 1e−16 |
| `g = 2`, `𝔽_3`, `f = 2 + 2x³ + x⁵` | 1, 13, 28, 73 | `1 − 3T + 6T² − 9T³ + 9T⁴` | FE ✓, ✓, 6e−15, 2e−15 |
| `g = 3`, `𝔽_3`, `f = 2 + x² + x³ + 2x⁵ + x⁷` | 3, 13, 15, 81, 248, 793 | `1 − T + 2T² − 6T³ + 6T⁴ − 9T⁵ + 27T⁶` | FE ✓, ✓, 1e−14, 0 |
| `g = 2`, `𝔽_5`, `f = 4 + 2x + 4x² + x³ + x⁵` | 9, 37, 108, 625 | `1 + 3T + 10T² + 15T³ + 25T⁴` | FE ✓, ✓, 6e−15, 4e−15 |
| `g = 3`, `𝔽_5`, `f = 4 + 4x + x² + 2x³ + 4x⁵ + x⁷` | 7, 23, 91, 651, 3157, 15875 | `1 + T − T² − 13T³ − 5T⁴ + 25T⁵ + 125T⁶` | FE ✓, ✓, 5e−15, 0 |

(`ν̂`: maximal difference between the counted `ν̂_X(m)`, `m ≤ 2g`, and `Σ_i e^{imθ_i}` from the zeros;
EF: zero side minus count side of the explicit formula on a random `h` supported in `[−g, g]`.)

Floors by support `D` (half-width `a = D log q`): even sector / odd sector; rank of `T_{2D+1}` in
parentheses. Exact zeros are reported by `eigvalsh` as `±0.0` (below `10⁻¹²` in every case).

| curve | `D = 0` | `D = 1` | `D = 2` | `D = 3` | `D = 4` | closes at |
|---|---|---|---|---|---|---|
| `g = 1`, `𝔽_3` | 2 (1) | **0** / 2.667 (2) | 0 / **0** (2) | 0 / 0 (2) | | `a = log 3 = 1.099` |
| `g = 2`, `𝔽_3` | 4 (1) | 1.000 / 5.000 (3) | **0** / 2.000 (4) | 0 / **0** (4) | 0 / 0 (4) | `a = 2 log 3 = 2.197` |
| `g = 3`, `𝔽_3` | 6 (1) | 4.543 / 7.000 (3) | 1.853 / 4.441 (5) | **0** / 4.441 (6) | 0 / **0** (6) | `a = 3 log 3 = 3.296` |
| `g = 2`, `𝔽_5` | 4 (1) | 0.707 / 6.200 (3) | **0** / 1.923 (4) | 0 / **0** (4) | 0 / 0 (4) | `a = 2 log 5 = 3.219` |
| `g = 3`, `𝔽_5` | 6 (1) | 5.600 / 5.400 (3) | 2.664 / 2.534 (5) | **0** / 2.228 (6) | 0 / **0** (6) | `a = 3 log 5 = 4.828` |

In every case: positive definite for `D < g`, rank `2D + 1`; singular from `D = g` on with rank frozen at
`2g`; the even floor reaches zero at `D = g` and the odd floor at `D = g + 1`, as the Proposition says.

**Reconstruction at the horizon.** For each curve the eigenvector of the zero eigenvalue of `T_{2g+1}`,
read as a polynomial, has all `2g` roots on the unit circle (moduli `1.000…`) at the angles `θ_i`, with
maximal angle error `≤ 9 · 10⁻¹⁶`. The counts inside the cone at `a = g log q` return the zeros.

**The impostor test.** Replacing the conjugate pair `α_1, α_2` by `√q · r, √q / r` with
`r ∈ {1.05, 1.1, 1.2, 1.5, 2}` (a pair off the circle, functional equation kept), the first support at which
`T_{2D+1}` has a negative eigenvalue is `D = 1, 2, 3, 2, 2` (`g = 1, 2, 3`, `𝔽_3`) and `D = 2, 2, 2, 2, 1`,
`3, 3, 2, 2, 2` (`g = 2, 3`, `𝔽_5`) for the five radii: always `≤ g`, as the Proposition requires, and
usually at the horizon itself for a small displacement. This is the function-field Davenport–Heilbronn:
the configuration that no closed point of degree `< 2D` can distinguish from a true zero measure is caught
by the degrees `≤ 2g`, never later.

## 4. The figure

`FUNCTION_FIELD_CONE.png` (script `rh_ff_cone_fig.py`). Panel A is the magenta cone drawn for the genus-3
curve over `𝔽_3`: the edge `log n = 2a`, the closed points of degree `d` entering at `a = ½ d log q` with
their counts `b_d` (3, 5, 4, 17, 49, 128, 338, 796), the floors `λ(D)` along the time axis, the dashed line
`a = ½ g log q` where `N_1, …, N_g` are inside (enough to fix the zeros when the genus is known), and the
horizon `a = g log q` where the floor is zero and the cone has closed; above it the shaded region is
`𝒫_X = {ν_X}`. Panel B plots the floors of the three `𝔽_3` curves against `a = D log q`, each reaching
zero exactly at `g log q` (even) and `(g + 1) log q` (odd). Panel C is the number-field floor of the odd Weil
form from `floor_grid_K40.csv` and the certified enclosures of `floor_grid_arb.csv`, on a log scale: it
decays through `10⁻⁹³` by `a = 1.5` and is positive at every support computed, with the entries of
`2, 3, 4, 5, 7, 8, 9, 11` marked.

## 5. With the genus known: the half horizon

The cone closes at `a = g log q` on positivity alone. If one is also told that the measure has exactly `2g`
atoms and is symmetric under `θ ↦ −θ` (the genus and the functional equation), then `N_1, …, N_g`, which
are inside the cone at `a = ½ g log q`, already determine `P(T)` through Newton's identities and
`P_{2g−k} = q^{g−k} P_k`. The dashed line in Panel A is this half horizon. In the number field the
functional equation is known but the "genus" is infinite, so there is no half horizon either.

## 6. The comparison with `ζ`, and what this does and does not say

1. **No finite horizon for `ζ`.** For `h ∈ C_c^∞(−a, a)`, `ĥ` is entire of exponential type `a`, with at
   most `(2a/π + o(1)) T` zeros of height `≤ T` (Jensen), while `ζ` has `(T/2π) log T` zeros there. So `ĥ`
   cannot vanish on all the zeros unless `h = 0`: on RH the Weil form is positive definite at every finite
   `a`, `λ(a) > 0` for all `a`, and `𝒫(a)` is infinite for every `a` (every nondegenerate truncated moment
   problem has infinitely many solutions). The collapse `∩_a 𝒫(a) = {ν_ζ}` of `atlas_potential_set.pdf`,
   Proposition "limit", happens at `a = ∞` only. The computed decay of `λ(a)` in Panel C is the number-field
   shadow of the function-field zero: the floor goes down because more and more of the zero measure is
   being pinned, and it never arrives because the zero measure has infinitely many atoms to pin.
2. **Where the positivity comes from over a curve.** Weil's proof that `|α_i| = √q` is the Hodge index
   theorem on the surface `X × X` applied to the graph of Frobenius (the Castelnuovo–Severi inequality):
   a statement about intersection numbers of divisors on a surface, not about the counts `N_m`. The counts
   are its consequences. In the cone's language the pattern `N_1, …, N_{2g}` is finite, and the sign that
   closes it is supplied by a geometry standing outside that pattern. This is the one setting in the
   programme where the thesis of `RH_GODEL_THEOREM.pdf`, Section 8, holds as a theorem: the whole pattern is
   explained by a framework above it, and because the pattern has a last stone (degree `2g`) the framework
   can be checked against all of it. For `ζ` no finite-dimensional `H¹` is known, the pattern has no last
   stone, and that is exactly the `Π⁰₁` form of RH.
3. **What is new here and what is not.** The Proposition is elementary and classical in substance (Weil;
   Carathéodory–Toeplitz; Newton's identities); the packaging as the closing of the cone, the identification
   of the null vector with the zeros, and the computed table are the only contributions of this note. The
   note proves nothing about `ζ`. It gives the magenta cone its calibration: the cone is the right picture,
   it closes where the genus is finite, and the distance it never covers over the integers is the infinite
   genus of `ζ`.

## 7. Files

- `rh_ff_cone.py`: the computation (brute-force point counts in `𝔽_{q^m}`, zeta polynomial, checks, floors,
  reconstruction, impostor test). `python3 rh_ff_cone.py 3 "2,0,1,1,0,2,0,1"` runs the genus-3 curve in
  0.3 s; the `𝔽_5` genus-3 curve takes 3 s.
- `data/ff_cone_g{1,2,3}_q3.json`, `data/ff_cone_g{2,3}_q5.json`: outputs.
- `rh_ff_cone_fig.py`, `FUNCTION_FIELD_CONE.png`: the figure.
