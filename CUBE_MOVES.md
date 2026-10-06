# The cube moves: fixed-ratio increments of the Möbius energy

2026-10-06. On the directive `received/Michaels_Arithmophysics_Rubiks_Cube_Continuation_Prompt.pdf` (8 pages): the
fixed-ratio step L = ⌊N/6⌋, the scale increment Δ₆R(N) = R(N) − R(L), the three exact moves (prefix cancellation, the
{2,3} packet diagonal, the separated family), the target (A), and the return (19)–(21). Script `rh_cube_moves.py`
(data `data/cube_moves.json`), figure `rh_cube_moves_fig.py` → `CUBE_MOVES.png`. Notation as in the directive:
R(N) = ‖c_N‖²_B with c_N(n) = μ(n)/n, B = min, M(x) = Σ_{n≤x} μ(n), S(X) = max_{n≤X} R(n), m(x) = Σ_{n≤x} μ(n)/n.

## 0. What was verified exactly

In rational arithmetic for every N ≤ 120 (the double sum to 60): the two forms of (1); the zero total charge of
r_L and its B-orthogonality to b_L and to u; both splits in (5); (6) = (7) = (8) = R(N) − R(L); the components of
(10) summing to I₂₃(N) = R(N) − D₂₃(N); (13); the bound (12) (maximum 0.800 on that range); the telescoping (19) along
every chain. In floating point to 10⁷: (6), (7), (8) agree to nine digits at N = 16384, 10⁵, 10⁶, 10⁷; (19) holds
for every N ≤ 10⁷ with error below 10⁻⁹; the separated family satisfies ⟨h_N, D_ℓ h_N⟩_B = 0 to 10⁻²⁰,
G_N = (1 + 1/ℓ_N)‖h_N‖²_B, and (16), (17) to nine digits at N = 46656, 10⁶, 4·10⁶, 10⁷. The calibration of the
directive reproduces to the digit: at N = 4,000,000, ℓ_N = 52919, K_N = 12, parents {1, 5, 7, 11},
R(N) = 1.142012922 − 0.319880236 + 0.888378414 = 1.710511100. The packet components at N = 16384 are those of
PARITY_DECAY.md (−0.537728, −0.024996, 0.112661).

## 1. The move sequence and the strongest proved inequality

No new arithmetic inequality for Δ₆R(N) is obtained. The sequence of moves is exact, and each step's complete
cost is accounted for below; the results that are proved are the following.

**Move 1 (prefix cancellation) and the local form of the increment.** For N ≥ 6, L = ⌊N/6⌋, u = c_N − c_L:

    Δ₆R(N) = ‖u‖²_B + 2M(L)(m(N) − m(L)),
    ‖u‖²_B = Σ_{k=L}^{N−1} (M(k) − M(L))²/(k(k+1)) + (M(N) − M(L))²/N.

Proof: B(m, n) = m for m ≤ n, so a vector of zero physical total supported in [1, L] is B-orthogonal to every
vector supported in [L, N]; r_L = c_L − (M(L)/L)e_L has zero total, which gives (5) and (6); expanding
‖b_L + u‖² gives (7), and Abel summation of the tail square gives (8) and the centred form above. The
increment depends on μ only through M on [L, N]; the prefix energy ‖r_L‖²_B is cancelled, not estimated.

**Proposition 1 (two-sided local bound).** For every N ≥ 6,

    M(N)²/N − M(L)²/L  ≤  Δ₆R(N)  ≤  (1 + log(N/L)) · max_{L≤k≤N} M(k)²/k,

with log(N/L) ≤ log 6 + log(1 + 5/(N − 5)). Proof: from (8), the difference between Δ₆R(N) and the left side is
the block sum Σ_{L≤k<N} M(k)²/(k(k+1)) ≥ 0; for the right side bound each M(k)²/k by the window maximum, use
Σ_{k=L}^{N−1} 1/(k+1) ≤ log(N/L), and drop −M(L)²/L. Consequences: (i) if M(x)² ≤ C′x for all x, then
Δ₆R(N) ≤ (1 + log 6 + o(1))C′ uniformly; (ii) if Δ₆R(N) ≤ C uniformly, then by (19) R(N) ≤ max_{n<6}R(n) +
C⌈log N/log 6⌉ and M(N)² ≤ N·R(N) = O(N log N); (iii) (A) for every ε > 0, with any 0 ≤ γ < 6, is equivalent to
the Goal R(N) ≪_ε N^ε: the direction (A) ⇒ Goal is (19)–(21) (Section 2), the converse is Δ₆R(N) ≤ R(N) with
γ = 0. The uniform constant budget sits between the Mertens-type bound M² ≤ C′x and the logarithmic bound on R.
Sharpness on the data (20,373 horizons N ≥ 100 to 10⁷): Δ₆R(N) is at most 0.389 of the upper bound (at
N = 7,108,798) and 0.110 of it on average; at the three largest increments above 10³ the lower bound carries
most of the value (0.177 of 0.228 at N = 42968; 0.180 of 0.227 at 300551; 0.151 of 0.203 at 1,065,673). Those
horizons are new maxima of |M|/√N (M = −88, 240, 422). So the peaks of the increment are the terminal
M(N)²/N at new maxima, and a bound on sup_N Δ₆R(N) is a bound on the block maxima of M(k)²/k.

**Proposition 2 (the cross term is the only piece with an unconditional saving).** With m(x) = Σ_{n≤x} μ(n)/n,
|m(x)| ≪ exp(−c√log x) (the prime number theorem with de la Vallée Poussin's error term, by partial summation
from M(x) ≪ x exp(−c√log x)); hence

    |2M(L)(m(N) − m(L))| ≤ 4|M(L)| max_{x≥L}|m(x)| ≪ √(L·S(L)) · exp(−c√log L),

using M(L)² ≤ L·R(L) ≤ L·S(L). Measured, the cross term is −0.052, −0.029, +0.0024, +0.0006 at N = 16384, 10⁵,
10⁶, 10⁷, against tail squares 0.110, 0.082, 0.088, 0.144: at 10⁷ the increment is the tail square to three
digits. The saving exp(−c√log L) does not reach the factor √L without a bound on S(L) of the kind the Goal
asserts, and the tail square ‖u‖²_B, the centred Mertens energy of one hexadic block, has no estimate below
Proposition 1.

**Move 2 (the {2,3} packet) relocates the increment into one signed band.** Proved: (a) the completed family
is g_N(n) = μ(n)/n · 1{n/gcd(n,6) ≤ N/6}, so ‖g_N‖²_B is itself a Möbius energy; (b) with Kern(N) :=
‖g_N‖²_B − D₂₃(N) = 2Σ_{m<n≤N/6, n<6m, (mn,6)=1} μ(m)μ(n)K₂₃(m,n), the kernel is horizon-free and the pairs with
both parents ≤ L/6 cancel in the difference, so

    Kern(N) − Kern(L) = 2 Σ_{L/6 < n ≤ N/6} Σ_{n/6 < m < n} μ(m)μ(n) K₂₃(m,n)   ((mn, 6) = 1),

verified exactly at N = 7776, 16384, 46656 (−0.363789, −0.334246, −0.324769); (c) the diagonal increment (12):
for L ≥ 6, with a = ⌊L/6⌋ ≥ 1, Σ_{a<m≤L} 1/m ≤ 1/(a+1) + log(L/(a+1)) < 1/2 + log 6, so D₂₃(N) − D₂₃(L) ≤
(2/3)(1/2 + log 6) = 1.528 (and ≤ (2/3)(1 + log L) < 1.862 for L < 6). Measured over 13 horizons from 216 to 10⁷
(Section 3 of the script, figure panel B): ΔD₂₃ between 0.362 and 0.378 for N ≥ 10³, the band ΔKern between −0.364 and
−0.264 (between −0.364 and −0.305 for N ≥ 7776), their sum Δ‖g_N‖² between −0.0001 and +0.058 for N ≥ 7776, Δcross
from +0.26 at 1296 down to |·| ≤ 0.013 for N ≥ 2.8·10⁵, Δunf between −0.27 and +0.085; ΔI₂₃ between −0.352 and −0.219, Δ₆R between +0.011 and
+0.144. The per-parent contributions to the band lie in [−0.021, +0.015] at 7776 and [−0.010, +0.009] at
46656 while the sum of their absolute values grows (1.25, 1.85, 2.95): the inner sums over m are already signed
and the sign-free bound over the newly completed parents diverges. The move pays the diagonal (constant) and
asks, in exchange, that the band be bounded *negatively*: Kern(N) − Kern(L) + Δcross + Δunf ≤ −ΔD₂₃ + allowance.
By (13) this is (A) with the constant shifted by at most 1.862; the directive says so, and the data confirm it:
the band tracks −ΔD₂₃ to within 0.06 for N ≥ 7776.

**Move 3 (the separated family) is proved with its cost, and its cost is loose by three orders.** The chain of
the directive is verified line by line: f = Σ_{d∈⟨2,3⟩, d≤K} D_{d,K} c_K (the Möbius inversion of the 6-part),
‖D_{d,K}c_K‖²_B = d⁻¹R(⌊K/d⌋) (no truncation occurs because every displayed descendant fits below N), hence
‖f‖_B ≤ √S(K)·Π_{p∈{2,3}}(1 − p^{−1/2})⁻¹, ‖(I − D₂)(I − D₃)f‖_B ≤ Π(1 + p^{−1/2})‖f‖_B, the zero-charge
separation ⟨h, D_ℓh⟩_B = 0, and (1 + 1/ℓ_N) ≤ 18/17 for N ≥ 216; so G_N ≤ A*·S(N^{1/6}) with
A* = 500.981190149. Measured: G_N = 0.867, 1.076, 1.142, 1.217 at N = 46656, 10⁶, 4·10⁶, 10⁷ against allowances
718.1, 718.1, 718.1, 861.3 (ratios 828, 667, 629, 708). On the whole range K_N ≤ N^{1/6} ≤ 14: the family has at
most five parents (1, 5, 7, 11, 13), forty vertices with their dilates, and W_N carries the rest of R(N). The
increments (18): G_N − G_L = +0.198, −0.0002, +0.066, +0.141 and W_N − W_L = −0.115, +0.091, −0.055, +0.003 at the
four horizons; ℓ_N runs through every prime from 17 to 113,557 as N runs to 10⁷ (10,701 values) and K_N changes with N, so neither piece
telescopes and neither has a sign.

## 2. The complete return

Proved (the directive's (19)–(21), reproduced): with N₀ = N, N_{j+1} = ⌊N_j/6⌋, stopping at N_J < 6, (19) is an
identity (verified for every N ≤ 10⁷). If Δ₆R(n) ≤ C for all n ≥ 6, then R(N) ≤ max_{n<6}R(n) + C⌈log N/log 6⌉ =
1.4333 + C⌈log N/log 6⌉. If Δ₆R(n) ≤ C_ε n^ε, then R(N) ≤ 1.4333 + C_ε N^ε/(1 − 6^{−ε}). Under (A): for N ≤ X,
R(N) ≤ 1.4333 + C_ε[1 + S(X^{1/6})]^γ Σ_j N_j^ε ≤ 1.4333 + C_ε X^ε[1 + S(X^{1/6})]^γ/(1 − 6^{−ε}), so
1 + S(X) ≤ C′_ε X^ε[1 + S(X^{1/6})]^γ; with ρ = limsup log(1 + S(X))/log X ≤ 1 this gives ρ ≤ ε + γρ/6, i.e.
ρ ≤ 6ε/(6 − γ) for every ε, so ρ = 0. The amortized form: Σ_{j<J}[Δ₆R(N_j)]₊ = R(N) − R(N_J) + Σ_j[Δ₆R(N_j)]₋
and, by Proposition 1, [Δ₆R(n)]₋ ≤ M(⌊n/6⌋)²/⌊n/6⌋ ≤ R(⌊n/6⌋), so R(N) − R(N_J) ≤ Σ_j[Δ]₊ ≤ R(N) + Σ_{j≥1}R(N_j):
the amortized budget and the chain's energies bound each other.

Finite return on [6, 10⁷] (a finite statement, not a theorem): the uniform measured budget is
max_{6≤N≤10⁷} Δ₆R(N) = 1.219 at N = 13, which returns R(10⁷) ≤ 1.433 + 8·1.219 = 11.19; the two-piece version
with S(100) = 1.738 and max_{100<N≤10⁷} Δ₆R(N) = 0.2453 (at N = 110) returns R(N) ≤ 1.738 + 0.2453·⌈log(N/100)/log 6⌉
for every N ≤ 10⁷, i.e. R(10⁷) ≤ 3.455 against the actual 1.8365 and S(10⁷) = 1.8904 (at N = 6,481,601). Chain
budgets: max_N Σ_j[Δ]₊ = 1.651 at N = 607,140 (R = 1.725, seven steps, negative parts 0.426); max_N Σ_j[Δ]₋ =
0.450 at N = 3,643,021; from 10⁷: positive parts 0.543 = 1.8365 − 1.4333 + 0.140 over eight steps.

## 3. The cost that matters

Proved allowances: the packet diagonal, constant (≤ 1.862); the separated family's square, G_N ≤ A*S(N^{1/6}),
a sixth-root factor with exponent γ = 1 (1/6 < 1), for a family of at most forty vertices. The remainder, in
every formulation (‖u‖²_B, or ΔKern + Δcross + Δunf, or W_N − W_L), has no proved bound below Proposition 1,
whose unconditional form is a fixed power: max_{L≤k≤N} M(k)²/k ≤ N·max(|M(k)|/k)² ≪ N exp(−2c(log N)^{3/5}
(log log N)^{−1/5}) (Vinogradov–Korobov), exponent 1 up to the sub-exponential saving. No smaller-envelope factor
with exponent below six is obtained for it. Finite observations, kept apart: on 100 < N ≤ 10⁷ the increment is
at most 0.2453 (declared range 100 < N ≤ 10⁶: 0.2453 at 110, 0.2276 at 42968; held-out range (10⁶, 10⁷]: 0.2034
at 1,065,673), its mean is 0.0548 (0.0306 per e-fold), its 99th percentile 0.15 to 0.18 per decade, and it is
negative at 9 percent of horizons with minimum −0.089 on the top decade. The 555,555 packet-completion horizons
N = 6m, (m, 6) = 1, have the same mean (0.0548) and maximum (0.226) as the rest.

## 4. The reproducible computation

`rh_cube_moves.py` (33 s to 10⁷, no dense matrices): Section 1 the rational checks; Section 2 Δ₆R(N) for every
N from the prefix sums with (6), (7), (8) compared, the decade table, declared and held-out maxima, the eight
largest separated peaks with their terminal, cross and tail parts, the completion horizons, the chain budgets
B±(N) for every N and the exact check of (19); Section 3 the fixed-packet difference (13) with every signed
component at thirteen horizons and the band formula checked at three; Section 4 the separated family with
(14), (16), (17), the calibration and the increments (18); Section 5 the controls, each a single realization used
at both horizons of every pair, zero positions preserved; Section 6 the finite return; Section 7 the sharpness
of Proposition 1. `rh_cube_moves_fig.py` draws `CUBE_MOVES.png`: panel A the increment at every N ≤ 10⁷ in 234
log bins (min, mean, max) against the bin maxima of the controls; panel B the components of (13).

Controls (max and mean of Δ₆ over all N ≤ 10⁷; R(10⁷)): Möbius 1.219 (at 13; 0.245 above 100) and 0.0548;
1.8365. All-positive charges: 1.67·10⁷ and 8.3·10⁶; 2·10⁷ (exact: Δ₆ = 2(N − L) − (H_N − H_L)). Squarefree
positive: 6.16·10⁶ and 3.1·10⁶; 7.39·10⁶. Independent signs on the squarefree support: 3.70 (at 2805) and 0.268
with negative fraction 0.40; 6.52 (expectation of one increment (6/π²) log 6 = 1.09; this realization runs low
on the top decade). Hexadic block shuffle, which keeps every M(6^j) and every zero: 2.26 (at 15) and 0.559;
3.80, decade maxima 0.74 to 1.45 above 100. The Möbius increments are ten times smaller than those of the shuffle that
preserves its block sums and five times smaller than this realization of independent signs (twenty times
smaller than their expectation); the increment of the actual source is
an arithmetic fact about μ inside each hexadic block, not a property of its block sums or of generic signs.

## 5. The surviving signed expression and the inequality still needed

The surviving expression, fully specified, is the first move's remainder:

    Δ₆R(N) = Σ_{k=L}^{N−1} (M(k) − M(L))²/(k(k+1)) + (M(N) − M(L))²/N + 2M(L) Σ_{L<n≤N} μ(n)/n,   L = ⌊N/6⌋,

or, after the second move, the band form Δ₆R(N) = [D₂₃(N) − D₂₃(L)] + 2Σ_{L/6<n≤N/6}Σ_{n/6<m<n} μ(m)μ(n)K₂₃(m,n)
+ Δcross + Δunf. The inequality still needed is (A): for every ε > 0, Δ₆R(N) ≤ C_ε N^ε[1 + S(N^{1/6})]^γ with
some γ < 6; by Proposition 1 it is implied by, and up to the sixth-root factor amounts to, max_{N/6≤k≤N} M(k)²/k ≤
C_ε N^ε[1 + S(N^{1/6})]^γ, and (Proposition 1(iii)) it is equivalent to the Goal. Which moves reduced the cost:
the first removed the prefix energy exactly, so the increment depends on μ only on (L, N]; the second paid the
diagonal at a constant and reduced the rest to one bilinear band of newly completed parents at ratios below 6
under a bounded positive-definite filter; the third pays its family at A*S(N^{1/6}) but covers forty vertices.
Which estimate lost the control: the tail square ‖u‖²_B, the centred Mertens energy of the block. Its only
available bounds are the window bound of Proposition 1 (loose by 2.6 at the peaks, by 9 on average) and, in
the band form, the sign-free sum over parents (1.25, 1.85, 2.95 at 7776, 16384, 46656, growing); the cross term,
the one piece with a proved saving, is 10⁻³ of the increment at 10⁷. The uniform statements that would close by
(19) are, in order of strength, M(x)² ≤ C′x, then sup_N Δ₆R(N) < ∞, then R(N) = O(log N); the data to 10⁷ are
consistent with all three (window maxima of M(k)²/k at most 0.214, increments at most 0.245 above 100), and
none is proved here.
