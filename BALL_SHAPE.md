# The shape of the final estimate, in ball arithmetic

2026-10-06. On `received/Michaels_Mobius_Energy_and_Light_Plane_Correspondence.pdf` (44 pages, "Möbius Energy and the
Light-Plane Correspondence"), whose single remaining input is Hypothesis 1.1 on the complete scale increment, and on the
request to determine the shape of that estimate with certified numbers, across the horizon in both directions.
Script `rh_ball_shape.py` (python-flint 0.9, arb/acb balls at 200 bits; 361 s to 10⁷), full output `data/ball_shape.log`,
data `data/ball_shape.json`, figure `rh_ball_shape_fig.py` → `BALL_SHAPE.png`. Every number below that is written as
[m ± r] is a ball: the true value lies within r of m. Where a maximum over 10⁷ horizons is reported, the position of the
maximum was found in a floating-point pass whose a-priori rounding error is below 5·10⁻⁹ (sums of nonnegative terms), and
the value at that position is then certified by the ball computation. Status labels as in the manuscript's convention:
**[source]** proved there, **[derived]** proved here, **[finite]** certified on a finite range, **[hypothesis]**.

## 1. The estimate whose shape is asked for

With R(N) = Σ_{k<N} M(k)²/(k(k+1)) + M(N)²/N, L_N = ⌊N/6⌋, δ_N = R(N) − R(L_N), S(X) = max{1, R(n): n ≤ X}, B = 1 + S:

    Hypothesis 1.1:   δ_N ≤ C (1 + log N)^a exp{A (log N)^ρ} B(N^{1/6})⁵      (N ≥ 6),

with the normalized coefficient c_N = [δ_N]₊/B(N^{1/6})⁵ and the three shapes: constant (a = A = 0), polylogarithmic
(a > 0), stretched-exponential (A > 0, 0 < ρ < 1). Theorem 1.2 of the manuscript turns any of them into R(N) = O(N^ε) and
the Mellin conclusion; Theorem 11.2 gives the constant-coefficient return explicitly. By (11.7) the same number δ_N can be
read in four planes: the field increment ‖t‖² + 2M(L)(h(N) − h(L)), the Hardy trace difference tr(Q_N − Q_L), the two-mode
signed trace B⁵·tr D_N, and the two-phase circle difference. The shape question is: which of the three forms do the
certified values of c_N across horizons require, and how do the four readings distribute the same increment.

## 2. Certified values

**Base values, exact rationals [derived].** R(1) = 1, R(2) = 1/2, R(3) = R(4) = 5/6, R(5) = 43/30, R(6) = 14/15,
R(7) = R(8) = R(9) = 143/105, R(10) = 223/210, R(11) = R(12) = 3083/2310, R(13) = 51629/30030, R(14) = 20452/15015. Hence
S(X) = 1 for X ≤ 4, 43/30 for 5 ≤ X ≤ 12, 51629/30030 for X = 13, 14, and B(N^{1/6})⁵ = 32 (N < 5⁶ = 15625), 85.312
(15625 ≤ N < 13⁶ = 4826809), 148.677 (4826809 ≤ N < 15⁶) on the whole computed range.

**The chain from 10⁷ [finite, certified].** B⁵ is the exact envelope factor; c_N = δ_N/B⁵; ν± are the two eigenvalues of
the rank-two readout (10.10), computed exactly from p₀ = R(L)/B⁵, q₀ = ‖t‖²/B⁵, z₀ = M(L)(h(N) − h(L))/B⁵ because the
calibrated lift is an isometry (10.6), so no lift has to be constructed; cos is z₀/√(p₀q₀), the angle between the carried
prefix and the tail in the lifted space.

| N | L | R(N) | δ_N | cross 2M(L)(h(N)−h(L)) | ‖t‖² | B⁵ | c_N | ν₊ | ν₋ | cos |
|---|---|---|---|---|---|---|---|---|---|---|
| 10⁷ | 1666666 | 1.83653342597 ± 6e-13 | 0.1441432521 ± 1e-11 | +0.000558 | 0.143585 | 148.68 | 0.00096951 | +0.003836 | −0.002866 | +0.0006 |
| 1666666 | 277777 | 1.69239017386 | 0.0723526683 | −0.003759 | 0.076112 | 85.31 | 0.00084810 | +0.004562 | −0.003714 | −0.0054 |
| 277777 | 46296 | 1.62003750554 | 0.0255878851 | −0.049414 | 0.075002 | 85.31 | 0.00029993 | +0.004196 | −0.003896 | −0.0714 |
| 46296 | 7716 | 1.59444962045 | 0.0744958225 | −0.003144 | 0.077639 | 85.31 | 0.00087322 | +0.004487 | −0.003614 | −0.0046 |
| 7716 | 1286 | 1.51995379794 | 0.0536959630 | −0.004241 | 0.057937 | 32 | 0.00167800 | +0.009986 | −0.008308 | −0.0073 |
| 1286 | 214 | 1.46625783491 | 0.0626801441 | 0 | 0.062680 | 32 | 0.00195875 | +0.010300 | −0.008341 | 0 |
| 214 | 35 | 1.40357769085 | 0.1102262802 | +0.035118 | 0.075108 | 32 | 0.00344457 | +0.011598 | −0.008153 | +0.0563 |
| 35 | 5 | 1.29335141062 | −0.1399819227 | −0.235546 | 0.095564 | 32 | −0.00437444 | +0.008993 | −0.013368 | −0.3182 |

All radii are below 5·10⁻⁹ (most below 10⁻¹¹); the log holds the full balls and the chains from 8388608, 4·10⁶, 6⁸,
10⁶, 10⁵, 16384 and 10⁴ (55 rows).

**Decade maxima of δ_N and of c_N [finite, certified].**

| range | max δ_N | at N | B⁵ | c_N there |
|---|---|---|---|---|
| [10, 100] | 1.21924741925 ± 3e-12 | 13 | 32 | 0.03810148185 ± 1.5e-12 |
| [100, 10³] | 0.245249781112 ± 5e-13 | 110 | 32 | 0.007664055660 ± 2.4e-13 |
| [10³, 10⁴] | 0.206881209106 ± 3e-14 | 2838 | 32 | 0.006465037785 ± 4.4e-13 |
| [10⁴, 10⁵] | 0.227638335587 ± 4e-14 | 42967 | 85.312 | 0.002668316701 ± 3.8e-13 |
| [10⁵, 10⁶] | 0.226518477568 ± 1e-13 | 300551 | 85.312 | 0.002655190020 ± 2.1e-13 |
| [10⁶, 10⁷] | 0.203372918230 ± 4e-13 | 1065673 | 85.312 | 0.002383883861 ± 1.2e-13 |

The global maximum of c_N over N ≥ 6 is at N = 13: c* = δ₁₃/32 = (R(13) − R(2))/32 = 18307/480480 = 0.038101482…, an
exact rational. The maximum over N ≥ 101 is 0.007664055660 (N = 110); over N ≥ 10⁴ it is 0.00629455534 (N = 11773,
certified separately); over N ≥ 10⁶ it is 0.002383883861 (N = 1065673). The envelope B(X) = 1 + S(X) at X = 10, …, 10⁷
is 2.433, 2.738, 2.738, 2.738, 2.796, 2.823, 2.890: log B(10⁷) = 1.0614.

**The explicit constant-coefficient return (Theorem 11.2) [derived, balls].** With c* = 18307/480480: α = log 5/log 6 =
0.8982444017, C* = 73/30 + c*/log 6 = 2.454598175, K(c*) = 5(5/4·log 2 + ¼ log C* + 5/16·log 6) = 8.254247889, and the
returned bound log B(X) ≤ K (log X)^α equals 100.263 at X = 10⁷ (actual 1.061), 138.127 at 10¹⁰, 1092.76 at 10¹⁰⁰. With
the coefficient above 100 instead (0.007664), K = 8.2455670: the constant is insensitive to c* because 73/30 dominates C*.

**Hardy coefficients across the horizon [finite, certified; (9.4)–(9.6)].** γ_{k,N} = Σ_{n≤N} μ(n)/n·P_k(log n) were
computed from the generating function exp(−x z/(1 − z)) by a ball series exponential per squarefree n (6.1·10⁶ terms in
81 s), which avoids the cancellation of the explicit Laguerre sums. At N = 10⁷: γ₀ = 0.00010152 (= h(N)), γ₁ = 0.99836095,
γ₂ = 0.43437582, γ₃ = 0.0595562, γ₄ = 0.0578294, γ₅ = −0.2762894, γ₆ = −0.0554286, γ₇ = −0.0678697, γ₈ = −0.1855482,
γ₉ = −0.1271233, γ₁₀ = −0.0015426, γ₁₁ = 0.0227590, γ₁₂ = −0.0337187 (radii ≤ 5·10⁻⁹). With K = K_N = ⌊(log N)^{3/4}⌋ at
the head of each pair and the same K at both ends, along the chain from 10⁷:

| N | L | K | H_K(N) | H_K(L) | T_K(N) − T_K(L) | δ_N | allowance (9.9) | allowance (9.10) |
|---|---|---|---|---|---|---|---|---|
| 10⁷ | 1666666 | 8 | 1.31074066 | 1.26004258 | 0.0934451775 | 0.1441432521 | 4.7·10²² | 4.4·10²³ |
| 1666666 | 277777 | 7 | 1.24784664 | 1.24354797 | 0.0680539947 | 0.0723526683 | 1.9·10¹⁹ | 1.4·10²¹ |
| 277777 | 46296 | 6 | 1.22607935 | 1.24342946 | 0.0429379961 | 0.0255878851 | 1.0·10¹⁶ | 4.4·10¹⁸ |
| 46296 | 7716 | 5 | 1.20843274 | 1.20833514 | 0.0743982234 | 0.0744958225 | 7.0·10¹² | 1.2·10¹⁶ |
| 7716 | 1286 | 5 | 1.20833514 | 1.20244736 | 0.0478081889 | 0.0536959630 | 1.3·10¹² | 6.7·10¹⁴ |
| 1286 | 214 | 4 | 1.19246326 | 1.18027414 | 0.0504910306 | 0.0626801441 | 1.4·10⁹ | 1.6·10¹² |
| 214 | 35 | 3 | 1.17100753 | 1.09419997 | 0.0334187173 | 0.1102262802 | 2.2·10⁶ | 3.2·10⁹ |

Theorem 9.1 at N = 1000 (K up to 800, 400-bit balls): H_K(1000) = 0.9408, 1.1860, 1.2031, 1.2622, 1.2651, 1.2709, 1.2734,
1.2799, 1.2999, 1.3367 at K = 1, 2, 5, 10, 20, 50, 100, 200, 400, 800, against R(1000) = 1.4594094175788 ± 9·10⁻¹⁵; the
tail is still 0.1227 at K = 800, and |γ_k(1000)| = 0.970, 0.0659, 0.00070, 0.0088, 0.0188 at k = 1, 10, 100, 400, 800.

**The three-step modifier (2.23)–(2.31) [finite, certified].** On the chain from 10⁷, δ_N − Q⁽³⁾_N = −1.06·10⁻⁷ (10⁷),
+9.5·10⁻⁷ (1666666), +2.6·10⁻⁶ (277777), +2.5·10⁻⁶ (46296), −2.4·10⁻⁵ (7716), +1.0·10⁻³ (1286), −1.7·10⁻³ (214),
−0.107 (35), against the proved allowance 2 log(N/L) = 3.58 at every step; ĉ⁽³⁾_N = [Q⁽³⁾]₊/B⁵ agrees with c_N to seven
digits at 10⁷ (0.000969507 against 0.000969506).

**The two-phase circle readout (8.2)–(8.4) [finite, certified with tails].** The series over k ∈ ℤ was summed to |k| ≤ K
with the omitted tail enclosed by A²ℓ/(π²(2K − 1)), A = Σ_{n≤N} μ(n)²n^{−1/2}. N = 5, p = 127, K = 2·10⁵: C_{p,0} =
1.44632 ± 5.7·10⁻⁶ against R + 2d₀Mh = 1.446316855; C_{p,π} = 1.42247 ± 6.3·10⁻⁶ against 1.422466214; w₀C₀ + w_πC_π =
1.43334 ± 9.2·10⁻⁶ against R(5) = 43/30. N = 12, p = 631, K = 1.5·10⁵: 1.3347 ± 3.6·10⁻⁵, 1.3346 ± 4.0·10⁻⁵, and the
two-phase combination 1.3346 ± 6.4·10⁻⁵ against R(12) = 1.334632035. |C_{p,π} − R| ≤ 1 holds (0.011 and < 10⁻⁴).

## 3. The shape, read from the certified values

1. **The coefficient is constant on the range, and its constant is attained at N = 13.** The running maximum of
   c_N is c* = 18307/480480 = 0.0381 and is never approached again: the decade maxima fall 0.0381, 0.00766, 0.00647,
   0.00267, 0.00266, 0.00238. On the data, the forms with a > 0 or A > 0 are not needed; the constant-coefficient
   version δ_N ≤ c* B(N^{1/6})⁵ holds for every 6 ≤ N ≤ 10⁷ with the exact rational c*. [finite]
2. **The fifth-power budget is far from binding.** Above 100 the increments themselves are at most 0.2453, while the
   envelope factor is 32, 85 or 149; the decline of the decade maxima of c_N is the growth of B(N^{1/6})⁵ under a flat
   numerator (figure, panel A: the certified maxima track 0.2453/B⁵). The envelope B(X) rises only from 2.74 to 2.89
   over 10² to 10⁷, so F(t) = log B(e^t) is flat (1.0 to 1.06) against the recursion F(t) ≤ 5F(t/6) + … of Lemma 11.1,
   whose right side is 5.0 at t = log 10⁷. [finite]
3. **The explicit return is loose by two orders of magnitude in the exponent.** Theorem 11.2 returns log B(X) ≤
   8.254 (log X)^{0.898}, which is 100 at 10⁷ against the actual 1.06. The looseness comes from the recursion, not
   from c*: the constant K barely moves when c* is replaced by the coefficient above 100. [derived]
4. **In the lifted plane the increment is a small difference of two eigenvalues of opposite sign.** Along the chain
   from 10⁷, ν₊ ≈ +0.004 and ν₋ ≈ −0.003 to −0.004 for N ≥ 10⁴, while their trace c_N is 0.0003 to 0.001; the carried
   prefix and the tail are nearly orthogonal in the calibrated space (|cos| ≤ 0.10 for N ≥ 100), so ω ≈ p₀q₀ and
   ν± ≈ ±√(p₀q₀) + (q₀ + 2z₀)/2. The sign of the increment is decided by the cross term 2M(L)(h(N) − h(L)) (negative
   at 277777, 8388608, 4·10⁶; positive at 10⁷, 214), and the tail square ‖t‖² is 0.06 to 0.14 at every step. [finite]
5. **The Hardy head is of size one and its proved allowance is 10¹² to 10²³.** With K = ⌊(log N)^{3/4}⌋ the head
   H_K(N) lies between 1.07 and 1.31 on the whole chain, while the allowances (9.9) and (9.10) are 10⁶ to 4·10²³; the
   stretched-exponential form of (11.8) with ρ = 7/8 is the shape of that allowance, not of the data. The fixed-K head
   converges as N grows: γ₁ → (1/ζ)′(1) = 1 and γ₂ → 1 − γ_E = 0.4228 (measured 0.99836 and 0.4344 at 10⁷), the jet of
   1/ζ at s = 1, a consequence of the prime number theorem [source: classical]. The increment is carried by the tail:
   T_K(N) − T_K(L) is 0.093 of δ_N = 0.144 at 10⁷, 0.068 of 0.072, 0.074 of 0.0745, and exceeds δ_N where the head
   decreases. At N = 1000 the Hardy energy is spread over very many modes: 800 modes hold 92 percent of R(1000) and the
   coefficients stop decreasing (|γ₈₀₀| > |γ₄₀₀| > |γ₁₀₀|), so no fixed finite jet captures the increment. [finite]
6. **The three-step compression is exact to six digits at large N.** The proved per-step cost 2 log(N/L) = 3.58 is
   realized as 10⁻⁷ to 10⁻⁵ for N ≥ 7716: the modified coefficient ĉ⁽³⁾ is c_N itself, and the transfer constant
   log 11/16 = 0.15 of (2.33) is loose by 10⁵. The reason is the weight 1/(k(k+1)) in (2.24): the interior residual sums
   are at most two in modulus and sit at k ≥ L. [finite; the mechanism is visible in (2.24)]
7. **The circle and field readings agree to the certified radii.** The two-phase isometry (8.4) returns R(N) exactly
   within 10⁻⁵, with the antiperiodic correction (8.6) at 0.011 for N = 5. The same increment is read identically in
   all four planes of (11.7); none of them changes its size. [finite]

## 4. What this says about the final estimate

On 6 ≤ N ≤ 10⁷, certified: the final estimate has the constant-coefficient shape, δ_N ≤ c* B(N^{1/6})⁵ with
c* = 18307/480480, and in fact δ_N ≤ 0.2453 for every N > 100, with the fifth-power envelope never used. The planes agree.
The explicit return of Theorem 11.2 is then available with K(c*) = 8.2542, and it is loose by a factor of about 100 in
log B at 10⁷. What the computation cannot do is extend the range: the hypothesis is a statement about every N, and the
shape above 10⁷ is extrapolation. The data are consistent with the stronger law δ_N = O(1), i.e. R(N) = O(log N), the
weak-Mertens form, which would give the return R(N) ≤ 43/30 + 0.2453·⌈log(N/100)/log 6⌉ for all N if the per-step
bound persisted; that is a conjecture, labelled as such, and the previous notes (CUBE_MOVES.md, GLOBAL_CLOSURE.md) record
that no proof of any per-step bound below the trivial one is available.

Going over the horizon in both directions: the increments were computed down every chain (N → ⌊N/6⌋ → …) and the
decade maxima up the range; the modified, lifted, Hardy and circle readings were taken at both ends of each step with the
same K, the same prime and the same weights, as the manuscript requires. If more is wanted, the certified pass extends to
10⁸ in a few minutes and to 10⁹ in about an hour with the same code; the Hardy pass to 10⁸ in about fifteen minutes.
