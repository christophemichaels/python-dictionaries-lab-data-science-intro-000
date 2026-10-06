# The boundary carry: compression of the source to its boundary charges, and what survives (2026-10-06)

Continuation on `received/Michaels_Arithmophysics_Boundary_Carry_Continuation_Prompt.pdf`, which builds on the
received *Signed Prefix–Tail Comparison* (not in this repository; its stated results are verified below where
they can be) and *Source Truncation and Growing Mode Bounds* (Section 6 of `HORIZON_ROUND_TRIP.md`). Script
`rh_boundary_carry.py`. Dependency statement, once: (G) ⟺ (T) ⟺ (BC)-with-(T) ⟺ RH; nothing below proves any
of them; every theorem is unconditional and every diagnostic uses no premise.

## 0. What the directive supplies, verified

With a_L = M(L), q = c_{≤L} − (a_L/L)e_L and w = c_{>L} + (a_L/L)e_L: q^T B_N q = I_L = R(L) − M(L)²/L,
q^T B_N w = 0, R(N) = I_L + w^T B_N w (verified to 10⁻⁹ and 10⁻¹⁶ at the four (N, J, L) of the supplied table);
r = V_N^T q has r_j = Σ_{k<L} M(k) D_j(k) (Abel, verified to 10⁻¹⁹); the table's R(L)/λ_J and I_L Ξ_{N,J}(L) are
reproduced to all seven digits (5.029876·10⁻⁵, 1.196549·10⁻⁷; 3.460066·10⁻⁶, 2.679424·10⁻⁹; 2.162977·10⁻⁷,
5.528069·10⁻¹¹; 1.060421·10⁻⁵, 2.057884·10⁻⁸); the inner products at (4096, 1, 1024) are 2.509743·10⁻⁹ and
2.491264·10⁻⁹ as reported. The bound (B), ‖Π_J r‖² ≤ C_∂ J²/N² at L = ⌊(N/J)^{5/6}⌋, holds with ratio
10⁻⁴ to 10⁻⁷ for the Möbius source at N = 4096 to 65536, J = 1, 4, 32. The interval identities (I), the formula
for I_r, and the B-orthogonality of disjoint zero-charge remainders are verified to 10⁻¹⁷ (N = 1024, dyadic).

## 1. The strongest new inequality proved: compression above L to √N·J^{5/2} boundary charges

**Theorem P (partition compression, coherent worst case).** Let 1 ≤ J ≤ N, L = ⌊(N/J)^{5/6}⌋, Λ = N^{3/2} J^{−5/2},
and partition (L, N] into consecutive intervals (a_r, b_r] of length ℓ_r = max(1, ⌊Λ/a_r⌋). For every real charge
s with |s_n| ≤ 1 and c_n = s_n/n, with Q_r = Σ_{a_r<n≤b_r} s_n and q^{(r)} = c·1_{(a_r,b_r]} − (Q_r/b_r) e_{b_r},

    ‖ Π_J V_N^T Σ_r q^{(r)} ‖²  ≤  C₂ J²/N²,        C₂ = 64π⁶/378 = 162.8,        m ≤ N^{1/2} J^{5/2}.

Consequently, with W^P_{N,J} the population of the compressed source (the boundary charge M(L)/L at L and the m
increment charges Q_r/b_r at the right endpoints),

    | √S_N(J) − √W^P_{N,J} |  ≤  (√C_∂ + √C₂) J/N,        C_∂ = 2π⁶/315.

Proof. For each interval, (V^T q^{(r)})_j = Σ_{a<k<b}[M_s(k) − M_s(a)] D_j(k) by (I), so by Cauchy–Schwarz with
the weights k(k+1), |(V^T q^{(r)})_j|² ≤ I_r Ξ_{r,j} with I_r = Σ_{a<k<b}[M_s(k) − M_s(a)]²/(k(k+1)) ≤
Σ_{i<ℓ} i²/a² ≤ ℓ³/(3a²) and Ξ_{r,j} = Σ_{a<k<b} k(k+1) D_j(k)². With g(x) = sin(xθ)/x, |g′(x)| ≤ θ³x/3 because
|t cos t − sin t| ≤ t³/3, so |D_j(k)| ≤ (2/√D)θ_j³(k+1)/3 and Ξ_{r,j} ≤ (4/(9D)) θ_j⁶ ℓ b⁴. Hence
√(I_rΞ_{r,j}) ≤ (2/√(27D)) θ_j³ ℓ² b²/a. Summing over the intervals by the triangle inequality (no sign is used),
and using ℓ_r ≤ Λ/a_r and b_r ≤ 2a_r (which holds since a_r ≥ L ≥ √Λ),

    Σ_r ℓ_r² b_r²/a_r ≤ Σ_r ℓ_r·(Λ/a_r)·4a_r = 4Λ(N − L) ≤ 4ΛN,

so |Σ_r (V^T q^{(r)})_j| ≤ (8ΛN/√(27D)) θ_j³, and with Σ_{j≤J}θ_j⁶ ≤ (π⁶/D⁶)(2J)⁷/14 and D ≥ 2N,

    ‖Π_J V^T Σ_r q^{(r)}‖² ≤ (64Λ²N²/(27D)) Σ_{j≤J} θ_j⁶ ≤ (64·128π⁶/(27·14)) Λ²N²J⁷/D⁷ ≤ (64π⁶/378) Λ²J⁷/N⁵ = C₂ J²/N².

The count: m ≤ Σ_r 1 ≤ 1 + Σ_r 2ℓ_r a_r/Λ ≤ 1 + 2N²/(2Λ)·… ≤ N^{1/2}J^{5/2} (each ⌊Λ/a⌋ ≥ Λ/(2a) when Λ ≥ a, and
singletons where Λ < a contribute zero error). The second display is the reverse triangle inequality with (B). ∎

Measured, with the Möbius source, random signs and the all-positive control on the same squarefree support
(N = 4096 … 65536; J = 1, 4, 8): the total error is 10⁻⁷ to 10⁻⁹ of the budget for Möbius and random signs, and
10⁻³ to 10⁻⁴ of it for the positive control, which is the coherent worst case the theorem is built for; the
number of intervals is 30, 61, 125 at J = 1 (bounds 64, 128, 256) and 4255 at (65536, 4). So the constant 163 is
loose by about 10³ even in the worst case, and the compression is exact in the sense that matters: the first J
modes of the source cannot see anything inside the intervals, only the increments of M across them.

This is a new proved statement relative to (B): (B) replaces the prefix below (N/J)^{5/6} by its single charge;
Theorem P replaces the whole range above it by √N·J^{5/2} charges, and the two together say that
**the low-mode population of the Möbius source is determined, to within (√C_∂ + √C₂)² J²/N² in the square, by
M(L) and the increments of M over intervals of length N^{3/2}J^{−5/2}/x.** At J = 1 the intervals at the top have
length √N, the random-walk scale of M.

## 2. The cost accounting

| step | replaced | retained signed data | cost (in the square) | proved for |
|---|---|---|---|---|
| boundary carry (B) | the prefix n ≤ L = (N/J)^{5/6} | M(L) at L | ≤ C_∂ J²/N², measured 10⁻⁴ to 10⁻⁷ of it | all bounded charges |
| Theorem P | the tail (L, N] | m ≤ √N J^{5/2} increments Q_r at b_r | ≤ C₂ J²/N², measured 10⁻³ of it at worst | all bounded charges |
| return (T′) | — | S_N(J) from W^P | \|√S − √W^P\| ≤ (√C_∂ + √C₂) J/N | exact |
| (Return) | — | R(N) from (BC) and (GM) | 9/4·A_η[1 + 2 log(3N)]N^η + C(N^η + ℓ²_N) | as in the directive |

Nothing in the compression costs a power of N: the errors are at the target scale J²/N² with no η, uniformly in
N and J, and they accumulate additively in the square root. The return chain is the directive's, unchanged, with
A_η = (√B_η + √C_∂ + √C₂)² if (BC) is established for W^P.

## 3. The first calculation: W split exactly, and where the saving is

With L = ⌊(N/J)^{5/6}⌋ and the exact split W = boundary² + tail² + mixed (`rh_boundary_carry.py`, Section 2):

| N | J | boundary² | tail² | mixed | W | W·N²/J² | max-\|M\| comparison, loss |
|---|---|---|---|---|---|---|---|
| 4096 | 1 | 1.09·10⁻⁹ | 5.69·10⁻⁹ | +4.98·10⁻⁹ | 1.18·10⁻⁸ | 0.197 | ×3 |
| 4096 | 4 | 7.99·10⁻⁸ | 9.81·10⁻⁸ | −1.23·10⁻⁷ | 5.53·10⁻⁸ | 0.058 | ×185 |
| 16384 | 1 | 3.14·10⁻¹⁰ | 1.63·10⁻⁹ | −1.43·10⁻⁹ | 5.10·10⁻¹⁰ | 0.137 | ×3 |
| 16384 | 4 | 1.34·10⁻⁹ | 3.42·10⁻⁹ | −2.38·10⁻⁹ | 2.39·10⁻⁹ | 0.040 | ×206 |
| 65536 | 1 | 1.55·10⁻¹¹ | 2.17·10⁻¹¹ | −3.66·10⁻¹¹ | 5.20·10⁻¹³ | 0.002 | ×80 |
| 65536 | 4 | 3.95·10⁻¹⁰ | 5.41·10⁻¹⁰ | −9.20·10⁻¹⁰ | 1.56·10⁻¹¹ | 0.004 | ×3217 |
| 65536 | 32 | 1.64·10⁻⁸ | 6.18·10⁻⁸ | −4.85·10⁻⁸ | 2.97·10⁻⁸ | 0.125 | ×2799 |

The answer to the directive's question is unambiguous: **the saving is joint.** Neither the boundary square nor
the tail square is small; the mixed term is of their size and of the opposite sign, and the sum is 1 to 50 % of
either (at (65536, 1): 1.4 %). In the representation (A) this is the statement that
M(L)·sinc(Lθ_j) + Σ_{L<n≤N} μ(n) sin(nθ_j)/n is the smoothed Möbius sum at the scale of the mode, M̃_j(N), which is
smaller than the boundary charge M(L) and than the tail's increment separately. Taking a maximum of |M| across
the interval, as the directive allows, loses the factors in the last column: ×3 at J = 1 (where the window is
monotone and the comparison is nearly sharp) and ×200 to ×3200 at J ≥ 4, where the window oscillates. So no
separate estimate of the boundary or of the tail can reach (BC); only the joint quantity can, and the joint
quantity is M at the mode scales.

## 4. The second calculation: the Gram form of the increments, and the controls

With Theorem P's partition the retained form is (P), W^P = Σ_{r,s} (Q_r Q_s/(x_r x_s)) K_{N,J}(x_r, x_s) plus the
boundary charge, with K_{N,J}(a, b) = (1/D)[H_J(π(a−b)/D) − H_J(π(a+b)/D)]. Measured (N = 65536; N = 10⁶ for the
increments alone):

| N, J | source | W^P·N²/J² | increment coherence (ΣQ_r)²/ΣQ_r² | ΣQ_r² / ((6/π²)(N − L)) |
|---|---|---|---|---|
| 65536, 1 | Möbius / random / positive | 0.0019 / 0.079 / 90 500 | 0.05 / 0.06 / 96 | — |
| 65536, 4 | Möbius / random / positive | 0.0040 / 2.26 / 28 100 | 0.00 / 1.58 / 2484 | — |
| 10⁶, 1 | Möbius / random / positive | 0.073 / 0.101 / 1.4·10⁶ | 0.15 / 0.05 / 30 | 0.83 / 0.69 / 18 000 |
| 10⁶, 4 | Möbius / random / positive | 0.030 / 8.82 / 4.2·10⁵ | 0.05 / 0.91 / 247 | 0.89 / 0.88 / 2 400 |

Reading. (a) For J small the kernel is nearly rank one, K(x_r, x_s)/(x_r x_s) ≈ (4/D)Σ_{j≤J}θ_j² sinc(x_rθ_j)
sinc(x_sθ_j), so W^P ≈ (4/D)(Σθ_j²)(M(L)/L·L + Σ_r Q_r·w_r)² with slowly varying weights w_r: the quadratic form of
the increments is the square of their weighted sum, and (BC) for it is |M(L) + Σ_r Q_r w_r| ≪ N^{1/2+η}, which is
M̃(N) ≪ N^{1/2+η} at the scale of the mode. (b) The increments' mean square is the squarefree density times the
length, ΣQ_r² ≈ (0.83 to 0.89)(6/π²)(N − L) for Möbius and random alike, so the candidate independent inequality
"square-root cancellation of the increments", (Σ_r Q_r w_r)² ≤ C N^η Σ_r Q_r², is equivalent in content to
M(N) − M(L) ≪ N^{1/2+η}: it is the target, not an input. (c) The controls discriminate: the Möbius increments are
anti-correlated at these spacings (coherence 0.00 to 0.26, below the random level 1), the random control sits at 1,
the positive control at m. **The arithmetic content of (BC) in this representation is exactly that the Möbius
increments over intervals of length N^{3/2}J^{−5/2}/x do not add up coherently, to within N^η of square-root
cancellation; no independent inequality for that was found, and the data say Möbius does better than random.**

A recurrence of the directive's form did not emerge: the compression produces a representation of S_N(J) by
≤ √N J^{5/2} numbers, and bounded-charge inputs control the compression error but not the retained form; the
retained form has no bounded-charge bound below W^P ≤ (positive control) ≍ (J²/N²)·N/J.

## 5. The arithmetic input, located once

The compression (B) and Theorem P use |μ| ≤ 1 only; they are as good for the positive control as for Möbius,
and the positive control shows they are sharp up to a constant. The actual Möbius source improves on the
bounded-charge controls at exactly one place: the value of the retained quadratic form W^P, i.e. the weighted sum
of the increments of M. Everything the horizon framework does before that point is lossless to within J²/N²; the
step after it is the Mertens function at the scales N^{3/2}J^{−5/2}/x ≤ x ≤ N.

## 6. Status of (BC), and the resume note

(BC) is not closed. The newly proved statement is Theorem P and its consequence: (BC) is equivalent, up to the
cost (√C_∂ + √C₂)²J²/N², to the same bound for the compressed form W^P on 1 ≤ J < J_η, a quadratic form in at
most √N·J^{5/2} increments of M. The surviving term is that quadratic form; its most informative diagnostic is the
coherence of the increments, 0.00 to 0.26 for Möbius against 1 for random signs; the most informative independent
inequality attempted is square-root cancellation of the weighted increments, which is the target itself.

The next exact calculation the compression points to: the pair correlation of the increments,
⟨Q_r Q_s⟩ for |x_r − x_s| ≳ N^{3/2}J^{−5/2}/x, is the quantity that decides W^P, and in the explicit formula it is
Σ_{ρ,ρ′} (x_r^ρ − x_{r−1}^ρ)(x_s^{ρ′} − x_{s−1}^{ρ′})/(ρζ′(ρ)ρ′ζ′(ρ′)), whose diagonal ρ = ρ′ gives the random-walk
level ΣQ_r² and whose off-diagonal is where the sub-random coherence of Möbius comes from. That is the explicit
form of what is being asked, and it has no premise-free estimate.
