# Colours inside the subpower bound: the boundary compression, the colour phase, the packets, the descent (2026-10-06)

Continuation on `received/Michaels_Arithmophysics_Color_and_Descent_Closure_Prompt.pdf`. The papers it carries forward
(*Orthogonal Boundary Compression*, *Prime Color Dual Kernel Bridge*, *Five-Level Architecture and Direct Energy
Descent*, the *Theory Atlas*) are not in this repository; what the directive states of them is verified below where
it can be and taken as stated otherwise. Scripts `rh_colors_descent.py`, `rh_color_phase.py`; figure `COLOR_PHASE.png`.
Dependency statement, once: (Goal) E_P(N) = O(N^ε) ⟺ (G) ⟺ RH by (3); nothing below proves any of them; every
theorem is unconditional and every diagnostic uses no premise.

## 0. The orthogonal boundary compression, verified and proved through the h-form

**Verified.** On the mesh x_{r+1} = x_r + max(1, ⌊N^{−1/6} x_r^{2/3}⌋) with the two-moment boundary coefficients (2):
m = 47, 94, 191, 382 at N = 256, 1024, 4096, 16384 (bounds 13√N + 2 = 210 … 1666); R(N) = E_P + E_res with
z^T B_N p = 0 to 10⁻¹⁵; E_P equals Σ_r [u(x_r) − u(x_{r−1})]²/(x_r − x_{r−1}) to all nine digits with
u(x) = M(x) + x(h(N) − h(x)); E_res = 0.0268, 0.0264, 0.0171, 0.0136 ≤ 2^{−4/3} = 0.397. At N = 16384, m = 382
agrees with the directive's baseline row; E_P = 1.589799 and E_res = 0.013641 here against 1.588343 and 0.015097
there, the sum 1.603440 agreeing to nine digits, so the split differs by 1.5·10⁻³ through a convention of the
projection at shared boundaries, not through the mesh. At large N by the u-form: E_res = 0.0091, 0.0056, 0.0034 at
N = 10⁵, 10⁶, 10⁷, with E_P = 1.6125, 1.7026, 1.8331.

**A one-line proof of the compression.** Since min(a, b) = ∫_0^N 1_{x<a} 1_{x<b} dx, the energy of any coefficient
vector is ∫_0^N (Σ_{n>x} c_n)² dx; for c_n = μ(n)/n this is R(N) = ∫_0^N (h(N) − h(x))² dx, the h-form of the
record, and h(N) − h(x) = u′(x). The two-moment projection replaces u on each mesh interval by its chord, whose
slope is the interval mean of u′ (that is (5)); the residual is u′ minus its interval means, so

    E_P(N) = Σ_r (Δu_r)²/Δx_r,        E_res(N) = Σ_r ∫_{x_{r−1}}^{x_r} ( h(x) − h̄_r )² dx,

the within-interval variance of h (verified to nine digits), orthogonal to the chords by construction. On (a, b] the
step function h oscillates by at most Σ_{a<n≤b} 1/n ≤ ℓ/a, so each term is at most ℓ(ℓ/a)²/4, and with
ℓ ≤ N^{−1/6}a^{2/3} every term is at most N^{−1/2}/4 and E_res ≤ (13√N + 2)/(4√N) ≤ 3.3 uniformly; the directive's
2^{−4/3} is the same bound with a sharper count. For any real |a_n| ≤ 1 the same holds. This is the theorem: the
whole energy is carried, up to an additive constant, by two numbers per mesh interval, the chord of u.

## 1. The colour Gram, and the exact relation that creates the gain

**The Gram.** With μ(n) = μ(n)²(−1)^{ω(n)}, the colour sources c^{[k]} = μ² 1_{ω=k}/n, their projections
p_k = T_N c^{[k]} and G_{kℓ} = p_k^T B_N p_ℓ (N = 16384, colours k = 0..5):

| | 0 | 1 | 2 | 3 | 4 | 5 |
|---|---|---|---|---|---|---|
| 0 | 1.000 | 2.535 | 2.319 | 0.931 | 0.152 | 0.007 |
| 1 | | 521.0 | 1064.9 | 714.8 | 170.7 | 10.2 |
| 2 | | | 2216.2 | 1508.5 | 364.9 | 22.1 |
| 3 | | | | 1042.6 | 256.3 | 15.9 |
| 4 | | | | | 64.3 | 4.1 |
| 5 | | | | | | 0.271 |

E_P = Σ(−1)^{k+ℓ}G_{kℓ} = 1.5898; the same-colour sum ΣG_kk = 3845.4 (0.235N); the mixed-colour signed part is
−3843.8. At N = 1024 and 4096: 276.5 and 1021.2 against 1.4447 and 1.5768. The unprojected Gram gives the same
to three digits. **The matrix is nearly rank one**: G₁₂² = 0.98·G₁₁G₂₂, G₂₃² = 0.98·G₂₂G₃₃. The entries that carry
the cancellation are the adjacent-colour covariances G_{k,k+1}, each of the size of the diagonal.

**The exact relation.** The colour generating function is the Euler product

    Σ_k z^k Σ_{ω(n)=k} μ(n)² n^{−s} = Π_p (1 + z p^{−s}),

z = −1 giving 1/ζ(s), the Möbius source, and z = 1 giving ζ(s)/ζ(2s), the positive control. Adjacent colours are
related by prime insertion, (k+1)·1_{ω=k+1}(n) = Σ_{p|n} 1_{ω=k}(n/p), which is the dilation structure of the
prime packets of Section 3, so the adjacent-colour covariances are the packet cross terms. The rank-one part of G
is the smooth part of each colour: by Selberg–Delange, Σ_{n≤x} μ(n)² z^{ω(n)} = (C(z)/Γ(z)) x (log x)^{z−1}(1 + O(1/log x))
with C(z) = Π_p (1 + z/p)(1 − 1/p)^z, so each colour's tail Σ_{n>x} μ²1_{ω=k}/n is a smooth function of log x of
size ≍ (log log N)^{k−1}/(k−1)!·(…), all colours share the profile, and the parity-signed combination of the smooth
profiles is the z = −1 value, where 1/Γ(z) vanishes. **The gain in the parity direction is the zero of 1/Γ(z) at
z = −1**: the smooth (rank-one) energy of size N is annihilated by parity because the Selberg–Delange main term of
Σ_{n≤x} μ(n) vanishes identically (that is the prime number theorem, M(x) = o(x/(log x)^A) for every A), and what
survives at z = −1 is the fluctuation energy R(N), the zeros.

**The energy as a function of the colour phase** (`rh_color_phase.py`, figure). With F_N(ζ) = Σ_k ζ^k p_k and
E(t) = ‖F_N(e^{it})‖²_B = Σ_{k,ℓ} G_{kℓ} cos((k−ℓ)t):

| N | E(0) measured / 2N(6/π²)² | mean over t = ΣG_kk | E(π/2) / SD | E(3π/4) / SD | E(0.95π) / SD | E(π) = R(N) |
|---|---|---|---|---|---|---|
| 16384 | 12122 / 12110 | 3858 | 1670 / 2875 | 87.5 / 73.9 | 2.26 / 0.37 | 1.6034 |
| 10⁵ | 73929 / 73915 | 22109 | 7712 / 11614 | 269 / 204 | 2.96 / 0.81 | 1.6216 |
| 10⁶ | 739165 / 739151 | 207635 | 56156 / 75696 | 1255 / 923 | 6.49 / 3.02 | 1.7082 |

Here SD is the Selberg–Delange main term |C(z)/(zΓ(z))|² ∫_1^N |(log N)^z − (log x)^z|² dx at z = e^{it}, which is
exact at t = 0 (E(0) = 2N(6/π²)² to 0.1 %: the positive control's energy is a theorem), tracks E(t) within a factor
0.6 to 1.4 over the phase (the O(1/log N) corrections), and vanishes at t = π. **The cost of passing from the
angular average to the parity phase is therefore the whole problem**: the average, 0.235N, is the pole of ζ seen at
t = 0; the parity phase is where every term of the Selberg–Delange expansion vanishes (1/Γ(−1 − j) = 0 for all
j ≥ 0), and E(π) is what is left after the main terms are gone. A nonnegative trigonometric polynomial of degree
K (the number of colours, K ≤ 7 here) can vanish at a point while its mean is anything, so no inequality from the
mean to the point exists; and E(t) near π behaves like (π − t)²·N/(log N)⁴ until the floor R(N) takes over at
π − t ≈ (log N)²·√(R/N). The candidates of the directive are answered by this: (1) the factorisation isolating the
coherent component is the Selberg–Delange main term, annihilated by parity exactly, by the prime number theorem;
(2) the mixed-colour correction against a diagonal budget is the passage from the mean of E(t) to E(π), which has
no inequality; (3) the finite difference across adjacent colours is prime insertion, Section 3.

**The small-prime refinement.** Grading by (2 | n, 3 | n) and the factor count outside {2, 3} gives 19 classes at
N = 16384; the parity energy is still 1.5898; the signed combination within each (2, 3)-signature has energy
14.93, 4.73, 7.17, 2.33 (sum 29.2, the rest is cross-signature). The (0, 0) class is the Möbius source coprime to
6, with energy 14.93 = E(f_{{2,3}}) of Section 3: the dependence on particular primes is exactly the packet
structure below, and the colour count alone does not separate it.

## 2. The prime packets, and what the complete packet does

The identity (9), c_N = Π_{p∈C}(I − D_{p,N}) f_{C,N}, holds to 10⁻¹⁸ at C = {2, 3} and {2, 3, 5} (N = 4096, 16384),
with every cutoff pm ≤ N inside D_p. Expanding the product over the subsets S ⊆ C gives the subset Gram of the
dilated coprime sources T_N D_S f_C, whose signed sum is E_P to six digits. At N = 16384:

| C | E(f_C) | diagonal of the subset Gram | sum of diagonal | signed sum = E_P |
|---|---|---|---|---|
| {2, 3} | 14.94 | 14.93, 7.17, 4.73, 2.33 | 29.2 | 1.5898 |
| {2, 3, 5} | 34.44 | 34.44, 16.21, 10.37, 5.79, 4.82, 2.70, 1.62, 0.66 | 76.6 | 1.5898 |

**Removing primes increases the energy.** Along the removal of the primes p ≤ y (N = 16384):

| y | 1 | 2 | 3 | 5 | 7 | 11 | 20 | 50 | 100 | 300 | 1000 | 3000 | 8192 | 16384 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| E(f_{C_y}) | 1.60 | 4.63 | 14.94 | 34.44 | 64.8 | 97.5 | 209 | 377 | 457 | 423 | 336 | 205 | 58.0 | 1.00 |

The per-prime factor E((I − D_p)v)/E(v) is 2.9, 3.2, 2.3, 1.9, 1.5 for p = 2, 3, 5, 7, 11: each complete packet step
(I − D_p), applied to the actual coprime source, *multiplies* the energy by about 1 + 1/p + 2/√p·(correlation), and
the energy of the coprime sources rises to 457 at y = 100 (the hump of Section 4 of `SUBPOWER_CLOSURE.md`, now in
the full energy) before the large primes bring it back to 1.60 and finally to 1.00 = E(e_1). So a "compulsory
prime-packet response with a uniform estimate" cannot be a contraction per prime: the contraction is at the end of
the packet, among the primes above N/2, and the base object f_C of any proper packet has energy larger than the
target by a factor that grows with |C|. The connection to Section 1 is exact: the colour Gram's cancellation is
the packet product's cancellation, and both are Π_p(1 − p^{−s}) at the parity phase.

## 3. The six-term reconstruction and the descent

The identity (14), μ(n) = Σ_{j=1}^6 (−1)^{j−1} C(6, j) (A^{*j} * 1^{*(j−1)})(n) for n ≤ N = K⁶ with A = μ·1_{[1,K]},
holds exactly (checked at N = 4096, 15625, 16384; it is 1/ζ = (1/ζ)[1 − (1 − Âζ)⁶] with 1 − Âζ supported above K).
The projected energies of the six signed terms at N = 16384 (K = 5) are 52, 371, 996, 1946, 1899, 289, with
E_P = 1.5898: the terms cancel by a factor 3500, and the triangle inequality over them gives a bound of order
10⁴, not of order S(K)^q with q < 6. The multiplicity structure of (14) is C(6, j) with 1^{*(j−1)} convolutions,
whose energies grow like powers of N far above S(K)^j; so the exponent in a descent (13) from (14) is not below
6 at scale N^{1/6}, and qβ = 1 is the critical case the directive names. The record-descent target (12) and the
reduction (10)–(11) belong to a paper not in this repository and are not evaluated here; what the six-term
calculation establishes is that (14) by itself supplies no contraction.

## 4. The signed work

The work identity R(N) − D_N = 2Σ_{n≥2} μ(n)M(n−1)/n holds exactly (verified to 10⁻¹² at N = 10³ … 10⁷), and the
trial inequality E_P ≤ D_N holds at every tested N with a margin that grows: D_N − R(N) = 3.78, 5.06, 6.42, 7.73,
9.01 at N = 10³ … 10⁷, with R/D_N falling from 0.28 to 0.17. The trial inequality says R(N) ≤ (6/π²) log N + O(1):
it is the weak Mertens conjecture in energy form (∫_1^X M²/x² dx ≪ log X), which is stronger than RH and implies
simple zeros; its truth in the data is the sub-random cancellation seen throughout. The work W_N = (R − D_N)/2 is
−1.9, −2.5, −3.2, −3.9, −4.5, about −0.30 log N. No correction converting W_N into a boundary term plus a
nonnegative remainder was found; the only decompositions available are the identities of the record, which do not
change the sign structure.

## 5. The controls (N = 16384, same mesh)

| source | E_P | same-colour | mixed-colour | M(N) |
|---|---|---|---|---|
| Möbius | 1.590 | 3845.4 | −3843.8 | −32 |
| all positive | 12122 | 3845.4 | +8276.8 | 9962 |
| random signs on squarefree support | 5.69 | 7.87 | −2.19 | 126 |
| Möbius signs shuffled within blocks of 64 | 4.97 | 4.46 | +0.51 | −32 |
| Möbius signs shuffled within blocks of 1024 | 6.68 | 5.17 | +1.51 | −32 |

The shuffles preserve the sign counts in each block (M(N) = −32 exactly) and destroy the arithmetic placement; the
energy rises from 1.59 to 5.0 and 6.7, the random level. So the factor three to four by which Möbius sits below
random is in where the signs are, not how many there are; and the colour decomposition shows why the colours
matter for Möbius and not for the controls: the Möbius colours are sign-pure, so the coherent energy (3845) is
concentrated in the same-colour terms and cancelled between colours, while the controls have no coherent part to
cancel (same-colour 4 to 8) and no sub-random structure either.

## 6. What the calculation establishes, and the resume point

Proved: the compression theorem of Section 0 with its h-form proof; the exact Euler-product relation of the colours
and the identification of the colour Gram's rank-one part with the Selberg–Delange main term; the positive control's
energy 2N(6/π²)²; the packet identity and the per-prime growth of the coprime sources; the six-term identity and the
absence of a contraction in its norms. Not proved: (Goal), or any of its equivalents. The exact colour relation that
creates the gain is the zero of 1/Γ(z) at z = −1 in Π_p(1 + zp^{−s}): the parity phase annihilates every main term,
and this is the prime number theorem; the subpower bound is the statement about what remains, which is the energy
of the zeros, the same quantity in every representation of this record.

The resume point the calculation justifies: the profile E(t) near t = π. Its main term is (π − t)²·N/(log N)⁴ up
to constants, its floor is R(N), and the crossover is at π − t ≈ (log N)²√(R(N)/N). The next exact inequality is a
bound on the derivative of the fluctuation part of E(t) at t = π in terms of the zeros: in the Selberg–Delange
expansion with the zeros included, Π_p(1 + zp^{−s}) = ζ(s)^z H(s, z) with H analytic for Re s > ½, and the zeros
enter as branch points of ζ(s)^z that become the poles of 1/ζ exactly at z = −1. The energy of the colour source at
phase t is the energy of ζ^{e^{it}} on the line, and E(π) = O(N^ε) is the statement that the branch points stay on
Re s = ½ as t → π. That is the colour form of the question, and it is the same question.
