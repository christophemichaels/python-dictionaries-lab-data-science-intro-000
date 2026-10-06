# Decay inside the growth: the parity defect, the packet remainder, and the finite descent (2026-10-06)

Continuation on `received/Michaels_Arithmophysics_Parity_Decay_and_Finite_Descent_Prompt.pdf`. Scripts
`rh_parity_decay.py` (every finite object of the directive) and `rh_parity_decay_fig.py` (figure `PARITY_DECAY.png`).
Dependency statement, once: E_P(N) = O(N^ε), δ_N = O(N^{−1+ε}), [I_C(N)]_+ = O(N^ε) and the retained-record target
(24) are each equivalent to (G), hence to RH; nothing below proves any of them; every uniform statement is marked as
such and every number is a measurement at the horizon named.

## 0. Reconciliation of the common definitions

With the step computed as the integer sixth root of ⌊x_r⁴/N⌋ and the projection (3), N = 16384 gives m = 382,
R = 1.6034402209313035, E_P = 1.5883431079039456, E_res = 0.0150971130273582, z^T B p = −2·10⁻¹⁷: the directive's
values to fourteen digits. The earlier split differing by 1.5·10⁻³ came from a floating-point mesh (x^{2/3} rounded
below at some cubes); the map and indexing were the same. The uniform cap 3.3 of `COLORS_AND_DESCENT.md` is a bound
on the same quantity from a cruder count; 2^{−4/3} is the sharper proved cap and is used here.

## 1. Strongest new proved inequality

No new arithmetic inequality for E_P, δ_N, I_C or the retained records is obtained. The strongest rigorous results
of this pass are the following, each with its proof or its exact check.

**(i) The relative defect is the energy over a linear background, exactly.** A_N/N = 0.3745, 0.3700, 0.3696, 0.3696
at N = 10³, 16384, 10⁵, 10⁶ … 8.4·10⁶, converging to the Selberg–Delange value for the even and odd colour sources:
each of ‖e_N‖², ‖o_N‖² is (1 + o(1))·N(6/π²)²/2, so A_N = (1 + o(1))·0.3695 N and δ_N = (2.706 + o(1))·E_P(N)/N.
The proof of (8) in the directive is verified (c_sf = 2 − π²/6 = 0.355; the measured constant is above c²_sf/2 =
0.063 and below 2); (9) holds, and the decay variable is E_P/N up to the constant 2.706.

**(ii) The Gram split (11)–(12) is exact.** At N = 16384, G − qqᵀ/N − C = 0 to 10⁻¹², C ⪰ 0 (smallest eigenvalue
1.4·10⁻³), and E_P = M(N)²/N + sᵀCs = 0.0625 + 1.5258. So the centred parity energy sᵀCs carries all of E_P but the
terminal.

**(iii) The packet remainder, its kernel, and its sign.** For C = {2, 3} at N = 16384 the four components of
(17)–(18) are 2.053503, −0.537728, −0.024996, 0.112661, summing to R(N) = 1.603440 (the directive's values), and
I_C(N) = R(N) − (2/3)Σ_{m≤N/6,(m,6)=1} μ(m)²/m = −0.450063. The kernel formula (19) is verified entry by entry
(K(5,7) = 0.057143, K(7,25) = −0.016190, K(5,29) = −0.001149, K(5,31) = 0 = K(1,6) = K(1,7)), k_C(0) = κ_C = 2/3, and
k_C vanishes outside [−log 6, log 6] (10⁻¹⁷ at log 6 and at 1.9). **The scan to 10⁷** (`rh_parity_decay.py`,
Section 4): I_{2,3}(N) > 0 at exactly 977 values of N, the last being N = 2837 (I = 0.0108); for every N from 2838
to 10⁷ it is negative, with values −0.0285, −0.372, −0.799, −1.179, −1.517 at 10³ … 10⁷ and I_C(N) + 0.19 log N
between 1.38 and 1.55 over 10⁴ … 10⁷. So [I_{2,3}(N)]_+ = 0 on (2837, 10⁷], at every energy record on that range
(−0.39 at 24137, −1.38 at 6481601), and the second target of (24) holds there with C_ε = 0. This is a finite
statement; its uniform version is the inequality in Section 5.

**(iv) The record accounting (22), verified at X = 10⁷, and the structure of the retained set.** Through 10⁷ there
are 213 strict energy records, every one a first passage of |M| to a new height (checked), with heights 2 … 1060 at
N = 5 … 6481601. The retained set (21) is **not empty**: 179 of the 212 noninitial records are retained, starting at
N = 24153 (the records through 16384 are 1, 5, 13, 31 and none is retained, as the directive reports). The reason is
structural: at a first passage N to height h, |M(N − 1)| = h − 1, so whenever h − 1 was itself first reached inside the
same climb, N − τ_{h−1} is 1, 2, 3, … (50 of the records have gap exactly 1), while the threshold √N/(1 + S(N^{1/6}))^{5/2}
is 16.8 at N = 24153 and 208.8 at N = 6.48·10⁶. Hence **every top of a climb that exceeds the previous maximum of |M|
by two or more is retained**, and the retained set contains essentially all new-maximum points of |M|. On these,
T_ret(10⁷) = max(1, max h²/N) = 1 (the maximum of h²/N over retained records is 0.2143 at N = 24185), S(10⁷) = 1.8904,
and the bound (22) evaluates to (4 log X + 4)·1 + 2(log X)²(1 + S(14))⁵ = 68.5 + 77250.5 with S(14) = R(13) = 1.7192:
valid and loose by 4·10⁴. The accounting is correct as proved in the directive; what it needs is (24) at the
retained records, which by the structure just shown is M(N)² ≤ C_ε N^{1+ε}(1 + S(N^{1/6}))⁵ at essentially every new
maximum of |M|, and M(x) = O(x^{1/2+ε}) holds for all x if and only if it holds at the new maxima. So the reduction
to retained records removes the non-record horizons and the slow climbs, and leaves the Mertens bound at its own
maxima, with a free factor (1 + S(N^{1/6}))⁵ that the induction itself must keep at N^{o(1)}. It does not weaken the
arithmetic input.

**(v) The profile at parity, exact finite values.** E′_N(π) = 0 within numerical error and the curvature (30) equals
its finite-difference value: E″_N(π) = 46.04, 90.58, 325.6 at N = 16384, 10⁵, 10⁶ (positive: parity is a local
minimum of the profile at all three, measured, not assumed); the width to 2E(π) is π − t = 0.228, 0.170, 0.099,
within 10 % of √(2E(π)/E″(π)) = 0.263, 0.188, 0.102.

**(vi) The conditional-shuffle expectation (31), exact.** At N = 16384: E_P(ā) = 1.473291335, V_N = 0.561934217,
expectation 2.035225552 (the directive's digits), Monte Carlo over 200 draws 2.047 ± 0.022, actual E_P = 1.588343108,
Δ_N = −0.447, V_N ≤ 2D_N = 13.9. The actual value sits 1.5 standard deviations of a single draw below the shuffle
mean: the dyadic conditional shuffle keeps the block sums of μ, hence M at every power of two, and with them most
of the energy; the remaining signed combination E_P(ā) + Δ_N is the energy of the dyadic coarse-graining of μ, which
is Σ_j (M(2^{j+1}) − M(2^j))²/2^j up to bounded factors, the dyadic form of R itself.

**(vii) The controls.** At N = 16384: all positive 32757.5 (2N − H_N); squarefree positive 12122.1 (2N(6/π²)²);
random signs on squarefree support 3.52; conditional shuffle 1.59 (one draw; expectation 2.035); a_{2,3}(n) =
μ²(−1)^{#{p∈{2,3}: p|n}} 336.85 against 2α²_C N = 336.4 (the linear law with α_C = (6/π²)(1/3)(1/2), verified);
Möbius 1.588.

## 2. Complete return

There is nothing new to carry. The chains (4), (8), (18), (22) are verified as transfers: E_P + E_res = R exactly;
δ_N·A_N = E_P exactly; κ_C D^{cop} + I_C = R exactly; S(10⁷) ≤ 77319 by (22). Each is an identity or a valid bound
whose arithmetic input is the target itself under another name: E_P subpower, δ_N N^{1−ε} bounded, [I_C]_+ subpower,
or M(N)² ≤ N^{1+ε}(1 + S(N^{1/6}))⁵ at the new maxima of |M|. The contraction parameter of (25) is γβ = 5/6 < 1 and
the finite iteration (27) is correct; the recurrence is not established because its input (24) is not. No step in
this pass loses a fixed power; every step that is proved costs an additive constant (compression), a logarithm
(packet diagonal κ_C(1 + log N), shuffle variance V_N ≤ 2D_N) or a polynomial in S(N^{1/6}) (the accounting), and the
step that is not proved is the one that would supply the power.

## 3. Growth-and-decay finding

Across 27 horizons from 1024 to 8388608 (ratio √2, with 10⁶ added), N·δ_N runs from 3.86 to 4.66 with local
exponents between −0.134 and +0.247 and mean +0.03 (figure, panel A). The data do not support a power law in
either direction: they are consistent with N·δ_N = 2.706·E_P(N) and E_P(N) ≈ 1.2 + 0.03 log N, the logarithmic
growth of the record with the slope of the order of Ng's constant, fluctuating by ±0.05 from horizon to horizon.
Separately: the relative parity decay is δ_N ≈ (3.2 + 0.09 log N)/N, a 1/N law with a logarithmic drift; the
absolute energy grows slowly and non-monotonically (E_P = 1.742 at 2097152 and 1.689 at 2965820); the phase width
shrinks like (log N)/√N·c with the curvature growing like N/(log N)^{2 to 4}, so the parity phase becomes a sharper
minimum of the profile while its floor rises slowly. The two terms of (10): the directional term is 93 to 99.99 %
of E_P at every horizon; the amplitude term (a − b)² is erratic between 4·10⁻⁶ and 0.17 and is not a monotone
function of M(N)²/N (0.098 against 0.016 at 1024; 0.015 against 0.045 at 10⁶). No contracting recurrence emerges:
E_P(N) − E_P(N/6) is of order 0.05, additive, not a fraction of E_P. A scaling law "δ_N ∝ N^{−1} log N" fits the data
with relative error under 5 % and is the weak Mertens law in energy form, a conjecture stronger than RH; it is not
a tool for descent because it does not relate scales multiplicatively.

## 4. Reproducible package

`rh_parity_decay.py` computes, in one run from a sieve to 10⁷: the exact mesh and projection with the reconciliation;
the parity table with both terms of (10) from the slope representation (no subtraction of large energies); the
Gram split; the packet components, the kernel checks and the sign scan of I_{2,3}(N); the records, first passages
and the retained set; the curvature and width; the shuffle expectation with its Monte Carlo; the controls.
`rh_parity_decay_fig.py` draws the figure from the printed rows and the remainder scan. Measured numbers are those
at a stated N; the uniform statements are (4), (8), (22) and the identities, all from the directive and verified;
nothing measured here is a uniform bound.

## 5. The surviving expression and the next attempt

The fully specified surviving expression is the {2, 3} packet remainder,

    I_{2,3}(N) = R(N) − (2/3) Σ_{m≤N/6, (m,6)=1} μ(m)²/m = 2 Σ_{m<n≤N/6, n<6m, (mn,6)=1} μ(m)μ(n) K_C(m,n) + 2⟨g_N, t_N⟩_B + ‖t_N‖²_B,

and the inequality still needed is [I_{2,3}(N)]_+ ≤ C_ε N^ε (1 + S(N^{1/6}))⁵ at the retained records, which the
structure of Section 1(iv) shows to be the Mertens bound at the new maxima of |M|. The estimate tested in this pass
is the sign scan: it establishes [I_{2,3}(N)]_+ = 0 on 2838 ≤ N ≤ 10⁷ and nothing beyond 10⁷. The parity alignment
identifies the required size exactly: I_C(N) must stay below the affordable diagonal by a margin, and the data say
it does so by about 0.19 log N + 1.4.

The next attempt the calculation justifies is the conjecture it suggests, stated as a target and not as a result:

    I_{2,3}(N) ≤ 0 for all N ≥ 2838,   i.e.   R(N) ≤ (2/3) Σ_{m≤N/6, (m,6)=1} μ(m)²/m ≈ 0.203 log N + O(1).

It is a weak-Mertens inequality with an explicit constant (the conjectural slope of R is β = 0.0288, well inside
0.203), true for every N up to 10⁷, and its proof would require exactly the control of the signed sums in (17) that
no independent input supplies. The route the directive names for it, a collective estimate of the completed
coprime packets against their unfinished parents, has its cross term measured at −0.025 and its unfinished square
at 0.113 at N = 16384: both are small, and the whole of the remainder's negativity is the local overlap −0.538,
the signed interaction of the coprime Möbius parents within the packet range n < 6m. That is the term to estimate,
with the kernel (19): it is a bilinear form in μ restricted to coprime-to-6 arguments at ratio below 6, with a
positive-definite kernel whose Fourier transform is the prime filter Π_{p∈{2,3}}|1 − p^{−1/2+iξ}|²/(ξ² + ¼). The
filter lies between fixed constants, so the estimate is a statement about the Möbius source coprime to 6 at bounded
ratios, and the arithmetic that would make it negative is the same arithmetic as everywhere in this record.
