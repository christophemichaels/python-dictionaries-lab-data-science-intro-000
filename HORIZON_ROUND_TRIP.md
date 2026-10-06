# The horizon round trip: what the transfer costs, and the signed correlations that remain (2026-10-06)

Proposed next research section, on the directive `received/Michaels_Arithmophysics_Horizon_Round_Trip_Prompt.pdf`.
It continues Sections 8–9 of `SUBPOWER_CLOSURE.md`. The paper *Arithmophysics: Prime Modes and Complete Packet
Cancellation* named there has not been received; its counterexamples to a uniform per-prime contraction and to a
generic positive log-weight anticommutator are retained by reference and not restated. Script `rh_round_trip.py`,
figure `HORIZON_ROUND_TRIP.png`, data `data/round_trip.json`. The convolution identities and their collapse are
treated as checkpoints and were not re-evaluated.

Dependency statement, once: (G) ⟺ (T) ⟺ (3) ⟺ RH. Nothing below proves any of them. Every estimate marked proved
is unconditional; the diagnostics use no premise.

## 1. The strongest inequality proved, and the one attempted

**No arithmetic improvement is proved.** The attempted inequality is (10) for the unweighted source,

    O^{(1)}_{N,J} = 2 Σ_{a<b≤N} μ(a)μ(b) K_{N,J}(a,b)/(ab) ≤ C_ε J² N^{−2+ε},        1 ≤ J < J₀(N) = ⌈N/log²(eN)⌉,

and the obstruction is stated at the end of this section. What is proved is the following.

**(i) The supplied lemmas hold, with the measured slack.** The two-sided bound (2): the true multiplier
1 + log(λ_1/λ_N) is 15.30, 16.69, 18.07 at N = 1000, 2000, 4000 against the certified 1 + 2 log(3N) = 17.01,
18.40, 19.79. The localization (6): on vectors supported on n ≥ 21 the inverse D_N^{−1} has energy norm 0.5321 at all
three N, against A_21 = 0.9676; the full inverse has 25.45, carried by n = 2, 3. The diagonal bounds (8) and (9) hold
with ratios 0.60–0.70 and 0.01–0.12. The kernel (13) agrees with the direct sum to all printed digits. The
high-mode bound holds with enormous slack: Σ_{j>J₀} λ_j F_j² is 1.31, 1.39, 1.56 against (9/4)C_log ℓ⁴ = 6677, 9345,
12741 (C_log = 0.7589).

**(ii) The exact cost accounting of the localization (proved; it is an identity).** For any L ≥ 2, since
1/max(a, b) = 1/b whenever a < L ≤ b,

    R(N) = R(L−1) + 2 M(L−1) (h(N) − h(L−1)) + R_tail(N),        R_tail(N) = Σ_{a,b≥L} μ(a)μ(b)/max(a,b),

with h(x) = Σ_{n≤x} μ(n)/n. At L = 21: R(20) = 1.601733, M(20) = −3, h(20) = −0.050716, so the head and the cross
term together tend to 1.601733 − 6·0.050716 = 1.2974 as N → ∞ (h(N) → 0 by the prime number theorem), and every bit
of the growth of R(N) is in the tail: R_tail = 0.188, 0.272, 0.321, 0.412, 0.540 at N = 10³ … 10⁷. On that tail the
return trip costs at most A_21 < 0.968 (measured 0.532). So the large inverse constant is a fixed-index effect of the
sources n ≤ 20, which carry a bounded, explicitly known energy, and the scale-dependent part of the problem lives
entirely where the transfer costs less than one.

**(iii) The unconditional baseline and the trivial range, for the record.** By partial summation against the
window g_j(x) = sin(xθ_j)/x, with |g_j′| ≤ θ_j² on x ≤ 1/θ_j and ≤ 2θ_j/x beyond, and |M(x)| ≤ B(x) = x e^{−c√log x}
from the zero-free region,

    |b_{j,N}| ≪ (2/√D) θ_j B(N) ≪ j N^{−1/2} e^{−c√log N},        ‖Π_J b_N‖² ≪ J³ N^{−1} e^{−2c√log N}.

(The directive quotes the paper's baseline as J² N^{−1} e^{−c√log N}; the J-exponent differs from what partial
summation gives here and is immaterial to the N-exponent.) By Parseval ‖Π_J b_N‖² ≤ Σ μ²/n² ≤ Z = 15/π², so (3)
holds trivially for J ≥ √Z·N^{1−ε/2}, which contains J₀ for large N. Between the trivial range and J = 1 nothing
unconditional reaches N^{−2+ε}: at J = 1 the baseline is N^{−1}e^{−2c√log N} against the target N^{−2+ε}, the gap
being the full factor N^{1−ε}.

**(iv) The obstruction.** With a, b both below N/J the kernel is K_{N,J}(a,b) ≈ (4/D)Σ_{j≤J}(aθ_j)(bθ_j) ∝ ab·J³/N³,
and with the Möbius weights 1/(ab) the off-diagonal is the square of the smoothed Möbius sum at scale N/J, up to the
diagonal; at every J the quantity to bound is the low-frequency energy of the source μ(n)/n in the band |θ| ≲ J/N.
The independent tools behave as follows. Parseval and the large sieve give the total energy Σ μ²/n² and cannot
resolve a frequency band below the average share, so they return O(1), no J-dependence, and cannot distinguish the
random control from the Möbius source. Bilinear (Type I/II) bounds give nothing for a kernel that is smooth at low
frequency. The zero-free region enters only through M(x) and gives (iii). The diagnostics of Section 3 show that
the arithmetic information needed is specifically that the Möbius source is more cancelled than random at low
frequency by a factor of about twenty, and no estimate independent of (G) was found that sees this.

## 2. What each transfer costs and what the round trip preserves

| step | representation | cost | nature |
|---|---|---|---|
| (2), K(z) → ‖z‖²_Λ | modes | factor ≤ 1 + 2 log(3N) (true 2 log N + 1.5) | scale-dependent, logarithmic |
| (5), the inverse (L_N − I)^{−1} | energy norm | 25.45 (finite-dimensional, N-independent to four decimals); certified < 216 | geometric, fixed indices n = 2, 3 |
| (6), the inverse on n ≥ 21 | energy norm | ≤ A_21 = 0.968 (measured 0.532) | geometric, no scale dependence |
| the localization at L = 21 | arithmetic | head R(20) = 1.6017, cross → −0.3043; both exact; the tail carries all growth | fixed-index, exact |
| (9), the log weight | source | factor ℓ²_N = log²(eN) on the diagonal | scale-dependent, logarithmic |
| the high modes j > J₀ | modes | ≤ (9/4)C_log ℓ⁴_N (measured 1.3–1.6) | scale-dependent, logarithmic |
| the completion at P ≥ √N (Section 8) | horizons | exact; sign-free recombination costs Σ_{p>√N} p^{−1/2} ≍ √N/log N | the only power-of-N loss, and it is sign-free |

The round trip preserves the energy identity exactly, the mode-scale correspondence, and the localization
identity. It loses nothing but logarithms and fixed constants on every step that keeps the signs. The one step
that loses a power of N is the recombination of child horizons without their signs (the triangle inequality in
the completion formula, Proposition 4 of `SUBPOWER_CLOSURE.md`), and that loss, √N/log N in the norm, is exactly
the gap between the baseline and the target. So the cost is not an arithmetic obstruction; the arithmetic
obstruction is that the signed recombination has no independent estimate.

On the recurrence (14): the transfers do not produce one. The completion at P = √N expresses b_N through
b_{smooth} and the complete children b_{⌊N/p⌋}, p > √N, with coefficients whose sign-free total
Σ_{p>√N} p^{−1/2}√(R(⌊N/p⌋)/R) is of order √N/log N, far from a contraction; any a_{N,k} summing to less than one
would have to come from the signs, which is the missing estimate itself. Inserting the target for the children is
excluded by the directive and would be circular.

## 3. Diagnostics (N = 10³ … 10⁷, all J < J₀(N), by one FFT per source; `data/round_trip.json`)

Normalised by N²/J², for the unweighted source (w = 1):

| N | J₀ | J | D·N²/J² | O·N²/J² | total·N²/J² |
|---|---|---|---|---|---|
| 10⁴ | 96 | 1 / 10 / 95 | 2.32 / 2.94 / 3.05 | −2.20 / −2.87 / −2.79 | 0.117 / 0.071 / 0.262 |
| 10⁵ | 639 | 1 / 10 / 100 / 638 | 2.32 / 2.94 / 3.00 / 3.04 | −2.27 / −2.87 / −2.86 / −2.74 | 0.055 / 0.065 / 0.138 / 0.296 |
| 10⁶ | 4556 | 1 / 10 / 100 / 1000 / 4555 | 2.32 / 2.94 / 2.99 / 3.01 / 3.03 | −2.24 / −2.86 / −2.87 / −2.86 / −2.69 | 0.079 / 0.075 / 0.123 / 0.144 / 0.343 |
| 10⁷ | 34127 | 1 / 10 / 100 / 1000 / 10⁴ / 34126 | 2.32 / 2.94 / 2.99 / 3.00 / 3.01 / 3.02 | −2.11 / −2.81 / −2.88 / −2.87 / −2.86 / −2.78 | 0.208 / 0.130 / 0.117 / 0.128 / 0.144 / 0.245 |

Controls at N = 10⁵, total·N²/J² at J = 1, 10, 100, 638: random signs 2.14, 3.16, 2.59, 1.37 (five draws; the
expectation is the diagonal); all-positive signs 1.4·10⁵, 1.8·10⁴, 1.8·10³, 2.9·10² (growing like N/J).

What the table says. (a) The diagonal is what (8) predicts, 2.32 at J = 1 rising to 3.0, below the bound π²/2.
(b) **The signed correlations are negative at every N and every J tested, and cancel the diagonal almost
completely**: O/D is between −0.90 and −0.97. The Möbius population is 3 to 40 times below the diagonal, which is
the random level; the random control sits at the diagonal, the positive control N/J above it. (c) The residual,
total·N²/J², stays between 0.03 and 0.35 over four decades of N and all J < J₀: (3) holds on the whole tested range
with ε = 0 and C = 0.35. (d) At large J the residual settles near 0.13–0.14 while the diagonal is 3.0, a ratio
≈ 0.045; this is β/(6/π²) = 0.0288/0.608 = 0.047, Ng's conditional constant over the squarefree density, which is
what the mean square of M(x)/√x against its random-walk level should be under the log law. The diagnostics are
consistent with that law and prove nothing.

For the log weight the picture is the same with logarithms: D^{(log)}·N²/J² is 54 to 453 (growing like ℓ²_N),
O^{(log)} ≈ −0.9 D^{(log)}, and the residual 1.2 to 48, growing like ℓ²_N times the same fluctuating factor.

**The dyadic Gram blocks (11) at N = 2000.** F_N = Σ_r q_r with q_r the block of dilations d ∈ [2^r, 2^{r+1}),
r = 1..10 (verified: the blocks sum to F_N to 10⁻¹⁶). At J = 10: ‖Π_J F‖² = 8.20·10⁻⁵, the diagonal blocks sum to
3.64·10⁻⁵ and the off-diagonal to +4.56·10⁻⁵; at J = 50: 2.02·10⁻³ = 1.04·10⁻³ + 0.97·10⁻³. The cross-block
interactions add energy on the whole (positive excess), with suppression concentrated between adjacent blocks
d ∈ [8,16), [16,32), [32,64) (the pairs (4,5) and (3,4) are the largest negative entries) and the first block
d ∈ {2, 3} slightly negative against everything. Panel B. No family of cross-block terms with a common sign or a
common bound emerged; the block structure organises the log-weighted energy, it does not estimate it.

**Mechanism discrimination.** A geometric estimate (Parseval, large sieve, the diagonal bound) returns the random
level or above; the positive control shows such an estimate is sharp for sign-free sources. The Möbius source is
twenty times below that level. The mechanism that would close must certify a negative O of size at least
(1 − N^{−1+ε}·C)·D: not a bound on the signed correlations but a near-equality, O = −D + O(J²N^{−2+ε}). No
independent inequality tested here produces a negative upper bound on O at all.

## 4. The return, with all terms, and the sharpest remaining question

If (10) held for w = 1 on 1 ≤ J < J₀ with exponent ε/2, then with (8), ‖Π_J b_N‖² ≤ (π²/2 + C) J² N^{−2+ε/2}, so
λ_J ‖Π_J b_N‖² ≤ (D²/(4(2J−1)²))·(π²/2 + C) J² N^{−2+ε/2}·(π²/4)… ≤ C′ N^{ε/2} on that range (using
λ_J ≤ D²/(4(2J−1)²)·(π²/4)... i.e. λ_J ≤ π²D²/(16(2J−1)²)); for J ≥ J₀, λ_J ≤ λ_{J₀} ≤ (9/4)ℓ⁴_N and ‖Π_J b_N‖² ≤ Z, so
λ_J‖Π_J b_N‖² ≤ (9/4)Zℓ⁴_N. Hence K(b_N) ≤ max(C′N^{ε/2}, (9/4)Zℓ⁴_N) and by (2), R(N) ≤ [1 + 2 log(3N)]·K(b_N) ≤
C_ε N^ε: this is (G), with the terminal term inside R(N) throughout, no omitted range (the trivial range J ≥ J₀ is
covered by Parseval and the high-mode bound), and the return through (2) only, without the inverse (5). The
log-weighted route returns through (9), the high-mode bound and (5) at the extra cost ℓ²_N·25.45², also subpower.

The sharpest remaining question is (10) for w = 1 on 1 ≤ J < N^{1−ε/2}, equivalently

    ‖Π_J b_N‖² ≤ C_ε J² N^{−2+ε}        (1 ≤ J < N^{1−ε/2}),

known trivially above that range, known nowhere below it, observed with C = 0.35 and ε = 0 for N ≤ 10⁷ and
J < J₀(N). The version with ε = 0 is presumably false in the limit (it would give M(x) = O(√x), against the expected
unboundedness of M(x)/√x), so the N^ε is not decorative. What would answer it is a negative upper bound on the
signed correlations O^{(1)}_{N,J}, of the form O ≤ −D + C J²N^{−2+ε}, proved from an input independent of (G); the
diagnostics say this is the true size, and this pass found no such input.
