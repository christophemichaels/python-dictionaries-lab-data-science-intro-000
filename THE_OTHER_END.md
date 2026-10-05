# The other end: what comes out when every zero is forced onto the point (2026-10-05)

Assume it exactly: every zero reflects onto itself, (1 − ρ̄)/ρ = 1 for all ρ, no oscillation around the
point. That is the Riemann Hypothesis taken as given. Push it through the explicit formula and something
does come out the other end, on the integers, and all of it is a theorem under that assumption. This note
lists what comes out, and checks the four most tangible items numerically (`rh_other_end.py`, 12 seconds
to x = 10⁷).

## 1. What comes out

Each line is a theorem of the form "if every zero is on the point, then …", or an equivalence "every zero
is on the point if and only if …", with its author.

| what comes out the other end | form | author |
|---|---|---|
| the primes are counted by li(x) to within √x log x / 8π, for every x ≥ 2657 | if RH then | Schoenfeld 1976 |
| the prime powers sum to x within √x log² x / 8π, for every x ≥ 73.2 | if RH then | Schoenfeld 1976 |
| π(x) = li(x) + O(√x log x) | if and only if | von Koch 1901 |
| M(x) = O(x^{1/2+ε}): the Möbius signs cancel at the square-root scale | if and only if | Littlewood |
| consecutive primes are never further apart than C √p log p | if RH then | Cramér 1920 |
| the sum of the divisors of every n > 5040 is below e^γ n log log n | if and only if | Robin 1984 |
| σ(n) ≤ H_n + e^{H_n} log H_n for every n ≥ 1, equality only at n = 1 | if and only if | Lagarias 2002 |
| the Farey fractions of order N deviate from equal spacing by O(N^{1/2+ε}) in total | if and only if | Franel–Landau 1924 |
| the determinant of the Redheffer 0–1 matrix is O(n^{1/2+ε}) | if and only if | Redheffer 1977 |
| Li's coefficients λₙ are nonnegative for every n, and grow like (n/2) log n | if and only if; growth if RH | Li 1997; Lagarias 2007 |
| the Weil floor λ(a) is positive at every support a, with no null vector | if and only if; strictness if RH | Weil; Prop. 2.11 of Across the Horizon |
| the potential set collapses to the single measure ν_ζ in the limit | if and only if | atlas, potential set |
| the Möbius state has only phase modes: no amplifying mode, subexponential power | if and only if | Crossing 3 |
| the horn has no missing ring, and its closed form is a cylinder with GUE-repelling rings | if and only if (rings); GUE conjectural | horn note §6 |
| the pair correlation of the zeros is the GUE law, in Montgomery's range | if RH then | Montgomery 1973 |
| the zeros are the spectrum of multiplication by γ on L²(ν_ζ), and the Weil form is that norm | if and only if | Connes 1999 (trace form) |

The first seven are statements about the integers alone, with no zero in them. That is what "comes out
the other end": force every zero onto the point and the primes, the prime powers, the Möbius signs, the
prime gaps and the divisor sums each receive a law with explicit constants.

## 2. Four of them checked

`python3 rh_other_end.py 1e7 1e6`:

| consequence | checked to | result | RH says |
|---|---|---|---|
| \|π(x) − li(x)\| ≤ √x log x / 8π for x ≥ 2657 | x ≤ 10⁷ | max ratio 0.977, attained near the start; over x ≥ 10⁵ the max ratio is 0.369; at x = 10⁷ the bound allows 2,028 and the error is 339.5 | < 1 |
| \|ψ(x) − x\| ≤ √x log² x / 8π for x ≥ 74 | x ≤ 10⁷ | max ratio 0.804 | < 1 |
| \|M(x)\| / √x | 100 ≤ x ≤ 10⁷ | max 0.567 at x = 199 | O(x^ε) |
| σ(n) / (e^γ n log log n) for n > 5040 | n ≤ 10⁶ | max 0.98582 at n = 10080 | < 1 |
| σ(n) − (H_n + e^{H_n} log H_n) | n ≤ 10⁶ | max 0, at n = 1 only | ≤ 0 |
| (p′ − p) / (√p log p) | p ≥ 100, p ≤ 10⁷ | max 0.279 at p = 113 (gap 14) | O(1) |

Every one of these is a finite check, hence inert: it proves nothing about the universal statement, by the
Gödel theorem. What it shows is the shape of the law on the other end and how much room is in it. Robin's
ratio at 10080 is within 1.5 percent of the wall; the Schoenfeld bound for π(x) is tight only where it starts
and has a factor of 6 to spare by x = 10⁷; the Möbius sum uses about half of the square-root allowance in this
range.

## 3. What this is

It is the Hypothesis read forward instead of backward. The programme spends its effort on the backward
direction, from the stones to the point: subpower, census, floor, positivity. This note is the forward
direction, from the point to the stones, which is where every theorem above lives. None of it is new; all
of it is what the one dot buys. The two directions meet nowhere short of infinity: the forward laws hold
exactly if the zeros are on the point, and no finite stretch of the stones obeying them puts a zero there.

## 4. Measuring the difference: forward from the dot against backward from the stones (added later on 2026-10-05)

The user asked for the two directions computed and their difference measured. Forward: every zero placed
literally on the point, β = ½, the first K of them, rebuilt into the prime-power sum by the explicit formula,

    ψ_K(x) = x − Σ_{n≤K} 2 Re( x^{½+iγₙ} / (½+iγₙ) ) − log 2π − ½ log(1 − x⁻²).

Backward: the actual ψ(x) from the primes, by sieve to 10⁶. The difference D_K(x) = ψ(x) − ψ_K(x) is
measured in the dot's own unit, √x, and averaged in log x over 20,000 non-integer points of [10², 10⁶]
(`rh_one_dot_difference.py`, figure `THE_DIFFERENCE.png`). The prediction, if every zero is on the point and
none is missing, is Parseval for the almost periodic normalized error: the mean square of D_K/√x equals the
sum over the rings not yet placed, Σ_{n>K} 2/|ρₙ|², with the part beyond the 6,000th zero estimated from the
counting law as (log(T/2π) + 1)/(πT) at T = γ₆₀₀₀ = 6365.85.

| K zeros on the dot | height | mean of D_K/√x | rms measured | rms predicted | max |
|---|---|---|---|---|---|
| 10 | 49.8 | −0.0001 | 0.1378 | 0.1383 | 0.578 |
| 30 | 101.3 | +0.0001 | 0.1082 | 0.1086 | 0.414 |
| 100 | 236.5 | +0.0002 | 0.0784 | 0.0789 | 0.292 |
| 300 | 541.9 | +0.0003 | 0.0540 | 0.0566 | 0.230 |
| 1,000 | 1,419.4 | +0.0002 | 0.0343 | 0.0379 | 0.145 |
| 3,000 | 3,533.3 | +0.0002 | 0.0219 | 0.0257 | 0.094 |
| 6,000 | 6,365.9 | +0.0002 | 0.0162 | 0.0199 | 0.073 |

Three things come out of the measurement.

1. **The difference is the rings not yet placed, and nothing else.** For K ≤ 100 the measured and predicted
   rms agree to half a percent. For larger K the prediction is dominated by the extrapolated tail beyond the
   last known zero and the measurement falls 10 to 20 percent below it, within what that extrapolation and
   the finite averaging length allow. There is no residual: no drift, no constant, no growth. A zero off the
   point at height γ₀ would add a term of size x^{β−½}/|ρ₀|, growing with x; a ring missing from the dot
   would raise the rms above the prediction by 2/|ρ₀|². Neither is seen.
2. **The mean is zero.** For every K the mean of the difference is within 3 × 10⁻⁴ of zero. The stones alone
   have mean −0.037 in this unit, and that number is the constant: −log 2π averaged over the grid in 1/√x
   gives −0.0395. The constants we have are exactly where the forward side puts them.
3. **The difference shrinks at the rate of the horn's own law.** Σ_{n>K} 2/|ρₙ|² is, by the counting law,
   (log(γ_K/2π) + 1)/(πγ_K): the rings beyond K contribute inversely to their height, times the density of
   the horn at that height. Placing more zeros on the dot removes the difference at exactly the rate the horn
   flares.

What this is: the forward and backward directions of the explicit formula compared on the integers, and
found to agree to the precision of the tail. It is a finite computation, inert, and it is the same computation
as Figure 27A of the One Dot paper done as a measurement rather than a picture. What it is not: a test that
can see a zero off the point at height beyond 6,366, since such a zero contributes to the difference exactly
what an on-point ring at that height would, up to a factor x^{β−½} that is invisible below x = 10⁶ unless
β − ½ is large. The one dot is exact to the horizon of the zeros used, and the horizon is where every measurement
in this programme ends.
