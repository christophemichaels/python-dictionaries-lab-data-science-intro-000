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
