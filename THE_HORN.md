# The horn of rings, reopened (2026-10-05)

The magenta horn is Figure 26 of *Primes, Folds, and the One Dot* (`received/Primes_Folds_One_Dot_paper.pdf`,
paper page 35). This note records how it was built, why it is hollow, what shape it actually is, and where π
enters. The checks are recomputed on the first 6,000 zeros (`rh_horn.py`, `data/zeros_6000.txt`,
`data/horn_rings.csv`, `data/horn_stones.json`, figure `THE_HORN.png`). Nothing here bears on the Riemann
Hypothesis; the horn is the counting function drawn, and every statement below is a statement about the
counting function.

## 1. How it was built

Section 10 of the paper, verbatim: *draw the n-th zero ½ + iγ_n as a circle of circumference γ_n at height n.
The resulting surface of revolution (a "horn") has profile given by the inverse of the counting function.*

So the rotation is not computed from anything; it is a drawing device. Ring n has

    circumference = γ_n,     radius r_n = γ_n / 2π,     height = n.

The first ring has circumference 14.13 and radius 2.25; ring 1,000 has circumference 1,419.4 and radius 225.9;
ring 6,000 has circumference 6,365.9 and radius 1,013.2. The paper's Computation 10.2 (first 1,000 rings) is
reproduced exactly: growth in radius per ring 0.2239 against 0.2238 from 1/log r. Over 6,000 rings the mean
growth is 0.1685 against 0.1685.

The profile is the inverse of the Riemann–von Mangoldt counting function. In the radius variable it reads

    n = r (log r − 1) + 7/8 + S,        S = π⁻¹ arg ζ(½ + iγ_n) = O(log γ_n),

and the figure's "rings counted minus the law" is this S, evaluated just above each zero. Over 6,000 rings it
stays in [−0.35, 1.36] with root mean square 0.58 (the paper's ±0.73 for 1,000 rings is the same quantity
centered on the half-count). Only zeros on the line can be drawn as rings, while the counting law counts every
zero in the strip; that is why, for the Davenport–Heilbronn function, the ring count drops by exactly 2 at each
off-line pair (paper, Figure 26C). The horn is the one picture in the paper in which an off-line zero shows up
as something missing rather than something added.

## 2. Why it is hollow

The horn is hollow because the construction gives it no interior. The heights γ_n are the only data; each
becomes a ring; the rings are stacked; the surface they trace is the horn. No quantity is assigned to the
inside, so the inside is empty by construction, not by discovery.

There is one thing that can be placed inside it without inventing anything, and it comes from the light cone.
The prime power m enters the cone of the odd Weil form at half-width a_m = ½ log m, and the horizon there is

    T*(a_m) = 2π e^{2 a_m} = 2π m,   i.e.   ring radius  r* = m.

So, in the horn's own unit, **the stone m sits at ring radius m.** When m switches on, the cone has resolved the
rings of radius ≤ m, which are the lowest N(2πm) rings of the horn:

| stone m | enters at a | horizon 2πm | rings resolved N(2πm) | law m(log m − 1) + 7/8 |
|---|---|---|---|---|
| 2 | 0.3466 | 12.57 | 0 | 0.26 |
| 3 | 0.5493 | 18.85 | 1 | 1.17 |
| 4 | 0.6931 | 25.13 | 3 | 2.42 |
| 5 | 0.8047 | 31.42 | 4 | 3.92 |
| 7 | 0.9730 | 43.98 | 8 | 7.50 |
| 8 | 1.0397 | 50.27 | 10 | 9.51 |
| 9 | 1.0986 | 56.55 | 12 | 11.65 |
| 11 | 1.1989 | 69.12 | 16 | 16.25 |
| 13 | 1.2825 | 81.68 | 21 | 21.22 |
| 16 | 1.3863 | 100.53 | 29 | 29.24 |
| 97 | 2.2874 | 609.47 | 348 | 347.62 |
| 997 | 3.4524 | 6264.34 | 5888 | 5887.91 |

The full table for every prime power to 1,000 is `data/horn_stones.json`. Two things in it are worth a look.
At the entry of 2 the horizon is 4π = 12.57, below the first zero at 14.13: when the first stone switches on,
no ring has been resolved yet, which is the identity phase of the potential set in the horn's language. And the
law column is the counting law with π removed: the stone m resolves m(log m − 1) + 7/8 rings, up to the wobble.

The axis of the horn is therefore the prime line in a new unit. The stones climb it, stone m at radius m, and
the horn around them is what each stone has seen. This uses the horizon T*(a) = 2πe^{2a}, which the paper
labels Heuristic 12.4 as a statement about resolution; the identity r*(a_m) = m is exact arithmetic given that
definition.

## 3. Where π enters

Five exact places, and one where it does not.

1. **Radius = height / 2π is wavenumber.** The zero ½ + iγ contributes to the prime line the wave x^{iγ}, which
   in u = log x has period 2π/γ. Its ring radius r = γ/2π = 1/period is the number of oscillations of that
   zero's wave per e-fold of x. The horn's radius axis is "oscillations per e-fold"; π is the conversion
   between a height and a period.
2. **The counting law loses its π in the radius variable.** N(T) = (T/2π) log(T/2πe) + 7/8 + S becomes
   n = r(log r − 1) + 7/8 + S: no π remains. Every 2π in the Riemann–von Mangoldt formula is the same
   2π, the one in r = γ/2π.
3. **One ring per π of phase.** The Riemann–Siegel phase is θ(t) = π[r(log r − 1) − 1/8] + O(1/t) with
   r = t/2π. Over the 6,000 zeros the mean advance of θ per ring is 1.00003 π. Gram points are the heights
   where θ is a multiple of π, and of the first 5,990 Gram intervals, 85.3 percent hold exactly one ring,
   7.4 percent none, 7.3 percent two (Gram's law, in the horn: one ring per turn of π).
4. **The horizon loses its π.** T*(a) = 2πe^{2a} is r* = e^{2a}; at the entry of stone m it is r* = m. The
   2π in the horizon is again the conversion of a height to a radius.
5. **The slope is the mean spacing.** The mean gap between zeros at height T is 2π/log(T/2π); in radius
   units it is 1/log r, which is exactly the growth per ring measured in Section 1. The local half-angle of
   the horn is the mean zero spacing, with its 2π absorbed.
6. **Where π does not come from.** Nothing in ζ makes the horn round. The circle is the drawing; a square of
   perimeter γ_n would give a square horn with the same profile. The one geometric content is the profile,
   the inverse counting function, and the π in it is the π of the Riemann–von Mangoldt density
   (1/2π) log(t/2π), which is the density of zeros per unit height and is in turn the Mellin dual of the
   density of prime powers per unit of log n. That is the π of the explicit formula, and it is the same π
   in every item above.

## 4. Horn, not cone

It looks like an ice-cream cone because it is drawn over 1,000 rings with the height axis compressed. It is not
a cone. A cone has r proportional to n; the horn has r ~ n/log n, and its slope dr/dn = 1/log r keeps
decreasing:

| ring | radius | local slope 1/log r |
|---|---|---|
| 1 | 2.25 | 1.233 |
| 10 | 7.92 | 0.483 |
| 100 | 37.64 | 0.276 |
| 1,000 | 225.91 | 0.184 |
| 6,000 | 1,013.16 | 0.144 |

The horn steepens forever, because the zeros get denser forever: it is the picture of "the genus of ζ is
infinite" from `FUNCTION_FIELD_CONE.md`. The same construction for a curve of genus g over 𝔽_q gives an exact
cone. The zeros of ζ_X are periodic in height with period 2π/log q, so the rings sit at radii
(θ_i/2π + k)/log q, exactly 2g of them per unit of radius · log q, and the profile is the straight line
r = n/(2g log q). For the genus-3 curve over 𝔽₃ the slope is 1/(6 log 3) = 0.152, between ζ's slope at ring
1,000 (0.184) and at ring 6,000 (0.144). Panel B of `THE_HORN.png` shows ζ's slope falling through the curve's
constant. An off-circle pair on a curve removes two rings in every period, so the curve's deficit grows
linearly where ζ's deficit from one off-line pair is a constant 2 (paper, Figure 26C): over a curve an impostor
is caught in every period, and the cone over a curve catches it by degree g (`FUNCTION_FIELD_CONE.md`).

## 5. What this does and does not say

The horn is N(T) drawn as a surface of revolution; its wobble is S(T); its slope is the mean spacing; its
interior can hold the stones at radius m by the horizon identity; and its shape is a horn because the zeros
have infinite genus. Each of these is a restatement of a classical fact in the paper's picture. None of it
constrains where the zeros are. What the horn does do, and the paper already used it for, is make an off-line
zero visible as a missing ring, which is the same event the cone over a curve catches at a finite support and
the cone over the integers never catches at any finite support.

## 6. Forcing the horn to close (added later on 2026-10-05)

The horn flares forever because the zeros get denser forever. There are exactly three ways to force it to
close, and only one of them keeps ζ.

1. **Cut it off** at ring N. The zero measure becomes finite, the plane at one half becomes
   N-dimensional, and the cone over it closes at every support at once, because a function of any compact
   support can vanish at N given points. The hollow image loses its meaning. Nothing is learned.
2. **Make it periodic.** This is what closes the horn of a curve: ζ_X is a function of q^{−s}, the zeros
   repeat with period 2π/log q in height, and the horn, read modulo the period, is 2g points on a circle.
   ζ has no period, and the Riemann–von Mangoldt law forbids one.
3. **Unfold it.** Measure height not in γ but in the Riemann–Siegel phase, φ = θ(γ)/π + 1. Then ring n
   sits at φₙ = n − S(γₙ): at integer height up to the wobble, and every ring has unit radius in the unit of
   the local mean spacing. The flare is gone; the horn is a cylinder. This is the unfolding of random-matrix
   theory, it loses nothing, and it is the one honest closure.

In the plane at one half the unfolded horn is the following object. The space is L²(ν_ζ); the unfolding
sends the atoms γₙ to the points n − S(γₙ) of the line; the measure becomes the counting measure on the
integers displaced by S. What remains after closing is therefore exactly S(T), and the law of the rings'
spacing.

Computed on the 6,000 zeros (`rh_horn_closed.py`, figure `THE_HORN_CLOSED.png`): the unfolded rings have
mean spacing 1.00003 and standard deviation 0.39; no two rings share a height (smallest spacing 0.046);
the displacement from integer height stays in [−1.36, 0.35]; and the spacing histogram sits on the GUE
Wigner surmise (32/π²) s² e^{−4s²/π} at L² distance 0.09 against 0.67 for the Poisson law e^{−s}. Only
0.9 percent of spacings are below a quarter of the mean, against 22 percent if the rings were independent.
The rings repel. That repulsion is Montgomery's pair correlation, proved in its restricted range, and it is
the quantum signature of the plane: a spectrum whose levels repel is the spectrum of a Hermitian operator,
not of independent events.

So forcing the horn to close does not remove the open question; it isolates it. After the closure the
horn is a cylinder of unit rings, and the entire content of the zeros is in two things: the wobble S, which
is bounded on average and is the subject of the census, and the repulsion law of the spacings, which is
the quantum theory. An off-line pair in the Davenport–Heilbronn function shows up in the closed horn as a
missing pair of rings, a step of two in the count; in ζ no such step has ever been seen.
