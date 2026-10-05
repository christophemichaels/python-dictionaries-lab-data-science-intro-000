# The zero hole: the plane inverted into an onion (2026-10-05)

The user asked to place all the zeros at one point, call it the zero hole, warp the whole space into an onion
whose top points to that one location, and see how the masses interact with the point. There is an exact
warp that does this, and the record already owns the quantity that comes out of it. Script `rh_zero_hole.py`,
figure `THE_ZERO_HOLE.png`, computed to N = 10⁷.

## 1. The warp

The gravity plane of `THE_GRAVITY_PLANE.md` is the half-line of radii with the masses m(n) on shells of
radius n and the potential 1/max(r, R). Invert it: r ↦ 1/r, the Kelvin transform, which preserves Newtonian
potentials up to the factor 1/r. The shell of radius n becomes a layer at radius 1/n; the layers nest inward,
1, ½, ⅓, …, accumulating at the centre; and the point at infinity, where the far field lives and where the
zeros are the modes, becomes the single point r = 0. That is the onion, and its centre is the zero hole. In
log coordinates the inversion is u ↦ −u, the reflection of the fold: the outside of the star folded onto the
inside, the point at infinity onto the centre, the zeros, which sit on the fixed line of the fold, onto
themselves. The one dot is the centre of the onion.

Panel A of the figure draws the first sixty layers of the Möbius onion, blue where μ(n) = +1, orange where
μ(n) = −1, faint where the layer carries no mass.

## 2. The interaction with the one point

Inside every shell the potential is constant, 1/n. So the potential of all the masses at the centre is

    Φ(0) = Σ_n m(n) / n,

the Dirichlet series of the masses at s = 1. That is how the space interacts with the one point: every layer
contributes its mass over its radius, and the sum is the series at the pole of ζ.

**The Möbius star.** Φ(0) = Σ μ(n)/n = 0: the potential of the Möbius star at the zero hole vanishes, and
that statement is equivalent to the prime number theorem (Landau). Its partial sums h(N) = Σ_{n≤N} μ(n)/n are
the mean obligation of `RESEARCH_H_RECORD.pdf`, Section 7, the quantity the route around h reduced the record
estimate to, and the rate at which they vanish is the Hypothesis: h(N) = O(N^{−1/2+ε}) if and only if RH.
Measured:

| N | h(N) | √N · h(N) |
|---|---|---|
| 10³ | +0.004412 | +0.140 |
| 10⁴ | −0.002083 | −0.208 |
| 10⁵ | −0.000487 | −0.154 |
| 10⁶ | +0.000201 | +0.201 |
| 10⁷ | +0.000102 | +0.321 |

Over 10² ≤ N ≤ 10⁷ the scaled potential √N · h(N) stays in [−0.44, +0.41] with root mean square 0.17
(Panel B). The potential of the Möbius star at the hole vanishes at the square-root rate throughout the
range, with the sign changing as the layers are added.

**The prime star.** Φ(0) = Σ_{n≤N} Λ(n)/n = log N − γ + o(1), Mertens' theorem: the uniform part of the star
gives the divergent log N, and the fluctuation gives exactly minus Euler's constant. Measured at 10⁷: the
potential is 15.540737 = log N − 0.577359, against γ = 0.577216 (Panel C). Under RH the error is
O(log² N / √N); its scaled size (log N − γ − Φ) √N / log² N stays within ±0.002 over the range.

So the constants we have appear once more, as potentials at the one point: 0 for the Möbius star, and −γ for
the prime star's fluctuation.

## 3. The energy in onion coordinates

With s = 1/t the field energy ∫_1^N M(t)² dt/t² becomes ∫_{1/N}^1 M(1/s)² ds: the energy per unit of onion
depth is the square of the mass enclosed outside that depth. The subpower target says the energy stored in
the core of the onion, between depth 1/N and the centre, grows slower than any power of N: the singularity
at the zero hole is integrable to within every ε. The modes of that singularity are the zeros, and an
amplifying mode, a zero off the point, would make the core's energy grow like a power.

## 4. What happens, then

Placing all the zeros at one point and folding the space onto it does not change the mathematics; it
changes which quantity sits in front. In the plane it was the field energy and its growth. At the hole it is
the potential, the Dirichlet series at s = 1, and its two values: zero for the Möbius star, which is the
prime number theorem, and −γ for the prime star, which is Mertens. The Hypothesis is the rate at which the
first one is reached, √N, and that rate is the mean obligation of the record estimate, the object the route
around h arrived at from the other direction. The onion and the record meet at the same number.

What it is not: a construction of anything new. The Kelvin transform is classical, Landau's and Mertens'
theorems are classical, and the measurements are finite and inert. The hole is the point at infinity with a
name, and the name is the one dot's.
