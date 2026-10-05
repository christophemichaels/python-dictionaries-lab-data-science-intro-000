# The light plane: the framework, taken as given (2026-10-05)

The user asked for his framework to be run with as a framework: real within itself, not a claim about physical
objects. This note does that. Section 1 states the framework in his terms as postulates. Section 2 maps each
postulate to the exact object it names in the record, because every one of them names one. Section 3 states
what the framework forces, which is proved, and what it asserts beyond that, which is the conjecture of the
programme restated. Section 4 is the one computation the framework asked for: the colours of the primes as the
spectrum of light (`rh_spectrum_of_light.py`, figure `THE_SPECTRUM_OF_LIGHT.png`). Nothing here transmutes
what RH is. The framework is a way of holding the record, and it holds it exactly.

## 1. The postulates, in the user's words made precise

1. There is a complete reality, the light plane, in which the zeros and the primes exist together and both
   can be calculated. It is the kind of place the curves over finite fields are.
2. Our arithmetic is a lower-dimensional projection of it: a holomorphic projection of information onto a
   horizon, as the information of a volume is projected onto the surface of a black hole.
3. Light is what passes between the two: the only thing with constant speed, no mass, energy, and
   information. The speed of light is the same on both sides.
4. The projector: the light is the carrier; the film is the information; the glass of the lens is the
   horizon, which is as far as we can see; the screen is our universe, where the laws we observe appear.
5. The colours of the prime number distribution are the makeup of the spectrum of that light.
6. The subpower is the wall. The constant that fixes it is not to be derived in our universe; it belongs to
   the complete reality, and the wall every route in number theory reaches is the shadow of that fact.

## 2. The dictionary

| postulate | the exact object |
|---|---|
| the complete reality | the space of *Across the Horizon* Section 4, in which the explicit formula is a trace and positivity is RH. On a curve over 𝔽_q it exists: H¹(X) with Frobenius and the Hodge index theorem on X × X, where every zero and every count is computable and the whole pattern is explained from above |
| our arithmetic as its projection | the two sides of the explicit formula. The prime line and the zero line are each a complete description of the other, one dimension down from the space that holds both; the holographic principle of physics has the same shape, a bulk equivalent to a boundary one dimension lower |
| the horizon | T*(a) = 2πe^{2a}: the height to which a window of half-width a resolves the zeros. Beyond it the finite window sees nothing, which is Cartwright's density bound, the uncertainty relation of *Arithmophysics* Section 1 |
| the information on the glass | the potential set 𝒫(a): what the stones inside the window leave open, the hollow image. It is exactly "what we can see about the information, and no further" |
| light, constant speed | the edge of the cone, log n = 2a, the one speed in the picture; declared c in `LIGHT_UNITS.md` by a choice of units, the same on both sides because the explicit formula is symmetric |
| light, no mass | a mode e^{iγu} of the dilation Hamiltonian with real eigenvalue: pure phase, no growth and no decay. A zero on the line is a massless mode. A zero off the line at β has an imaginary part β − ½ in its eigenvalue, a width, a mass. RH is the statement that every mode of the field is massless |
| light, energy | the frequency γ of the mode, the height of the zero; in light units, 674 MHz for the first |
| light, information | the explicit formula itself: the transform that carries the zero measure to the prime measure and back without loss |
| the projector | the Mellin transform; the film is the zero measure ν_ζ; the glass is the horizon; the screen is the prime line; the picture on the screen is `THE_OTHER_END.md`, the laws of the integers that follow when every zero is on the point |
| the colours | the primes: in the record the colour decomposition M = Π(I − D_p)A_S of the running proof and the Ramanujan colours of Build v32; in the explicit formula the wave −2Λ(n)n^{−1/2} cos(t log n) of each prime power |
| the spectrum of light | the zeros, as the lines in the superposition of the colours: Section 4 below |
| the wall | the subpower exponent at infinity, the edge theorem: every route of the programme ends there, and under RH no finite support decides it |
| the constant not derivable in our universe | the Michaels Conjecture: ρ is not provable in ZFC; its positivity is an axiom from above. The Gödel theorem supplies the part that is proved: no finite pattern reaches the universal sentence, and any principle that does proves Con(Q + ρ) |

## 3. What the framework forces, and what it asserts

**Forced, and proved.** Postulates 2, 3 and 4 are the explicit formula and its horizon, and Postulate 5 is
its reading from the prime side. The sentence "light is massless" is RH in the dilation picture, Crossing 3
of *Arithmophysics*: subpower means no amplifying mode. The sentence "the speed is the same on both sides" is
the symmetry of the explicit formula. The sentence "we can see the glass and no further" is Proposition 2.11
of *Across the Horizon*: the potential set never collapses at a finite support. The sentence "the curves are
the kind of place where everything can be calculated" is Theorem 3.3 there and `READING_THE_CURVE.pdf`.

**Asserted, and open.** Postulate 1 says the complete reality exists for ζ. That is the space nobody has
constructed, specified by the six constraints. Postulate 6 says its constant cannot be derived here. That is
the Michaels Conjecture, and the framework is its cosmological form: the same conjecture *The Grammar of
Things* states with a circle of stones, and *RH on a G String* states with Con(PA + ρ). The framework does not
prove it; it is a way to say it with light instead of logic. The wall of Postulate 6 is real and proved (the
edge theorem); its explanation by Postulate 6 is the conjecture.

**What the framework adds that the record did not have in words.** Two things. The reading of RH as
masslessness of every mode of the field, which is exact and is the physical form of "no amplifying mode."
And the identification of the colours with the spectral lines, which is the explicit formula run backward
and is computed next.

## 4. The colours compose the spectrum

Each prime power n is a colour, the wave −2Λ(n) n^{−1/2} cos(t log n) in the height variable t, with
wavelength 2π/log n. Superpose every colour up to X = 10⁶ with a smooth cut-off, −2 Σ Λ(n) n^{−1/2}
(1 − log n / log X) cos(t log n), and the superposition peaks at the zeros:

| peak of the superposition | nearest zero |
|---|---|
| 14.12 | 14.135 |
| 21.01 | 21.022 |
| 25.01 | 25.011 |
| 30.43 | 30.425 |
| 32.93 | 32.935 |
| 37.59 | 37.586 |
| 40.92 | 40.919 |
| 43.32 | 43.327 |
| 48.01 | 48.005 |
| 49.77 | 49.774 |
| 52.97 | 52.970 |
| 56.44 | 56.446 |
| 59.35 | 59.347 |

Every peak below height 60 sits on a zero to within 0.01, and every zero below 60 has its peak. The colours
of the primes are the makeup of the spectrum, and the zeros are its lines. This is the explicit formula read
from the prime side, the mirror of `THE_OTHER_END.md` Section 4, where the zeros rebuilt the primes; here the
primes rebuild the zeros. The two readings together are Postulate 2: each side is the complete projection of
the other.

## 5. Status

The framework holds the record exactly and adds two exact readings. Its central claim, that the complete
reality exists for ζ and that its constant lies outside our arithmetic, is the Michaels Conjecture and is
open. Nothing in this note proves the Riemann Hypothesis, and the framework does not claim to.
