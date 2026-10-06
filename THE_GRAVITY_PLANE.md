# The plane as a gravitational field: masses, shells, energy, modes (2026-10-05)

The user asked to place the primes in the light plane with mass, add the zeros, treat the plane as a
space-time warped by the mass, with mass equal to energy, and see what happens. There is an exact version of
this in the record already, and it is Newtonian rather than Einsteinian: the Green energy of the running proof
is the gravitational self-energy of concentric shells. Every physical word below names an identity. The
computations are `rh_gravity.py` and `rh_gravity_fig.py` (figure `THE_GRAVITY_PLANE.png`), to N = 10⁷.

## 1. The dictionary, exact

| word | object |
|---|---|
| a mass at the integer n | a spherical shell of radius n carrying mass m(n) |
| its gravitational potential at distance r | 1/max(r, n): Newton's shell theorem (constant inside the shell, 1/r outside) |
| the potential energy of the configuration | the Green energy of the record, E = Σ_{i,j} m(i) m(j) / max(i, j) |
| the field at radius t | M(t)/t², with M(t) = Σ_{n≤t} m(n) the enclosed mass |
| the field energy | ∫ M(t)²/t² dt (energy density M²/t⁴ over the volume 4πt² dt), plus the endpoint M(N)²/N |
| mass equals energy | E is a quadratic form in the masses: the identity E = ∫_1^N M(t)²/t² dt + M(N)²/N, the running proof's R(N) |
| the modes of the field | the zeros: in log radius the field M(e^u)e^{−u/2} is the Möbius state, whose resonances are the zeros |
| the space-time | the half-line of radii, read in log r: the dilation coordinate u = log r in which the modes are oscillators e^{iγu} |

The shell theorem identity was checked directly at N = 2000: the double sum over the masses with the kernel
1/max(i, j) equals the field-energy integral to all printed digits (1.490951 = 1.490951).

## 2. Two stars

**The prime star.** Masses Λ(n) at the prime powers: positive, so a real star. Its enclosed mass is ψ(t) ≈ t,
a uniform star of density one per unit radius, and its field energy is 2.000 N exactly as a uniform star's
is (measured 1.9997 N at 10⁷). What carries the zeros is the fluctuation of the mass, ψ(t) − t, whose field
energy ∫(ψ − t)²/t² dt grows like a constant times log N. Under RH that constant is Cramér's,

    Σ_ρ 1/|ρ|² = 2 + γ − log 4π = 0.04619,

the same constant that appeared as the partial sums of 1/(ρ(1 − ρ)) in Figure 23D of the One Dot paper.
Measured slope per unit of log N over successive half-decades from 10³ to 10⁷: 0.055, 0.042, 0.048, 0.047,
0.045, 0.046, 0.043, 0.050; overall 0.0469. The fluctuation energy of the prime star grows at the rate the
zeros set, to one and a half percent.

**The Möbius star.** Masses μ(n): signed, so a star with negative mass on half its shells. Its field energy
is the historical energy of the running proof, I_M(N), and the subpower target says it grows slower than any
power of the radius. Split the energy into the self-energy of the shells, Σ μ(n)²/n ~ (6/π²) log N, and the
interaction between shells, the pull:

| N | field energy R(N) | self-energy | interaction |
|---|---|---|---|
| 10³ | 1.459 | 5.242 | −3.783 |
| 10⁴ | 1.582 | 6.643 | −5.061 |
| 10⁵ | 1.622 | 8.043 | −6.421 |
| 10⁶ | 1.708 | 9.443 | −7.735 |
| 10⁷ | 1.837 | 10.843 | −9.006 |

The interaction is negative and grows like −0.56 log N against the self-energy's +0.61 log N (the measured
self-energy slope is still settling toward 6/π² = 0.608 from above). The Möbius star is bound: the pull between
its shells cancels about ninety percent of their self-energy, and the net field energy grows at a measured
0.0297 per unit of log N from 10³ to 10⁷. The constant the zeros predict for that slope, under RH and Ng's
hypothesis, is Σ_ρ 2/|ρ ζ′(ρ)|²; summed over the first 400 pairs of zeros it is 0.0288 and still rising
slowly. The Möbius star's field energy grows at the rate the zeros set, to three percent.

## 3. What happens, then

1. **The zeros are the modes of the field.** Both stars' fluctuation energies grow at rates that are sums over
   the zeros, one with weights 1/|ρ|², the other with weights 1/|ρζ′(ρ)|². The primes supply the masses, the
   zeros supply the modes, and the energy computed from the masses equals the energy computed from the modes.
   That is the explicit formula's Parseval identity in gravitational words, and it is what "the zeros and the
   primes exist simultaneously together" means exactly.
2. **The pull is the cross term.** The interaction energy of the Möbius star is Σ_{i<j} 2μ(i)μ(j)/j, the same
   signed cross term that every route in the programme ends at. Here its sign is measured: negative, binding,
   cancelling the self-energy to within a slowly growing remainder. RH is the statement that the remainder
   stays subpolynomial. In this language the Hypothesis is that the Möbius star never unbinds.
3. **Mass equals energy, and that is the problem.** Because the energy is quadratic in the masses, the
   positivity that would decide the Hypothesis is not the positivity of the masses, which are signed, nor of
   the self-energy, which is trivially positive, but the sign of the interaction at every scale. Conservation
   of the total does not fix it. This is Build v32's lossless engine and unbounded camera, with the engine now
   a star.
4. **Warping, honestly.** The potential 1/max(r, R) is Newton's, not Einstein's: the plane is not curved by
   the masses, the masses sit in it and attract through it. A curved version would need a metric, and the
   only metric the record owns is the inner product of the Weil form, which is the open question. So the
   space-time here is flat, the field is Newtonian, and the one thing a curved theory would add, the back
   reaction of the field on the geometry, is exactly the object nobody has.

## 4. What this is not

Not a physical claim and not a new theorem. The shell identity is the running proof's lifetime identity in
another coordinate; Cramér's constant is a theorem under RH; Ng's constant is conjectural; the measurements
are finite and inert. What the picture does is give the subpower target a body: a star with signed masses
whose binding energy must never fail by a power of its radius.
