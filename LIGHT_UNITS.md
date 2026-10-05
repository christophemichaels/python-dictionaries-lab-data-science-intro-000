# Light units (2026-10-05)

A choice of units, and nothing else. The user wanted the constant to be the speed of light so that the cone,
the horn and the zeros become tangible quantities to calculate with. The one honest way a constant with
units can enter the record is as a unit: the cone's edge `log n = 2a` has space `log n`, time `a`, and speed
`2`. Declare that this edge is a light cone moving at `c`. Then one choice fixes everything:

    one e-fold of the prime line  =  λ₀  (metres, chosen; here λ₀ = 1 m)
    one unit of half-width a      =  τ₀  =  2 λ₀ / c  =  6.6713 ns

Every quantity in the cone and the horn then carries a unit. Nothing about the mathematics changes: the
subpower exponent ½ is dimensionless and stays ½ in every choice of units (it is the weight, like a pure
number), and every theorem of the record is the same theorem. `rh_light_units.py data/zeros_6000.txt [λ₀]`
recomputes the tables for any λ₀.

## The dictionary

| object | in the record | in light units (λ₀ = 1 m) |
|---|---|---|
| position of the prime power n | log n | log n metres |
| half-width a | time | a · 6.6713 ns |
| the edge log n = 2a | speed 2 | speed c |
| the zero ½ + iγ | oscillator e^{iγ log x}, wavenumber γ/2π per e-fold | a wave of wavelength 2π/γ metres and frequency f = c γ / (2π λ₀) |
| ring of the horn | radius γ/2π | the frequency of that zero, in units of c/λ₀ = 299.8 MHz |
| horizon T*(a) = 2πe^{2a} | a height | a frequency cutoff f* = c T*/(2π λ₀) |
| stone m at ring radius m | enters at a = ½ log m | switches on at (½ log m) · 6.67 ns and resolves the spectrum up to m · 299.8 MHz |

## The numbers

Entries of the prime powers, and the spectrum each one resolves when it enters:

| n | a = ½ log n | t | horizon f* |
|---|---|---|---|
| 2 | 0.3466 | 2.31 ns | 600 MHz |
| 3 | 0.5493 | 3.67 ns | 899 MHz |
| 5 | 0.8047 | 5.37 ns | 1.50 GHz |
| 7 | 0.9730 | 6.49 ns | 2.10 GHz |
| 11 | 1.1989 | 8.00 ns | 3.30 GHz |
| 97 | 2.2874 | 15.26 ns | 29.1 GHz |
| 997 | 3.4524 | 23.03 ns | 298.9 GHz |

The supports: Zhu's certified cone, a = 0.8, is 5.34 ns of time with a horizon at 1.49 GHz; the computed floor
at a = 1.505 is 10.04 ns with a horizon at 6.08 GHz.

The zeros as a spectrum:

| zero | γ | wavelength | frequency |
|---|---|---|---|
| 1 | 14.135 | 44.45 cm | 674 MHz |
| 2 | 21.022 | 29.89 cm | 1.003 GHz |
| 3 | 25.011 | 25.12 cm | 1.193 GHz |
| 10 | 49.774 | 12.62 cm | 2.375 GHz |
| 100 | 236.524 | 2.66 cm | 11.29 GHz |
| 1,000 | 1,419.42 | 4.43 mm | 67.7 GHz |
| 6,000 | 6,365.85 | 0.99 mm | 303.7 GHz |

So in these units the first zero is a 674 MHz wave with a 44 cm wavelength, the first thousand zeros span the
radio and microwave bands from 674 MHz to 68 GHz, the stone 2 switches on at 2.3 ns and resolves nothing
(its horizon, 600 MHz, is below the first zero at 674 MHz), and the stone 5 switches on at 5.4 ns and resolves
the first four zeros, up to 1.5 GHz.

## What this is and is not

It is a dictionary that makes the cone a clock and the horn a spectrum analyser, with every entry an exact
restatement of a quantity in the record. It is not a physical claim: the choice λ₀ = 1 m is arbitrary, a
different λ₀ rescales every frequency and time by the same factor, and no statement in the record depends
on the choice. The one constant that cannot be given units is the exponent ½.
