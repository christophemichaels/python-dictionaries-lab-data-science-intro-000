# The arithmetic field: the Green energy as a Yukawa field of mass ½ (2026-10-05)

The user relayed a proposal: in logarithmic coordinates the shell kernel becomes the Green function of a
definite operator, so the gravitational picture has an exact field equation. The proposal is correct, and
checked here; it also yields two things it did not state. Script `rh_field.py`, figure `THE_FIELD.png`.

## 1. The field equation

Put u = log x. The shell kernel of `THE_GRAVITY_PLANE.md` factors as

    1 / max(a, b)  =  (ab)^{−1/2} · e^{−|log a − log b| / 2},

and e^{−|u−v|/2} is the Green function of

    H₀ = −d²/du² + ¼,        H₀ e^{−|u−v|/2} = δ(u − v).

So with the Möbius source and its response,

    J_N(u) = Σ_{n≤N} μ(n) n^{−1/2} δ(u − log n),        H₀ φ_N = J_N,

the field is φ_N(u) = Σ_{n≤N} μ(n) n^{−1/2} e^{−|u − log n|/2} and its energy is the Green energy of the record:

    ∫_ℝ ( φ_N′(u)² + ¼ φ_N(u)² ) du  =  ⟨J_N, H₀⁻¹ J_N⟩  =  Σ_{m,n≤N} μ(m)μ(n) / max(m,n)  =  R(N)  =  I_M(N) + M(N)²/N.

Checked to six decimals at N = 10³, 10⁴, 10⁵, 10⁶ (1.459409, 1.582361, 1.621567, 1.708249), and from the energy
density on a grid at N = 10⁴ (1.5816 against 1.5824, the difference being the grid at the pulses).

H₀ is the Klein–Gordon operator of a field of mass ½ in one dimension. **The exponent ½ of the whole programme
is the mass of the field.** Its spectrum is [¼, ∞), and on the mode e^{iγu} it takes the value γ² + ¼ = |ρ|²:
the quantity 1/|ρ|² that weights every Parseval sum in the record is the inverse energy of the field at the
frequency of the zero. The Yukawa range of the kernel is two e-folds.

## 2. The two solutions, and the equipartition

Between pulses H₀φ = 0, whose solutions are e^{−u/2} and e^{u/2}, that is x^{−1/2} and x^{1/2}: the two sides
of the critical point. The field is exactly

    φ_N(log x)  =  x^{−1/2} M(x)  +  x^{1/2} ( h(N) − h(x) ),        h(x) = Σ_{n≤x} μ(n)/n,

the inward solution carried by the enclosed mass, which is the Möbius state m(u) of *Across the Horizon*, and
the outward solution carried by the mean obligation h of `RESEARCH_H_RECORD.pdf`, Section 7. The energy
density between the pulses at k and k+1 is ½ M(k)²/x + ½ x (h(N) − h(k))², and integrating:

    R(N)  =  ½ [ I_M(N) + M(N)²/N ]  +  ½ [ h(N)² + Σ_{k<N} (h(N) − h(k))² ].

Since the left side is I_M(N) + M(N)²/N, the two brackets are equal. The energy is shared exactly in half
between the inward and the outward solution, and the outward half gives a new exact form of the energy:

    R(N)  =  h(N)²  +  Σ_{k<N} ( h(N) − h(k) )²,

the Green energy is the sum of the squared mean obligations. (Proof in one line: Σ_{k<N} Σ_{m,n>k}
μ(m)μ(n)/(mn) counts, for each pair, the k below min(m,n), and min(m,n)/(mn) = 1/max(m,n).) Checked to six
decimals at the same four N. So the subpower target reads, equivalently: the mean square of the tails of
Σ μ(n)/n is subpolynomial. Under RH each tail is O(k^{−1/2+ε}) and the sum is O(N^{2ε}); the route around h
of the record estimate is this identity read at the records.

## 3. The spectral form

Fourier transforming in u, the energy is

    R(N)  =  (1/2π) ∫_ℝ | Ĵ_N(k) |² / (k² + ¼) dk,        Ĵ_N(k) = Σ_{n≤N} μ(n) n^{−1/2 − ik},

the Möbius Dirichlet polynomial on the critical line, weighted by the inverse energy of the field. As N grows
the polynomial tends, in the sense of the record's Mellin step, toward 1/ζ(½ + ik), whose poles are the zeros:
the energy grows because the source's spectrum sharpens onto the lines of the spectrum, and the subpower says
it sharpens slower than any power. This is the same statement as Crossing 3 of *Arithmophysics*, now with the
weight 1/(k² + ¼) identified as the field's propagator.

## 4. What can be done with it, honestly

| proposal | what it produces | assessment |
|---|---|---|
| develop the field equation | Sections 1–3: the equation, the mass, the equipartition, the h-form, the spectral form | done; exact |
| build a simulation | the field and its energy density for any N; Figure `THE_FIELD.png` for N = 10⁴ | done; cheap to extend |
| quantum dynamics with H₀ | free propagation of a massive field in one dimension: dispersion of the pulses, nothing arithmetical, because H₀ knows nothing of ζ; the arithmetic is entirely in the source | not worth a project |
| the zero spectrum | the operator whose spectrum is the zeros: H₀ has continuous spectrum [¼, ∞) and no arithmetic; what is needed is the identification that makes the source's resonances into eigenvalues, the inner product of Section 4 of *Across the Horizon* | the open problem, restated |

So the field equation adds two exact facts to the record, the mass-½ reading and the h-form of the energy, and
it locates the open problem in the same place as before: not in the operator, which is free, but in the
positivity that would turn the source's spectral resonances into the eigenvalues of something.
