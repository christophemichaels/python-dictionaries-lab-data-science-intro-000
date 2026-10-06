# The arithmetic field: the Green energy as a Yukawa field of mass ½ (2026-10-05)

The user relayed a proposal: in logarithmic coordinates the shell kernel becomes the Green function of a
definite operator, so the gravitational picture has an exact field equation. The proposal is correct, and
checked here; it also yields two things it did not state. Script `rh_field.py`, figure `THE_FIELD.png`.
Sections 5–7 (2026-10-06) add the exact recovery law, the successor as a field event, and the criterion at
infinite scale stated correctly (subpower growth of the energy, not boundedness). Script `rh_field_recovery.py`,
figure `THE_FIELD_RECOVERY.png`.

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

the Möbius Dirichlet polynomial on the critical line, weighted by the inverse energy of the field. At finite N
the numerator is an entire function of k; the only poles of the integrand are the propagator's, at k = ±i/2.
The zeros of ζ are not in this formula at any finite N. They enter only in the passage N → ∞: the Mellin
transform of the inward channel is ∫_1^∞ M(x) x^{−s−1} dx = 1/(sζ(s)) for Re s > 1, and the growth of the
energy as N → ∞ is governed by how far that meromorphic function continues to the left (Section 7). So the
energy grows because the source's spectrum sharpens, with N, onto the lines of 1/ζ, and the subpower says it
sharpens slower than any power. This is the same statement as Crossing 3 of *Arithmophysics*, now with the
weight 1/(k² + ¼) identified as the field's propagator, and with the finite-to-infinite passage kept explicit.

## 4. What can be done with it, honestly

| proposal | what it produces | assessment |
|---|---|---|
| develop the field equation | Sections 1–3: the equation, the mass, the equipartition, the h-form, the spectral form | done; exact |
| build a simulation | the field and its energy density for any N; Figure `THE_FIELD.png` for N = 10⁴ | done; cheap to extend |
| quantum dynamics with H₀ | free propagation of a massive field in one dimension: dispersion of the pulses, nothing arithmetical, because H₀ knows nothing of ζ; the arithmetic is entirely in the source | not worth a project |
| exact recoverability | Section 5: the field's value and slope return M and the tail of h; its endpoint values are h(N) and M(N)/√N; its derivative jumps are the source | done; exact |
| the zero spectrum | the operator whose spectrum is the zeros: H₀ has continuous spectrum [¼, ∞) and no arithmetic; what is needed is the identification that makes the source's resonances into eigenvalues, the inner product of Section 4 of *Across the Horizon* | the open problem, restated |

So the field equation adds two exact facts to the record, the mass-½ reading and the h-form of the energy, and
it locates the open problem in the same place as before: not in the operator, which is free, but in the
positivity that would turn the source's spectral resonances into the eigenvalues of something.

## 5. Exact recoverability

The field was built from the arithmetic; the arithmetic comes back out of the field. Away from the pulses write

    A_N(u) = e^{−u/2} M(e^u),        B_N(u) = e^{u/2} ( h(N) − h(e^u) ),        0 < u < log N.

Between pulses M and h are constant, so A′ = −A/2 and B′ = B/2, and from φ_N = A_N + B_N,

    φ_N′ = ( B_N − A_N ) / 2,

which inverts to the **recovery law**

    M(e^u)          =  e^{u/2}  ( φ_N/2 − φ_N′ ),
    h(N) − h(e^u)   =  e^{−u/2} ( φ_N/2 + φ_N′ ).

The field's value and slope at any one point jointly encode the enclosed mass and the remaining obligation.
At the two ends of the source,

    φ_N(0) = h(N),        φ_N(log N) = M(N) / √N,

since M(1) = h(1) = 1; and the jump of φ_N′ at u = log n is −μ(n) n^{−1/2}, the source itself. So value, slope and
derivative jumps return μ, M and h entirely.

Checked with the field evaluated directly as the sum of its pulses, never through the closed form: at N = 10⁴
and the 9,999 midpoints log(k + ½), the recovered M agrees with the sieve to 5·10⁻¹⁴ and the recovered tail to
10⁻¹⁶; the endpoint values agree to 10⁻¹¹ (`rh_field_recovery.py`, figure `THE_FIELD_RECOVERY.png`, panels A, B).

The equipartition of Section 2 is the factorisation of the operator. H₀ = (−∂_u + ½)(∂_u + ½), and the two
channels are the kernels of the two factors (A is killed by ∂_u + ½, B by −∂_u + ½), which gives

    φ_N′² + ¼ φ_N²  =  ½ A_N² + ½ B_N²

identically. Extended over the whole line with the truncated sums (A = 0 for u < 0, A = e^{−u/2} M(N) for
u > log N; B = e^{u/2} h(N) for u < 0, B = 0 for u > log N), each channel's squared integral is R(N) by itself:

    ∫_ℝ A_N² du  =  I_M(N) + M(N)²/N  =  R(N),        ∫_ℝ B_N² du  =  h(N)² + Σ_{k<N} (h(N) − h(k))²  =  R(N),

so each supplies exactly half of the completed field's energy. Checked at N = 10³ … 10⁷ (1.459409, 1.582361,
1.621567, 1.708249, 1.836533 on both channels).

## 6. The successor as a field event

The arithmetic evolution is the arrival of the next pulse. Adding the integer n changes the field by

    φ_n(u) − φ_{n−1}(u)  =  μ(n) n^{−1/2} e^{−|u − log n|/2},

and the energy by the self-energy of the pulse plus its interference with the field already present,
2⟨φ_{n−1}, δφ⟩_{H₀} = 2 μ(n) n^{−1/2} φ_{n−1}(log n):

    R(n) − R(n−1)  =  ( μ(n)² + 2 μ(n) M(n−1) ) / n,        because  φ_{n−1}(log n) = M(n−1) / √n.

This is the successor identity M(n)² − M(n−1)² = μ(n)² + 2μ(n)M(n−1) of the record, read geometrically: the
entire earlier source enters the new energy through the one number the existing field has at the point where
the pulse lands. That is the field form of the reach-back interaction of `THE_GRAVITY_PLANE.md`. Checked: the
increments summed from n = 1 reproduce R(N) to 10⁻¹³ at N = 10³ … 10⁷, and φ_{n−1}(log n) computed from the
pulses equals M(n−1)/√n at n = 30, 31, 1009, 9973 (panel C of the figure shows n = 31).

## 7. The criterion at infinite scale, stated correctly

The sentence "the energy is bounded iff RH", said in conversation, is wrong, and this section replaces it.

The energy is never bounded. The undamped inward integral ∫_0^∞ |m(u)|² du = ∫_1^∞ M(t)²/t² dt diverges
unconditionally (Proposition 3.2 (iii) of *Across the Horizon*: m ∉ L²), so R(N) → ∞ for every arithmetic.
The precise statement is **subpower growth**:

    RH  ⟺  R(N) = O_ε(N^ε) for every ε > 0.

One direction: M(N)²/N ≤ R(N), so subpower energy gives M(N) = O(N^{1/2+ε}), which is RH. The other: under RH,
M(t) = O_δ(t^{1/2+δ}), and inserting this in R(N) = ∫_1^N M(t)²/t² dt + M(N)²/N gives O(N^{2δ}). Equivalently,
in the field's own language, **arbitrarily weak exponential damping in logarithmic distance makes the
accumulated energy finite**:

    RH  ⟺  ∫_0^∞ e^{−2δu} |m(u)|² du < ∞ for every δ > 0,

the Mellin transform of e^{−δu}m being 1/(sζ(s)) on the line Re s = ½ + δ, and square integrability on every
such line giving the continuation of 1/ζ to Re s > ½. The damped energies at δ = 0.05, 0.1, 0.25 beside the
undamped I_M(N):

| N | δ = 0 | δ = 0.05 | δ = 0.1 | δ = 0.25 |
|---|---|---|---|---|
| 10³ | 1.455409 | 1.217483 | 1.040113 | 0.716688 |
| 10⁴ | 1.529461 | 1.250738 | 1.055115 | 0.718103 |
| 10⁵ | 1.598527 | 1.275178 | 1.063802 | 0.718504 |
| 10⁶ | 1.663305 | 1.293383 | 1.068940 | 0.718622 |
| 10⁷ | 1.728997 | 1.308117 | 1.072259 | 0.718661 |

Undamped square integrability is a different assertion, and false. Even the logarithmic bound on the undamped
integral, ∫_1^X M(t)²/t² dt = O(log X), is the weak Mertens conjecture, which is stronger than the subpower
statement and implies the simplicity of the zeros as well.

What the field adds, then, is a single object on which all of this is read: the arithmetic constructs φ_N; the
value and slope of φ_N give the arithmetic back; the next integer is a pulse whose energy cost is set by the
field already present where it lands; and the growth of the total energy in N is the criterion.
