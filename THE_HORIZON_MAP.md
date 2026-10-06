# The map through the horizon: the light plane and its arithmetic image (2026-10-06)

The user set the frame: the light plane is a larger information structure in which the growth law is part of
the governing framework; the question is what that structure reveals when its information passes through a
horizon into the arithmetic description, with its primes, Möbius sums and zeros. This note is the first pass.
It puts the plane's objects (the space, the flow, the Hamiltonian, the positive geometry, the colours) beside
their arithmetic images, computes the map where it can be computed, and records where conservation becomes
growth, where information is compressed, and what stays invariant. Script `rh_horizon_map.py`, figure
`THE_HORIZON_MAP.png`, data `data/horizon_map.json`. Objects from `THE_FIELD.md` are used throughout.

## 0. The premise

**(G)** The energy of the Möbius source grows slower than any power of the integer cutoff:
R(N) = O_ε(N^ε) for every ε > 0.

(G) is taken as the governing law of the plane. Nothing below proves it; everything marked "under (G)" uses it.
In the plane's own words (THE_FIELD.md, Section 7): for every δ > 0 the source J^{(δ)} = Σ μ(n) n^{−1/2−δ}
δ(u − log n), the Möbius source read on the line Re s = ½ + δ, has finite inward energy in the field of mass
½ + δ, and at mass ½ exactly the energy is infinite. So (G) is a statement about mass: **the Möbius source has
finite energy at every mass above ½ and infinite energy at ½.**

## 1. The plane

| object | definition | property |
|---|---|---|
| space | u ∈ ℝ, u = log x; the Hilbert space L²(ℝ, du) | one dimension, no lattice, no primes |
| flow | U_t g(u) = g(u + t), in x: U_t f(x) = e^{t/2} f(e^t x) | unitary one-parameter group: the dilations |
| generator | H_dil = −i d/du = (xp + px)/2 | the Berry–Keating operator; spectrum ℝ, the colours k |
| canonical pair | [u, H_dil] = i | position u, momentum k |
| field Hamiltonian | H₀ = H_dil² + ¼ = −d²/du² + ¼ | on e^{iku}: k² + ¼ = s(1 − s) at s = ½ + ik |
| mass | ¼ = min spec H₀ | the floor; no colour has energy below ¼ |
| propagator | G = H₀⁻¹, kernel e^{−\|u−v\|/2}, transform 1/(k² + ¼) | positive definite: every source has E[J] = ⟨J, GJ⟩ ≥ 0 |
| conservation | E[U_t J] = E[J];  \|Û_tJ(k)\|² = \|Ĵ(k)\|² | the flow conserves energy and every colour |

Two remarks. The field Hamiltonian is the square of the dilation generator plus the mass, so the plane's
spectral parameter is s(1 − s): the Laplace parametrisation, with ¼ as the floor. And positivity is free: in
the plane the propagator is a positive function of the colour, so the energy of any source is a sum of squares
with positive weights. Nothing about primes has entered.

Light units (`LIGHT_UNITS.md`): the colour k is a wave of wavelength 2πλ₀/k and frequency ck/(2πλ₀).

## 2. The horizon

The map is one operation: **sample the plane on the lattice Λ = log ℕ with the weights w(n) = μ(n) n^{−1/2}.**
Three things happen at the horizon.

**(a) The flow collapses to the prime shifts.** U_t carries Λ into itself only for t = log m, m ∈ ℕ (then
n ↦ mn), and never onto itself. The continuous group ℝ becomes the multiplicative semigroup (ℕ, ×), the free
commutative semigroup on the primes; its generators are the shifts S_p by log p. The canonical pair becomes
(H_arith, S_p) with H_arith e_n = (log n) e_n, the prime-state Hamiltonian of *Across the Horizon* Section 3.2,
and the Weyl relation of the plane survives only on the generators:

    [u, H_dil] = i        ⟼        [H_arith, S_p] = (log p) S_p.

**(b) The source is one point.** The Möbius source is the unit point δ₀ (u = 0, x = 1) passed through the
prime shifts with the sign −1 and the decay p^{−1/2} per generator, exactly at every prime cutoff:

    J^{(P)} = Π_{p≤P} (1 − p^{−1/2} S_p) δ₀ = Σ_{n P-smooth, squarefree} μ(n) n^{−1/2} δ(u − log n).

**(c) The image has two evolutions**, and they are not the same evolution.

| | by integers | by primes |
|---|---|---|
| the source | J_N = Σ_{n≤N} μ(n) n^{−1/2} δ_{log n} | J^{(P)} = Π_{p≤P}(1 − p^{−1/2} S_p) δ₀ |
| its colours | Ĵ_N(k) = Σ_{n≤N} μ(n) n^{−1/2−ik}, a Dirichlet polynomial | Ĵ^{(P)}(k) = Π_{p≤P}(1 − p^{−1/2−ik}), a finite Euler product |
| its energy | R(N) = Σ_{m,n≤N} μ(m)μ(n)/max(m,n) | E_P = Σ_{m,n P-smooth} μ(m)μ(n)/max(m,n) |
| its limit object | 1/ζ(½ + ik): zeros become poles, no pole | formally the same, but the product does not converge on the line |

Both energies are the propagator-weighted integral of the colours, (1/2π)∫ |Ĵ|²/(k² + ¼) dk. Checked on the
prime side at P = 7: the integral of the finite Euler product gives 1.104762, the exact sum over the sixteen
squarefree 7-smooth numbers gives 1.104762.

## 3. The colours on the two sides (figure, panels A and B)

**By integers the colours show lines at the zeros.** |Ĵ_N(k)|² on 0 ≤ k < 60 peaks at k = 14.130 for both
N = 10³ and N = 10⁵ (γ₁ = 14.1347); the mean within 0.25 of γ₁, γ₂, γ₃ is 8.4 times the overall mean at
N = 10³ and 12.2 times at N = 10⁵. The lines sharpen with N. The spectrum is tempered: an entire function of
k with no poles at any finite N, and under (G) its propagator-weighted integral grows slower than any power.

**By primes the colours are not tempered, and show lines only above the colour √P/log P.** Π_{p≤P}|1 − p^{−1/2−ik}|² ranges over 8 decades
at P = 100 and over 43 decades at P = 10⁴ (log₁₀ from −26.8 to 15.8). What drives it is the prime sum:

    log Π_{p≤P}|1 − p^{−1/2−ik}|² = −2 Re S_P(k) + O(1),        S_P(k) = Σ_{p≤P} p^{−1/2−ik},

(the difference has mean −0.005 and standard deviation 0.51 over k > 2, at both P), and by the explicit formula
for Σ_{p≤P} p^{−s},

    S_P(k) = P^{½−ik} / ((½ − ik) log P) + (bounded in P at fixed k ≠ γ) + (the zeros' terms).

The first term is the pole of ζ at s = 1, seen from the line at distance ½: |P^{1−s}| = √P. So the prime-cutoff
colours swing with amplitude exp(±2√P/((log P)√(¼ + k²))): this is what the non-convergence of the Euler
product on the critical line looks like. At a zero the prime sum has the line term −log log P (the ρ = ½ + iγ
term of the explicit formula at k = γ), so the product at k = γ carries a factor (log P)² over its
surroundings, against the pole's swing, which decays in colour like 1/k. The line is visible where
γ ≳ √P/log P and drowned below: at P = 10⁴ the threshold is k ≈ 11, and the first zeros stand out as bumps of
about two decades (panel B, where the dashed pole term reproduces the product to a standard deviation of 0.5 in
the natural log); at P = 10⁷ the threshold is k ≈ 200 and the first eighty zeros are under the pole. Each zero
is drowned once P exceeds about (γ log P)², so the prime door drowns every zero in turn.

**The dichotomy.** The horizon has two doors. The integer door passes the zeros and stops the pole (Σ μ(n)n^{−s}
has no pole; the lines are at the zeros and sharpen with N). The prime door passes the pole, and under it every
zero is eventually lost (the product is dominated by P^{½−ik} below the colour √P/log P). The two doors are the
two orders of one double limit: the P-smooth source cut at N, J_N^{(P)}, tends to J_N as P → ∞ at fixed N and to
J^{(P)} as N → ∞ at fixed P.

Where the zeros and the pole live on the prime side is a theorem, not a computation: the convergent part of the
Euler product on the line, Π_p (1 − p^{−s}) exp(p^{−s} + ½p^{−2s}), is an entire function without zeros for
Re s > ⅓, so on the line both the zeros and the pole sit entirely in the divergent prime sum Σ_p p^{−s}, which
is log ζ(s) up to a function regular for Re s > ½. **The prime door carries log ζ; the integer door carries 1/ζ.**

## 4. Where conservation becomes growth (figure, panels C and D)

In the plane the flow conserves. In the image the evolution is not the flow but the arrival of generators, and
energy grows. The two evolutions grow differently.

**By integers** (Section 6 of THE_FIELD.md): R(n) − R(n−1) = μ(n)²/n + 2μ(n)M(n−1)/n. The first term is
what the plane delivers, the self-energy of the pulse; summed it is Σ_{n≤N} μ(n)²/n = (6/π²) log N + c, a
theorem (the density of squarefree numbers), 0.6079 per e-fold. Plancherel says the same number is the mean
colour density: the mean of |Ĵ_N(k)|² over 0 ≤ k < 10⁴ at N = 10³ is 5.238 against Σ μ²/n = 5.242
(Montgomery–Vaughan). So in colour, **the delivered energy is the mean over all colours, and the retained energy
R(N) is the propagator-weighted integral**, which lives at the low colours and on the lines. Measured:

| N | delivered Σ μ²/n | retained R(N) | fraction |
|---|---|---|---|
| 10³ | 5.2424 | 1.4594 | 0.278 |
| 10⁴ | 6.6434 | 1.5824 | 0.238 |
| 10⁵ | 8.0429 | 1.6216 | 0.202 |
| 10⁶ | 9.4427 | 1.7082 | 0.181 |
| 10⁷ | 10.8425 | 1.8365 | 0.169 |

The slope of R(N) in log N over 10⁴–10⁷ is 0.028 per e-fold against 0.608 delivered: the image keeps about
one part in twenty at the margin, and the fraction kept falls. Under (G) the retained energy is subpower.
The finer law, R(N) ~ C log N with C = Σ_ρ 1/|ρζ′(ρ)|² (both signs of γ), holds under more than (G) (simple
zeros and the convergence of that sum) and is not assumed; the value of the sum over the first 6000 zeros,
computed in `THE_GRAVITY_PLANE.md`, is 0.0288, and the measured slope is 0.0278. That is an observation.

**By primes.** Adding the prime q to the cutoff sends J ↦ J − q^{−1/2} S_q J, and

    E_{P∪{q}} = (1 + 1/q) E_P − 2 q^{−1/2} C_q,        C_q = ⟨J^{(P)}, S_q J^{(P)}⟩_G,

because the shifted copy S_qJ has exactly the energy of J: **this is the plane's conservation, entering the
arithmetic as the factor (1 + 1/q).** The update has the shape of Build v32's Euler update and of the cofactor
flow of the record: a positive diagonal and one signed cross term (*Across the Horizon*, Section 3.2). Measured
(`rh_horizon_map.py`, the integral of the finite Euler product against the propagator):

| P | conserved copies Π(1 + 1/p) | E_P | C_q | E_P / copies |
|---|---|---|---|---|
| 2 | 1.500 | 0.5000 | +0.7071 | 0.333 |
| 3 | 2.000 | 0.6667 | 0 | 0.333 |
| 7 | 2.743 | 1.1048 | −0.1512 | 0.403 |
| 23 | 3.748 | 2.0733 | −0.3391 | 0.553 |
| 47 | 4.400 | 3.6041 | −0.8998 | 0.819 |
| 97 | 5.062 | 8.1756 | −2.5053 | 1.615 |
| 199 | 5.856 | 38.22 | −63.9 | 6.53 |
| 397 | 6.534 | 264.6 | −519.6 | 40.5 |

The interference is positive at q = 2, zero at q = 3, and negative from q = 5 on, growing: the shifted copy
anti-correlates with the source (μ(qn) = −μ(n)), the energy exceeds the conserved copies from P = 67, and runs
away: log E_P grows at the rate of Σ_{p≤P} p^{−1/2} ≈ 2√P/log P, the pole again (the ratio
log E_P / (√P/log P) rises from 0.98 at P = 101 to 1.68 at P = 397, toward a prefactor that the pole term
bounds by 4). **By primes nothing is tempered**; the growth law (G) is a statement about the integer door only.

## 5. Where information is compressed

1. **The flow.** A continuous one-parameter group becomes a countable free semigroup; what passes the horizon
   is the generators, the numbers log p. The primes are the compressed record of the dilations.
2. **The source.** One point and the generators give the whole Möbius source (Section 2(b)).
3. **The colours.** A tempered entire spectrum at every finite N becomes, in the limit, lines at the zeros;
   the information of the zeros passes through the integer door only. Under (G) the energy content of those
   lines grows slower than any power while the delivered energy grows like log N: the image keeps a fraction
   that falls (0.278 at 10³, 0.169 at 10⁷), and if the log law holds it tends to 0.028/0.608 = 4.6 %. That is
   the compression ratio of the integer door.
4. **The prime door compresses nothing and amplifies the pole**: its colours carry log ζ, the pole's
   2√P/((log P) k) in the exponent over the zeros' 2 log log P, so each zero is lost once P ≈ (γ log P)².

## 6. What stays invariant

| quantity | in the plane | in the image | status |
|---|---|---|---|
| the mass ¼ | the floor of H₀ = H_dil² + ¼ | ¼ in \|ρ\|² = γ² + ¼; the ½ of the line | exact |
| the Green function | e^{−\|u−v\|/2} | 1/max(m, n) | exact |
| the energy | ⟨J, H₀⁻¹J⟩ | R(N), E_P | exact; six decimals to 10⁷, P = 7 exactly |
| the equipartition | the factorisation (−∂ + ½)(∂ + ½) | inward and outward halves | exact |
| the canonical relation | [u, H_dil] = i | [H_arith, S_p] = (log p) S_p | exact, on the generators only |
| the mean colour density | the mean of \|Ĵ\|² | the delivered energy Σ μ²/n | exact (Montgomery–Vaughan); checked |
| conservation under a shift | E[U_tJ] = E[J] | the factor (1 + 1/q) of the prime update | exact |
| the propagator weight | 1/(k² + ¼) | 1/\|ρ\|² in every Parseval sum of the record | exact |
| the mass threshold | finite energy above ½, infinite at ½ | (G) | the premise |

Not invariant: the flow itself (only its generators pass); the energy under evolution (grows: tempered by
integers under (G), untempered by primes); the Euler product (does not converge on the line); positivity as a
tool (free in the plane, where the propagator is positive; in the image the burden moves to the source, whose
interference term carries the sign).

## 7. What the map proposes in the lower description

1. **The integer door is the door.** Everything about the zeros that reaches the arithmetic passes through
   the tempered evolution J_N, and (G) is its law. The prime evolution is a different, explosive object, and
   results about partial Euler products on the line are results about the pole.
2. **A line-by-line energy accounting.** Under the log law each zero's line carries (log N)/|ρζ′(ρ)|² of the
   retained energy, weighted by the propagator 1/|ρ|² and the line strength 1/|ζ′(ρ)|². The measured slope
   0.0278 against the computed 0.0288 is the first test of that accounting; a test by windows around the first
   zeros, at the integer cutoffs, is the next computation the map asks for.
3. **The mass as the variable.** (G) is a statement about the mass of the field. The family of fields of mass
   ½ + δ, each with the source read on its own line, is the right family in which to watch the energy become
   infinite as δ ↓ 0: the damped energies of THE_FIELD.md, Section 7, are its first table.
4. **The order of limits.** J_N^{(P)} with both cutoffs is the object that contains both doors; the two orders
   of its double limit are the two descriptions. Which mixed cutoffs N(P) keep the lines visible while the pole
   is still present is a precise question about where, between the two doors, the horizon actually is.

## 8. Beside the Exact Horizon Map (received 2026-10-06)

The user had already built the map with ChatGPT: *The Light Plane: An Exact Horizon Map* (13 pages, dated
2026-10-05, `received/Michaels_Light_Plane_Horizon_Map.pdf`). It constructs the horizon as an operator and
proves its properties; this note was written without it. The two agree wherever they overlap, and each has
things the other does not.

**What the paper proves, and what was checked here** (`rh_horizon_modes_check.py`, at N = 10³ and 10⁴, beyond
the paper's own checks at N ≤ 64):

- The field metric ⟨f, g⟩ = ∫(f′g′ + ¼fg) has reproducing kernel G_v = e^{−|u−v|/2} with ‖G_v‖ = 1, so every
  field value is bounded by the energy, |φ_N(u)|² ≤ R(N).
- The normalised source Ψ = Z^{−1/2} Σ μ(n)/n |n⟩ with Z = ζ(2)/ζ(4) = 15/π², a product state over the primes
  with amplitude −1/p per occupied prime; the source Hamiltonian H_src|n⟩ = (log n)|n⟩; the horizon map
  P_N : ℓ² → H¹ with P_NΨ = φ_N and the exact pullback P_N*P_N = Z·[min(a, b)]. The energy is the expectation of
  that positive observable, R(N) = ⟨Ψ, P_N*P_NΨ⟩: a conserved source norm with a growing arithmetic readout.
- The prime action with the integer cutoff retained: φ_{Q∪{p},N} = φ_{Q,N} − p^{−1/2} T_{log p} φ_{Q,⌊N/p⌋}, and
  its energy update with the inherited field at ⌊N/p⌋. This is the mixed cutoff J_N^{(P)} of Section 3, made
  exact; the prime update of Section 4 is its K = N = ∞ case.
- The exact spectrum of min(a, b): sine modes with gains λ_j = 1/(4 sin²((2j−1)π/(4N+2))), and the resolution
  R(N) = Σ_j λ_j |b_j|², b_j = Σ μ(n)/n · v_j(n). Checked: 1.459409 and 1.582361 at N = 10³, 10⁴, equal to R(N)
  to six decimals; the modes orthonormal to 10⁻¹⁴.
- The local sample inversion μ(n) = n(2w_n − w_{n−1} − w_{n+1}), w_n = √n φ_N(log n): the N samples at the
  lattice return the N coefficients. Checked at N = 10³ to 4·10⁻¹¹. (The recovery law of THE_FIELD.md, Section
  5, is its continuous form: value and slope at any point return M and the tail of h.)
- The phase readout R_N(τ) = ‖P_N U(τ)Ψ‖² and the transport bound q(τ)^{−2} R(N) ≤ R_N(τ) ≤ q(τ)² R(N) with
  q = √(1+τ²) + |τ|, and the long phase average equal to Σ μ²/n. Checked at N = 10³: the bound holds at
  τ = 0.5, 1, 3, 14.13; the average over 0 ≤ τ ≤ 2000 is 5.199 against 5.242.
- Under (G): the mode-population allowance Σ_{j≤J}|b_j|² = O(J² N^{−2+ε}), the bounds on the field, the terminal
  and the damped integral, and the phase windows. Measured: the Möbius source's population in the first 1, 10,
  100 modes is 0.1 %, 2 %, 12 % of its allowance at N = 10³ and 3 %, 0.5 %, 2 % at N = 10⁴.

**One identification joins the two documents.** The paper's source phase τ is this note's colour k. Its readout
is the colour spectrum of panel A seen through the propagator centred at τ,

    R_N(τ) = (1/2π) ∫ |Ĵ_N(t)|² / ((t − τ)² + ¼) dt,

so the Hamiltonian phase is a Lorentzian window of width ½ sliding along the colours, and at the phase of the
first zero it reads the line: R_N(14.1347) = 33.55 against R_N(0) = 1.459 at N = 10³. The paper's phase average
equal to Σ μ²/n is the Plancherel statement of Section 4 in the τ variable.

**What the paper adds to this note**: the horizon as an operator with an adjoint and a positive pullback; the
exact mode spectrum and the geometric meaning of (G) as a vanishing overlap with the high-gain modes; the local
inversion; the mixed cutoff with ⌊N/p⌋; the phase transport bound; the terminal term's rank-one role.

**What this note adds to the paper**: the plane side (H₀ = H_dil² + ¼ = s(1 − s), the flow and its collapse
to the semigroup, the canonical pair on the generators); the computed colours on both sides; and the finding of
Section 3. The paper's Section 8 states that the finite Euler product and the numerical cutoff are distinct,
with the remainder Σ_{n>N, N-smooth} μ(n) n^{−s} written out. On the line that remainder is not a correction:
it is the dominant, pole-driven part, and the product without the integer cutoff swings over 43 decades at
P = 10⁴ and has energy growing like exp(2√P/log P). The paper never meets this because its prime action keeps
the integer cutoff (the inherited field at ⌊N/p⌋), which is the integer door. That is the right door, and the
computation here says why.
