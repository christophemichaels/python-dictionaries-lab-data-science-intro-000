# Global closure and Connes geometry: the light plane, the scaling site, and the arithmetic return

2026-10-06. On the directive `received/Michaels_Light_Plane_Global_Closure_and_Connes_Geometry_Prompt.pdf` (19 pages).
Script `rh_global_closure.py` (data `data/global_closure.json`), figure `rh_global_closure_fig.py` → `GLOBAL_CLOSURE.png`.
Sources read for this note: [CC-S] A. Connes, C. Consani, *The scaling site*, arXiv:1507.05818v2 (fetched; the JHU-hosted
copy returned 403), and [CC-A] *Geometry of the arithmetic site*, arXiv:1502.05580 (fetched, used for the definition of the
site and its points); [C-T] was not fetched and no statement from it is used. In the arXiv version of [CC-S] the
Riemann–Roch theorem is Theorem 6.7 and Definition 6.6 is the norm; the directive's "Theorem 6.6" is read as that theorem.
Status labels, as required: **[source]** proved in a supplied or cited source, **[derived]** proved in this investigation,
**[finite]** finite diagnostic, **[hypothesis]** research hypothesis. Notation as in the directive: R(N) = ‖φ_N‖²_E, M, h,
S(X), c_N(n) = μ(n)/n, V_N(x) = Σ_{n≤N} μ(n)/n·min(x, n), f_N(u) = V_N(eᵘ), Δ_sc = ∂_u² − ∂_u.

## 1. Strongest result

No inequality of the forms (27), (28), (29), (31), (32)–(33) or (SP) is obtained from the geometry. The sharpest completed
results are four theorems, each with its proof status, and a finite return.

**Theorem 1 (the ambient object; the Möbius family is the inverse of the orbit sum). [derived]** Let G₁(x) = min(x, 1) on
[0, ∞) and, for n ∈ ℕ^×, let ρ_n act on potentials by (ρ_n V)(x) = V(x/n), with the horizon rule ρ_n: V_N ↦ V_{⌊N/n⌋}(·/n)
(ρ_n is the inverse direction of the action γ_n(ξ)(λ) = ξ(nλ) of [CC-S], eq. (2)). Then for every N ≥ 1:

    (i)   V_N = Σ_{n≤N} μ(n) ρ_n G₁;
    (ii)  Σ_{d≥1} ρ_d V_{⌊N/d⌋} = G₁                                   (Möbius inversion, horizon compatible);
    (iii) Σ_{k≥0} ρ_{p^k} V_{⌊N/p^k⌋} = V_N^{(p)}  for every prime p   (the p-free potential), and
          V_{Q∪{p},N} = V_{Q,N} − ρ_p V_{Q,⌊N/p⌋}                      (the family action (46) in potential coordinates);
    (iv)  −Δ_sc f_N = Σ_{n≤N} μ(n) δ_{log n}; in the sense of [CC-S] Definition 6.2 the order of V_N at λ = n is
          n·(V_N′(n⁺) − V_N′(n⁻)) = −μ(n), and the degree of the divisor is Σ orders = −M(N);
    (v)   R(N) = ∫₀^∞ |V_N′(x)|² dx = h(N)² + Σ_{k<N} (h(N) − h(k))², and ∫₀^∞ |G₁′|² dx = 1.

Proof. (i) is (18). (ii): the coefficient of ρ_n G₁ in Σ_d ρ_d V_{⌊N/d⌋} = Σ_d Σ_{m≤N/d} μ(m) ρ_{dm} G₁ is Σ_{d|n} μ(n/d) =
[n = 1]. (iii): the same computation with d = p^k gives Σ_{k≤a} μ(p^{a−k}m) = μ(m)[a = 0] for n = p^a m, p ∤ m; the
family action is the term n = pm of the definition, μ(pm)/(pm)·min(x, pm) = −μ(m)/m·min(x/p, m). (iv): V_N is affine
between integers with slope h(N) − h(⌊x⌋), the slope changes by −μ(n)/n at x = n, and Δ_sc(V(eᵘ)) = x²V″(x) with
x²δ(x − n) = n δ_{log n}. (v): ∫|V′|² = Σ_k (h(N) − h(k))²·1 + h(N)²·1 (cells [k, k+1) and [0, 1)), and Abel summation
turns it into (1). Verified: (ii) and (iii) exactly in rationals at N = 60 on all half-integers and in floating point at
N = 10⁵ (1.8·10⁻⁹); (iv) at N = 50 (jumps −μ(n) to 3·10⁻¹⁶); (v) exactly at N = 50 (R(50) = 1.479629309286). ∎

Consequence. The exact ambient invariant is the unit potential: energy 1 and degree 1 at every horizon, conserved by (ii).
The readout R(N) is the energy of the Möbius inverse of the orbit-sum operator Σ: (V_{⌊N/d⌋})_d ↦ G₁ on horizon-compatible
families. The inequality (27), R(N) ≤ A_ε N^ε·L_ε(Ξ), with L_ε(Ξ) = 1 the invariant, is the statement that this inverse has
norm O(N^{ε/2}) on the Möbius family, which is the Goal (G). The geometry supplies the invariant and the exact inversion
and no bound on the inverse; Section 7 displays this as the resume point.

**Theorem 2 (the fold onto Cₚ: the circle carries the degree ledger, not the energy). [derived, on [CC-S] Lemma 6.3,
Lemma 6.4, Theorem 6.5]** Let p be prime and N ≥ 1 with M(N) = 0. The two-sided fold F_{p,N}(x) = Σ_{k∈ℤ} V_N(p^k x)
converges for every x > 0, satisfies F(px) = F(x), is continuous and piecewise affine, and (p − 1)·c_p(N)·F_{p,N} is a
global section of the sheaf of quotients K_p on C_p = ℝ₊^*/p^ℤ of [CC-S] Lemma 6.4 (piecewise affine, slopes in
H_p = ℤ[1/p]), where c_p(N) is the prime-to-p part of lcm(1, …, N). Its divisor on C_p is: order −(p − 1)c_p(N)μ(m) at the
class [m] = m·p^ℤ of every unfinished parent (m ≤ N < pm, p ∤ m), and order 0 at the class of every complete p-packet
{m, pm} ⊂ [1, N], whose fold is the constant μ(m); so F = M^{(p)}(⌊N/p⌋) + (regularized folds of the unfinished parents),
with M^{(p)} the Mertens function of the p-free integers. The degree is −(p − 1)c_p(N)·Σ_{unfinished} μ(m) =
−(p − 1)c_p(N)M(N) = 0, as [CC-S] Lemma 6.4(ii) requires of a global section. If M(N) ≠ 0 the partial folds grow by M(N)
per unit of k and there is no global section: this is the balance condition (PB) of the directive, and (PC) is the
statement that complete packets fold to constants. The factor p − 1 is needed (N = 2, p = 5: the slope on (1, 2) is
−3/8, and 2·(−3/8) ∉ ℤ[1/5]). The archimedean energy ∫ e^{−u}|F′|² du over one period is multiplied by p^{−1} when the
period is translated by log p, so it is not a function on C_p; the norm of [CC-S] Definition 6.6 is max |h(λ)|_p/λ, the
p-adic size of the slope, and the Riemann–Roch dimension of [CC-S] Theorem 6.7 equals the real degree. Neither involves the
L² norm of the slope in x, which is R(N).

Proof. Convergence: V_N(y) = y·h(N) for y ≤ 1 and V_N(y) = M(N) = 0 for y ≥ N, so the terms with p^k x < 1 form a geometric
series and the terms with p^k x ≥ N vanish. Periodicity is the index shift. The slope on a cell is F′(x) = Σ_k p^k
V_N′(p^k x); for x ∈ [1, p) this is h(N)/(p − 1) + Σ_{k≥0, p^k x<N} p^k(h(N) − h(⌊p^k x⌋)), whose denominators divide
(p − 1)·lcm(1..N), hence (p − 1)c_p(N)F′ ∈ ℤ[1/p]. Orders: the breakpoint of V_N(p^k x) at x = n/p^k has slope jump
p^k·(−μ(n)/n) and order (n/p^k)·p^k·(−μ(n)/n) = −μ(n); the class [m] contains the breakpoints n = m and n = pm, and
−μ(m) − μ(pm)·[pm ≤ N] = −μ(m)[pm > N]. The fold of a complete packet, Σ_k [G₁(p^k x/m) − G₁(p^k x/(pm))], telescopes to
lim_{k→∞} G₁(p^k x/m) − lim_{k→−∞} = 1. The energy scaling is e^{−(u+log p)} = p^{−1}e^{−u} with F′ periodic. ∎
Verified at N = 39, p = 5 (unfinished parents 11, 13, 14, 17, 19, 21, 22, 23, 26, 29, 31, 33, 34, 37, 38, 39, charge 0;
complete parents 1, 2, 3, 6, 7, constant −1 = M^{(5)}(7)): periodicity 9·10⁻¹⁶, orders at every class to 3·10⁻⁹, degree
0, slopes with denominators 1 after the rescaling, energy ratio between consecutive periods 0.200000; at N = 38
(M = −1) the partial folds −5.59, −10.59, −20.59 at x = 2 for k < 5, 10, 20. Figure, panels A and B.

**Theorem 3 (the local bridge and the domain facts). [derived]** (17)–(22) hold as displayed in the directive: f_N =
e^{u/2}φ_N, −Δ_sc f_N = Σμ(n)δ_{log n}, and on every finite interval ∫_a^b (|φ′|² + φ²/4) du = ∫_a^b e^{−u}|f′|² du −
½[e^{−u}f²]_a^b, the boundary term telescoping across the breakpoints because f is continuous. On the twisted domain
ψ(u + log p) = p^{−1/2}ψ(u) of (23), the Wronskian form (24) survives with the factor q² − 1 = p^{−1} − 1, and the modes
e^{−u/2}e^{iκu} have H₀-eigenvalue κ² + iκ: H₀ is not symmetric there. A continuous function with periods log 2 and
log 3 is constant (ℤ log 2 + ℤ log 3 is dense: log 3/log 2 is irrational since 2^a ≠ 3^b), so imposing (23) for two primes
on one scalar ψ makes e^{u/2}ψ constant. The compatible all-prime object is the commuting family ρ_p of Theorem 1 with
the inherited horizons ⌊N/(pq)⌋: it keeps every prime's information because the orbit sums over different primes are
distinct operators composing to the full orbit sum (ii). The Möbius potential has slopes h(N) − h(k) ∈ ℚ: it is not a
section of K_p on any C_p (slopes not in ℤ[1/p], not periodic), and lcm(1..N)·V_N is a section of the sheaf of quotients K of
[CC-S] Proposition 6.1 over every bounded open interval of [0, ∞), with integer slopes; the germ of V_N at the rational point p_ℚ (the
stalk of K at a rank-one subgroup H has slopes in H, [CC-S] Theorem 4.2(i) and Proposition 6.1) is a section without
rescaling. Verified: (19) to 3·10⁻¹⁵ and (7) to 5·10⁻¹⁶ at N = 50; (21) on three intervals to ten digits with the left side
computed by Gauss–Legendre quadrature of the Green-function sum and the right side from the potential, and on [−30, 30]
reproducing R(50); (23), (24), (25) numerically; the translation scaling 0.200000 = 1/5 and the half-line sum cell·p/(p−1);
(46) in potential coordinates to 2·10⁻¹⁶ at Q = {2, 3}, p = 5, N = 1000; the compactified metric (26) at N = 10
(1.061904762 = R(10)); the spectral identity (6) at N = 20 (1.601733 with the tail estimate, R(20) = 1.601733).

**Theorem 4 (the exact return with a positive storage term). [derived; statements (37), (42) of the directive]** Let
N₀ = N, N_{j+1} = ⌊N_j/6⌋, N_J < 6, and let the mesh of N be refined by the chain points. With the two-moment projection
on this mesh, ε_j := E_res(N_j) − E_res(N_{j+1}) ≥ 0 and Σ_j ε_j = E_res(N) − E_res(N_J) ≤ C₀ = 2^{−4/3}, because the
residual is a sum over mesh intervals of within-interval energies and the intervals of N_{j+1} are those of N_j below
N_{j+1}. Hence, with Δ̂_j = E_P(N_j) − E_P(N_{j+1}) and d₂₃(N) = D₂₃(N) − D₂₃(⌊N/6⌋),

    R(N) ≤ 43/30 + D₂₃(N) + B̂₂₃(N) + C₀,     B̂₂₃(N) = Σ_j [Δ̂_j − d₂₃(N_j)]₊,     D₂₃(N) = (2/π²) log N + O(1).

Moreover 0 ≤ d₂₃(N) ≤ 4/5 for every N ≥ 6, with equality exactly for 30 ≤ N ≤ 35: for a = ⌊L/6⌋ ≥ 11 (N ≥ 396),
Σ_{a<m≤L, (m,6)=1} 1/m ≤ (1/3)log((L+1)/(a−5)) + 1/(a+1) ≤ (1/3)log(6(a+1)/(a−5)) + 1/(a+1) ≤ 0.912, so d₂₃ ≤ 0.608
there; below 396 the maximum is computed. Finite: B̂₂₃(N) = 0 on the chains from 10⁶ and 10⁷ (Δ̂_j between −0.14 and
+0.14 against d₂₃(N_j) ≥ 0.35 at every step; ε_j ≤ 0.0021, Σε_j = 0.0056 and 0.0034), so (42) reads 1.708 ≤ 4.717 and
1.837 ≤ 5.184. The unprojected all-horizon budget B₂₃(N) = Σ_j[Δ₆R(N_j) − d₂₃(N_j)]₊ ≥ B̂₂₃(N) is at most 0.5526 for every
N ≤ 10⁷ (at N = 13), vanishes for 73 % of the horizons, and takes its positive parts only from chain elements ≤ 65; hence
R(N) ≤ D₂₃(N) + 1.053 for every N ≤ 10⁷ (the maximum of R − D₂₃ is 1.0526 at N = 13). Figure, panel C.

Arithmetic consequence of the four theorems: none beyond the exact identities. (G) and the Mellin conclusion remain
conditional on the single input displayed in Section 7.

## 2. Overall geometry

| Role (directive, Section 7) | Object identified | Status |
|---|---|---|
| Overall light plane L | The scaling site [0, ∞) ⋊ ℕ^× of [CC-S] Definition 2.2: sheaves on it are the ℕ^×-equivariant sheaves on [0, ∞) (Prop. 2.1); its points are the adele classes mod Ẑ^* (Theorems 3.1, 5.1): support {0} ↔ points of N̂^×, support ≠ {0} ↔ rank-one subgroups λH of ℝ. The structure sheaf O: convex piecewise affine functions with integer slopes; its quotient sheaf K: piecewise affine functions (Prop. 6.1). | [source] |
| Admissible arithmetic state | The unit potential G₁ = min(λ, 1), with −G₁ = max(−λ, −1) a global section of O; its ℕ^×-orbit {ρ_n G₁}; energy 1, degree 1. | [derived] |
| Arithmetic view A | The horizon family {V_N}: the Möbius inverse of the orbit of G₁ (Theorem 1(i)–(ii)). It recovers μ(n) as the orders at λ = n, M(N) as the value for x ≥ N (the degree), h(N) as the slope at 0, and the explicit-formula side as the Mellin transform ∫₀^∞ V_N′(x)x^{s−1}dx = D_N(1 − s)/s, whose Plancherel form is (6). | [derived] |
| Projection / observation | V_N ↦ the pair (M(k), h(N) − h(k))_{k≤N} ↦ φ_N = e^{−u/2}f_N = P_NΨ, with ‖φ_N‖²_E = R(N) (Theorem 1(v), (19), (22)). It is injective at every finite horizon (the potential determines μ(n), n ≤ N); there are no fibres. The "additional direction" of the drawing is the orbit coordinate d of ρ_d in (ii) (one coordinate encoding infinitely many modes: a finite set of displayed coordinates, infinitely many arithmetic modes); the "2D projection" is the potential at one horizon; the structure is base and coordinate, not boundary and bulk nor a quotient. | [derived] |
| Connes-related geometry | The periodic orbits C_p = ℝ₊^*/p^ℤ with the sheaves O_p (convex, slopes in ℤ[1/p], Lemma 6.3) and K_p (Lemma 6.4), degree-zero law for global sections (Lemma 6.4(ii)), Jacobian ℤ/(p−1)ℤ (Theorem 6.5), p-adic slope norm (Definition 6.6), real Riemann–Roch (Theorem 6.7), and the operator Δ′ = ∂_u² − ∂_u = λ²∂_λ² (last page). The arithmetic site [CC-A] is the same ℕ^×-structure over the Boolean semifield; [CC-S] is its scalar extension to ℝ^max₊. The map used: η_p: ℝ₊^*/p^ℤ → C_p (Lemma 6.3) and the fold of Theorem 2. No theorem of [C-T] is transferred. | [source] for the objects, [derived] for the fold |
| Marked points / two ends | u → −∞ carries the slope h(N) (exterior energy h(N)²/2; h(N) → 0 by the prime number theorem); u → +∞ carries the value M(N) (exterior energy M(N)²/2N). Under v = tanh(u/2) they are v = ∓1 with the metric (26), verified at N = 10. They are spatial infinities of ℝ, not boundary components of a compact curve. The "bounded-looking central region" is [0, log N], unbounded as N grows. | [derived] |
| Global control | The exact invariant: Σ_d ρ_d V_{⌊N/d⌋} = G₁ with energy 1 (Theorem 1(ii)); the degree-zero law on C_p (Theorem 2); the storage ledger (Theorem 4). No positive form, compactness or comparison theorem bounding the inverse of the orbit sum is available in [CC-S] or [CC-A]. | [derived]; the bound is the [hypothesis] |
| Return | Theorem 4 and Section 5. | [derived] |

Answers to the eight questions of Section 7 of the directive. (1) The orbit coordinate d carries the dilation by d with
the inherited horizon ⌊N/d⌋; nothing else resides in it. (2) No global relation excluding high-gain source directions was
found: the invariant constrains the orbit sum, which mixes every horizon, not the state at one horizon. (3) The
alternating signs are the Möbius inversion (i); the squarefree rule is μ(p²m) = 0; the cutoffs are the horizon rule of
ρ_n. (4) The constant quantity is the energy and the degree of the unit potential, both 1. (5) Theorem 1(ii) relates it to
every horizon, as an identity. (6) The circle transitions act on the fold, where the signed cross terms between a
complete packet and anything else vanish (the packet folds to a constant); the cross terms of the unfinished part are
preserved on the circle but their energy is not a circle quantity (Theorem 2). (7) The infinite prime family is a
compatible system of commuting operators ρ_p with the uniform identity (iii); its completion over all primes is (ii). The
complete prime products of Section 10 of the directive have the exact energies 1/2, 2/3, 13/15, 116/105, 1433/1155,
21436/15015 for the first six primes, and R_{Q∪{p}}(N) = R_Q(N) + p^{−1}R_Q(⌊N/p⌋) + cross with cross = 0, 1/15, 0.1143,
0.0355, 0.0915, 0.0915, 0.164 for p = 3, …, 19 at N = 10⁷ (all complete there); the growth log E_P ~ 4√P/log P is the
audit's result [source]. (8) At each marked point the exterior energy is determined by the source (h(N) and M(N));
repeated passage through the horizon is the chain of Theorem 4, with the storage ε_j ≥ 0 paid once.

Ledger of the distinct quantities (directive, Section 4). Source norm ‖Ψ‖_{ℓ²} = 1, conserved [source]. Source
evolution U(τ): unitary, with the phase theorem (14) at constant q(τ)² for fixed τ [source]. Field energy R(N) = ‖P_NΨ‖²_E,
the quantity to bound; R(10⁷) = 1.836533426 [finite]. Horizon gain ‖P_N‖² ~ 4ZN²/π² on arbitrary directions [source].
Change of coordinates ψ ↔ f: the exact identity (21) with the boundary term; v = tanh(u/2): the metric (26) [derived].
New ambient invariant: the unit potential, energy 1 and degree 1 [derived]. "Constant growth": the measured derivative
of R with respect to log N is 0.029–0.031 per unit of log N over 10³ to 10⁷ [finite]; a constant derivative with respect
to log N means R(N) = O(log N), power 0 in N, with the consequence M(N)² ≤ N·R(N) = O(N log N); it is not proved.

## 3. Operator and domain correspondence

The local identities. For f = e^{u/2}ψ: f′ = e^{u/2}(ψ′ + ψ/2), f″ = e^{u/2}(ψ″ + ψ′ + ψ/4), so (∂_u² − ∂_u)f =
e^{u/2}(ψ″ − ψ/4) and e^{−u/2}(−Δ_sc)e^{u/2}ψ = −ψ″ + ψ/4 = H₀ψ: this is (17), an equality of differential expressions
[derived]. For V piecewise affine, Δ_sc(V(eᵘ)) = x²V″(x) away from breakpoints and, at a breakpoint n with slope change
σ, the distribution x²·σδ(x − n) = nσ δ_{log n}; with σ = −μ(n)/n this is (20). The energy identity (21): e^{−u}|f′|² =
|ψ′ + ψ/2|² = |ψ′|² + |ψ|²/4 + ½(|ψ|²)′, and ∫_a^b ½(|ψ|²)′ = ½[|ψ|²]_a^b = ½[e^{−u}|f|²]_a^b; on a cell where V = α + βx
both sides equal (α²/2)(e^{−a} − e^{−b}) + (β²/2)(e^b − e^a), and the boundary terms telescope across breakpoints since f
is continuous there (ψ′ jumps, |ψ|² does not) [derived]. The full-line boundary term vanishes because e^{−u}f² → 0 at
both ends (f ~ e^u h(N) at −∞, f → M(N) at +∞), giving (22) [derived]. All verified in Section 1 (Theorem 3).

The circle domain. On the smooth twisted domain (23) with q = p^{−1/2}: ψ(ℓ) = qψ(0), ψ′(ℓ) = qψ′(0) for ψ = e^{−u/2}f with f
periodic, so [ψχ′ − ψ′χ]_0^ℓ = (q² − 1)(ψ(0)χ′(0) − ψ′(0)χ(0)), nonzero in general: H₀ is not symmetric on this domain,
consistently with the complex eigenvalues κ² + iκ of (25) [derived; verified]. A Hilbert-space spectral use of Δ′ on C_p
therefore needs its own domain; [CC-S] uses Δ′ only to characterize piecewise affinity (its last page), and makes no
spectral statement. The raw ψ-energy of a periodic f over consecutive periods decreases by the factor p^{−1}; the
half-line total is the cell energy times p/(p − 1); the energy over ℝ of a nonconstant periodic f diverges at −∞. So a
scalar circle energy does not exist for the ψ-metric; the quantity that lives on C_p is the divisor (order, degree,
χ-invariant) of the fold [derived; verified]. The compensation the directive asks for in 6.1 is therefore not a metric
transition: the half-density e^{−u/2} is exactly what converts the unit-modulus action ρ_p on potentials (coefficient
−1, pure dilation, Theorem 1(iii)) into the p^{−1/2}T_{log p} of (46) and the twist of (23); both factors have this one
origin, and the minus sign is μ(p) = −1, the smaller horizon is pm ≤ N [derived; (46) verified].

Prime information. The two-period scalar model collapses (Theorem 3), so a single scalar function with all the twists
(23) is constant. The object that keeps the primes apart is the ℕ^×-action itself: distinct commuting operators ρ_p,
one orbit sum per prime (iii), composing to (ii). Its topology is that of the scaling site (equivariant sheaves on
[0, ∞)); its representation on the arithmetic fields is ρ_n on potentials, equivalently n^{−1/2}T_{log n} with horizon
⌊N/n⌋ on the fields. The slope groups: V_N has slopes in (1/lcm(1..N))ℤ, a section of K over (0, ∞) after the integer
rescaling, and a section at the rational point without it; on C_p only the fold of a degree-zero potential is a section,
after the rescaling by (p − 1)c_p(N) (Theorem 2). This is the global section result: **the Möbius potential is a section
of the quotient sheaf of the scaling site over every bounded interval (after an integer rescaling), and its fold is a
global section of K_p on the orbit C_p exactly when M(N) = 0, with divisor carried by the unfinished p-packets.**

## 4. The bounding mechanism

What the geometry supplies, with its constants: the invariant energy 1 and degree 1 of G₁ (Theorem 1); the degree-zero
law on every C_p, which the complete packets satisfy and the horizon truncation breaks by exactly M(N) (Theorem 2); the
storage ε_j ≥ 0 with total allowance C₀ = 2^{−4/3} = 0.3969 and the diagonal allowance 4/5 per step, 2/π² per unit of
log N (Theorem 4). Where the actual Möbius signs enter: only in the inversion (i)–(ii). Every other identity of this note
holds verbatim for the squarefree-positive control μ²: the bridge (17)–(22), the fold (Theorem 2 with μ replaced by μ²,
where the condition becomes Q(N) = 0, never satisfied, so no fold converges), and the ledger of Theorem 4 with its own
diagonal. The control's orbit sum is Σ_d ρ_d V⁺_{⌊N/d⌋} = Σ_{n≤N} 2^{ω(n)}ρ_n G₁ instead of G₁, and its energy is
2(6/π²)²N(1 + o(1)) (ratio 1.000005 at 10⁷) [finite; the asymptotics derived in COLORS_AND_DESCENT.md]. The fixed-sign
control a_C has orbit sum Σ_n (Σ_{d|n} a_C(n/d))ρ_n G₁ and sixth-scale increments (5/3)α_C²N (ratio 1.0000 at 10⁶ and
10⁷) [finite]. So the structural hypothesis that distinguishes the Möbius state is the inversion identity (ii), and the
theorem that would close is a bound on the inverse of the orbit sum restricted to horizon-compatible families, in the
metric (3). It is not in [CC-S]: the Riemann–Roch theorem is an equality of real dimensions and degrees, the Jacobian is
finite, and the norm is p-adic; no positive intersection form on C_p or on the square of the arithmetic site ([CC-A],
Section 6) is constructed there, and the finite-field calibration (Hodge index on the surface C × C) has no transferred
analogue with matched spaces. The attempts, each specific:

- (29), low-mode population from ambient constraints: the invariant constrains Σ_d ρ_d V_{⌊N/d⌋}, which mixes all horizons;
  no constraint on b_{j,N} at one horizon follows. Not obtained.
- (31), bounded damped energies: equivalent to (G) (THE_FIELD.md, §7; HORIZON_ROUND_TRIP.md). Not independent.
- (32)–(33), storage and increments: Theorem 4 is the storage ledger with Q_N = E_res(N) ≥ 0 and allowance C₀; the
  increment estimate reduces to B̂₂₃(N) ≤ C_ε N^ε[1 + S(N^{1/6})]^γ, the arithmetic input of CUBE_MOVES.md.
- (SP) and Weil positivity: H₀ on H¹(ℝ) and the diagonal operator log n have positive realizations and no zero modes;
  the zeros enter only through D_N(s) → 1/ζ(s), i.e. through the inverse of the orbit sum (the Mellin side of (ii) is
  Σ_d d^{s−1}D_{⌊N/d⌋}(1 − s) = 1). No faithful positive realization of the zeros is constructed here, and no theorem of
  [C-T] is invoked.

## 5. Complete return

Conditional theorem [derived from (49)–(52) of the directive and Theorem 4]. If B̂₂₃(N) ≤ C_ε N^ε for every ε > 0 (or
the amortized positive-increment allowance of CUBE_MOVES.md §2), then R(N) ≤ 43/30 + (2/π²)log N + O(1) + C_εN^ε + C₀,
hence (G); then |M(N)|² ≤ N·R(N) gives M(x) = O_η(x^{1/2+η}); the Mellin representation 1/ζ(s) = s∫₁^∞ M(x)x^{−s−1}dx,
absolutely convergent for ℜs > 1 and holomorphic on ℜs > 1/2 by the Mertens bound, continues 1/ζ there, so ζ has no zero
with ℜs > 1/2, and the functional equation excludes ℜs < 1/2. Unconditional: R(N) ≤ 43/30 + D₂₃(N) + B̂₂₃(N) + C₀ exactly,
and the elementary S(X) ≤ 2X. Finite [finite]: R(N) ≤ D₂₃(N) + 1.053 for every N ≤ 10⁷; along the chains from 10⁶ and 10⁷
the retained budget B̂₂₃ is zero. No partial saving in R(N) or in the mode populations beyond the previous baseline is
obtained from the geometry.

## 6. Verification

`rh_global_closure.py` (49 s to 10⁷): checkpoints R(10³) = 1.4594094176, R(10⁴) = 1.5823613269, R(10⁵) = 1.6215665260,
R(10⁶) = 1.7082491339, R(10⁷) = 1.8365334260, R(16384) = 1.6034402209, original-mesh E_P(16384) = 1.5883431079039456,
E_res = 0.0150971130273582 (the directive's values, the projection pair to fourteen digits; the two are the integer sixth-root mesh with 382
intervals); the separated-family calibration (53) is in CUBE_MOVES.md. Section 2: (19), (7), (20), (22), (21) as in
Theorem 3. Section 3: (23), (24), (25), the translation scaling, the two-period convergents (|1054 log 2 − 665 log 3| =
4.4·10⁻⁵). Section 4: (46) in potential coordinates; the two-prime and all-prime assembly with every cross term at
N = 10⁴ and 10⁷. Section 5: the orbit sums exactly at N = 60 and in floats at 10⁵; the fold at N = 39, p = 5 with its
orders, degree, slopes and energy ratio, and the divergence at N = 38. Section 6: (6) at N = 20 and (26) at N = 10.
Section 7: the chain-refined ledger at 10⁶ and 10⁷, the all-horizon budget B₂₃(N). Section 8: d₂₃(N) ≤ 4/5 at every
N ≤ 10⁷, equality on 30–35, the analytic bound above 396. Section 9: the controls. Controls that identify the arithmetic
content: μ² passes every identity and fails the inversion (orbit sum 2^{ω}); a_C likewise; the scalar circle model shows
the metric and domain facts; the two-period model shows what a scalar gluing loses; the random and shuffle controls of
CUBE_MOVES.md remain the comparison for the increments.

## 7. The precise resume point

The remaining equation, fully specified, is Theorem 1(i)–(ii):

    V_N = Σ_{n≤N} μ(n) ρ_n G₁,     Σ_{d≥1} ρ_d V_{⌊N/d⌋} = G₁,     R(N) = ∫₀^∞ |V_N′|² dx,     ∫₀^∞ |G₁′|² dx = 1,

and the estimate still needed is ‖V_N‖ ≤ A_ε N^{ε/2}‖G₁‖ for every ε > 0, a bound on the inverse of the orbit-sum operator
on horizon-compatible families in the metric (3); equivalently (G). Classification, as the directive asks: it is not a
definition still required (every object above is defined and verified); it is not a domain mismatch (Theorem 2 settles the
circle: C_p carries the degree ledger of the fold and no archimedean energy, which is a proved fact about scalar circle
energies, and the ℕ^×-object of Theorem 1 is well defined without any circle); the unproved geometric theorem would be a
positive form on the scaling site whose restriction to Möbius inverses is bounded, which [CC-S] and [CC-A] do not contain;
and the unproved arithmetic inequality is the displayed bound, which is the Goal itself. The geometry has supplied the
invariant and the exact ledger; it has not supplied the inequality.
