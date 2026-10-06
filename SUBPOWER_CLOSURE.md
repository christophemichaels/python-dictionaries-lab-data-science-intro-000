# The subpower closure attempt through the horizon modes (2026-10-06)

The research directive (`received/Michaels_Arithmophysics_Subpower_Closure_Prompt.pdf`) asks for a direct
analytical attempt at the closing target (T) through the exact horizon map, carried out rather than planned:
verify (A)–(D), evaluate the transfer matrix, investigate the cancellation in complete prime packets, and
either prove (T) or return the proved identities with one explicit residual. This note is the attempt. Script
`rh_subpower_closure.py` (N = 2000 and 4000; data `data/subpower_closure_{2000,4000}.json`), figure
`SUBPOWER_CLOSURE.png`.

**Dependency statement, once.** (G) ⟺ (T) ⟺ (T′) by (A) and the comparability λ_j ≍ N²/j², and (T) ⟺ RH
through M(N)² ≤ N R(N) and the classical criterion. Nothing below proves any of them. Everything marked
"theorem" or "proposition" is unconditional; the computations use no premise; (G) is never used as an input.
The result of the attempt is stated in Section 6: no closing mechanism was found, two proposed mechanisms are
shown to fail for a precise reason, and the residual is stated in Section 7.

## 1. (A)–(D), resolved and verified

All vectors are real, so conjugation is trivial. Indices: n, m are lattice points, j, ℓ are modes; V_N has
columns v_{j,N}; B_N = V_N Λ_N V_N^T with B_N = [min(a, b)].

**(A).** K_N ≤ R(N) since R = Σ_j λ_j |b_j|² ≥ λ_J Σ_{j≤J}|b_j|² for each J (gains decreasing). The upper
bound by Abel summation, R = λ_N S(N) + Σ_{j<N}(λ_j − λ_{j+1})S(j), S(j) ≤ K_N/λ_j, and 1 − r ≤ −log r. The
multiplier 1 + log(λ_1/λ_N) is 16.69 at N = 2000 and 18.07 at N = 4000, about 2 log N + 1.5. Verified: at
N = 2000, 0.394460 ≤ 1.490951 ≤ 6.5817; at N = 4000, 0.394448 ≤ 1.524962 ≤ 7.1283.

**(B).** For p ∉ Q and K = ⌊N/p⌋: c_{Q∪{p},N} = c_{Q,N} − p^{−1} E_p c_{Q,K}, because the Q∪{p}-smooth
squarefree n ≤ N not Q-smooth are exactly n = pm with m Q-smooth squarefree, m ≤ K, and μ(pm)/(pm) =
−p^{−1}·μ_Q(m)/m. Applying V_N^T gives (B) with T_{p;N,K} = V_N^T E_p V_K. In the computation this is used in
the form T b_{Q,K} = V_N^T E_p c_{Q,K} (since V_K V_K^T = I), so no small-horizon modes are needed. Verified to
10⁻¹⁶ at every prime step at both N.

**(C).** T^T T = V_K^T E^T V_N V_N^T E V_K = V_K^T E^T E V_K = I_K. And E^T B_N E has entries min(pm, pm′) =
p min(m, m′), so T^T Λ_N T = V_K^T E^T B_N E V_K = p Λ_K. Verified at p = 2, 3, 7, 101 (N = 2000) to 10⁻¹³ and
10⁻⁸ (the second involves λ_1 ≈ 1.6·10⁶).

**(D).** Expanding ‖Π_J(b_{Q,N} − p^{−1} T b_{Q,K})‖². Verified to 4·10⁻¹⁶ at every step and every J tracked.

**The slices.** One identity makes the induction transparent: p^{−1} E_p c_{Q_p,K} = −c·1_{P⁺(n) = p}, where
Q_p is the set of primes below p and P⁺(n) is the largest prime factor (verified to 10⁻¹⁸). So the prime
induction in increasing order adds, at the prime p, exactly the slice of the Möbius source whose largest
prime factor is p; the squared term of (D) is the low-mode population of that slice, and the cross term is the
interference of the slice with everything built before it.

## 2. The transfer matrix (Theorem 2)

**Closed form.** With α^± = p θ_{j,N} ± θ_{ℓ,K} and the Dirichlet kernel D_K(α) = Σ_{m=1}^K cos(mα) =
sin(Kα/2) cos((K+1)α/2) / sin(α/2), D_K(0) = K,

    T_{jℓ} = (2 / √((2N+1)(2K+1))) · [ D_K(α⁻) − D_K(α⁺) ].

Hence |T_{jℓ}| ≤ (2/√((2N+1)(2K+1))) · [min(K, 1/|sin(α⁻/2)|) + min(K, 1/|sin(α⁺/2)|)], uniformly in N, p, j, ℓ.

**Near-diagonal.** Write N = pK + r, 0 ≤ r < p. Then α⁻ at j = ℓ equals (2j−1)π (p − 1 − 2r)/((2N+1)(2K+1)),
so Kα⁻ = O(j/K) and D_K(α⁻) = K(1 + O(j²/K²)); the second kernel is O(K/j). Since 2K/√((2N+1)(2K+1)) =
p^{−1/2}(1 + O(1/K)),

    T_{jj} = p^{−1/2} ( 1 + O(1/j + j²/K² + 1/K) ).

Measured: T_{11} = 0.70728 (p = 2, p^{−1/2} = 0.70711), 0.57706 (p = 3; 0.57735), 0.37759 (p = 7; 0.37796),
0.09795 (p = 101; 0.09950). The floor-dependent error is the factor (p − 1 − 2r) in α⁻.

**Aliases.** D_K(α^±) is of size K exactly when α^± ≡ 0 (mod 2π), i.e. at the rows

    j = ½ + (2N+1)(2m ∓ (2ℓ−1)/(2K+1)) / (2p),        m = 0, 1, …, ⌊p/2⌋,

m = 0 (sign −) being the diagonal. For odd p this gives the diagonal and p − 1 aliases; each carries squared
amplitude ≈ 1/p, so by T^T T = I the mode ℓ of the small horizon is transported with **equal squared
amplitude 1/p to p modes of the large horizon** (up to the rounding of the alias positions, which splits an
alias over two neighbouring rows). Verified at p = 7, ℓ = 1: the seven largest entries of the column sit at
rows 1, 571/573, 1143/1144, 1715/1716 with p·T² = 1.00, 0.55/0.89, 1.18/0.65, 0.27/1.13; predicted rows
572/573, 1143/1144, 1715/1716. At p = 101 the 101 largest entries all have p·T² between 0.3 and 1.2 and sit
on the predicted rows.

**Gain-weighted energy.** By T^T Λ_N T = pΛ_K, Σ_j λ_{j,N} T_{jℓ}² = p λ_{ℓ,K}, and the diagonal alone gives
λ_{ℓ,N}/p = p λ_{ℓ,K} (1 + O(1/K)): **the gain-weighted energy goes entirely to the diagonal to leading order,
while the amplitude is shared equally with the p − 1 aliases**, which sit at j ≥ (2N+1)/p − ℓ where the gains
are O(p²/m²). High modes of the small horizon (ℓ near K, m = 1, sign −) alias into low modes of the large
horizon: this is the channel by which the top of the Parseval mass at horizon K reaches the target modes at
horizon N.

## 3. Mode equals scale (Proposition 3)

    b_{j,N} = (2/√(2N+1)) Σ_{n≤N} μ(n) sin(nθ_{j,N})/n = (2θ_{j,N}/√(2N+1)) Σ_{n≤N} μ(n) sinc(nθ_{j,N}),

so the mode j reads the Möbius sum through the window sinc(nθ_j) at the scale 1/θ_j = (2N+1)/((2j−1)π). For
j = 1 the window is sinc on [0, π/2], monotone from 1 to 2/π. Partial summation (the window has total
variation ≤ θ_1) gives, unconditionally,

    |b_{1,N}| ≤ (4π/(2N+1)^{3/2}) · max_{x≤N} |M(x)|,

and (T′) at J = 1 reads |Σ_{n≤N} μ(n) sinc(πn/(2N+1))| ≤ C_ε N^{1/2+ε}: the sinc-smoothed Möbius criterion,
implied by RH. (T′) at general J is the same criterion through the windows sin(nθ_j)/(nθ_j), which oscillate
over (2j−1)/2 half-periods on [1, N]; mode j at horizon N is mode 1 at the horizon (2N+1)/(2j−1) continued
over the longer range.

**The random benchmark.** For a source with the same moduli μ(n)²/n and independent random signs, E S(J) =
(4/(2N+1)) Σ_n (μ(n)²/n²) Σ_{j≤J} sin²(nθ_j) ≍ J²/N², so (T′) holds for a random-sign source with ε = 0 and
R ≍ log N. Measured (40 draws): S(J)N²/J² between 2.2 and 5.4 for 1 ≤ J ≤ 1000 at N = 2000, R = 6.51 (and
6.25 at N = 4000), K_N = 1.53 (1.51). For the Möbius source S(J)N²/J² is 0.013, 0.24, 0.78, 1.20, 1.85 at
J = 1, 10, 50, 100, 1000 (N = 2000) and 0.049, 0.15, 0.47 at J = 1, 10, 50 (N = 4000): **the Möbius source
populates the low modes 3 to 200 times less than a random-sign source at these horizons**, and its R(N) is
1.49 and 1.52 against 6.5 and 6.3. (T′) says: the Möbius source is at most N^ε worse than random in every low-mode
population, uniformly in J and N.

## 4. Two mechanisms that fail, and why (Propositions 4 and 5)

**Proposition 4 (sign-free bounds give the trivial Möbius bound).** For every source with |c_n| ≤ μ(n)²/n,

    ‖Π_J V_N^T c‖ ≤ Σ_n (μ(n)²/n) ‖Π_J V_N^T e_n‖ ≤ C √(J/N) (1 + log J),

with ‖Π_J V^T e_n‖² = (4/(2N+1)) Σ_{j≤J} sin²(nθ_j) ≤ (4/(2N+1)) min(J, (4π²/3) J³ n²/(2N+1)²) and the sum
over n split at N/J. The bound is attained up to constants by the squarefree indicator c_n = μ(n)²/n, whose
first mode is (2θ_1/√(2N+1))·(6/π²)N(1+o(1)) ≍ N^{−1/2}. So S(J) ≤ C(J/N)(1 + log J)² is sharp for sign-free
arguments, it exceeds the target J²N^{−2} by the factor N/J, and it is exactly the trivial bound M(x) ≤ x at the
scale x = N/J. The Cauchy–Schwarz bound on the cross term of (D) at every prime step is a sign-free argument,
hence cannot do better. Measured: the median of |cross|/(Cauchy–Schwarz) over the primes is 1.00 at J = 1
(one mode: the bound is an equality), 0.24 at J = 10, 0.05 at J = 50 (N = 2000). Every closing argument must
use the signs of μ across primes, and (D) locates where: in the cross terms between the slices.

**Proposition 5 (the increasing-prime induction is not an invariant region).** Let Q = {p ≤ N/2}. Then
c_Q = c + Σ_{N/2<p≤N} e_p/p, since the only n ≤ N that are not Q-smooth are the primes in (N/2, N]. For such
p, pθ_{1,N} ∈ (π/4 − o(1), π/2), so sin(pθ_1) ≥ 0.7 for N ≥ 10, and

    b_{Q,1} = b_1 + (2/√(2N+1)) Σ_{N/2<p≤N} sin(pθ_1)/p  ≥  b_1 + 1.4 (π(N) − π(N/2)) / (N √(2N+1)).

By the prime number theorem π(N) − π(N/2) ~ N/(2 log N), and by Proposition 3 with M(x) = o(x/log x)
(de la Vallée Poussin), |b_1| = o(1/(√N log N)). Hence |b_{Q,1}| ≥ (0.7 + o(1))/(√(2N) log N) and, with
λ_1 ~ (2N+1)²/π²,

    K_Q(N) ≥ λ_1 |b_{Q,1}|² ≥ (c + o(1)) · N / log² N,        c = 0.49·2/π² ≈ 0.099,

unconditionally. Measured: λ_1 S_Q(1) = 11.63 = 0.336·N/log²N at N = 2000 and 17.43 = 0.300·N/log²N at
N = 4000; the maximum of K_Q along the induction is larger still (13.68 after p = 661; 20.69 after p = 1327),
35 and 52 times the final K_N = 0.394. The same holds in decreasing order (Q = {p > N/2} gives b_{Q,1} =
v_1(1) − the same sum). **Any argument that propagates a bound K_Q(N) ≤ C N^ε through the induction over primes
must pass through states where K_Q ≍ N/log²N, a power of N above the target.** The subpower cancellation is a
property of the complete set of primes up to N; every initial segment {p ≤ y}, y ≤ N/2, misses it by the
Möbius sum over the y-rough numbers, which is of order N/log y (sieve heuristic; proved here only for y = N/2).
Panel B of the figure shows the hump: K_Q rises from 0.3 to 13.7 as the primes up to N/2 are added and falls
back to 0.394 only when the last packet, the primes in (N/2, N], arrives one pulse at a time.

## 5. What the computation shows about the cancellation

**Where the target lives.** K_N is attained at J* = 1629 (N = 2000) and 3255 (N = 4000), in the high modes,
where λ_J ≈ ¼ and S(J) ≈ Σ μ²/n² ≈ Z: there K_N ≈ 0.394 is the Parseval mass and bounded for every source
(λ_J S(J) ≤ ½ Z for J ≥ N/2). The content of (T) is entirely in the low modes J ≤ N/2, where λ_J S_N(J) is
small at these N (0.005 to 0.2) because the Möbius populations are far below random. Panel A.

**The induction at low J** (N = 2000; the same pattern at 4000):

| J | base S_∅(J) | Σ_p squared terms | Σ_p cross terms | S_N(J) | J²/N² |
|---|---|---|---|---|---|
| 1 | 6.2·10⁻¹⁰ | 2.43·10⁻⁷ | +2.41·10⁻⁷ | 3.17·10⁻⁹ | 2.5·10⁻⁷ |
| 10 | 8.2·10⁻⁷ | 2.87·10⁻⁵ | +2.35·10⁻⁵ | 6.03·10⁻⁶ | 2.5·10⁻⁵ |
| 50 | 1.0·10⁻⁴ | 7.37·10⁻⁴ | +3.53·10⁻⁴ | 4.87·10⁻⁴ | 6.3·10⁻⁴ |

Three facts. (i) The sum of the squared terms is itself of the target size: 0.97, 1.15, 1.18 times J²/N² at
N = 2000 and 1.06, 1.17, 1.22 at N = 4000. (ii) The total interference is positive at every J and both N: the
slices cancel each other, and the final population is below the incoherent sum (by a factor 77 at J = 1, 4.8
at J = 10, 1.5 at J = 50). (iii) The sign of the cross term is not uniform over primes (negative at 102, 53,
118 of the 303 primes at J = 1, 10, 50), and the dyadic prime packets alternate in sign too; the cancellation
is between packets, not within them.

**The squared terms are the sliced historical integral.** With M_p(x) = Σ_{m≤x, P⁺(m)<p} μ(m),
sq_p(1) = (4θ_1²/(2N+1)) (Σ_m μ_{Q_p}(m) sinc(pmθ_1))² ≈ (4π²/(2N+1)³) M_p(N/p)², and for p > √N every
m ≤ N/p is p-smooth, so M_p(N/p) = M(N/p). Measured: Σ_p sq_p(1) = 2.43·10⁻⁷ against (4π²/(2N+1)³)Σ_p M_p(N/p)²
= 3.49·10⁻⁷ and R(N)/N² = 3.73·10⁻⁷ at N = 2000; 6.62·10⁻⁸, 9.57·10⁻⁸, 9.53·10⁻⁸ at N = 4000. So the
incoherent slice sum is comparable to R(N)/N², the historical integral over the scales N/p: bounding it is (G)
again, not an independent input.

**Dyadic mode bands against dyadic scales.** The band energies Σ_{j∈[2^k,2^{k+1})} λ_j|b_j|² are 0.005,
0.013, 0.028, 0.050, 0.074, 0.139, 0.225, 0.225, 0.093, 0.291, 0.349 for k = 0..10 (N = 2000), against the
historical integral over x ∈ (N/2^{k+1}, N/2^k]: 0.023, 0.020, 0.028, 0.033, 0.058, 0.079, 0.139, 0.238,
0.279, 0.083, 0.500. Comparable band by band within factors of about 3, equal in total: the mode bands are a
scale decomposition of the historical integral, which is what Proposition 3 says.

## 6. Result of the attempt

**Strongest theorems obtained** (all unconditional): Theorem 2 (the transfer matrix in closed form; the
near-diagonal law p^{−1/2}(1 + O(1/j + j²/K² + 1/K)); the alias rows; equal amplitude sharing among p modes
with the gain-weighted energy on the diagonal), Proposition 3 (mode j is the sinc-windowed Möbius sum at scale
(2N+1)/((2j−1)π); the unconditional bound on b_1), Proposition 4 (sign-free arguments give exactly the trivial
Möbius bound, sharply), Proposition 5 (the prime induction passes through K_Q ≍ N/log²N, so no invariant-region
argument over primes can carry (T)). The identities (A)–(D) and the slice identity are verified.

**No closing mechanism was found.** The attempt stops at the cross terms of (D): their total is observed to be
positive and of the size needed, but the only tools available for them without the signs of μ are sign-free and
give the trivial bound (Proposition 4), and the natural induction that would organise the signs is not
monotone (Proposition 5). The strongest unconditional statement about the target remains
K_N ≤ R(N) ≪ N exp(−c (log N)^{3/5} (log log N)^{−1/5}), from the zero-free region through M(x); nothing here
improves it.

## 7. The residual

The induction gives, exactly, for every J,

    S_N(J) = S_∅(J) + Σ_{p≤N} ‖Π_J V_N^T c^{(p)}‖² − Σ_{p≤N} (2/p) ⟨Π_J b_{Q_p,N}, Π_J T_p b_{Q_p,⌊N/p⌋}⟩,

with c^{(p)} = c·1_{P⁺(n)=p} the slice of the Möbius source at largest prime factor p, S_∅(J) ≤ (16π²/3)
J³/(2N+1)³, and the squared sum equal to the incoherent slice population. The residual inequality is

    (R)   for every ε > 0 there is C_ε such that for all N and all 1 ≤ J ≤ N/2:
          Σ_{p≤N} (2/p) ⟨Π_J b_{Q_p,N}, Π_J T_p b_{Q_p,⌊N/p⌋}⟩  ≥  Σ_{p≤N} ‖Π_J V_N^T c^{(p)}‖²  −  C_ε J² N^{−2+ε}.

Source class: the Möbius source only, sliced by largest prime factor, with the inherited cutoffs ⌊N/p⌋ kept.
Quantifiers: uniform in N and J ≤ N/2. Sign requirement: the total interference between the slices must cancel
the incoherent slice population to within the target; it need not be positive termwise, and it is not
(Section 5(iii)). Dependencies: (R) ⟺ (T′) exactly, by the identity above and S_∅ ≤ CJ²/N²; so (R) is (T′)
written with its two halves separated, and it is not weaker. What (R) adds to (T′): it names the two quantities
that must balance, the incoherent slice sum (observed 0.97 to 1.22 times J²/N², and comparable to R(N)/N²) and
the total cross term (observed positive, 0.4 to 1.0 times the slice sum), and it says where the attempt stops:
the incoherent sum is (G) in sliced form, and the cross term is where the signs of μ across the primes act. The
strongest refinement tested is the observation that the balance holds at N = 2000 and 4000 at every J tracked,
with the cross term never negative in total; the strongest refinement proved is Proposition 5, which says the
balance cannot be reached by adding the primes one at a time.

## 8. The completion of the prime construction, and the identity that governs it (2026-10-06, later)

The user, with ChatGPT, read panel B exactly: at the stage Q = {p ≤ N/2} every composite's coefficient is
already present, and the remaining additions are the primes in (N/2, N], each with coefficient −1, so in every
mode b_{j,N} = b_{Q,N}(j) − d_{j,N} with the explicit packet d_{j,N} = Σ_{N/2<p≤N} v_{j,N}(p)/p. Their first-mode
energies λ_1|b_1|², before and after the packet, are reproduced here to the last digit (`rh_completion.py`):

| N | before, λ_1 b_{Q,1}² | after, λ_1 b_1² | the packet, λ_1 d_1² | packet / (N/log²N) |
|---|---|---|---|---|
| 1000 | 6.674594 | 0.001316 | 6.488490 | 0.3096 |
| 2000 | 11.628940 | 0.005139 | 11.145157 | 0.3220 |
| 4000 | 17.432057 | 0.019781 | 18.626284 | 0.3203 |

Their asymptotic for the packet, from the prime number theorem through the sine window,
d_{1,N} ~ (√2/(√N log N)) ∫_{1/2}^1 sin(πx/2)/x dx and λ_1 d_1² ~ (8I²/π²) N/log²N = 0.30338 N/log²N, is
confirmed (I = 0.61178629); the measured ratios 0.31–0.32 approach it slowly, as the packet's own
corrections of relative size 1/log N dictate. This sharpens Proposition 5's constant from the crude 0.099
to the true 0.30338.

**The completion principle at every P ≥ √N.** Every n ≤ N that is not P-smooth has exactly one prime factor
above P, and its cofactor m ≤ N/p < √N ≤ P is P-smooth automatically, so

    c_N = c_{Q_P,N} − Σ_{P<p≤N} p^{−1} E_p c_{⌊N/p⌋},        φ_N = φ_{Q_P,N} − Σ_{P<p≤N} p^{−1/2} T_{log p} φ_{⌊N/p⌋},

with the complete Möbius sources at the small horizons. Verified exactly at P = ⌈√N⌉, N/4, N/2. The populations
of the unfinished state, of the completion and of the survivor (λ_J times the population of the first J modes,
N = 2000):

| P | J = 1: unfinished, completion, survivor, cosine | J = 10 | J = 50 |
|---|---|---|---|
| 45 | 0.381, 0.474, 0.0051, 1.000000 | 0.357, 0.536, 0.027, 0.990 | 0.162, 0.082, 0.081, 0.708 |
| 500 | 10.25, 9.80, 0.0051, 1.000000 | 0.112, 0.098, 0.027, 0.873 | 0.083, 0.004, 0.081, 0.186 |
| 1000 | 11.63, 11.15, 0.0051, 1.000000 | 0.066, 0.040, 0.027, 0.767 | 0.082, 0.002, 0.081, 0.123 |

Two readings. In the first mode the unfinished state and its completion are parallel to six decimals at every P,
because M_{Q_P}(N) = M(N) + Σ_{P<p≤N} M(⌊N/p⌋) (for P ≥ N/2 every inherited cutoff is 1 and the sum is
π(N) − π(P); at N = 4000, P = 1000 the inherited values give 189 = the actual M_{Q_P}(N), while the naive
π(N) − π(P) would give 373) and the second term dominates: **the hump is the pole**, the uniform density of the
primes seen through the first mode, and it cancels identically. In the higher low modes the
completion is negligible (at J = 50 it is 1 to 3 % of the unfinished state) because the window sin(pθ_j)
oscillates over the packet: integrating by parts against the smooth prime density gives the main term
d_{j,N} = O(1/(j √N log N)), so the packet's mode-j energy is O(N/(j⁴ log²N)) and exceeds N^ε only for
j ≲ N^{1/4}; as a bound on the actual prime packet uniformly in j this needs the prime discrepancy π(x) − li(x)
controlled against the window, which the prime number theorem with error term supplies for j ≤ N^{1/2−δ} and not
beyond. The survivor at J = 50 is simply
the unfinished state: there the last packet neither adds nor cancels anything.

**The identity that controls the cancellation.** Two exact identities, valid at every horizon and in every
mode, with the inherited cutoffs ⌊N/d⌋ kept:

    (U)   Σ_{d≤N} d^{−1} E_d c_{⌊N/d⌋} = e_1,

the uniform completion of all horizons is the unit source (this is Σ_{m|n} μ(m) = [n = 1]); and, from
μ·log = −(Λ * μ), i.e. (1/ζ)′ = (ζ′/ζ)(1/ζ),

    (Λ)   b^{log}_{N} := V_N^T (c · log n) = −Σ_{d≤N} (Λ(d)/d) T_d b_{⌊N/d⌋}.

Subtracting (U) from (Λ):

    b^{log}_{N} + V_N^T e_1  =  −Σ_{d≤N} ((Λ(d) − 1)/d) T_d b_{⌊N/d⌋}.

Verified to 2·10⁻¹⁶ at N = 1000, 2000, 4000. This is the answer to the question asked: **the cancellation between
an unfinished source and its exact completion is governed by Λ(d) − 1, the fluctuation of the primes about
their uniform density.** The uniform part of every completion, the pole, cancels identically through (U); what
survives in every mode is the bilinear form of the prime fluctuation against the transported complete
amplitudes at all smaller horizons. Measured (λ_J times populations, N = 4000): at J = 1 the log-weighted
amplitude has 0.782, the unit term 0.0005, the fluctuation term 0.743; at J = 10: 0.650, 0.002, 0.641; at
J = 50: 0.849, 0.009, 0.803. The survivor is the fluctuation term, nothing else.

**What the identity gives and what it cannot.** It is the mode form of the classical identity
M(x) log x = ∫_1^x M(t) dt/t − 1 − Σ_{d≤x} (Λ(d) − 1) M(x/d), through which the prime number theorem transfers to
M(x) = o(x) (Axer's theorem handles the bilinear form). Its quantitative content is exactly the error term of the
prime number theorem, since Λ(d) − 1 is small only on average, ψ(x) − x = −Σ_ρ x^ρ/ρ + O(1). So the identity
locates the survivor in the zeros, uniformly across horizons and modes, and it bounds the survivor only as well
as ψ(x) − x is bounded: x e^{−c√log x} unconditionally, x^{1/2+ε} under RH. No elementary identity does better,
because x^{1/2+ε} for ψ is RH. The productive question is answered in that sense: the control is Λ − 1, and
Λ − 1 is the zeros.

## 9. The coupled fluctuation equation, and what its inversion gives (2026-10-06, later)

The user, with ChatGPT, refined the fluctuation identity by separating its d = 1 term (Λ(1) − 1 = −1, T_1 = I),
so that every horizon on the right is strictly smaller than N:

    (L_N − I) b_N = −u_N − F_N,        L_N = V_N^T diag(log n) V_N,   u_N = V_N^T e_1,   F_N = Σ_{d=2}^N ((Λ(d) − 1)/d) T_d b_{⌊N/d⌋},

and proposed the inversion bound ‖(L_N − I)^{−1}‖ ≤ 216 in the energy norm ‖z‖²_Λ = Σ λ_j z_j², hence
√R(N) ≤ 216 (1 + ‖F_N‖_Λ), and the sufficient estimate ‖F_N‖²_Λ ≤ C_ε N^ε. Everything here is verified and
then evaluated (`rh_fluctuation.py`).

**The equation** holds to 2·10⁻¹⁶ at N = 1000, 2000, 4000, and ‖u_N‖_Λ = 1 exactly (it is B_{11}).

**The inversion is bounded uniformly in N**, and the mechanism is the one stated: in coefficient space the
operator is multiplication by g(n) = 1/(log n − 1), the energy norm of a coefficient vector is the H^{−1}-norm
of its source Σ √n c_n δ_{log n} with respect to H₀, and multiplication by a function G on the log line with
G(log n) = g(n) is bounded on H^{−1} by the H¹-multiplier norm of any bounded Lipschitz extension G. The
sampled values are −1, −3.259, 10.141, 2.589, 1.640, …; the steepest pair is (log 2, log 3), with slope 33.1
and sup 10.14, which gives the multiplier bound √(8·33.1² + 2·10.14²) ≈ 95 by the elementary estimate
‖Gf‖²_E ≤ (8‖G′‖²_∞ + 2‖G‖²_∞)‖f‖²_E. The exact value, the 2-norm of C^T diag(g) C^{−T} (C the cumulative-sum
matrix, B = CC^T), is

    ‖(L_N − I)^{−1}‖_{Λ→Λ} = 25.4536        at N = 1000, 2000 and 4000 alike,

N-independent to four decimals because the extremal vectors live at n = 2, 3 where g jumps. So the chain is
√R(N) ≤ 25.46 (1 + ‖F_N‖_Λ); measured 1.235 ≤ 25.46 (1 + 2.165) = 80.6 at N = 4000. The bound 216 is valid and
loose by a factor eight.

**The coupled expression collapses exactly.** In coefficient space F_N = V_N^T f with

    f_n = (1/n) Σ_{d | n, d ≥ 2} (Λ(d) − 1) μ(n/d) = (1/n) [ (Λ * μ)(n) − (1 * μ)(n) + μ(n) ] = −(log n − 1) μ(n)/n    (n ≥ 2),   f_1 = 0,

by Λ * μ = −μ·log and 1 * μ = δ. Verified to 6·10⁻¹⁷. Therefore

    ‖F_N‖²_Λ = R^{log}(N) := Σ_{a,b=2}^N μ(a)μ(b)(log a − 1)(log b − 1)/max(a, b) = Σ_{k=2}^{N−1} G(k)²/(k(k+1)) + G(N)²/N,

with G(x) = Σ_{2≤n≤x} μ(n)(log n − 1) = M(x)(log x − 1) − ∫_1^x M(t) dt/t + 1. Values:

| N | R^{log}(N) | R(N) | R^{log}/(R log²N) | G(N) | G(N)²/N |
|---|---|---|---|---|---|
| 10³ | 2.500 | 1.459 | 0.036 | 22.5 | 0.51 |
| 10⁴ | 8.472 | 1.582 | 0.063 | −173.2 | 3.00 |
| 10⁵ | 13.890 | 1.622 | 0.065 | −481.3 | 2.32 |
| 10⁶ | 28.072 | 1.708 | 0.086 | 2754.7 | 7.59 |
| 10⁷ | 58.071 | 1.837 | 0.122 | 15727.7 | 24.74 |

So **the signed contributions from the smaller horizons combine exactly, by the convolution identity, into the
log-weighted Möbius vector −(log n − 1)μ(n)/n.** There is no cancellation left in F_N to exploit: all of it has
been used in forming μ·log from Λ * μ. The equation (L_N − I)b_N = −u_N − F_N is, once F_N is evaluated, the
identity (L_N − I)b_N = −u_N + (L_N − I)b_N + u_N, and the inversion bound reads √R(N) ≤ 25.46 (1 + √R^{log}(N)):
true, and circular, since R^{log} is the same energy with logarithmic weights.

**The sufficient estimate is the target with logarithmic weights, equivalent to it.** R^{log}(N) = O_ε(N^ε) ⟺ RH:
(⟸) under RH, M(x) = O(x^{1/2+ε}) gives G(x) = O(x^{1/2+ε} log x) and R^{log}(N) = O(N^{3ε}). (⟹) G(N)² ≤ N·R^{log}(N)
gives G(x) = O(x^{1/2+ε}); by partial summation the Dirichlet series Σ_{n≥2} μ(n)(log n − 1) n^{−s} =
ζ′(s)/ζ(s)² − 1/ζ(s) + 1 then converges, hence is analytic, on Re s > ½; at a zero ρ of ζ of order m with Re ρ > ½
the first term has a pole of order m + 1 and the second a pole of order m, so the sum has a pole; contradiction.
Hence R^{log} subpower ⟺ RH ⟺ (G) ⟺ (T). Also two-sidedly, without any hypothesis: R(N) ≤ 25.46² (1 + √R^{log}(N))²
and R^{log}(N) ≤ (log N − 1)²·R(N)·(1 + o(1)) + …, so √R and √R^{log} are equivalent up to the factor log N and the
constant 25.46.

**Where the identity is not circular.** Only when the horizons are kept separate and an independent input on
Λ(d) − 1 is fed in. Sign-free in d, (C) gives ‖T_d b_{⌊N/d⌋}‖_Λ = √(d R(⌊N/d⌋)) and the bound
√R(N) ≤ 25.46 (1 + Σ_{d≥2} |Λ(d) − 1| d^{−1/2} √R(⌊N/d⌋)) ≤ C √N max R, the trivial bound again. With the signs
of Λ(d) − 1 but the vectors T_d b_{⌊N/d⌋} treated as given, the available input is ψ(x) − x, and the transfer
(Axer's theorem) yields M(x) = o(x) from ψ(x) ~ x and the zero-free-region bound from the zero-free region; it
yields x^{1/2+ε} only from ψ(x) − x = O(x^{1/2+ε}), which is RH. Whether the coupled form has cancellation beyond
this is answered by the collapse: as a vector, F_N is μ·(1 − log n)/n, and any cancellation it carries is the
cancellation of μ.
