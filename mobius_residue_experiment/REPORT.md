# Results report: Michaels prime–zero shape experiment

Run `results/run`, 8 October 2026. Software: qiskit 2.5.2, qiskit-aer 0.17.2, qiskit-ibm-runtime 0.50.0, mpmath 1.4.1,
python-flint 0.9.0, numpy 2.4.6, Python 3.11.15. Wall times: verify 0.8 s, arithmetic-shape 107 s (65,531 cutoffs,
seven complete block decompositions), prime-zero-shape 114 s, simulate 88 s (80 cases + extras), prepare-hardware 0.4 s.

## 1. Deterministic gates (verify.json): 47 of 47 PASS

Largest discrepancies: fixture state-preparation errors ≤ 5.4·10⁻¹⁵; statevector P₀ against the exact rational
≤ 10⁻¹²; Weyl matrix and commutation errors ≤ 3.8·10⁻¹⁵; Hadamard-test expectations against cos θ, sin θ ≤ 10⁻¹⁰;
feature Gram identity, residue coarsening and the exact fixtures hold with zero error in rational arithmetic;
arithmetic identities δ = κ + V hold exactly for N = 6 … 256 and the selected-family prefix formula equals its
defining pair sum to 10⁻¹²; block decompositions at 64 and 128 reproduce V − A to 10⁻¹⁴; the prime–zero return
identity holds to 10⁻⁸ at X = 4 with 10 zeros and the mode Gram energy equals its quadrature to 10⁻⁸.

| Fixture | R | Q | P₀ | expected | statevector P₀ |
|---|---|---|---|---|---|
| A Möbius N=4, d=4 | 5/6 | 11/6 | 5/44 | 5/44 | 0.113636363636 |
| B squarefree-positive N=4, d=4 | 25/6 | 11/6 | 25/44 | 25/44 | 0.568181818182 |
| C packet N=6, d=2 | 2/3 | 1 | 1/3 | 1/3 | 0.333333333333 |
| C packet N=6, d=4 | 2/3 | 5/3 | 1/10 | 1/10 | 0.100000000000 |
| C packet N=6, d=8 | 2/3 | 2 | 1/24 | 1/24 | 0.041666666667 |
| D all-positive N=8, d=8 | 3719/280 = 2N − H_N | 761/280 = H_N | 3719/6088 | (2N−H_N)/(dH_N) | — |

## 2. Section 9A: the signed arithmetic remainder, 6 ≤ N ≤ 65536

Definitions as in the prompt (L = ⌊N/6⌋, K_red = ⌊N^{1/6}⌋, V_N = Σ_{k=L}^{N−1}(h(N)−h(k))², κ_N = 2M(L)H + LH²,
δ_N = R(N) − R(L) = κ_N + V_N, A_N the selected family, C_N = V_N − A_N, Γ_N = κ_N + C_N). Numerical precision of the
sweep: float64 with the direct positive sum for V (the prefix formula is kept as a cross-check) and prefix arrays
F_{q,j} for A; measured identity errors max |δ − (κ+V)| = 9.6·10⁻¹⁵, max |V_direct − V_prefix| = 4.3·10⁻¹⁴. Exact
rationals at N = 6, 7, 8, 12, 36, 64, 128, 256, 512 (`arith_exact_small.csv`); 200-bit ball enclosures at 14 cutoffs
including every sixth-power threshold ±1 (`arith_summary.json`, "balls"), with the float values inside the balls to
2·10⁻¹¹ at worst.

**Finding 1: the complement and Γ are negative at every cutoff above 63.** C_N > 0 only at N = 6, 7 and a few cutoffs
below 64 (10 % of 6 ≤ N ≤ 63); Γ_N > 0 only at N = 12, 13. The global maxima are C₇ = 0.4683 (exact 59/126) and
Γ₁₃ = 0.5555; the global maximum of δ_N is 1.2192 at N = 13. For every N > 100: max C_N = −0.178 (N = 221), max
Γ_N = −0.178 (N = 110). The records of [C]₊/N^η and [Γ]₊/N^η for η = 0.05, 0.1, 0.2 are therefore all set at N ≤ 13 and
never move. Declared range [6, 16384] against held-out (16384, 65536]: max C 0.468 (at 7) against −0.306; max Γ 0.556
(at 13) against −0.314. No growth trend exists to survive or fail: above 63 the positive parts are identically zero.

**Finding 2: the complement is a small negative difference of two large spectral energies.** With the masked
dyadic blocks of Section 13 (ordered block pairs, all squarefree g ≤ N/(K_red+1), SVD split S = E⁺ − E⁻):

| N | blocks | E⁺ total | E⁻ total | E⁺ − E⁻ = C_N | κ_N | Γ_N | status |
|---|---|---|---|---|---|---|---|
| 64 | 65 | 1.6749 | 2.0241 | −0.349238 | +0.2045 | −0.1448 | complete |
| 128 | 139 | 2.1465 | 2.5081 | −0.361595 | −0.0100 | −0.3716 | complete |
| 256 | 261 | 2.4188 | 2.7771 | −0.358347 | −0.0089 | −0.3672 | complete |
| 512 | 555 | 3.3309 | 3.6796 | −0.348759 | −0.0462 | −0.3949 | complete |
| 1024 | 801 | 3.9971 | 4.4186 | −0.421446 | +0.0111 | -0.4103 | complete |
| 2048 | 1570 | 4.5860 | 4.9848 | −0.398789 | -0.0259 | -0.4247 | complete |
| 4096 | 3107 | 5.6058 | 5.9748 | −0.369032 | −0.1430 | −0.5120 | complete |

Reconstruction errors of the aggregate against V − A are ≤ 2.1·10⁻¹²; the largest positive single blocks are the
g = 1 top dyadic blocks (U = V = N/2 scale: 0.12 at 64, 0.19 at 4096). Both E⁺ and E⁻ grow with N while their
difference stays between −0.35 and −0.42. Every decomposition listed is complete; no partial decomposition was needed
because no positive excursion exists above N = 63.

**Finding 3: the selected family changes by a jump at each sixth-power threshold, the complement by the opposite
jump.** A_N: 0.607 → 0.416 at 64, 0.410 → 0.480 at 729, no change at 4096 (K_red → 4 adds no new coprime pair below
the squarefree condition), 0.478 → 0.599 at 15625, 0.599 → 0.531 at 46656; C_N jumps by the negatives; δ_N and c_N
are continuous there. Band means (K_red = 2, 3, 4, 5, 6): A 0.410, 0.478, 0.478, 0.599, 0.531; V 0.081, 0.078, 0.073,
0.073, 0.093; κ −0.028, −0.027, −0.018, −0.024, −0.033; δ +0.053, +0.051, +0.056, +0.050, +0.060.

**Finding 4: the original normalized allowance.** c_N = [δ_N]₊/B(K_red)⁵ has its maximum 0.038101 at N = 13 (exact
18307/480480), 0.00847 on 64 ≤ N < 729, 0.00647 on the next band, 0.00629, 0.00267, 0.00254 on the last three
bands; the signed diagnostic Γ_N/B(K_red)⁵ is negative at every N above 13 (maximum above 63: −0.0037).

## 3. Section 9B: prime discrepancy against the listed zeros

Zeros: the first 200 critical-line zeros from `mpmath.zetazero(n, info=True)` at 30 digits (maximum ordinate
396.3819); provenance in `zeros_mpmath.json` (Gram brackets, Rosser-block patterns); completeness of the census is
not certified by this computation. Windows X = 2, 4, …, 4096, ratio a = 6, nested lists of 10, 25, 50, 100, 200 zeros
plus conjugates. P is computed exactly in the staircase weights; Z, I, D by Gauss–Legendre on every jump interval,
subdivided in log x so that γ_max·Δ(log x) ≤ 0.5 with 12 nodes; doubling the nodes and halving the pieces changes Z and
D by ≤ 3·10⁻¹⁵, and the quadrature of (s−x)² reproduces the exact P to 10⁻¹⁴. The return identity P = Z + I + D holds
to 1.6·10⁻¹⁴ on all 60 rows and both norm bounds hold.

| X | zeros | P | Z | I | D | D/P | E_Z | off-diag | λ_Z | utilization |
|---|---|---|---|---|---|---|---|---|---|---|
| 2 | 200 | 8.9012 | 8.8781 | −0.0018 | 0.0250 | 0.003 | 0.614 | −0.122 | 23.3 | 0.17 |
| 16 | 200 | 1.7865 | 1.7360 | +0.0009 | 0.0496 | 0.028 | 0.680 | −0.056 | 23.3 | 0.18 |
| 128 | 200 | 0.8110 | 0.7305 | −0.0019 | 0.0824 | 0.102 | 0.639 | −0.097 | 23.3 | 0.17 |
| 1024 | 200 | 0.7484 | 0.6708 | +0.0006 | 0.0770 | 0.103 | 0.647 | −0.089 | 23.3 | 0.18 |
| 4096 | 200 | 0.9741 | 0.8970 | +0.0018 | 0.0753 | 0.077 | 0.901 | +0.164 | 23.3 | 0.24 |
| 4096 | 10 | 0.9741 | 0.5745 | +0.0276 | 0.3720 | 0.382 | 0.576 | +0.102 | 13.9 | 0.48 |
| 4096 | 50 | 0.9741 | 0.8210 | +0.0083 | 0.1448 | 0.149 | 0.825 | +0.176 | 19.6 | 0.32 |

The retained list reconstructs 90–97 % of the prime-discrepancy energy for X ≥ 32 (D/P between 0.05 and 0.10 at 200
zeros) and 99.7 % at X = 2; D falls roughly like the inverse of the zero count at fixed X (0.372, 0.221, 0.145, 0.107,
0.075 at X = 4096 for 10 … 200 zeros). The signed cross term I is small, |I| ≤ 0.035, and changes sign with the zero
count (at X = 4096: +0.028, +0.035, +0.008, −0.001, +0.002 for 10, 25, 50, 100, 200 zeros): the residual is nearly
orthogonal to the reconstruction in the window norm. The finite mode bound E_Z ≤ λ_Z W_Z holds on every row with
utilization 0.17–0.57; λ_Z grows slowly with the list (13.9, 17.3, 19.6, 21.6, 23.3) and the signed off-diagonal mode
energy is negative for X ≤ 1024 and positive at X = 4096 with 200 zeros. λ_Z is a finite-mode constant; nothing here
bounds the omitted spectrum.

## 4. Quantum calibration (Aer, ideal statevector sampling, 8192 shots, seed 20261008)

80 cases (4 sources × N ∈ {4, 8, 16, 32, 64} × r ∈ {2, 3, 4, 5}). Max |P_statevector − P_classical| = 1.7·10⁻¹²;
max normalization error 4.4·10⁻¹⁶; max state-preparation error 1.1·10⁻¹¹. Wilson 95 % intervals cover the exact P₀ in
75 of 80 cases (6.25 % misses against the nominal 5 %; the five misses are listed in `cases.csv`, `covered = False`,
and are sampling variation, not implementation failures: the statevector agrees with the exact value in every case).

Möbius main results (exact P₀ = R/(dQ); sampled estimate; reconstructed R̂ with its interval):

| N | r | R | Q | R−Q | P₀ | P̂₀ [Wilson] | R̂ [interval] | prep depth / CX | data qubits |
|---|---|---|---|---|---|---|---|---|---|
| 4 | 2 | 0.8333 | 1.8333 | −1.000 | 0.113636 | 0.1133 [0.107, 0.120] | 0.831 [0.78, 0.88] | 42 / 11 | 4 |
| 16 | 2 | 1.1621 | 2.4077 | −1.246 | 0.120667 | 0.1260 [0.119, 0.133] | 1.213 [1.15, 1.28] | 327 / 57 | 6 |
| 16 | 4 | 1.1621 | 2.7488 | −1.587 | 0.026423 | 0.0266 [0.023, 0.030] | 1.170 [1.03, 1.33] | 1145 / 247 | 8 |
| 32 | 3 | 1.7382 | 3.2982 | −1.560 | 0.065878 | 0.0658 [0.061, 0.071] | 1.736 [1.60, 1.88] | 1282 / 247 | 8 |
| 64 | 2 | 1.3329 | 4.4578 | −3.125 | 0.074750 | 0.0717 [0.066, 0.077] | 1.278 [1.18, 1.38] | 1060 / 247 | 8 |
| 64 | 5 | 1.3329 | 3.4795 | −2.147 | 0.011971 | 0.0121 [0.010, 0.015] | 1.346 [1.11, 1.64] | 12060 / 2036 | 11 |

Controls at N = 64, r = 2: R/D = 0.375 (Möbius), 0.511 (random signs on the squarefree support), 14.3
(squarefree-positive), 26.0 (all-positive). At d ≥ N the saturation Q = D holds exactly and P₀ ∝ 1/d while R̂ stays
constant (e.g. N = 4: P₀ = 5/44, 5/88, 5/176, 5/352 for r = 2…5): this decline is the known normalization effect.
Phase-modified readout (V₁ then H, N = 16, r = 2): formula P₁ = 0.243610, statevector 0.243610, sampled 0.2433
[0.234, 0.253]. Negative control (raw residue masses β_b): Möbius N = 16, d = 4 gives β = (0, −1, 2, −2) and
probability 1/36 = 0.0278, not R/(dQ) = 0.1207; all-positive gives 1.
Character spectra: exact P_t match R_{a,t}/(dQ) to 10⁻¹⁴ in all eight cases, sum to one, and P_t = P_{d−t}; sampled
Möbius N = 16, r = 2: (0.1227, 0.2444, 0.3861, 0.2468) against (0.1207, 0.2436, 0.3921, 0.2436).

## 5. Ancilla Weyl tests

W = V_t T_h V_t† T_h† = e^{iθ}I, θ = 2πth/d, verified as matrices for r = 2…5 (t = h = 1) and with zero controls.
Hadamard test (t = h = 1): r = 2, θ = π/2: ⟨X⟩ exact 0, sampled −0.0095 ± 0.022; ⟨Y⟩ exact 1, sampled 1.000 ± 0 (the
post-rotation state is a Y eigenstate, so the zero variance is genuine); r = 3, θ = π/4: ⟨X⟩ 0.7071 / 0.7063 ± 0.016,
⟨Y⟩ 0.7071 / 0.7034 ± 0.016; r = 4 exact 0.9239 / 0.3827. Joint-state case (N = 16 Möbius state plus ancilla, Y)
exact 1. The unconditional orders V T and T V give identical readout (0.243610135266 both) and a shift alone leaves
P₀ unchanged, as the global-phase argument requires. Transpiler report: at optimization level 3 the controlled loop
keeps 12 CX (r = 2) and 116 CX (r = 3); it is not collapsed to a hard-coded ancilla phase.

## 6. Resources

State preparation by `StatePreparation` (isometry synthesis), basis rz/sx/x/cx, optimization level 1, seed 1.
Depth and CX count double with each added qubit: 4 qubits 42 / 11; 6 qubits 327 / 57; 8 qubits 1060–1282 / 247;
10 qubits 5107–6080 / 1013; 11 qubits 10696–12218 / 2036 (the largest implemented circuit: Möbius/random N = 64,
r = 5: 11 data qubits, prep depth 12060, full depth 12062, 10067 one-qubit and 2036 two-qubit gates, transpile 0.6 s,
sampling under a second). The synthesis size is ≈ 2^{r+s}, i.e. linear in N·d: general amplitude loading is exponential
in the qubit count and no efficient quantum algorithm is claimed. Classical reference cost: O(N) for the amplitude
array and O(N) for R by the cumulative formula. Modeled noise (NoiseModel from FakeManilaV2, Möbius N = 4, r = 2,
depth 50 after routing): P₀ = 0.166 [0.158, 0.174] against the exact 0.1136, a bias toward the uniform outcome; this is
a simulator noise model, not device data.

## 7. Hardware stage: pending

No saved IBM Quantum account exists in this environment, so nothing was submitted and no job IDs exist. The candidate
batch was transpiled against the fake Brisbane target (`results/run/hardware/isa_batch.qpy`, `batch.json`): Möbius
N = 4, r = 2 (4 qubits, ISA depth 97, 19 two-qubit gates, P₀ = 5/44); squarefree-positive N = 4, r = 2 (depth 54, 11);
packet N = 6, r = 1 (depth 96, 19, P₀ = 1/3); Weyl loop r = 2 X and Y (3 qubits, depth 86, 14 each). Recommended
first batch: Möbius N = 4 r = 2, packet N = 6 r = 1 and the two Weyl circuits at 4000 shots (16,000 shots total).
Ready command: `python mobius_residue_experiment.py --results results/run submit-hardware --backend <name> --shots 4000 --confirm`,
which runs only with `--confirm` and a saved account and records job IDs, counts and metadata separately from any mitigation.

## 8. Interpretation

* **Forced by normalization:** the 1/d decline of P₀ at fixed N once d ≥ N (Q = D there); the proportionality
  R̂ = dQP̂₀; the growth of Q with N for the positive sources.
* **Differences between sign patterns:** R/D at N = 64 is 0.37 for Möbius and 0.51 for random signs against 14 and 26
  for the positive sources; the cross-residue contribution R − Q is negative for Möbius at every N and r (−1.0 to −3.1)
  and the character spectrum concentrates on t = d/2 for Möbius (0.39 at N = 16, r = 2; 0.27 at N = 32, r = 3).
* **Shot noise:** every sampled deviation is within or at the edge of its Wilson interval; five of eighty intervals miss,
  consistent with 95 % coverage; the modeled-noise case is the only one with a bias.
* **Mixed terms and terminal:** the feature vectors realize the full kernel 1/max(m,n) exactly (Gram identity exact),
  the joint array keeps every cross-residue term in Q and the terminal A_b(N)²/N in each residue channel; R = dQP₀
  holds for every case to 10⁻¹² in the statevector.
* **Practical benefit:** none after preparation and readout: the state is loaded classically at exponential cost and the
  readout estimates one number with shot-noise error.
* **Where the positive excursions of C and Γ occur:** only at N ≤ 13 (C: 6, 7; Γ: 12, 13); above 63 both are negative
  everywhere up to 65536; the complement is a difference of two growing block energies (E⁺ ≈ 5.6, E⁻ ≈ 6.0 at 4096);
  the jumps at the sixth-power thresholds come from the selected family A_N and cancel in δ_N.
* **How much prime-discrepancy energy the listed zeros reconstruct:** 90–97 % for X ≥ 32 and 99.7 % at X = 2 with 200
  zeros, the residual falling like the inverse zero count; the cross term is at most 4 % of P and of either sign.
* **Not measured:** any bound on the omitted zeros; any uniform bound on C_N beyond 65536; any identification of the
  residue precision r with a zero height or of C_N with D; the quantum success probability is not a prime–zero
  inequality.

## 9. Deviations from the specification

1. Hardware not run (no account); fake Brisbane used only as a transpilation target. 2. Block decompositions were
complete at all computed cutoffs (64 … 4096 and the small record cutoffs 6, 7, 12, 13); the partial-decomposition
branch was never exercised because no positive excursion exists above 63. 3. Character distributions were computed
exactly for eight cases and sampled for one. 4. One modeled-noise case only. 5. The sweep is float64 with reported
identity errors; exactness is established on the listed small cutoffs and by 200-bit balls at 14 cutoffs. 6. The
simulator seed per case is seed + 7N + r (recorded). 7. The arithmetic sweep used K_red thresholds 64, 729, 4096, 15625,
46656 as specified; A_N does not change at 4096 because K_red = 4 adds no new squarefree cofactor.

## 10. Status statement

Verified identities (exact or to the stated tolerance): the feature Gram identity, the residue lift and its two Q
formulas, R = dQP₀, the coarsening relation, the Weyl relations and the controlled-loop phase, the QFT character
formula, δ = κ + V, the prefix form of A, the block SVD split and its aggregate, the prime–zero return identity and the
finite mode bound. Finite empirical observations: C_N and Γ_N negative for all 63 < N ≤ 65536 with the listed maxima;
the band means; the reconstruction residual D and its decay with the zero count; the sampled probabilities. Unproved
arithmetic bounds, untouched by this computation: the all-scale aggregate hypothesis C_N ≤ D_η N^η for every N, the
normalized allowance for every N, and any statement about zeros not in the list.
