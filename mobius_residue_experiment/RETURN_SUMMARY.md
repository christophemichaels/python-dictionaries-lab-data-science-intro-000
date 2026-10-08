# Return summary for ChatGPT review (run results/run, 8 October 2026)

* **Deterministic gates: 47/47 PASS.** Largest discrepancies: state preparation 5.4e-15; statevector P0 vs exact
  rational 1.7e-12; Weyl matrices/commutation 3.8e-15; Hadamard-test expectations vs cos/sin 1e-10; exact fixtures,
  feature Gram and coarsening: zero error in rational arithmetic; block aggregates vs V−A 2.1e-12; prime–zero identity 1.6e-14.
* **Arithmetic cutoff coverage:** every integer 6 ≤ N ≤ 65536 (float64, identity errors ≤ 4.3e-14; exact rationals at
  N ≤ 512 listed cutoffs; 200-bit balls at 14 cutoffs incl. all sixth-power thresholds ±1). Largest positive records:
  C₇ = 59/126 = 0.46825 (N = 7), Γ₁₃ = 0.55552 (N = 13), δ₁₃ = 1.21925. **C_N < 0 and Γ_N < 0 for every N > 63** up to
  65536 (max above 100: C = −0.178 at 221, Γ = −0.178 at 110). Observed envelopes: c_N max 0.038101 (= 18307/480480 at
  N = 13), ≤ 0.00847 above 63, ≤ 0.00254 on 46656–65536; Γ/B(K_red)⁵ ≤ −0.0037 above 63. Held-out (16384, 65536]: max C −0.306, max Γ −0.314.
* **Spectral block decompositions:** COMPLETE at N = 64, 128, 256, 512, 1024, 2048, 4096 (and 6, 7, 12, 13); E⁺/E⁻ totals
  1.67/2.02 (64) … 5.61/5.97 (4096), differences −0.35 … −0.42; no partial decomposition was needed.
* **Prime–zero reconstruction:** 200 critical-line zeros from mpmath.zetazero (30 digits, max height 396.38; Gram/Rosser
  provenance saved; completeness not certified here). D/P with 200 zeros: 0.003 (X=2), 0.028 (16), 0.10 (128), 0.10
  (1024), 0.077 (4096); D at X = 4096 falls 0.372 → 0.075 from 10 to 200 zeros; |I| ≤ 0.035 with sign changes;
  P = Z + I + D to 1.6e-14; finite Gram bound holds everywhere, utilization 0.17–0.57, λ_Z = 13.9 … 23.3.
* **Fixtures (exact = measured):** A 5/44 → sampled 0.1133 [0.107, 0.120]; B 25/44; C 1/3, 1/10, 1/24; D 2N−H_N.
* **Main Möbius / controls:** 80 cases, Wilson 95 % coverage 75/80; e.g. N=64 r=5 (11 qubits) P0 = 0.011971, sampled
  0.0121 [0.0099, 0.0147], R̂ = 1.346 [1.11, 1.64] vs R = 1.3329. R/D at N=64: Möbius 0.37, random signs 0.51,
  squarefree-positive 14.3, all-positive 26.0. Saturation d ≥ N: Q = D, P0 ∝ 1/d, R̂ constant (normalization effect).
* **Ancilla Weyl phase:** θ = 2πth/d; r=2 (θ=π/2): ⟨X⟩ 0 / −0.010±0.022, ⟨Y⟩ 1 / 1.000; r=3 (θ=π/4): 0.7071 / 0.706±0.016
  and 0.7071 / 0.703±0.016. Sign convention: ancilla |+⟩, loop controlled on 1 in order T†, V†, T, V; X: H then measure; Y: S†, H, measure.
* **Largest circuit:** N=64, r=5: 11 data qubits (+1 ancilla in the Weyl variant), prep depth 12060, 2036 CX, 10067 one-qubit
  gates, transpile 0.6 s; synthesis size ≈ 2^(r+s).
* **Hardware:** NOT RUN (no saved IBM account); batch of 5 circuits transpiled against fake_brisbane (ISA depth 54–97,
  11–19 two-qubit); no job IDs exist; ready command in results/run/hardware/RUN_COMMAND.txt.
* **Deviations:** hardware pending; all block decompositions complete (no partial needed); characters exact for 8 cases,
  sampled for 1; one modeled-noise case; sweep in float64 with reported errors.
* **Status:** verified identities (Gram, lift, readout law, Weyl, QFT, δ = κ + V, prefix A, block split, prime–zero return,
  finite Gram bound) are exact; the finite observations are: C_N, Γ_N negative above 63 to 65536, D/P ≈ 0.03–0.10 with 200
  zeros; the unproved bounds (all-scale C_N ≤ D_η N^η, the normalized allowance for all N, the omitted spectrum) are untouched.
