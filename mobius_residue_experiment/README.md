# Michaels prime–zero shape experiment

Signed arithmetic remainder (Section 9A of the master prompt), direct prime–zero comparison from computed zeta zeros
(Section 9B), and an exact Qiskit realization of the original Green kernel 1/max(m, n) with a residue register, the
coherent energy readout, Weyl operators and an ancilla Hadamard test (Sections 3–8, 10–11). Everything here is
reproducible from the command line; every stage writes CSV/JSON under a results directory and `plot` draws all figures
from those files.

## Layout

```
mobius_residue_experiment.py   command-line entry point (verify, arithmetic-shape, prime-zero-shape, simulate, plot,
                               prepare-hardware, submit-hardware, all)
mre/arith.py                   Section 9A: R, V, A, C, kappa, delta, Gamma, c_N; exact rationals, float sweep, balls, blocks
mre/primezero.py               Section 9B: psi staircase, zero reconstruction, P/Z/I/D, finite mode Gram bound, zero provenance
mre/quantum.py                 Sections 3–8: feature register, residue lift, readout, Weyl/Hadamard tests, QFT characters, resources
mre/plots.py                   all figures
results/run/                   the delivered run: config_*.json, verify.json, fixtures.csv, arithmetic/ (arith_sweep.csv.gz, one row per cutoff), primezero/, quantum/, hardware/, plots/
logs/                          stage logs of the delivered run
REPORT.md                      results report;  RETURN_SUMMARY.md: the concise return summary
```

## Tested environment

Python 3.11.15, qiskit 2.5.2, qiskit-aer 0.17.2, qiskit-ibm-runtime 0.50.0, numpy 2.4.6, scipy 1.17.1, mpmath 1.4.1,
python-flint 0.9.0 (ball arithmetic), matplotlib 3.11.2, Linux x86_64. Installation into a fresh virtual environment:

```
python3 -m venv venv && . venv/bin/activate
pip install "qiskit>=2.5,<3" "qiskit-aer>=0.17" "qiskit-ibm-runtime>=0.50" numpy scipy mpmath python-flint matplotlib
```

(The system pip on the delivery machine refused to upgrade a distribution-owned PyJWT, which is why a venv was used;
no account environment was modified.)

## Commands

```
python mobius_residue_experiment.py --results results/run verify
python mobius_residue_experiment.py --results results/run arithmetic-shape --max-N 65536
python mobius_residue_experiment.py --results results/run prime-zero-shape --max-X 4096 --positive-zeros 200
python mobius_residue_experiment.py --results results/run simulate --shots 8192 --seed 20261008
python mobius_residue_experiment.py --results results/run plot
python mobius_residue_experiment.py --results results/run prepare-hardware --backend fake_brisbane      # or a real backend name
python mobius_residue_experiment.py --results results/run submit-hardware --backend <name> --shots 4000 --confirm
python mobius_residue_experiment.py --results results/run all
```

`submit-hardware` does nothing without `--confirm` and a saved IBM Quantum account; `prepare-hardware` never submits.
The delivered run found no saved account and transpiled the candidate batch against the fake Brisbane target; the
hardware stage is pending.

## Conventions

* Residue qubits are 0..r−1 (low-order bits), feature qubits r..r+s−1; flat amplitude index (k−1)·d + b. Verified against
  the installed Qiskit: `StatePreparation` places amplitude index Σ bit_j 2^j with qubit 0 least significant, and a counts
  string lists qubit 0 as its rightmost character. A measured residue string `b_{r−1}…b_0` therefore reads the classical
  residue b = Σ b_j 2^j; the success outcome is `0…0`.
* The feature vectors g_n are not normalized (squared norm 1/n); the joint array is normalized once by √Q.
* Energy readout: H on every residue qubit, measure the residue register, success = all-zero; R̂ = d·Q·P̂₀.
* Weyl loop W = V_t T_h V_t† T_h† is applied in the chronological order T_h†, V_t†, T_h, V_t (rightmost factor first),
  controlled on an ancilla prepared in |+⟩; X readout = H then measure, Y readout = S† then H then measure;
  ⟨X⟩ = cos θ, ⟨Y⟩ = sin θ, θ = 2π t h/d.
* Cyclic characters: `QFTGate` with QFT|b⟩ = d^{−1/2} Σ_t e^{2πi bt/d}|t⟩ (standard output order); P_t = R_{a,t}/(dQ).
* Sampled probabilities carry Wilson 95 % intervals; R̂ intervals are the P intervals times the known d·Q.
* Zeros: `mpmath.zetazero(n, info=True)` at 30 digits; the provenance record (Gram brackets, Rosser-block patterns) is in
  `results/run/primezero/zeros_mpmath.json`. Completeness of the zero census below the listed height is not certified by this
  computation; results are labelled "reconstruction from the listed zeros".
