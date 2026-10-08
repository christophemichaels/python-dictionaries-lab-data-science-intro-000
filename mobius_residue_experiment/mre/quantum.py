"""Sections 3-8, 10-11: the exact Green-feature register, residue lift, coherent energy readout, Weyl operators and the
ancilla Hadamard test, character-resolved readout, fixtures, resources, hardware preparation.  Qiskit >= 2.x.
Conventions (verified against the installed Qiskit): qubit 0 is the least significant bit of the statevector index and the
rightmost character of a counts string; residue qubits are 0..r-1, feature qubits r..r+s-1; flat index = (k-1)*d + b."""
import math, json, time, hashlib, io
from fractions import Fraction as Fr
import numpy as np
from qiskit import QuantumCircuit, ClassicalRegister, transpile, qpy
from qiskit.circuit.library import StatePreparation, QFTGate, UnitaryGate
from qiskit.quantum_info import Statevector, Operator
from qiskit_aer import AerSimulator

# ------------------------------------------------------------------ sources
def make_sources(mu, Nmax=64, seed=20261008):
    rng = np.random.default_rng(seed); signs = rng.choice([-1, 1], size=Nmax + 1)
    srcs = {
        "mobius": [0] + [int(mu[n]) for n in range(1, Nmax + 1)],
        "squarefree_positive": [0] + [int(mu[n] != 0) for n in range(1, Nmax + 1)],
        "random_sign": [0] + [int(signs[n]) if mu[n] != 0 else 0 for n in range(1, Nmax + 1)],
        "all_positive": [0] + [1]*Nmax,
    }
    return srcs, seed

def source_hash(a): return hashlib.sha256(",".join(str(int(v)) for v in a[1:]).encode()).hexdigest()[:16]

# ------------------------------------------------------------------ exact quantities
def R_pairwise(a):
    N = len(a) - 1; return sum(Fr(a[m]*a[n], max(m, n)) for m in range(1, N + 1) for n in range(1, N + 1))
def R_cumulative(a):
    N = len(a) - 1; A = 0; tot = Fr(0)
    for k in range(1, N + 1):
        A += a[k]; tot += Fr(A*A, k*(k + 1)) if k < N else Fr(A*A, N)
    return tot
def D_diag(a):
    return sum(Fr(a[n]*a[n], n) for n in range(1, len(a)))
def residue_cumulatives(a, d):
    """A_b(k) for k = 0..N, b = 0..d-1 (ints)."""
    N = len(a) - 1; A = np.zeros((N + 1, d), dtype=np.int64)
    for k in range(1, N + 1):
        A[k] = A[k-1]; A[k, k % d] += a[k]
    return A
def Q_pairwise(a, d):
    N = len(a) - 1; return sum(Fr(a[m]*a[n], max(m, n)) for m in range(1, N + 1) for n in range(1, N + 1) if (m - n) % d == 0)
def Q_feature_array(a, d):
    N = len(a) - 1; A = residue_cumulatives(a, d); tot = Fr(0)
    for b in range(d):
        for k in range(1, N + 1):
            tot += Fr(int(A[k, b])**2, k*(k + 1)) if k < N else Fr(int(A[N, b])**2, N)
    return tot
def exact_quantities(a, d):
    R1, R2 = R_pairwise(a), R_cumulative(a); Q1, Q2 = Q_pairwise(a, d), Q_feature_array(a, d)
    assert R1 == R2 and Q1 == Q2, "exact cross-checks disagree"
    if Q1 == 0: raise ValueError("zero source: Q = 0, state undefined")
    return dict(R=R1, Q=Q1, D=D_diag(a), cross=R1 - Q1, P0=R1/(d*Q1))
def R_character(a, d, t):
    """R_{a,t} = sum_{k<N} |sum_{n<=k} a(n) e^{2 pi i t n/d}|^2/(k(k+1)) + |...|^2/N (float)."""
    N = len(a) - 1; acc = 0j; tot = 0.0
    for k in range(1, N + 1):
        acc += a[k]*np.exp(2j*np.pi*t*k/d); tot += abs(acc)**2/(k*(k + 1)) if k < N else abs(acc)**2/N
    return tot
def negative_control(a, d):
    N = len(a) - 1; beta = np.zeros(d)
    for n in range(1, N + 1): beta[n % d] += a[n]
    den = d*float(np.sum(beta**2))
    return dict(beta=beta.tolist(), probability=(float(sum(a[1:]))**2/den) if den > 0 else None, defined=bool(den > 0))

# ------------------------------------------------------------------ feature register
def feature_gram_exact(N):
    """max |<g_m|g_n> - 1/max(m,n)| as an exact Fraction (should be 0)."""
    worst = Fr(0)
    for m in range(1, N + 1):
        for n in range(1, N + 1):
            ip = sum(Fr(1, k*(k + 1)) for k in range(max(m, n), N)) + Fr(1, N)
            worst = max(worst, abs(ip - Fr(1, max(m, n))))
    return worst
def feature_matrix(N):
    G = np.zeros((N, N))
    for n in range(1, N + 1):
        for k in range(n, N): G[n-1, k-1] = 1/math.sqrt(k*(k + 1))
        G[n-1, N-1] = 1/math.sqrt(N)
    return G
def amplitude_array(a, r):
    """normalized joint amplitudes, flat index (k-1)*d + b, auxiliary dimension padded to K = 2^s."""
    N = len(a) - 1; d = 2**r; s = max(1, math.ceil(math.log2(N))) if N > 1 else 1; K = 2**s
    A = residue_cumulatives(a, d); F = np.zeros((K, d))
    for k in range(1, N + 1):
        F[k-1, :] = A[k, :]/math.sqrt(k*(k + 1)) if k < N else A[N, :]/math.sqrt(N)
    Q = float(np.sum(F*F)); 
    if Q == 0: raise ValueError("zero source")
    return F.reshape(-1)/math.sqrt(Q), Q, s, K
def residue_coarsening_check(a, r):
    """q_b^(r) = q_b^(r+1) + q_{b+2^r}^(r+1) in the fixed-N feature space."""
    N = len(a) - 1; d = 2**r
    A1 = residue_cumulatives(a, d); A2 = residue_cumulatives(a, 2*d)
    return int(np.max(np.abs(A1 - (A2[:, :d] + A2[:, d:]))))

# ------------------------------------------------------------------ circuits
def prep_circuit(amps, r, s, label="prep"):
    qc = QuantumCircuit(r + s, name=label); qc.append(StatePreparation(amps), list(range(r + s))); return qc
def state_check(qc, amps):
    return float(np.max(np.abs(Statevector(qc).data - amps)))
def readout_circuit(prep, r):
    qc = prep.copy(); qc.add_register(ClassicalRegister(r, "res"))
    for j in range(r): qc.h(j)
    qc.measure(list(range(r)), list(range(r))); return qc
def P0_statevector(prep, r):
    qc = prep.copy()
    for j in range(r): qc.h(j)
    return float(Statevector(qc).probabilities(qargs=list(range(r)))[0])
def sample_counts(qc_meas, shots, seed):
    sim = AerSimulator(seed_simulator=seed); tq = transpile(qc_meas, sim, seed_transpiler=seed)
    return sim.run(tq, shots=shots).result().get_counts()
def wilson(k, n, z=1.959963984540054):
    if n == 0: return (float('nan'), float('nan'))
    p = k/n; den = 1 + z*z/n; c = (p + z*z/(2*n))/den; h = z*math.sqrt(p*(1 - p)/n + z*z/(4*n*n))/den
    return (max(0.0, c - h), min(1.0, c + h))
def resources(qc, basis=("rz", "sx", "x", "cx"), opt=1, seed=1):
    t0 = time.time(); tq = transpile(qc, basis_gates=list(basis), optimization_level=opt, seed_transpiler=seed); el = time.time() - t0
    ops = dict(tq.count_ops()); one = sum(v for k, v in ops.items() if k in ("rz", "sx", "x")); two = ops.get("cx", 0)
    return dict(depth=tq.depth(), ops=ops, one_qubit=int(one), two_qubit=int(two), transpile_s=round(el, 3)), tq

# ------------------------------------------------------------------ Weyl operators
def V_circuit(r, t):
    d = 2**r; qc = QuantumCircuit(r, name=f"V_{t}")
    for j in range(r): qc.p(2*math.pi*t*(2**j)/d, j)
    return qc
def T_circuit(r, h):
    """T_h |b> = |b + h mod 2^r> by h increments; each increment is the MCX cascade (high targets first), then X on qubit 0."""
    qc = QuantumCircuit(r, name=f"T_{h}")
    for _ in range(h % (2**r)):
        for j in range(r - 1, 0, -1): qc.mcx(list(range(j)), j)
        qc.x(0)
    return qc
def T_dense(r, h):
    d = 2**r; P = np.zeros((d, d))
    for b in range(d): P[(b + h) % d, b] = 1.0
    qc = QuantumCircuit(r, name=f"Tdense_{h}"); qc.append(UnitaryGate(P), list(range(r))); return qc
def weyl_checks(r, t, h):
    d = 2**r; V = Operator(V_circuit(r, t)).data; T = Operator(T_circuit(r, h)).data; Td = Operator(T_dense(r, h)).data
    Vexp = np.diag([np.exp(2j*np.pi*t*b/d) for b in range(d)]); Texp = np.zeros((d, d), dtype=complex)
    for b in range(d): Texp[(b + h) % d, b] = 1
    phase = np.exp(2j*np.pi*t*h/d)
    return dict(V_matrix_error=float(np.max(np.abs(V - Vexp))), T_matrix_error=float(np.max(np.abs(T - Texp))), T_dense_error=float(np.max(np.abs(Td - Texp))),
                commutation_error=float(np.max(np.abs(V @ T - phase*(T @ V)))), theta=2*math.pi*t*h/d)
def hadamard_test(r, t, h, which, input_circuit=None, extra_qubits=0):
    """ancilla = last qubit; controlled loop in chronological order T_h^dag, V_t^dag, T_h, V_t; X or Y measurement."""
    n = r + extra_qubits; qc = QuantumCircuit(n + 1, 1); anc = n
    if input_circuit is not None: qc.compose(input_circuit, qubits=list(range(input_circuit.num_qubits)), inplace=True)
    else: qc.x(0)                                                                     # simple residue input |b = 1>
    qc.h(anc)
    T = T_circuit(r, h).to_gate(); V = V_circuit(r, t).to_gate(); res = list(range(r))
    for g in (T.inverse(), V.inverse(), T, V):
        qc.append(g.control(1), [anc] + res)
    if which == "Y": qc.sdg(anc)
    qc.h(anc); qc.measure(anc, 0); return qc
def expectation_from_counts(counts):
    n0 = sum(v for k, v in counts.items() if k[-1] == "0"); n1 = sum(v for k, v in counts.items() if k[-1] == "1")
    return (n0 - n1)/(n0 + n1) if n0 + n1 else float('nan'), n0, n1
def expectation_exact(qc_meas):
    qc = qc_meas.remove_final_measurements(inplace=False); p = Statevector(qc).probabilities(qargs=[qc.num_qubits - 1]); return float(p[0] - p[1])
def orders_test(prep, r, t, h):
    """V_t T_h versus T_h V_t applied unconditionally: identical readout probabilities (global phase)."""
    out = {}
    for name, seq in (("VT", (T_circuit(r, h), V_circuit(r, t))), ("TV", (V_circuit(r, t), T_circuit(r, h)))):
        qc = prep.copy()
        for g in seq: qc.compose(g, qubits=list(range(r)), inplace=True)
        out[name] = P0_statevector(qc, r)
    qc = prep.copy(); qc.compose(T_circuit(r, h), qubits=list(range(r)), inplace=True); out["shift_only"] = P0_statevector(qc, r)
    out["base"] = P0_statevector(prep, r); return out
def character_probabilities(prep, r):
    """cyclic QFT on the residue register (QFT|b> = d^-1/2 sum_t e^{2 pi i bt/d}|t>, standard output order), marginal P_t."""
    qc = prep.copy(); qc.append(QFTGate(r), list(range(r)))
    return Statevector(qc).probabilities(qargs=list(range(r)))
def character_circuit(prep, r):
    qc = prep.copy(); qc.append(QFTGate(r), list(range(r))); qc.add_register(ClassicalRegister(r, "res")); qc.measure(list(range(r)), list(range(r))); return qc
def save_qpy(circuits, path):
    with open(path, "wb") as f: qpy.dump(circuits, f)

if __name__ == "__main__":
    from .arith import sieve_mu
    mu = sieve_mu(64)
    # fixtures
    fx = {"A_mobius_N4": ([0, 1, -1, -1, 0], 2, Fr(5, 6), Fr(11, 6), Fr(5, 44)), "B_sqfree_N4": ([0, 1, 1, 1, 0], 2, Fr(25, 6), Fr(11, 6), Fr(25, 44)),
          "C_packet_N6_d2": ([0, 1, -1, -1, 0, 0, 1], 1, Fr(2, 3), Fr(1), Fr(1, 3)), "C_packet_N6_d4": ([0, 1, -1, -1, 0, 0, 1], 2, Fr(2, 3), Fr(5, 3), Fr(1, 10)), "C_packet_N6_d8": ([0, 1, -1, -1, 0, 0, 1], 3, Fr(2, 3), Fr(2), Fr(1, 24))}
    for name, (a, r, Rx, Qx, Px) in fx.items():
        ex = exact_quantities(a, 2**r); amps, Q, s, K = amplitude_array(a, r); qc = prep_circuit(amps, r, s)
        print(f"{name}: R={ex['R']} ({Rx}) Q={ex['Q']} ({Qx}) P0={ex['P0']} ({Px}) exact match {ex['R']==Rx and ex['Q']==Qx and ex['P0']==Px}; state err {state_check(qc, amps):.1e}; P0 statevector {P0_statevector(qc, r):.12f}")
    N = 8; a = [0] + [1]*N; ex = exact_quantities(a, 8); print(f"D all-positive N=8 d=8: R={ex['R']} = 2N-H_N {2*N - sum(Fr(1,k) for k in range(1,N+1))}, Q={ex['Q']} = H_N, P0={ex['P0']}")
    print("feature Gram exact error N=16:", feature_gram_exact(16), "; coarsening check:", residue_coarsening_check([0]+[int(mu[n]) for n in range(1, 33)], 2))
    print("weyl r=2 t=1 h=1:", weyl_checks(2, 1, 1)); print("weyl r=3 t=1 h=1:", {k: (round(v, 12) if isinstance(v, float) else v) for k, v in weyl_checks(3, 1, 1).items()})
    for r in (2, 3):
        for which in ("X", "Y"):
            qc = hadamard_test(r, 1, 1, which); print(f"Hadamard test r={r} {which}: exact <{which}> = {expectation_exact(qc):.6f}  expected {math.cos(2*math.pi/2**r) if which=='X' else math.sin(2*math.pi/2**r):.6f}")
    a = [0] + [int(mu[n]) for n in range(1, 17)]; amps, Q, s, K = amplitude_array(a, 2); qc = prep_circuit(amps, 2, s)
    print("orders test N=16 r=2:", orders_test(qc, 2, 1, 1))
    pt = character_probabilities(qc, 2); ex = exact_quantities(a, 4); print("character probs:", pt, "formula:", [R_character(a, 4, t)/(4*float(ex['Q'])) for t in range(4)], "sum", pt.sum())
    counts = sample_counts(readout_circuit(qc, 2), 8192, 20261008); k = counts.get("00", 0); print("sampled N=16 r=2:", counts, "P0 hat", k/8192, "wilson", wilson(k, 8192), "exact", float(ex['P0']))
    t0 = time.time(); a = [0] + [int(mu[n]) for n in range(1, 65)]; amps, Q, s, K = amplitude_array(a, 5); qc = prep_circuit(amps, 5, s); res, tq = resources(qc); print("N=64 r=5 (11 qubits): prep", res, "state err", state_check(qc, amps), "total", round(time.time()-t0, 1), "s")
