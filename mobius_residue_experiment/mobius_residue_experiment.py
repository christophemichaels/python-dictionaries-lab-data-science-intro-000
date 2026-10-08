#!/usr/bin/env python
"""Michaels prime-zero shape experiment: command-line entry point.

  python mobius_residue_experiment.py verify                    exact fixtures and deterministic gates
  python mobius_residue_experiment.py arithmetic-shape --max-N 65536
  python mobius_residue_experiment.py prime-zero-shape --max-X 4096 --positive-zeros 200
  python mobius_residue_experiment.py simulate --shots 8192 --seed 20261008
  python mobius_residue_experiment.py plot --results <run_directory>
  python mobius_residue_experiment.py prepare-hardware --backend <name|fake_brisbane>
  python mobius_residue_experiment.py submit-hardware --backend <name> --confirm   (never runs without --confirm and a saved account)
  python mobius_residue_experiment.py all --results results/run

Every stage writes machine-readable CSV/JSON under the results directory; `plot` draws every figure from those files."""
import argparse, csv, json, math, os, sys, time, platform
from fractions import Fraction as Fr
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mre import arith, primezero

def fr(x):  # Fraction -> "num/den" string
    return f"{x.numerator}/{x.denominator}" if isinstance(x, Fr) else str(x)
def versions():
    out = dict(python=platform.python_version(), numpy=np.__version__)
    for m in ("qiskit", "qiskit_aer", "qiskit_ibm_runtime", "scipy", "mpmath", "flint", "matplotlib"):
        try: out[m] = __import__(m).__version__
        except Exception as e: out[m] = f"missing ({type(e).__name__})"
    return out
def write_csv(path, rows):
    if not rows: return
    keys = list(rows[0].keys())
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys); w.writeheader()
        for r in rows: w.writerow({k: (fr(v) if isinstance(v, Fr) else v) for k, v in r.items()})
def ensure(d): os.makedirs(d, exist_ok=True); return d

# ------------------------------------------------------------------ verify
def cmd_verify(args):
    from mre import quantum as qm
    out = ensure(args.results); gates = []; t0 = time.time()
    def gate(name, ok, detail): gates.append(dict(gate=name, status="PASS" if ok else "FAIL", detail=detail)); print(f"[{'PASS' if ok else 'FAIL'}] {name}: {detail}")
    mu = arith.sieve_mu(4096); gate("mobius sieve initial values n<=30", arith.check_mu(mu), str(arith.MU_FIRST_30))
    fixtures = {"A mobius N=4 d=4": ([0, 1, -1, -1, 0], 2, Fr(5, 6), Fr(11, 6), Fr(5, 44)), "B squarefree-positive N=4 d=4": ([0, 1, 1, 1, 0], 2, Fr(25, 6), Fr(11, 6), Fr(25, 44)),
                "C packet N=6 d=2": ([0, 1, -1, -1, 0, 0, 1], 1, Fr(2, 3), Fr(1), Fr(1, 3)), "C packet N=6 d=4": ([0, 1, -1, -1, 0, 0, 1], 2, Fr(2, 3), Fr(5, 3), Fr(1, 10)), "C packet N=6 d=8": ([0, 1, -1, -1, 0, 0, 1], 3, Fr(2, 3), Fr(2), Fr(1, 24))}
    fixture_rows = []
    for name, (a, r, Rx, Qx, Px) in fixtures.items():
        ex = qm.exact_quantities(a, 2**r); amps, Q, s, K = qm.amplitude_array(a, r); qc = qm.prep_circuit(amps, r, s)
        se = qm.state_check(qc, amps); psv = qm.P0_statevector(qc, r)
        ok = ex['R'] == Rx and ex['Q'] == Qx and ex['P0'] == Px and se < 1e-10 and abs(psv - float(Px)) < 1e-10
        gate(f"fixture {name}", ok, f"R={fr(ex['R'])} Q={fr(ex['Q'])} P0={fr(ex['P0'])}; state error {se:.1e}; statevector P0 {psv:.12f}")
        fixture_rows.append(dict(fixture=name, N=len(a) - 1, r=r, d=2**r, R_exact=fr(ex['R']), R_expected=fr(Rx), Q_exact=fr(ex['Q']), Q_expected=fr(Qx), P0_exact=fr(ex['P0']), P0_expected=fr(Px), state_error=se, P0_statevector=psv, status="PASS" if ok else "FAIL"))
    for N in (4, 8, 16):
        a = [0] + [1]*N; ex = qm.exact_quantities(a, 2**math.ceil(math.log2(N))); HN = sum(Fr(1, k) for k in range(1, N + 1)); d = 2**math.ceil(math.log2(N))
        ok = ex['R'] == 2*N - HN and ex['Q'] == HN and ex['P0'] == (2*N - HN)/(d*HN)
        gate(f"fixture D all-positive N={N} d={d}", ok, f"R={fr(ex['R'])} = 2N-H_N, Q={fr(ex['Q'])} = H_N, P0={fr(ex['P0'])}")
    for N in (4, 8, 16, 32): gate(f"feature Gram identity N={N} (exact)", qm.feature_gram_exact(N) == 0, "max |<g_m|g_n> - 1/max| = 0")
    a16 = [0] + [int(mu[n]) for n in range(1, 17)]
    for r in (1, 2, 3): gate(f"residue coarsening r={r}->{r+1} (N=16)", qm.residue_coarsening_check(a16, r) == 0, "max |q_b^(r) - q_b^(r+1) - q_(b+2^r)^(r+1)| = 0")
    for r in (2, 3, 4, 5):
        w = qm.weyl_checks(r, 1, 1); gate(f"Weyl matrices and commutation r={r} (t=h=1)", max(w['V_matrix_error'], w['T_matrix_error'], w['T_dense_error'], w['commutation_error']) < 1e-12, f"errors {w['V_matrix_error']:.1e} {w['T_matrix_error']:.1e} {w['T_dense_error']:.1e} {w['commutation_error']:.1e}, theta={w['theta']:.6f}")
        w0 = qm.weyl_checks(r, 0, 1); w1 = qm.weyl_checks(r, 1, 0); gate(f"Weyl zero controls r={r}", w0['commutation_error'] < 1e-12 and w1['commutation_error'] < 1e-12 and w0['theta'] == 0 and w1['theta'] == 0, "t=0 and h=0 give theta=0 and commuting matrices")
    for r in (2, 3, 4):
        th = 2*math.pi/2**r
        ex_ = {w: qm.expectation_exact(qm.hadamard_test(r, 1, 1, w)) for w in ("X", "Y")}
        gate(f"ancilla Hadamard test r={r} (controlled loop, chronological T^dag V^dag T V)", abs(ex_['X'] - math.cos(th)) < 1e-10 and abs(ex_['Y'] - math.sin(th)) < 1e-10, f"<X>={ex_['X']:.6f} cos={math.cos(th):.6f}; <Y>={ex_['Y']:.6f} sin={math.sin(th):.6f}")
    amps, Q, s, K = qm.amplitude_array(a16, 2); qc16 = qm.prep_circuit(amps, 2, s); ot = qm.orders_test(qc16, 2, 1, 1)
    gate("two Weyl orders identical probabilities (N=16 r=2)", abs(ot['VT'] - ot['TV']) < 1e-12, f"P0 after VT {ot['VT']:.12f}, after TV {ot['TV']:.12f}")
    gate("shift alone leaves P0 unchanged", abs(ot['shift_only'] - ot['base']) < 1e-12, f"{ot['shift_only']:.12f} vs {ot['base']:.12f}")
    # joint-state Hadamard test (residue + feature + ancilla)
    qcj = qm.hadamard_test(2, 1, 1, "Y", input_circuit=qc16, extra_qubits=s); gate("joint-state ancilla test (N=16 r=2, Y)", abs(qm.expectation_exact(qcj) - 1.0) < 1e-10, f"<Y>={qm.expectation_exact(qcj):.6f} expected 1")
    pt = qm.character_probabilities(qc16, 2); exq = qm.exact_quantities(a16, 4); form = [qm.R_character(a16, 4, t)/(4*float(exq['Q'])) for t in range(4)]
    gate("cyclic QFT character readout (N=16 r=2)", max(abs(pt[t] - form[t]) for t in range(4)) < 1e-10 and abs(pt.sum() - 1) < 1e-10 and abs(pt[1] - pt[3]) < 1e-12, f"P_t={np.round(pt, 8).tolist()} formula={np.round(form, 8).tolist()} sum={pt.sum():.12f}")
    gate("uniform projection equals R/(dQ) (statevector, N=16 r=2)", abs(qm.P0_statevector(qc16, 2) - float(exq['P0'])) < 1e-10, f"{qm.P0_statevector(qc16, 2):.12f} vs {float(exq['P0']):.12f}")
    for r in (4, 5):
        exs = qm.exact_quantities(a16, 2**r); gate(f"saturation d={2**r} >= N=16: Q = D", exs['Q'] == exs['D'], f"Q={fr(exs['Q'])} D={fr(exs['D'])} P0={fr(exs['P0'])}")
    # arithmetic exact identities
    for N in (6, 7, 12, 36, 64, 100, 128, 256):
        ex = arith.exact_case(mu, N); sw = arith.ArithSweep(300, mu=mu[:301]); row = sw.row(N)
        ok = ex['delta'] == ex['delta_via'] and abs(float(ex['A']) - row['A_selected']) < 1e-12 and abs(float(ex['V']) - row['V']) < 1e-12
        gate(f"arithmetic identities N={N}: delta = kappa + V exact; A pair sum = prefix formula; V direct = prefix", ok, f"delta={fr(ex['delta'])[:40]}{'...' if len(fr(ex['delta']))>40 else ''}, C={float(ex['C']):.9f}, Gamma={float(ex['Gamma']):.9f}")
    for N in (64, 128):
        rows_b, agg, comp, st = arith.block_decomposition(mu, N); ex = arith.exact_case(mu, N)
        gate(f"block decomposition N={N} aggregates to V - A", abs(agg - float(ex['C'])) < 1e-10 and max(r['svd_error'] for r in rows_b) < 1e-12 and comp, f"{len(rows_b)} blocks, aggregate {agg:.12f} vs C {float(ex['C']):.12f}, max SVD split error {max(r['svd_error'] for r in rows_b):.1e}")
    # prime-zero identity on a small case
    st = primezero.Staircase(200); ords, zh, cache = primezero.zeros_cached(10, os.path.join(out, "zeros_mpmath.json"))
    rh = primezero.zero_multiset(ords); P = st.P_exact_weights(4, 6); q = primezero.window_integrals(st, 4, 6, rh); g = primezero.mode_gram(4, 6, rh); zq = primezero.zero_part_quadrature(4, 6, rh)
    gate("prime-zero return identity P = Z + I + D (X=4, 10 zeros)", abs(P - (q['Z'] + q['I'] + q['D'])) < 1e-8 and abs(P - q['P_quadrature']) < 1e-8, f"P={P:.9f} Z={q['Z']:.9f} I={q['I']:+.9f} D={q['D']:.9f}")
    gate("mode Gram energy equals quadrature (X=4, 10 zeros)", abs(g['E_modes'] - zq) < 1e-8, f"Gram {g['E_modes']:.9f} quadrature {zq:.9f}; bound utilization {g['utilization']:.4f}")
    gate("norm bounds |sqrt P - sqrt Z| <= sqrt D and |I| <= 2 sqrt(ZD)", abs(math.sqrt(P) - math.sqrt(q['Z'])) <= math.sqrt(q['D']) + 1e-12 and abs(q['I']) <= 2*math.sqrt(q['Z']*q['D']) + 1e-12, "consistency of the norm structure")
    res = dict(gates=gates, fixtures=fixture_rows, versions=versions(), elapsed_s=round(time.time() - t0, 2), all_pass=all(g['status'] == "PASS" for g in gates))
    json.dump(res, open(os.path.join(out, "verify.json"), "w"), indent=1); write_csv(os.path.join(out, "fixtures.csv"), fixture_rows)
    print(f"verify: {sum(g['status']=='PASS' for g in gates)}/{len(gates)} gates pass; {res['elapsed_s']} s")
    return res

# ------------------------------------------------------------------ arithmetic shape
def cmd_arith(args):
    out = ensure(os.path.join(args.results, "arithmetic")); t0 = time.time(); Nmax = args.max_N
    sw = arith.ArithSweep(Nmax); rows = sw.sweep(6, Nmax)
    Ns = np.array([r['N'] for r in rows]); C = np.array([r['C_complement'] for r in rows]); G = np.array([r['Gamma'] for r in rows]); dl = np.array([r['delta'] for r in rows])
    rec = {"C": arith.running_records(Ns, C), "Gamma": arith.running_records(Ns, G), "delta": arith.running_records(Ns, dl), "c_normalized": arith.running_records(Ns, [r['c_normalized'] for r in rows])}
    for eta in (0.05, 0.1, 0.2):
        rec[f"Cpos_over_N^{eta}"] = arith.running_records(Ns, np.maximum(C, 0)/Ns**eta); rec[f"Gammapos_over_N^{eta}"] = arith.running_records(Ns, np.maximum(G, 0)/Ns**eta)
    recset = set()
    for k, v in rec.items():
        for N, _ in v: recset.add(N)
    for r in rows: r['record_flags'] = ",".join(k for k, v in rec.items() if any(N == r['N'] for N, _ in v)); r['numerical_precision'] = "float64 (direct positive sum for V; prefix arrays for A)"; r['identity_errors'] = max(r['err_delta'], r['err_V'])
    write_csv(os.path.join(out, "arith_sweep.csv"), rows)
    thresholds = [k**6 for k in range(2, 20) if k**6 <= Nmax]
    summary = dict(Nmax=Nmax, rows=len(rows), max_err_delta=float(max(r['err_delta'] for r in rows)), max_err_V=float(max(r['err_V'] for r in rows)), Kred_thresholds=thresholds,
                   max_C=[int(Ns[np.argmax(C)]), float(C.max())], max_Gamma=[int(Ns[np.argmax(G)]), float(G.max())], max_delta=[int(Ns[np.argmax(dl)]), float(dl.max())],
                   min_C=[int(Ns[np.argmin(C)]), float(C.min())], frac_C_positive=float(np.mean(C > 0)), frac_Gamma_positive=float(np.mean(G > 0)),
                   max_c_normalized=[int(Ns[np.argmax([r['c_normalized'] for r in rows])]), float(max(r['c_normalized'] for r in rows))],
                   held_out={"declared_range": [6, Nmax//4], "held_out_range": [Nmax//4 + 1, Nmax], "max_C_declared": float(C[Ns <= Nmax//4].max()), "max_C_held_out": float(C[Ns > Nmax//4].max()), "max_Gamma_declared": float(G[Ns <= Nmax//4].max()), "max_Gamma_held_out": float(G[Ns > Nmax//4].max())},
                   records=rec)
    # certified balls at selected cutoffs: the records of C and Gamma, the sixth-power thresholds, and a few fixed ones
    sel = sorted(set([int(N) for N, _ in rec['C'][-6:]] + [int(N) for N, _ in rec['Gamma'][-6:]] + thresholds + [t - 1 for t in thresholds] + [64, 512, 4096, Nmax]))
    balls = []
    for N in sel:
        if N < 6 or N > Nmax: continue
        b = arith.ball_case(sw.mu, N); row = sw.row(N)
        balls.append(dict(N=N, L=b['L'], K_red=b['K_red'], R=arith.ball_str(b['R']), delta=arith.ball_str(b['delta']), delta_via=arith.ball_str(b['delta_via']), V=arith.ball_str(b['V']), A=arith.ball_str(b['A']), C=arith.ball_str(b['C']), Gamma=arith.ball_str(b['Gamma']), kappa=arith.ball_str(b['kappa']),
                          float_C=row['C_complement'], float_Gamma=row['Gamma'], float_minus_ball_C=float(row['C_complement'] - float(b['C'].mid())), float_minus_ball_Gamma=float(row['Gamma'] - float(b['Gamma'].mid()))))
    json.dump(dict(summary=summary, balls=balls), open(os.path.join(out, "arith_summary.json"), "w"), indent=1)
    # exact rationals for small cutoffs
    exact_rows = []
    for N in (6, 7, 8, 12, 36, 64, 128, 256, 512):
        if N > Nmax: continue
        ex = arith.exact_case(sw.mu, N); exact_rows.append({k: (fr(v) if isinstance(v, Fr) else v) for k, v in ex.items()} | dict(C_float=float(ex['C']), Gamma_float=float(ex['Gamma']), delta_float=float(ex['delta'])))
    write_csv(os.path.join(out, "arith_exact_small.csv"), exact_rows)
    # block decompositions: complete where feasible, partial at the record excursions
    brows = []; bsum = []
    full_list = [N for N in (64, 128, 256, 512, 1024, 2048, 4096) if N <= Nmax]
    partial_list = sorted(set([int(N) for N, _ in rec['C'][-3:]] + [int(N) for N, _ in rec['Gamma'][-3:]]) - set(full_list))
    for N in full_list + partial_list:
        t1 = time.time(); rws, agg, comp, st = arith.block_decomposition(sw.mu, N, entries_budget=args.block_entries_budget); row = sw.row(N)
        Ep = sum(r['E_plus'] for r in rws); Em = sum(r['E_minus'] for r in rws)
        for r in rws: r['aggregate_coverage'] = "complete" if comp else "partial"; r['reconstruction_error'] = abs(agg - row['C_complement']) if comp else float('nan')
        brows += rws
        bsum.append(dict(N=N, blocks=len(rws), complete=comp, E_plus_total=Ep, E_minus_total=Em, signed_total=agg, C=row['C_complement'], carry=row['carry'], Gamma=row['Gamma'], reconstruction_error=abs(agg - row['C_complement']) if comp else None, entries=st['entries'], elapsed_s=round(time.time() - t1, 2),
                         top_blocks=sorted([(r['g'], r['U'], r['V'], r['signed_block']) for r in rws], key=lambda x: -x[3])[:5], max_svd_error=max((r['svd_error'] for r in rws), default=0.0)))
        print(f"blocks N={N}: {len(rws)} blocks, {'complete' if comp else 'PARTIAL'}, E+={Ep:.6f} E-={Em:.6f} signed={agg:.9f} C={row['C_complement']:.9f} ({time.time()-t1:.1f} s)")
    write_csv(os.path.join(out, "block_spectrum.csv"), brows); json.dump(bsum, open(os.path.join(out, "block_summary.json"), "w"), indent=1)
    print(f"arithmetic-shape: {len(rows)} cutoffs, max identity error {summary['max_err_delta']:.1e}/{summary['max_err_V']:.1e}; max C {summary['max_C']}, max Gamma {summary['max_Gamma']}; {round(time.time()-t0, 1)} s")
    return summary

# ------------------------------------------------------------------ prime-zero shape
def cmd_primezero(args):
    out = ensure(os.path.join(args.results, "primezero")); t0 = time.time()
    ords, zh, cache = primezero.zeros_cached(args.positive_zeros, os.path.join(out, "zeros_mpmath.json"), dps=args.zero_dps)
    Xs = [2**k for k in range(1, int(math.log2(args.max_X)) + 1)]; counts = [c for c in (10, 25, 50, 100, 200, 400, 800) if c <= args.positive_zeros]
    st = primezero.Staircase(6*args.max_X + 64)
    rows = primezero.campaign(st, ords, Xs, counts, a=6, zero_hash=zh)
    write_csv(os.path.join(out, "primezero.csv"), rows)
    # window samples for plots
    samples = {}
    for X in (64, 1024):
        if X > args.max_X: continue
        x = np.linspace(X, 6*X, 3000); rh = primezero.zero_multiset(ords); samples[str(X)] = dict(x=x.tolist(), e=(st.psi0(x) - x).tolist(), eZ=primezero.e_Z(x, rh).tolist(), zeros=len(ords))
    json.dump(dict(zero_hash=zh, provenance={k: v for k, v in cache.items() if k != 'ordinates' and k != 'gram_info'}, gram_info_first=cache['gram_info'][:5], samples=samples, elapsed_s=round(time.time() - t0, 1)), open(os.path.join(out, "primezero_meta.json"), "w"), indent=1)
    print(f"prime-zero-shape: {len(rows)} rows; zeros {len(ords)} (max height {ords[-1]:.3f}); max identity residual {max(abs(r['identity_residual']) for r in rows):.1e}; all bounds ok {all(r['bound1_ok'] and r['bound2_ok'] and r['gram_bound_ok'] for r in rows)}; {round(time.time()-t0,1)} s")
    return rows

# ------------------------------------------------------------------ simulate
def cmd_simulate(args):
    from mre import quantum as qm
    out = ensure(os.path.join(args.results, "quantum")); ensure(os.path.join(out, "circuits")); t0 = time.time()
    mu = arith.sieve_mu(max(args.N_list) + 1); srcs, seed_src = qm.make_sources(mu, Nmax=max(args.N_list), seed=args.seed)
    vers = versions(); rows = []; saved = []
    for sname, afull in srcs.items():
        for N in args.N_list:
            a = afull[:N + 1]; sh = qm.source_hash(a)
            for r in args.r_list:
                d = 2**r; t1 = time.time(); ex = qm.exact_quantities(a, d); amps, Q, s, K = qm.amplitude_array(a, r)
                qc = qm.prep_circuit(amps, r, s, label=f"{sname}_N{N}_r{r}"); se = qm.state_check(qc, amps); psv = qm.P0_statevector(qc, r)
                qcm = qm.readout_circuit(qc, r); counts = qm.sample_counts(qcm, args.shots, args.seed + N*7 + r); k = counts.get("0"*r, 0); lo, hi = qm.wilson(k, args.shots)
                prep_res, _ = qm.resources(qc); full_res, tq = qm.resources(qcm)
                Rhat = d*float(ex['Q'])*k/args.shots
                row = dict(source_name=sname, source_seed=seed_src if sname == "random_sign" else None, source_hash=sh, N=N, r=r, d=d, feature_qubits=s, data_qubits=r + s,
                           R_exact=fr(ex['R']), R_exact_float=float(ex['R']), Q_exact=fr(ex['Q']), Q_exact_float=float(ex['Q']), D_exact=fr(ex['D']), D_exact_float=float(ex['D']), cross_channel_exact=fr(ex['cross']), cross_channel_float=float(ex['cross']),
                           P_classical=float(ex['P0']), P_classical_exact=fr(ex['P0']), P_statevector=psv, successes=int(k), shots=args.shots, P_sampled=k/args.shots, P_ci_low=lo, P_ci_high=hi,
                           R_reconstructed=Rhat, R_ci_low=d*float(ex['Q'])*lo, R_ci_high=d*float(ex['Q'])*hi, dP0_classical=d*float(ex['P0']), dP0_sampled=d*k/args.shots, R_over_D=float(ex['R']/ex['D']),
                           normalization_error=abs(float(np.sum(amps**2)) - 1.0), probability_error=abs(psv - float(ex['P0'])), state_error=se, covered=bool(lo <= float(ex['P0']) <= hi),
                           prep_depth=prep_res['depth'], full_depth=full_res['depth'], one_qubit_gates=full_res['one_qubit'], two_qubit_gates=full_res['two_qubit'], prep_two_qubit=prep_res['two_qubit'],
                           backend_or_simulator="AerSimulator (ideal, statevector method)", software_versions=json.dumps({k: vers[k] for k in ('qiskit', 'qiskit_aer')}), elapsed_times=json.dumps(dict(case_s=round(time.time() - t1, 3), transpile_s=full_res['transpile_s'])), job_id_if_any=None)
                rows.append(row)
                if (N, r) in ((4, 2), (16, 2), (64, 5)) and sname in ("mobius",): saved.append((f"{sname}_N{N}_r{r}", qcm))
                print(f"{sname:20s} N={N:3d} r={r}: P0 exact {float(ex['P0']):.6f} sv {psv:.6f} sampled {k/args.shots:.6f} [{lo:.4f},{hi:.4f}] covered={lo <= float(ex['P0']) <= hi}; R={float(ex['R']):.4f} Rhat={Rhat:.4f}; prep depth {prep_res['depth']} cx {prep_res['two_qubit']}")
    write_csv(os.path.join(out, "cases.csv"), rows); json.dump(rows, open(os.path.join(out, "cases.json"), "w"), indent=1)
    for name, qc in saved: qm.save_qpy([qc], os.path.join(out, "circuits", f"{name}.qpy"))
    # subset: phase-modified readout (V_1 then H) sampled, characters, negative control, Weyl sampled, transpiler report, noise model
    extra = {}
    a16 = srcs["mobius"][:17]; amps, Q, s, K = qm.amplitude_array(a16, 2); qc16 = qm.prep_circuit(amps, 2, s); ex16 = qm.exact_quantities(a16, 4)
    qv = qc16.copy(); qv.compose(qm.V_circuit(2, 1), qubits=[0, 1], inplace=True); qvm = qm.readout_circuit(qv, 2); cnt = qm.sample_counts(qvm, args.shots, args.seed + 99); k = cnt.get("00", 0)
    extra["phase_modified_readout_N16_r2_t1"] = dict(P_t_formula=qm.R_character(a16, 4, 1)/(4*float(ex16['Q'])), P_statevector=qm.P0_statevector(qv, 2), sampled=k/args.shots, ci=qm.wilson(k, args.shots), counts=cnt)
    chars = {}
    for sname in srcs:
        for N, r in ((16, 2), (32, 3)):
            a = srcs[sname][:N + 1]; amps, Q, s, K = qm.amplitude_array(a, r); qc = qm.prep_circuit(amps, r, s); ex = qm.exact_quantities(a, 2**r)
            pt = qm.character_probabilities(qc, r); form = [qm.R_character(a, 2**r, t)/(2**r*float(ex['Q'])) for t in range(2**r)]
            chars[f"{sname}_N{N}_r{r}"] = dict(P_t=pt.tolist(), formula=form, sum=float(pt.sum()), max_error=float(max(abs(pt[t] - form[t]) for t in range(2**r))), symmetric=bool(max(abs(pt[t] - pt[(-t) % 2**r]) for t in range(2**r)) < 1e-12))
    qcc = qm.character_circuit(qc16, 2); cntc = qm.sample_counts(qcc, args.shots, args.seed + 7); chars["sampled_mobius_N16_r2"] = dict(counts=cntc, P_t_sampled=[cntc.get(format(t, "02b")[::-1] if False else format(t, "02b"), 0)/args.shots for t in range(4)], note="counts strings are little-endian: bitstring b1b0, t = b0 + 2 b1")
    extra["characters"] = chars
    extra["negative_control"] = {f"{sname}_N{N}_d{d}": qm.negative_control(srcs[sname][:N + 1], d) for sname in srcs for N in (4, 16) for d in (4, 8)}
    weyl = []
    for r in (2, 3):
        for which in ("X", "Y"):
            qc = qm.hadamard_test(r, 1, 1, which); exv = qm.expectation_exact(qc); cnt = qm.sample_counts(qc, args.shots, args.seed + 11*r + (which == "Y")); ev, n0, n1 = qm.expectation_from_counts(cnt)
            th = 2*math.pi/2**r; se_ = 2*math.sqrt(max(1e-12, (1 - ev*ev))/args.shots)
            from qiskit import transpile as _tr
            t3 = _tr(qc, basis_gates=["rz", "sx", "x", "cx"], optimization_level=3, seed_transpiler=3); t1_ = _tr(qc, basis_gates=["rz", "sx", "x", "cx"], optimization_level=1, seed_transpiler=3)
            weyl.append(dict(r=r, d=2**r, t=1, h=1, theta=th, which=which, expected=math.cos(th) if which == "X" else math.sin(th), exact=exv, sampled=ev, sampled_2sigma=se_, counts=cnt, shots=args.shots,
                             transpiled_ops_opt1=dict(t1_.count_ops()), depth_opt1=t1_.depth(), transpiled_ops_opt3=dict(t3.count_ops()), depth_opt3=t3.depth(), sign_convention="ancilla |+>, loop controlled on ancilla=1 in order T^dag, V^dag, T, V (rightmost of W = V T V^dag T^dag first); X: H then Z-measure; Y: S^dag, H, Z-measure; <X> = cos theta, <Y> = sin theta with theta = 2 pi t h/d"))
            print(f"Weyl r={r} {which}: exact {exv:.6f} expected {math.cos(th) if which=='X' else math.sin(th):.6f} sampled {ev:.4f} +- {se_:.4f}; opt3 ops {dict(t3.count_ops())}")
    extra["weyl"] = weyl
    try:
        from qiskit_ibm_runtime.fake_provider import FakeManilaV2
        from qiskit_aer.noise import NoiseModel
        fb = FakeManilaV2(); nm = NoiseModel.from_backend(fb); sim = qm.AerSimulator(noise_model=nm, seed_simulator=args.seed)
        a4 = srcs["mobius"][:5]; amps, Q, s, K = qm.amplitude_array(a4, 2); qc4 = qm.readout_circuit(qm.prep_circuit(amps, 2, s), 2)
        from qiskit import transpile as _tr
        tq = _tr(qc4, backend=fb, optimization_level=1, seed_transpiler=1); cnt = sim.run(tq, shots=args.shots).result().get_counts(); k = cnt.get("00", 0)
        extra["noise_model"] = dict(model=f"NoiseModel.from_backend(FakeManilaV2) [modeled noise, not device data]", case="mobius N=4 r=2", P0_exact=5/44, P0_noisy=k/args.shots, ci=qm.wilson(k, args.shots), counts=cnt, depth=tq.depth(), ops=dict(tq.count_ops()))
    except Exception as e:
        extra["noise_model"] = dict(error=f"{type(e).__name__}: {e}")
    extra["versions"] = vers; extra["elapsed_s"] = round(time.time() - t0, 1); extra["shots"] = args.shots; extra["seed"] = args.seed
    json.dump(extra, open(os.path.join(out, "extra.json"), "w"), indent=1, default=str)
    cov = sum(r['covered'] for r in rows); print(f"simulate: {len(rows)} cases, Wilson 95% coverage {cov}/{len(rows)}, max probability error {max(r['probability_error'] for r in rows):.1e}, max state error {max(r['state_error'] for r in rows):.1e}; {round(time.time()-t0,1)} s")
    return rows

# ------------------------------------------------------------------ hardware
def cmd_prepare_hardware(args):
    from mre import quantum as qm
    from qiskit.transpiler import generate_preset_pass_manager
    out = ensure(os.path.join(args.results, "hardware")); t0 = time.time()
    backend = None; account = False; note = ""
    try:
        from qiskit_ibm_runtime import QiskitRuntimeService
        svc = QiskitRuntimeService(); account = True
        backend = svc.backend(args.backend) if args.backend and not args.backend.startswith("fake") else svc.least_busy(operational=True, simulator=False); note = f"real backend target {backend.name} (no job submitted)"
    except Exception as e:
        from qiskit_ibm_runtime.fake_provider import FakeBrisbane, FakeManilaV2
        backend = FakeBrisbane() if (args.backend in (None, "", "fake_brisbane")) else FakeManilaV2(); note = f"no saved IBM account ({type(e).__name__}); transpiled against the fake backend {backend.name} as a stand-in target"
    pm = generate_preset_pass_manager(optimization_level=3, backend=backend, seed_transpiler=args.seed)
    mu = arith.sieve_mu(8)
    cands = {"mobius_N4_r2": ([0, 1, -1, -1, 0], 2, 5/44), "squarefree_positive_N4_r2": ([0, 1, 1, 1, 0], 2, 25/44), "packet_N6_r1": ([0, 1, -1, -1, 0, 0, 1], 1, 1/3)}
    batch = []; circuits = []
    for name, (a, r, p0) in cands.items():
        amps, Q, s, K = qm.amplitude_array(a, r); qc = qm.readout_circuit(qm.prep_circuit(amps, r, s), r); isa = pm.run(qc); ops = dict(isa.count_ops())
        batch.append(dict(name=name, data_qubits=r + s, P0_exact=p0, d=2**r, Q=Q, isa_depth=isa.depth(), isa_ops=ops, isa_two_qubit=sum(v for k, v in ops.items() if k in ("cz", "cx", "ecr")), success_bitstring="0"*r, shots=args.shots)); circuits.append(isa)
    for which in ("X", "Y"):
        qc = qm.hadamard_test(2, 1, 1, which); isa = pm.run(qc); ops = dict(isa.count_ops())
        batch.append(dict(name=f"weyl_loop_r2_{which}", data_qubits=3, expected=(0.0 if which == "X" else 1.0), theta=math.pi/2, isa_depth=isa.depth(), isa_ops=ops, isa_two_qubit=sum(v for k, v in ops.items() if k in ("cz", "cx", "ecr")), shots=args.shots)); circuits.append(isa)
    qm.save_qpy(circuits, os.path.join(out, "isa_batch.qpy"))
    est = dict(circuits=len(batch), shots_per_circuit=args.shots, total_shots=len(batch)*args.shots, backend=getattr(backend, "name", str(backend)), account_found=account, note=note,
               authorization="NOT SUBMITTED: hardware execution requires explicit user authorization; run submit-hardware --confirm with a saved account", recommended_subset=["mobius_N4_r2", "packet_N6_r1", "weyl_loop_r2_X", "weyl_loop_r2_Y"], elapsed_s=round(time.time() - t0, 1))
    json.dump(dict(batch=batch, estimate=est), open(os.path.join(out, "batch.json"), "w"), indent=1)
    open(os.path.join(out, "RUN_COMMAND.txt"), "w").write(f"# pending hardware stage; nothing submitted\npython mobius_residue_experiment.py submit-hardware --results {args.results} --backend <backend_name> --shots {args.shots} --confirm\n")
    print(f"prepare-hardware: {note}; batch of {len(batch)} circuits saved (ISA depths {[b['isa_depth'] for b in batch]}, two-qubit {[b['isa_two_qubit'] for b in batch]}); status PENDING, not submitted")
    return batch

def cmd_submit_hardware(args):
    if not args.confirm:
        print("submit-hardware: --confirm not given; nothing submitted (hardware stage pending)"); return
    try:
        from qiskit_ibm_runtime import QiskitRuntimeService, SamplerV2
        from qiskit import qpy
        svc = QiskitRuntimeService(); backend = svc.backend(args.backend)
        with open(os.path.join(args.results, "hardware", "isa_batch.qpy"), "rb") as f: circuits = qpy.load(f)
        sampler = SamplerV2(mode=backend); job = sampler.run(circuits, shots=args.shots)
        meta = dict(job_id=job.job_id(), backend=backend.name, shots=args.shots, circuits=len(circuits), submitted=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()))
        json.dump(meta, open(os.path.join(args.results, "hardware", "job.json"), "w"), indent=1); print("submitted:", meta)
    except Exception as e:
        print(f"submit-hardware: not submitted ({type(e).__name__}: {e})")

# ------------------------------------------------------------------ plot
def cmd_plot(args):
    from mre import plots
    plots.make_all(args.results); print("plots written to", os.path.join(args.results, "plots"))

def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--results", default="results/run")
    sub = p.add_subparsers(dest="cmd", required=True)
    sub.add_parser("verify")
    a = sub.add_parser("arithmetic-shape"); a.add_argument("--max-N", type=int, default=65536); a.add_argument("--block-entries-budget", type=int, default=2*10**8)
    z = sub.add_parser("prime-zero-shape"); z.add_argument("--max-X", type=int, default=4096); z.add_argument("--positive-zeros", type=int, default=200); z.add_argument("--zero-dps", type=int, default=30)
    s = sub.add_parser("simulate"); s.add_argument("--shots", type=int, default=8192); s.add_argument("--seed", type=int, default=20261008); s.add_argument("--N-list", type=int, nargs="+", default=[4, 8, 16, 32, 64]); s.add_argument("--r-list", type=int, nargs="+", default=[2, 3, 4, 5])
    sub.add_parser("plot")
    h = sub.add_parser("prepare-hardware"); h.add_argument("--backend", default="fake_brisbane"); h.add_argument("--shots", type=int, default=4000); h.add_argument("--seed", type=int, default=20261008)
    sh = sub.add_parser("submit-hardware"); sh.add_argument("--backend", required=True); sh.add_argument("--shots", type=int, default=4000); sh.add_argument("--confirm", action="store_true")
    al = sub.add_parser("all"); al.add_argument("--max-N", type=int, default=65536); al.add_argument("--max-X", type=int, default=4096); al.add_argument("--positive-zeros", type=int, default=200); al.add_argument("--shots", type=int, default=8192); al.add_argument("--seed", type=int, default=20261008)
    args = p.parse_args(); ensure(args.results)
    cfg = dict(command=args.cmd, args={k: v for k, v in vars(args).items()}, versions=versions(), time=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()))
    json.dump(cfg, open(os.path.join(args.results, f"config_{args.cmd}.json"), "w"), indent=1)
    if args.cmd == "verify": cmd_verify(args)
    elif args.cmd == "arithmetic-shape": cmd_arith(args)
    elif args.cmd == "prime-zero-shape": cmd_primezero(args)
    elif args.cmd == "simulate": cmd_simulate(args)
    elif args.cmd == "plot": cmd_plot(args)
    elif args.cmd == "prepare-hardware": cmd_prepare_hardware(args)
    elif args.cmd == "submit-hardware": cmd_submit_hardware(args)
    elif args.cmd == "all":
        class A: pass
        for cmd, extra in (("verify", {}), ("arithmetic-shape", dict(max_N=args.max_N, block_entries_budget=2*10**8)), ("prime-zero-shape", dict(max_X=args.max_X, positive_zeros=args.positive_zeros, zero_dps=30)),
                           ("simulate", dict(shots=args.shots, seed=args.seed, N_list=[4, 8, 16, 32, 64], r_list=[2, 3, 4, 5])), ("prepare-hardware", dict(backend="fake_brisbane", shots=4000, seed=args.seed)), ("plot", {})):
            ns = argparse.Namespace(results=args.results, cmd=cmd, **extra); print(f"==== {cmd} ===="); globals()["cmd_" + cmd.replace("-", "_").replace("prime_zero_shape", "primezero").replace("arithmetic_shape", "arith")](ns)

if __name__ == "__main__":
    main()
