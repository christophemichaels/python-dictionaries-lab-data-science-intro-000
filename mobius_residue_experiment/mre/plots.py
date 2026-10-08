"""All figures, drawn from the machine-readable results of a run directory."""
import csv, json, math, os
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
BLUE, ORANGE, INK, MUTED, GRID, SURF = "#2a78d6", "#eb6834", "#0b0b0b", "#52514e", "#e6e5e1", "#fcfcfb"
plt.rcParams.update({"font.size": 9.5, "axes.edgecolor": MUTED, "axes.labelcolor": INK, "xtick.color": MUTED, "ytick.color": MUTED, "text.color": INK, "axes.titlecolor": INK, "figure.facecolor": SURF, "axes.facecolor": SURF})
def style(ax):
    ax.grid(True, color=GRID, lw=0.8); ax.set_axisbelow(True)
    for s in ("top", "right"): ax.spines[s].set_visible(False)
def read_csv(path):
    import gzip
    if not os.path.exists(path) and os.path.exists(path + ".gz"): path = path + ".gz"
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as f: return list(csv.DictReader(f))
def num(rows, k): return np.array([float(r[k]) for r in rows])
def save(fig, out, name): fig.tight_layout(); fig.savefig(os.path.join(out, name), dpi=150, facecolor=SURF); plt.close(fig)
def mark_thresholds(ax, thr, ymin=None):
    for t in thr: ax.axvline(t, color=MUTED, lw=0.7, ls=":")

def arithmetic_plots(res, out):
    d = os.path.join(res, "arithmetic")
    if not (os.path.exists(os.path.join(d, "arith_sweep.csv")) or os.path.exists(os.path.join(d, "arith_sweep.csv.gz"))): return
    rows = read_csv(os.path.join(d, "arith_sweep.csv")); summ = json.load(open(os.path.join(d, "arith_summary.json")))['summary']
    N = num(rows, "N"); R = num(rows, "R"); S = num(rows, "S"); dl = num(rows, "delta"); A = num(rows, "A_selected"); C = num(rows, "C_complement"); kap = num(rows, "carry"); G = num(rows, "Gamma"); cN = num(rows, "c_normalized"); GB = num(rows, "Gamma_over_B5")
    thr = summ['Kred_thresholds']
    # 1. R and S
    fig, ax = plt.subplots(figsize=(7.5, 3.6)); style(ax); ax.plot(N, R, color=BLUE, lw=0.8, label="R(N) (terminal included)"); ax.plot(N, S, color=ORANGE, lw=1.4, label="running envelope S(N)"); ax.set_xscale("log"); ax.set_xlabel("N"); ax.set_ylabel("energy"); ax.set_title("Arithmetic view 1: Möbius energy and its envelope", loc="left"); ax.legend(frameon=False); mark_thresholds(ax, thr); save(fig, out, "arith_1_R_S.png")
    # 2. signed components
    fig, ax = plt.subplots(figsize=(9, 4)); style(ax)
    for y, lab, col, lw in ((dl, "δ_N", INK, 0.8), (A, "A_N selected", BLUE, 0.8), (C, "C_N = V − A", ORANGE, 0.9), (kap, "κ_N carry", MUTED, 0.7), (G, "Γ_N = κ + C", "#6a3d9a", 0.9)):
        ax.plot(N, y, color=col, lw=lw, label=lab)
    ax.axhline(0, color=INK, lw=0.8); ax.set_xscale("log"); ax.set_xlabel("N"); ax.set_ylabel("signed value (linear scale)"); ax.set_title("Arithmetic view 2: signed increment, selected family, complement, carry, Γ", loc="left"); ax.legend(frameon=False, ncol=5, fontsize=8); mark_thresholds(ax, thr); save(fig, out, "arith_2_signed.png")
    # 3. scaled by N^eta with running maxima of positive parts
    fig, axs = plt.subplots(2, 1, figsize=(9, 6.5), sharex=True)
    for ax, y, name in ((axs[0], C, "C_N"), (axs[1], G, "Γ_N")):
        style(ax); ax.plot(N, y, color=MUTED, lw=0.6, label=f"{name} (unscaled)")
        for eta, col in ((0.05, BLUE), (0.1, ORANGE), (0.2, "#6a3d9a")):
            v = np.maximum(y, 0)/N**eta; ax.plot(N, np.maximum.accumulate(v), color=col, lw=1.3, label=f"running max of [{name}]₊/N^{eta}")
        ax.axhline(0, color=INK, lw=0.8); ax.set_xscale("log"); ax.set_ylabel(name); ax.legend(frameon=False, fontsize=8, ncol=2); mark_thresholds(ax, thr)
    axs[0].set_title("Arithmetic view 3: positive parts over N^η and their running maxima", loc="left"); axs[1].set_xlabel("N"); save(fig, out, "arith_3_scaled_records.png")
    # 4. c_N and Gamma/B^5
    fig, ax = plt.subplots(figsize=(9, 3.8)); style(ax); ax.plot(N, cN, color=BLUE, lw=0.8, label="c_N = [δ_N]₊ / B(K_red)⁵"); ax.plot(N, GB, color=ORANGE, lw=0.8, label="Γ_N / B(K_red)⁵ (signed)"); ax.axhline(0, color=INK, lw=0.8); ax.set_xscale("log"); ax.set_xlabel("N"); ax.set_ylabel("normalized"); ax.set_title("Arithmetic view 4: the original normalized allowance; dotted lines = sixth-power thresholds of K_red", loc="left"); ax.legend(frameon=False); mark_thresholds(ax, thr); save(fig, out, "arith_4_normalized.png")
    # 5. block decomposition totals
    if os.path.exists(os.path.join(d, "block_summary.json")):
        bs = json.load(open(os.path.join(d, "block_summary.json")))
        fig, ax = plt.subplots(figsize=(8, 4)); style(ax); x = np.arange(len(bs)); w = 0.2
        ax.bar(x - 1.5*w, [b['E_plus_total'] for b in bs], w, color=BLUE, label="E⁺ total"); ax.bar(x - 0.5*w, [-b['E_minus_total'] for b in bs], w, color=ORANGE, label="−E⁻ total")
        ax.bar(x + 0.5*w, [b['signed_total'] for b in bs], w, color=INK, label="E⁺ − E⁻ = C_N"); ax.bar(x + 1.5*w, [b['Gamma'] for b in bs], w, color="#6a3d9a", label="Γ_N = C_N + κ_N")
        ax.set_xticks(x); ax.set_xticklabels([f"{b['N']}\n{'complete' if b['complete'] else 'partial'}" for b in bs], fontsize=8); ax.axhline(0, color=INK, lw=0.8); ax.set_ylabel("energy"); ax.set_title("Arithmetic view 5: masked-block spectral split of the complement (ordered block pairs, all g)", loc="left"); ax.legend(frameon=False, fontsize=8); save(fig, out, "arith_5_blocks.png")
        # block ratio diagnostic
        br = read_csv(os.path.join(d, "block_spectrum.csv"))
        if br:
            fig, ax = plt.subplots(figsize=(8, 3.8)); style(ax)
            for Nv, col in zip(sorted(set(int(float(r['N'])) for r in br))[:4], (BLUE, ORANGE, "#6a3d9a", INK)):
                sel = [r for r in br if int(float(r['N'])) == Nv]; ax.scatter([float(r['g'])*max(float(r['U']), float(r['V'])) for r in sel], [float(r['normalized_block']) for r in sel], s=8, color=col, label=f"N={Nv}")
            ax.axhline(0, color=INK, lw=0.8); ax.set_xscale("log"); ax.set_xlabel("g·max(U,V)"); ax.set_ylabel("g·max(U,V)/√(UV) · S_{g,U,V}"); ax.set_title("Optional diagnostic: signed block ratio for the stronger dyadic square-root target (not the aggregate criterion)", loc="left", fontsize=9); ax.legend(frameon=False, fontsize=8); save(fig, out, "arith_5b_block_ratio.png")
    # 6. zooms around the largest positive excursions
    fig, axs = plt.subplots(1, 2, figsize=(10, 3.8))
    for ax, y, name in ((axs[0], C, "C_N"), (axs[1], G, "Γ_N")):
        style(ax); k = int(np.argmax(y)); lo, hi = max(0, k - 60), min(len(N), k + 60); ax.plot(N[lo:hi], y[lo:hi], color=BLUE, lw=1.0, marker=".", ms=3); ax.axhline(0, color=INK, lw=0.8)
        ax.set_title(f"View 6: {name} around its maximum N = {int(N[k])} ({y[k]:+.5f})", loc="left", fontsize=9); ax.set_xlabel("N")
        for t in thr:
            if N[lo] <= t <= N[hi - 1]: ax.axvline(t, color=ORANGE, lw=1.0, ls="--", label=f"K_red threshold {t}"); ax.legend(frameon=False, fontsize=8)
    save(fig, out, "arith_6_zooms.png")

def primezero_plots(res, out):
    d = os.path.join(res, "primezero")
    if not os.path.exists(os.path.join(d, "primezero.csv")): return
    rows = read_csv(os.path.join(d, "primezero.csv")); meta = json.load(open(os.path.join(d, "primezero_meta.json")))
    Xs = sorted(set(int(float(r['X'])) for r in rows)); cnts = sorted(set(int(float(r['positive_zero_count'])) for r in rows))
    # 1. windows
    fig, axs = plt.subplots(2, 2, figsize=(11, 6.5), gridspec_kw=dict(height_ratios=[2, 1]))
    for j, (X, smp) in enumerate(meta['samples'].items()):
        if j > 1: break
        x = np.array(smp['x']); e = np.array(smp['e']); ez = np.array(smp['eZ'])
        ax = axs[0, j]; style(ax); ax.plot(x, e, color=INK, lw=0.7, label="e(x) = ψ(x) − x (prime powers, midpoint at jumps)"); ax.plot(x, ez, color=ORANGE, lw=0.9, label=f"e_Z(x), {smp['zeros']} listed zeros + conjugates"); ax.set_title(f"Prime–zero view 1: window X = {X}, a = 6", loc="left"); ax.legend(frameon=False, fontsize=8)
        ax = axs[1, j]; style(ax); ax.plot(x, e - ez, color=BLUE, lw=0.7); ax.axhline(0, color=INK, lw=0.6); ax.set_xlabel("x"); ax.set_ylabel("residual r_Z")
    save(fig, out, "primezero_1_windows.png")
    # 2. heatmap of D and D/P
    fig, axs = plt.subplots(1, 2, figsize=(11, 4))
    for ax, key, title in ((axs[0], "D_residual", "D(X; Z) = ‖e − e_Z‖²_X"), (axs[1], "D_over_P", "D/P")):
        M = np.full((len(cnts), len(Xs)), np.nan)
        for r in rows: M[cnts.index(int(float(r['positive_zero_count']))), Xs.index(int(float(r['X'])))] = float(r[key])
        im = ax.imshow(np.log10(M), aspect="auto", origin="lower", cmap="viridis"); ax.set_xticks(range(len(Xs))); ax.set_xticklabels(Xs, fontsize=8); ax.set_yticks(range(len(cnts))); ax.set_yticklabels([f"{c} (T≤{float([r for r in rows if int(float(r['positive_zero_count']))==c][0]['actual_max_height']):.0f})" for c in cnts], fontsize=8)
        ax.set_xlabel("X"); ax.set_ylabel("positive zeros listed (max height)"); ax.set_title(f"Prime–zero view 2: log10 {title}", loc="left"); fig.colorbar(im, ax=ax, shrink=0.8)
    save(fig, out, "primezero_2_heatmap.png")
    # 3. P, Z, I, D vs X for the largest zero list
    cmax = cnts[-1]; sel = [r for r in rows if int(float(r['positive_zero_count'])) == cmax]
    fig, ax = plt.subplots(figsize=(8, 4)); style(ax); xs = [float(r['X']) for r in sel]
    for key, lab, col, mk in (("P_prime", "P (prime discrepancy energy)", INK, "o"), ("Z_reconstruction", "Z (full reconstruction energy)", ORANGE, "s"), ("D_residual", "D (residual)", BLUE, "^"), ("I_cross", "I (signed cross term)", "#6a3d9a", "v")):
        ax.plot(xs, [float(r[key]) for r in sel], marker=mk, ms=4, color=col, lw=1.0, label=lab)
    ax.axhline(0, color=INK, lw=0.6); ax.set_xscale("log"); ax.set_xlabel("X"); ax.set_ylabel("window energy (X⁻² ∫_X^{6X})"); ax.set_title(f"Prime–zero view 3: P = Z + I + D with {cmax} listed zeros", loc="left"); ax.legend(frameon=False, fontsize=8); save(fig, out, "primezero_3_PZID.png")
    # 4. Gram bound
    fig, axs = plt.subplots(1, 2, figsize=(11, 4))
    ax = axs[0]; style(ax)
    for c, col in zip(cnts, plt.cm.viridis(np.linspace(0.1, 0.9, len(cnts)))):
        sl = [r for r in rows if int(float(r['positive_zero_count'])) == c]; ax.plot([float(r['X']) for r in sl], [float(r['bound_utilization']) for r in sl], marker="o", ms=3, color=col, lw=1.0, label=f"{c} zeros, λ={float(sl[0]['lambda_finite']):.1f}")
    ax.set_xscale("log"); ax.set_xlabel("X"); ax.set_ylabel("E_Z / (λ_Z W_Z)"); ax.set_title("Prime–zero view 4: utilization of the finite mode bound", loc="left"); ax.legend(frameon=False, fontsize=7)
    ax = axs[1]; style(ax); ax.plot(xs, [float(r['E_modes']) for r in sel], marker="o", ms=3, color=INK, label="E_Z = c*Gc"); ax.plot(xs, [float(r['E_diag']) for r in sel], marker="s", ms=3, color=BLUE, label="diagonal part"); ax.plot(xs, [float(r['E_offdiag']) for r in sel], marker="^", ms=3, color=ORANGE, label="signed off-diagonal part"); ax.plot(xs, [float(r['lambda_finite'])*float(r['W_weighted_diagonal']) for r in sel], ls="--", color=MUTED, label="λ_Z W_Z (bound)")
    ax.axhline(0, color=INK, lw=0.6); ax.set_xscale("log"); ax.set_yscale("symlog", linthresh=1e-3); ax.set_xlabel("X"); ax.set_title(f"Mode energy and bound, {cmax} zeros", loc="left"); ax.legend(frameon=False, fontsize=8); save(fig, out, "primezero_4_gram.png")
    # 5. convergence
    fig, ax = plt.subplots(figsize=(8, 3.6)); style(ax)
    ax.plot(xs, [float(r['quadrature_Z_change']) for r in sel], marker="o", ms=3, color=BLUE, label="|ΔZ| when nodes doubled and pieces halved"); ax.plot(xs, [float(r['quadrature_D_change']) for r in sel], marker="s", ms=3, color=ORANGE, label="|ΔD|"); ax.plot(xs, [abs(float(r['identity_residual'])) for r in sel], marker="^", ms=3, color=INK, label="|P − (Z+I+D)|"); ax.plot(xs, [float(r['D_residual']) for r in sel], ls="--", color=MUTED, label="D itself (reconstruction residual, for scale)")
    ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlabel("X"); ax.set_title("Prime–zero view 5: quadrature convergence, separate from the reconstruction residual", loc="left"); ax.legend(frameon=False, fontsize=8); save(fig, out, "primezero_5_convergence.png")

def quantum_plots(res, out):
    d = os.path.join(res, "quantum")
    if not os.path.exists(os.path.join(d, "cases.csv")): return
    rows = read_csv(os.path.join(d, "cases.csv")); extra = json.load(open(os.path.join(d, "extra.json"))) if os.path.exists(os.path.join(d, "extra.json")) else {}
    srcs = ["mobius", "squarefree_positive", "random_sign", "all_positive"]; rs = sorted(set(int(float(r['r'])) for r in rows)); Ns = sorted(set(int(float(r['N'])) for r in rows))
    cols = {"mobius": BLUE, "squarefree_positive": ORANGE, "random_sign": "#6a3d9a", "all_positive": MUTED}
    # 1. P0 vs N per r (mobius and controls)
    fig, axs = plt.subplots(1, 4, figsize=(15, 3.8), sharey=False)
    for ax, src in zip(axs, srcs):
        style(ax)
        for r, col in zip(rs, plt.cm.viridis(np.linspace(0.1, 0.9, len(rs)))):
            sel = [x for x in rows if x['source_name'] == src and int(float(x['r'])) == r]; xs = [float(x['N']) for x in sel]
            ax.plot(xs, [float(x['P_classical']) for x in sel], color=col, lw=1.0, label=f"r={r} classical"); ax.plot(xs, [float(x['P_statevector']) for x in sel], "x", color=col, ms=5)
            ax.errorbar(xs, [float(x['P_sampled']) for x in sel], yerr=[[float(x['P_sampled']) - float(x['P_ci_low']) for x in sel], [float(x['P_ci_high']) - float(x['P_sampled']) for x in sel]], fmt="o", ms=3, color=col, capsize=2)
        ax.set_xscale("log", base=2); ax.set_yscale("log"); ax.set_xlabel("N"); ax.set_title(f"P₀: {src}", loc="left", fontsize=9); ax.legend(frameon=False, fontsize=6.5)
    axs[0].set_ylabel("success probability"); fig.suptitle("Coherent Green-energy return: P₀ = R/(dQ); lines classical, × statevector, dots sampled with Wilson 95 %", x=0.01, ha="left", fontsize=10); save(fig, out, "quantum_1_P0.png")
    # 2. d P0 = R/Q
    fig, ax = plt.subplots(figsize=(8, 3.8)); style(ax)
    for src in srcs:
        for r, ls in zip(rs, ("-", "--", ":", "-.")):
            sel = [x for x in rows if x['source_name'] == src and int(float(x['r'])) == r]; ax.plot([float(x['N']) for x in sel], [float(x['dP0_classical']) for x in sel], ls=ls, color=cols[src], lw=1.0, label=f"{src} r={r}" if r == rs[0] else None)
    ax.set_xscale("log", base=2); ax.set_xlabel("N"); ax.set_ylabel("d·P₀ = R/Q"); ax.set_title("Finite residue precision: coherent-to-residue energy ratio (line styles: r = 2, 3, 4, 5)", loc="left"); ax.legend(frameon=False, fontsize=8); save(fig, out, "quantum_2_dP0.png")
    # 3. reconstructed R vs exact
    fig, ax = plt.subplots(figsize=(7, 5)); style(ax)
    for src in srcs:
        sel = [x for x in rows if x['source_name'] == src]; ax.errorbar([float(x['R_exact_float']) for x in sel], [float(x['R_reconstructed']) for x in sel], yerr=[[float(x['R_reconstructed']) - float(x['R_ci_low']) for x in sel], [float(x['R_ci_high']) - float(x['R_reconstructed']) for x in sel]], fmt="o", ms=3, color=cols[src], capsize=2, label=src)
    lim = [0.5, max(float(x['R_exact_float']) for x in rows)*1.2]; ax.plot(lim, lim, color=INK, lw=0.8); ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlabel("R exact"); ax.set_ylabel("R̂ = d Q P̂₀"); ax.set_title("Coherent Green-energy return: reconstructed energy against the exact value", loc="left"); ax.legend(frameon=False, fontsize=8); save(fig, out, "quantum_3_Rhat.png")
    # 4. R - Q and R/D
    fig, axs = plt.subplots(1, 2, figsize=(11, 3.8))
    for ax, key, lab in ((axs[0], "cross_channel_float", "R − Q (cross-residue contribution)"), (axs[1], "R_over_D", "R / D (energy over the fully separated diagonal)")):
        style(ax)
        for src in srcs:
            for r, ls in zip(rs, ("-", "--", ":", "-.")):
                sel = [x for x in rows if x['source_name'] == src and int(float(x['r'])) == r]; ax.plot([float(x['N']) for x in sel], [float(x[key]) for x in sel], ls=ls, color=cols[src], lw=1.0, label=f"{src}" if r == rs[0] else None)
        ax.axhline(0, color=INK, lw=0.6); ax.set_xscale("log", base=2); ax.set_yscale("symlog", linthresh=0.1); ax.set_xlabel("N"); ax.set_title(lab, loc="left", fontsize=9); ax.legend(frameon=False, fontsize=7)
    save(fig, out, "quantum_4_cross_and_ratio.png")
    # 5. Q vs r with saturation
    fig, ax = plt.subplots(figsize=(8, 3.8)); style(ax)
    for src in srcs:
        for N, mk in zip(Ns, ("o", "s", "^", "v", "D")):
            sel = [x for x in rows if x['source_name'] == src and int(float(x['N'])) == N]; ax.plot([float(x['r']) for x in sel], [float(x['Q_exact_float']) for x in sel], marker=mk, ms=4, color=cols[src], lw=0.8, label=f"{src} N={N}" if src == "mobius" else None)
            sat = [float(x['r']) for x in sel if 2**int(float(x['r'])) >= N]
            if sat: ax.scatter(sat, [float(x['Q_exact_float']) for x in sel if 2**int(float(x['r'])) >= N], s=60, facecolors="none", edgecolors=INK)
    ax.set_yscale("log"); ax.set_xlabel("residue precision r"); ax.set_ylabel("Q_{N,r}"); ax.set_title("Finite residue precision: Q versus r; circled = saturation d ≥ N where Q = D", loc="left"); ax.legend(frameon=False, fontsize=7); save(fig, out, "quantum_5_Q_vs_r.png")
    # 6. depth and 2q gates
    fig, axs = plt.subplots(1, 2, figsize=(11, 3.8))
    for ax, key, lab in ((axs[0], "prep_depth", "state-preparation depth (after decomposition)"), (axs[1], "two_qubit_gates", "two-qubit (CX) gates, full readout circuit")):
        style(ax); sel = [x for x in rows if x['source_name'] == "mobius"]; ax.scatter([float(x['data_qubits']) for x in sel], [float(x[key]) for x in sel], s=14, color=BLUE)
        q = np.array([float(x['data_qubits']) for x in sel]); ax.plot(sorted(set(q)), [2**k for k in sorted(set(q))], ls="--", color=MUTED, label="2^qubits (reference)"); ax.set_yscale("log"); ax.set_xlabel("data qubits r + s"); ax.set_title(lab, loc="left", fontsize=9); ax.legend(frameon=False, fontsize=8)
    save(fig, out, "quantum_6_resources.png")
    # 7. ancilla Weyl
    if extra.get("weyl"):
        fig, ax = plt.subplots(figsize=(7, 3.8)); style(ax); w = extra['weyl']; x = np.arange(len(w))
        ax.bar(x - 0.2, [e['expected'] for e in w], 0.4, color=MUTED, label="exact prediction cos θ / sin θ"); ax.errorbar(x + 0.2, [e['sampled'] for e in w], yerr=[e['sampled_2sigma'] for e in w], fmt="o", color=BLUE, capsize=3, label="sampled (±2σ)")
        ax.set_xticks(x); ax.set_xticklabels([f"r={e['r']} {e['which']}\nθ={e['theta']:.3f}" for e in w], fontsize=8); ax.axhline(0, color=INK, lw=0.6); ax.set_title("Ancilla Hadamard test of W = V_t T_h V_t† T_h† = e^{iθ} I (t = h = 1)", loc="left"); ax.legend(frameon=False, fontsize=8); save(fig, out, "quantum_7_weyl.png")
    # 8. characters
    if extra.get("characters"):
        ch = {k: v for k, v in extra['characters'].items() if isinstance(v, dict) and 'P_t' in v}; keys = list(ch)[:8]
        fig, axs = plt.subplots(2, 4, figsize=(14, 5.5)); axs = axs.ravel()
        for ax, k in zip(axs, keys):
            style(ax); pt = ch[k]['P_t']; ax.bar(range(len(pt)), pt, color=BLUE, label="QFT measurement (statevector)"); ax.plot(range(len(pt)), ch[k]['formula'], "x", color=ORANGE, ms=6, label="formula R_{a,t}/(dQ)"); ax.set_title(k, fontsize=8, loc="left"); ax.set_xlabel("character t")
        axs[0].legend(frameon=False, fontsize=7); fig.suptitle("Residue-character spectrum (cyclic QFT), not the zeta zeros", x=0.01, ha="left", fontsize=10); save(fig, out, "quantum_8_characters.png")

def make_all(res):
    out = os.path.join(res, "plots"); os.makedirs(out, exist_ok=True)
    arithmetic_plots(res, out); primezero_plots(res, out); quantum_plots(res, out)
