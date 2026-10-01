"""
The constants of the edge from the Wiener-Hopf representation (paper Proposition 8.4, Computation 8.6).

For the prime-free form, Corollary 7.24 gives the edge of the transform as A_0/(t sigma_-(t)), and the inverse transform of
that factor is the edge law with its constants fixed:
    f(a - delta) = C (log(a/delta) + beta)^{-1/2} (1 + O(log^{-3})),   C = |A_0|,   beta = -gamma - log(2 pi a) - lambda,
i.e. f(a - delta) = C (sigma_inf(1/delta) - gamma)^{-1/2}: the inverse square root of the symbol at frequency 1/delta, shifted by
Euler's constant.  With the boundary law (Proposition 8.3) C^2 = -lambda'/2.  This script tests both constants, with nothing
fitted, on the prime-free edge-FEM critical points: lambda from the FEM, lambda' from the arb dilation engine (data/dilation_*.log).

Usage: python3 rh_edge_constants.py
"""
import json, math, re, numpy as np
gamma = 0.5772156649015329
cases = [(0.3, "data/fem_a0.3_primefree.json"), (0.4, "data/fem_a0.4_primefree.json"), (0.5, "data/fem_a0.5_primefree.json"),
         (0.5, "data/fem_a0.5_primefree_471nodes.json"), (0.6, "data/fem_a0.6_primefree.json")]
for a, fn in cases:
    d = json.load(open(fn)); lam = d["lambda"]; dl = np.array(d["deltas"]); f = np.array(d["f"])
    txt = open(f"data/dilation_a{a}_primefree.log").read(); lp = float(re.search(r"a lambda' = \[([-0-9.e+]+)", txt).group(1))/a
    C = math.sqrt(-lp/2); beta = -gamma - math.log(2*math.pi*a) - lam
    print(f"{fn}: a = {a}, lambda = {lam:+.6f}, lambda' = {lp:+.5f};  predicted C = (-lambda'/2)^(1/2) = {C:.4f}, beta = -gamma - log(2 pi a) - lambda = {beta:+.4f}")
    o = np.argsort(dl); dl, f = dl[o], f[o]
    print("      delta      L = log(a/delta)   f(a-delta)   f / [C (L + beta)^(-1/2)]")
    for dd in (1e-12, 1e-10, 1e-8, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2):
        i = np.argmin(np.abs(np.log(dl[1:]) - math.log(dd))) + 1
        L = math.log(a/dl[i]); r = f[i]/(C*(L + beta)**-0.5)
        print(f"      {dl[i]:.2e}   {L:7.3f}           {f[i]:.5f}      {r:.4f}")
