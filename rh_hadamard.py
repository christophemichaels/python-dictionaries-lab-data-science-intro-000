"""
The boundary law lambda'(a) = -2 C(a)^2 (memo section 2.9; paper section 8).

C(a) is the edge amplitude of the minimizer, f(a - delta) = C (log(a/delta) + beta)^(-1/2), fitted on the output
of rh_edge_fem.py (1/f^2 linear in log(a/delta) on 14 < t < 20); lambda'(a) is the exact Hellmann-Feynman
derivative from rh_dilation.py.  This script reproduces the kappa table from those outputs.

Usage: python3 rh_hadamard.py fem1.json dil1.log [fem2.json dil2.log ...]
"""
import sys, json, math, re
import numpy as np

def edge_fit(fn, lo=14, hi=20):
    J = json.load(open(fn)); a = J["a"]; d = np.array(J["deltas"]); f = np.array(J["f"])
    m = d > 0; t = np.log(a/d[m]); ff = f[m]; w = (t > lo) & (t < hi)
    sl, ic = np.polyfit(t[w], 1/ff[w]**2, 1)
    return J["a"], J["primes"], J["lambda"], 1/math.sqrt(sl), ic/sl

def exact_derivative(fn):
    t = open(fn).read()
    a = float(re.search(r"^a = \[?([0-9.]+)", t, re.M).group(1))
    return float(re.search(r"sum = a lambda' = \[([0-9.e+-]+)", t).group(1))/a

if __name__ == "__main__":
    print("  a     primes  lambda(FEM)      C        beta    lambda'(exact)   kappa = -lambda'/(2 C^2)")
    for fj, dl in zip(sys.argv[1::2], sys.argv[2::2]):
        a, p, lam, C, beta = edge_fit(fj); lp = exact_derivative(dl)
        print(f"{a:6.3f}   {p}     {lam:+.6e}  {C:.5f}  {beta:+.3f}   {lp:+.6e}   {-lp/(2*C*C):.4f}")
