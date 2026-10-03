"""
Zeros of zeta to high precision with arb (acb.zeta_zeros): python3 rh_zeros_arb.py n0 n1 bits out.txt
Lines "n gamma_n" with 0.3*bits digits.  The plunge region of the tail-law computation (zeros below ~2.5 T*) needs
|delta gamma| << sqrt(lambda_K): 600 bits for a = 2 (lambda ~ 1e-276), 220 bits above the horizon.
"""
import sys, time
from flint import acb, ctx
n0, n1, bits, out = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), sys.argv[4]
ctx.prec = bits; digits = int(bits*0.30103) - 2
t0 = time.time()
with open(out, "w") as fh:
    n = n0
    while n <= n1:
        cnt = min(50, n1 - n + 1)
        zs = acb.zeta_zeros(n, cnt)
        for k, z in enumerate(zs):
            fh.write(f"{n + k} {z.imag.mid().str(digits, radius=False)}\n")
        fh.flush(); n += cnt
        print(f"  zeros to {n-1} ({time.time()-t0:.0f}s), last radius {zs[-1].imag.rad().str(3)}", flush=True)
