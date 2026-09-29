import mpmath as mp
from rh_weil_odd import OddWeil
W = OddWeil(K=24)
a = mp.mpf('0.8')
lam0, fa, _ = W.floor(a)
print(f"a=0.8 unperturbed: lambda = {mp.nstr(lam0,8)}")
for eps in ['3e-14', '-3e-14', '1e-13', '-1e-13']:
    lam, _, _ = W.floor(a, shift={2: eps})
    print(f"  log 2 -> log 2 + ({eps}):  lambda = {mp.nstr(lam,8)}   change/eps = {mp.nstr((lam-lam0)/mp.mpf(eps),5)}")
from rh_weil_odd import OddWeil

# Beurling-style: move the prime 3 (and 9) by a relative 1e-3; where does the floor fail?
for eps in ['1e-3', '-1e-3']:
    print(f"log 3 -> log 3 + ({eps})")
    for a in ['0.55', '0.60', '0.65', '0.70']:
        lam, _, _ = W.floor(mp.mpf(a), shift={3: eps})
        lam0, _, _ = W.floor(mp.mpf(a))
        print(f"   a={a}: lambda_pert = {mp.nstr(lam,6)}   (unperturbed {mp.nstr(lam0,6)})")
