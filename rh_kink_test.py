import mpmath as mp
from rh_weil_odd import OddWeil
W = OddWeil(K=24)
for n in [2, 3]:
    an = mp.log(n)/2
    h = mp.mpf('0.0005')
    lam = {}
    for k in [-2, -1, 0, 1, 2]:
        lam[k], fa, _ = W.floor(an + k*h)
        if k == 0: fa0 = fa
    dl_minus = (3*lam[0] - 4*lam[-1] + lam[-2])/(2*h)
    dl_plus  = (-3*lam[0] + 4*lam[1] - lam[2])/(2*h)
    pred = 4*mp.log(n)/mp.sqrt(n)*fa0**2
    print(f"n={n}  a_n={mp.nstr(an,8)}  lambda(a_n)={mp.nstr(lam[0],8)}  f(a_n)={mp.nstr(fa0,8)}")
    print(f"   lambda'(a_n^-) = {mp.nstr(dl_minus,8)}   lambda'(a_n^+) = {mp.nstr(dl_plus,8)}")
    print(f"   jump = {mp.nstr(dl_plus-dl_minus,8)}   predicted 4*Lambda(n)/sqrt(n)*f(a_n)^2 = {mp.nstr(pred,8)}   ratio={mp.nstr((dl_plus-dl_minus)/pred,6)}")
    print(f"   check monotone: |lambda'(a_n^-)| >= jump ? {abs(dl_minus) >= (dl_plus-dl_minus)}")
