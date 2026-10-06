import mpmath as mp
from rh_weil_odd import OddWeil
W = OddWeil(K=24)
# Table 1 (paper): n=2 entry 0.346574 fails 0.371601 ; n=3 entry 0.549306 fails 0.557323
for n, grid in [(2, [0.0, 0.010, 0.020, 0.025, 0.030]), (3, [0.0, 0.003, 0.006, 0.008, 0.010])]:
    an = mp.log(n)/2
    print(f"n={n}: deleted form lambda_del(a_n + d)")
    rows = []
    for d in grid:
        lam, fa, _ = W.floor(an + d, delete=n)
        rows.append((d, lam))
        print(f"   d={d:.3f}  lambda_del={mp.nstr(lam,6)}")
    # linear interpolation for crossing
    for (d0,l0),(d1,l1) in zip(rows, rows[1:]):
        if l0 > 0 and l1 < 0:
            dc = d0 + (d1-d0)*l0/(l0-l1)
            print(f"   crossing near d = {mp.nstr(dc,4)}  ->  a = {mp.nstr(an+dc,6)}   (paper: {0.371601 if n==2 else 0.557323})")
    slope = (rows[1][1]-rows[0][1])/(rows[1][0]-rows[0][0])
    print(f"   lambda_del(a_n)/|slope| = {mp.nstr(rows[0][1]/abs(slope),4)}   vs paper's gap {0.0250 if n==2 else 0.0080}")
