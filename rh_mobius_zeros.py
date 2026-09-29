# Zero-side constants sum 2/|rho|^2 and sum 2/|rho zeta(rho)|^2 from the first 4000 zeros (Section 14 of mobius_modifier_v2.tex).
# Zero-side constants: Cramer's sum 1/|rho|^2 (psi mean square) and Ng's sum 1/|rho zeta'(rho)|^2 (M mean square).
import mpmath as mp, time, json
mp.mp.dps = 15
t0 = time.time(); S1 = 0; S2 = 0; rows = []; zs = []
K = 4000
for n in range(1, K+1):
    rho = mp.zetazero(n)
    dz = mp.zeta(rho, derivative=1)
    zs.append((float(mp.im(rho)), float(mp.re(dz)), float(mp.im(dz))))
    S1 += 2/abs(rho)**2
    S2 += 2/abs(rho*dz)**2
    if n in (1, 10, 100, 500, 1000, 2000, 3000, 4000):
        rows.append((n, float(mp.im(rho)), float(S1), float(S2)))
        print(n, float(mp.im(rho)), float(S1), float(S2), f"{time.time()-t0:.0f}s", flush=True)
gam = float(mp.im(rho))
# tails: sum over gamma>T of 2/gamma^2 ~ (1/pi) log(T/2pi)/T ; for S2 use typical |zeta'|^2 ~ (log gamma)^2 /? -- report raw partial sums only
print("Cramer constant 2+gamma-log(4 pi) =", float(2+mp.euler-mp.log(4*mp.pi)))
print("partial S1 =", float(S1), " tail estimate (1/pi) log(T/2pi)/T =", float(mp.log(gam/(2*mp.pi))/(mp.pi*gam)))
print("partial S2 =", float(S2))
import numpy as np; np.save("zeros_dz.npy", np.array(zs))
json.dump({"rows": rows, "S1": float(S1), "S2": float(S2), "T": gam}, open("ngconst.json","w"))
