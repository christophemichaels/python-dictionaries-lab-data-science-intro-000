"""Light units for the cone and the horn (LIGHT_UNITS.md). A choice of units and nothing else:
one e-fold of the prime line = lambda0 metres; the cone's edge log n = 2a moves at c, so one unit of half-width a is tau0 = 2 lambda0 / c.
Zeros become frequencies f = c * (gamma / 2 pi) / lambda0, since gamma / 2 pi is the wavenumber of the zero's wave per e-fold.
Usage: python3 rh_light_units.py data/zeros_6000.txt [lambda0_m]"""
import sys, math
c = 299_792_458.0
g = [float(l.split()[1]) for l in open(sys.argv[1])]
lam = float(sys.argv[2]) if len(sys.argv) > 2 else 1.0
tau = 2 * lam / c
f = lambda T: c * (T / (2 * math.pi)) / lam
print(f"units: one e-fold of the prime line = {lam} m; edge speed c; one unit of half-width a = tau0 = {tau*1e9:.4f} ns")
print(f"the subpower exponent 1/2 is dimensionless and stays 1/2 in every choice of units")
print("\nentries of the prime powers (a_n = log n / 2):")
for n in (2, 3, 4, 5, 7, 8, 9, 11, 13):
    a = math.log(n) / 2; print(f"  {n:3d}: a = {a:.4f}  ->  t = {a*tau*1e9:6.3f} ns;  horizon T* = 2 pi n = {2*math.pi*n:7.2f}  ->  f* = {f(2*math.pi*n)/1e6:8.1f} MHz")
print("\nsupports: Zhu's certified cone and the computed floor")
for a, label in ((0.3466, "entry of 2"), (0.8, "Zhu certified"), (1.2, "floor 7e-37"), (1.505, "floor ~1e-93")):
    T = 2 * math.pi * math.exp(2 * a); print(f"  a = {a:6.4f} ({label:14s}): t = {a*tau*1e9:6.3f} ns; horizon T* = {T:8.2f} -> f* = {f(T)/1e9:7.3f} GHz")
print("\nzeros as frequencies (wave of the zero per e-fold, travelling at c):")
for n in (1, 2, 3, 10, 100, 1000, 6000):
    gm = g[n-1]; print(f"  zero {n:5d}: gamma = {gm:9.3f}; wavelength 2 pi/gamma = {2*math.pi/gm*lam*100:8.3f} cm; f = {f(gm)/1e9:8.3f} GHz")
print("\nstones at ring radius m -> frequency c m / lambda0:")
for m in (2, 3, 5, 7, 11, 97, 997):
    print(f"  stone {m:4d}: switches on at {math.log(m)/2*tau*1e9:6.3f} ns, resolves the spectrum to {c*m/lam/1e9:8.3f} GHz")
