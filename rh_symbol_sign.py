import numpy as np
def vonmangoldt(n):
    for p in range(2, n+1):
        if n % p == 0:
            m = n
            while m % p == 0: m //= p
            return np.log(p) if m == 1 else 0.0
    return 0.0
def redigamma(t):
    z = 0.25 + 0.5j*t
    # shift up for accuracy then asymptotic expansion
    s = 0
    for k in range(8):
        s += 1/z; z += 1
    return (np.log(z) - 1/(2*z) - 1/(12*z**2) + 1/(120*z**4) - 1/(252*z**6) - s).real
print(" a     T*     RMS(comb)  2a    frac Psi_a<0 on [0,T*]  on [T*,3T*]  on [3T*,10T*]   min Psi_a on [0,T*]")
for a in [0.6, 0.8, 1.0, 1.1, 1.5, 2.0]:
    T = 2*np.pi*np.exp(2*a)
    t = np.linspace(0.5, 10*T, 400000)
    comb = np.zeros_like(t)
    n = 2
    while np.log(n) < 2*a:
        L = vonmangoldt(n)
        if L: comb += 2*L/np.sqrt(n)*np.cos(t*np.log(n))
        n += 1
    psi = redigamma(t) - np.log(np.pi) - comb
    rms = np.sqrt(np.mean(comb[t<10*T]**2))
    f1 = np.mean(psi[t<T] < 0); f2 = np.mean(psi[(t>=T)&(t<3*T)] < 0); f3 = np.mean(psi[(t>=3*T)] < 0)
    print(f"{a:.1f} {T:7.1f} {rms:8.3f} {2*a:5.1f} {f1:12.3f} {f2:14.3f} {f3:14.3f} {psi[t<T].min():14.2f}")
