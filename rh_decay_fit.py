# Fit decay models for the odd Weil floor lambda(a) (values from "Primes, Folds, and the One Dot", Table 8, and Zhu a=1.2).
# Requires numpy. Run: python3 rh_decay_fit.py
import numpy as np
a  = np.array([0.4, 0.5, 0.6, 0.8, 1.0, 1.1, 1.2])
lam= np.array([1.47e-2, 1.94e-4, 5.96e-7, 1.56e-14, 1.49e-26, 7.68e-35, 4.76e-45])
T  = 2*np.pi*np.exp(2*a)
y  = -np.log(lam)
print(" a     T*      -log lam   ratio=-loglam/T*   ratio/a   ratio/log(T*)")
for ai,Ti,yi in zip(a,T,y):
    print(f"{ai:.1f} {Ti:7.2f} {yi:10.3f} {yi/Ti:10.4f} {yi/Ti/ai:12.4f} {yi/Ti/np.log(Ti):10.4f}")
# fits on a>=0.6 (4 dof-ish)
m = a>=0.6
def fit(X, name):
    A = np.vstack([X[m], np.ones(m.sum())]).T
    c,res,_,_ = np.linalg.lstsq(A, y[m], rcond=None)
    pred = A@c
    rms = np.sqrt(np.mean((pred-y[m])**2))
    print(f"{name:22s} coef={c[0]:.4f} const={c[1]:8.3f} rms={rms:.3f}")
    return c
fit(T, "-log lam ~ c*T*")
fit(T*a, "-log lam ~ c*a*T*")
fit(T*np.log(T), "-log lam ~ c*T*logT*")
# fit excluding 1.2, predict 1.2
m2 = (a>=0.6)&(a<1.2)
for X,name in [(T,"c*T*"),(T*a,"c*a*T*"),(T*np.log(T),"c*T*logT*")]:
    A=np.vstack([X[m2],np.ones(m2.sum())]).T; c=np.linalg.lstsq(A,y[m2],rcond=None)[0]
    p=c[0]*X[-1]+c[1]; print(f"predict a=1.2 from {name:10s}: lam={np.exp(-p):.2e} (Zhu 4.76e-45)")
# local slopes
for i in range(len(a)-1):
    print(f"d(-loglam)/dT* on [{a[i]},{a[i+1]}]: {(y[i+1]-y[i])/(T[i+1]-T[i]):.3f}")
# calibration: a needed to reach Platt-Trudgian height
H=3e12; aH=0.5*np.log(H/(2*np.pi)); print("a with T*=3e12:", aH, " digits of 1/lam at 1.8*T*:", 1.8*H/np.log(10))
