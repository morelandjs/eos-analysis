#!/usr/bin/env python
import numpy as np
import scipy 
import fortranformat as ff
from fortranformat import FortranRecordWriter
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy import integrate
from scipy.integrate import quad

# constants
hbarc = 0.19733

# temperature limits in GeV
T0 = 0.2
Tmin = 0.001
Tmax = 1.0
dT = 0.005
nT = int((Tmax-Tmin)/dT)
Tvec = np.linspace(Tmin,Tmax,nT)
t=Tvec/T0

# eos parameters
h0 = 0.1396
h1 = -0.1800
h2 = 0.0350 
f0 = 1.05 
f1 = 6.39 
f2 = -4.72 
g1 = -0.92 
g2 = 0.57

# trace anomaly: I/T^4
tr = np.exp(-h1/t-h2/t**2.)*(h0 + f0*(np.tanh(f1*t+f2)+1.)/(1.+g1*t+g2*t**2.)) 
tr = np.asarray(tr)

# pressure: p/T^4
p = []
for T in Tvec:
    p.append(integrate.quad(lambda x: np.exp(-h1/x-h2/x**2.)*(h0 + f0*(np.tanh(f1*x+f2)+1.)/(1.+g1*x+g2*x**2.))/x, 0.0, T/T0)[0])
p = np.asarray(p)
fp = interp1d(Tvec, p, kind='cubic')

# energy density: e/T^4
e = []
for iT in range(nT):
    e.append(tr[iT] + 3*p[iT])
e = np.asarray(e)
fe = interp1d(Tvec, e, kind='cubic')

# entropy density: s/T^3
s = []
for iT in range(nT):
    T = Tmin + iT*dT;
    s.append(e[iT]+p[iT])
s = np.asarray(s)
fs = interp1d(Tvec, s, kind='cubic')

# convert energy density to temperature
def e2T(e):
    counts = 0
    error = 1.0
    Tlb = Tmin
    Tub = Tmax
    Tmid = (Tub+Tlb)/2.

    while error > 1e-8: 
        e0 = fe(Tmid)*Tmid**4./hbarc**3;
        if counts > 50:
            print "error excessive bisections: ",counts
        if e0 > e:
            Tub = Tmid
        if e0 <= e:
            Tlb = Tmid
        Tmid = (Tub+Tlb)/2.
        error = np.abs((e0-e)/e)
        counts += 1
    return Tmid    

# e output mesh: GeV/fm^3
emin = 0.1e-2
emax = 0.310999e3
ne = 155500
de = (emax-emin)/float(ne)
evec = np.linspace(emin,emax,ne)

# T output mesh: GeV
Tinv = []
printed = False
for en in evec:
    Tinv.append(e2T(en))
    prog = int(en/emax*100)
    if prog % 5 == 0:
        if printed == False:
            print prog,"%"
        printed = True
    if prog % 5 != 0:
        printed = False
Tinv = np.asarray(Tinv)

# plot curves
plt.plot(Tvec,e,'o')
plt.plot(Tinv,fe(Tinv),'-')
plt.plot(Tvec,s,'o')
plt.plot(Tinv,fs(Tinv),'-')

# output arrays
eEOS = fe(Tinv)*Tinv**4./(hbarc**3)
pEOS = fp(Tinv)*Tinv**4./(hbarc**3)
sEOS = fs(Tinv)*Tinv**3./(hbarc**3)
TEOS = Tinv

line = FortranRecordWriter('(4E15.6)')

# open output file
with open('WB-EOS.dat','w') as wf:
    for i in range(len(eEOS)):
        wf.write(line.write([eEOS[i],pEOS[i],sEOS[i],TEOS[i]])+"\n")

#plt.legend(loc='upper right')
#plt.xlim(0.0, 0.8)
#plt.ylim(0.0, 5.0)
#plt.xlabel('$T$ MeV')
#plt.ylabel('$(e-3p)/T^4$')
#plt.title('Trace Anomaly')
#plt.savefig("trace.pdf")
plt.show()
