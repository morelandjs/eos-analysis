#!/usr/bin/env python
import numpy as np
import scipy 
import fortranformat as ff
from fortranformat import FortranRecordWriter
import matplotlib.pyplot as plt
from scipy import interpolate
from scipy.interpolate import interp1d
from scipy import integrate
from scipy.integrate import quad

# constants
hbarc = 0.19733

# temperature limits in GeV
Tmin = 0.05
Tmax = 0.7
dT = 0.001
nT = int((Tmax-Tmin)/dT)
T = np.linspace(Tmin,Tmax,nT)

# eos parameters
Tc = 0.154
ct = 3.8706
an = -8.7704
bn = 3.9200
cn = 0.0
dn = 0.3419
t0 = 0.9761
ad = -1.2600
bd = 0.8425
cd = 0.0
dd = -0.0475
tb = T/Tc
pid = 47.5*np.pi**2./90.

# pressure over T^4
p = 0.5*(1.+np.tanh(ct*(tb-t0)))*(pid+an/tb+bn/tb**2.+cn/tb**3.+dn/tb**4.)/(1.+ad/tb+bd/tb**2.+cd/tb**3.+dd/tb**4.)
dfp = interp1d(T, interpolate.splev(T,interpolate.splrep(T,p,s=0),der=1), kind='cubic')
fp = interp1d(T, p, kind='cubic')

# trace anomaly over T^4
I = []
for iT in T:
    I.append(dfp(iT)*iT)
fI = interp1d(T, I, kind='cubic')

# energy density over T^4 
e = []
for iT in range(nT):
    e.append(I[iT] + 3.*p[iT])
fe = interp1d(T, e, kind='cubic')

# entropy density over T^3
s = []
for iT in range(nT):
    s.append(e[iT] + p[iT])
fs = interp1d(T, s, kind='cubic')

# convert e in GeV to temperature in GeV
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

# output arrays
eEOS = fe(Tinv)*Tinv**4./(hbarc**3.)
pEOS = fp(Tinv)*Tinv**4./(hbarc**3.)
sEOS = fs(Tinv)*Tinv**3./(hbarc**3.)
TEOS = Tinv

line = FortranRecordWriter('(4E15.6)')

# open output file
with open('HotQCD-EOS.dat','w') as wf:
    for i in range(len(eEOS)):
        wf.write(line.write([eEOS[i],pEOS[i],sEOS[i],TEOS[i]])+"\n")

# plot curves
plt.plot(Tinv,eEOS,"-")
plt.plot(Tinv,sEOS,"-")
plt.plot(Tinv,pEOS,"-")

# plot properties
#plt.legend(loc='upper right')
#plt.xlim(0.0, 0.8)
#plt.ylim(0.0, 5.0)
#plt.xlabel('$T$ MeV')
#plt.ylabel('$(e-3p)/T^4$')
#plt.xlim(0,800)
#plt.ylim(0,20)
plt.show()
