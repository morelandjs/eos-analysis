#!/usr/bin/env python
import numpy as np
import scipy 
import fortranformat as ff
from fortranformat import FortranRecordWriter
import matplotlib.pyplot as plt
from scipy import interpolate
from scipy.interpolate import interp1d
from scipy.optimize import fminbound
from scipy import integrate
from scipy.integrate import quad

####################################################################################

# construct low temperature interaction measure 

#T,QS,PN,EN,PEN,SN,CV,SNCV = np.loadtxt("hrg-eos/OUT_2.5.DAT140_2014_bulk",dtype='float',unpack=True,skiprows=2)
#fIHRG = interp1d(T, QS, kind='cubic')

T, I, e, p, Cs = np.loadtxt("hrg-eos/hrg-urqmd-eos.dat",dtype='float',unpack=True,skiprows=1)
fI_lo = interp1d(T/1000., I, kind='cubic')

##################################################################################

# construct high temperature interaction measure

# temperature limits in GeV
Tmin = 0.065
Tmax = 0.800
dT = 0.001
nT = int((Tmax-Tmin)/dT)
T = np.linspace(Tmin,Tmax,nT)

# eos parameters
hbarc = 0.19733
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

# pressure over T**4
p = 0.5*(1.+np.tanh(ct*(tb-t0)))*(pid+an/tb+bn/tb**2.+cn/tb**3.+dn/tb**4.)/(1.+ad/tb+bd/tb**2.+cd/tb**3.+dd/tb**4.)
dfp = interp1d(T, interpolate.splev(T,interpolate.splrep(T,p,s=0),der=1), kind='cubic')
fp = interp1d(T, p, kind='cubic')

# trace anomaly over T**4
I = []
for iT in T:
    I.append(dfp(iT)*iT)
fI_hi = interp1d(T, I, kind='cubic')

###############################################################################

# blend UrQMD HRG interaction measure into HotQCD interaction measure

# temperature limits in GeV
Tmin = 0.005
Tmax = 0.800
nT = 10000
dT = (Tmax-Tmin)/float(nT)
T = np.linspace(Tmin,Tmax,nT)

I = np.zeros(nT)

for iT, T_ in enumerate(T):
    if T_ < 0.120:
        I[iT] = fI_lo(T_)
    elif T_ >= 0.120 and T_ < 0.180:
        w1 = lambda x: 0.005/(1+np.exp((x-0.176)/0.002))
        u1 = lambda x: (1+np.tanh((T_-0.150)/w1(T_)))/2
        I[iT] = fI_lo(T_) + u1(T_)*(fI_hi(T_)-fI_lo(T_))
    else:
        I[iT] = fI_hi(T_) 

fI = interp1d(T,I)   	         

# calculate p/T**4 from interaction measure
p = []
integral = 0.0
for iT,T_ in enumerate(T):
    p.append(integral)
    integral += I[iT]/T_*dT
fp = interp1d(T,p)
  
# calculate e/T**4 
e = []
for iT,T_ in enumerate(T):
    e.append(I[iT] + 3.*p[iT])
fe = interp1d(T,e)

# calculate s/T**3
s = []
for iT,T_ in enumerate(T):
    s.append(e[iT] + p[iT])
fs = interp1d(T,s)

# e output mesh: GeV/fm**3
evec = np.arange(1,650000,2)*1e-3
print(evec.size)

Tinv = []
for ie0, e0 in enumerate(evec):
    Tinv.append(fminbound(lambda x: abs(fe(x)*(x**4.)/(hbarc**3.)-e0),Tmin,Tmax))
Tinv = np.asarray(Tinv)

# output arrays
eEOS = evec
pEOS = fp(Tinv)*Tinv**4./(hbarc**3.)
sEOS = fs(Tinv)*Tinv**3./(hbarc**3.)
TEOS = Tinv

line = FortranRecordWriter('(4E15.6)')

# open output file
with open('HotQCD-EOS.dat','w') as wf:
    for i in range(len(eEOS)):
        wf.write(line.write([eEOS[i],pEOS[i],sEOS[i],TEOS[i]])+"\n")

# plot curves
#plt.plot(T, e, 'o')
#plt.plot(Tinv,fe(Tinv),"-")
#plt.plot(T, s, 's')
#plt.plot(Tinv,fs(Tinv),"-")
#plt.plot(Tinv,pEOS,"-")

# plot properties
#plt.legend(loc='upper right')
#plt.xlim(0.0, 0.8)
#plt.ylim(0.0, 5.0)
#plt.xlabel('$T$ MeV')
#plt.ylabel('$(e-3p)/T^4$')
#plt.xlim(0,800)
#plt.ylim(0,20)
#plt.show()
