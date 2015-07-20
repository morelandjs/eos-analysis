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

hbarc = 0.19733

####################################################################################

# construct low temperature interaction measure 

#T,QS,PN,EN,PEN,SN,CV,SNCV = np.loadtxt("hrg-eos/OUT_2.5.DAT140_2014_bulk",dtype='float',unpack=True,skiprows=2)
#fIHRG = interp1d(T, QS, kind='cubic')

T, I, e, p, Cs = np.loadtxt("hrg-eos/hrg-urqmd-eos.dat",dtype='float',unpack=True,skiprows=1)
fI_lo = interp1d(T/1000., I, kind='cubic')
#plt.plot(T/1000., I)
##################################################################################

# load s95_PCE EOS from the VISHNU table

e, p, s, T = np.loadtxt("old/s95-PCE-EOS-old.dat", unpack=True)
I = (e - 3.*p)*(hbarc**3)/(T**4)
slope = interp1d(T, interpolate.splev(T,interpolate.splrep(T,I,s=0),der=1))
fI_hi = interp1d(T,I)

# tame weird PCE divergence
for iT, T_ in enumerate(T):
    if T_ < 0.17:
        I[iT] = fI_hi(0.17) + slope(0.17)*(T_ - 0.17)
fI_hi = interp1d(T,I)

# extrapolate a little bit further out
I = np.zeros(10000)
Tvec = np.linspace(T[0], 0.8, 10000)
for iT, T_ in enumerate(Tvec):
    if T_ > 0.631:
        I[iT] = fI_hi(0.631) + slope(0.631)*(T_ - 0.631)
    else:
        I[iT] = fI_hi(T_)
fI_hi = interp1d(Tvec, I)
#plt.plot(Tvec, I)
#################################################################################

# blend UrQMD HRG interaction measure into s95 interaction measure

# temperature limits in GeV
Tmin = 0.005
Tmax = 0.7
nT = 10000
dT = (Tmax-Tmin)/float(nT)
T = np.linspace(Tmin,Tmax,nT)
I = np.zeros(nT)

# piecewise blending
for iT, T_ in enumerate(T):
    if T_ < 0.120:
        I[iT] = fI_lo(T_)
    elif T_ >= 0.120 and T_ < 0.180:
        u1 = lambda x: (1+np.tanh((T_-0.159)/0.005))/2
        I[iT] = fI_lo(T_) + u1(T_)*(fI_hi(T_)-fI_lo(T_))
    else:
        I[iT] = fI_hi(T_) 

fI = interp1d(T,I)   	         
#plt.plot(T, I)

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
evec = np.arange(1,311000,2)*1e-3
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
with open('s95-PCE-EOS.dat','w') as wf:
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
