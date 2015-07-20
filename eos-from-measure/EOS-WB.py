#!/usr/bin/env python
import numpy as np
import scipy 
import fortranformat as ff
from fortranformat import FortranRecordWriter
import matplotlib.pyplot as plt
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
T0 = 0.2
Tmin = 0.065
Tmax = 0.800
dT = 0.005
nT = int((Tmax-Tmin)/dT)
Tvec = np.linspace(Tmin,Tmax,nT)
t=Tvec/T0

# eos parameters
hbarc=0.19733
h0 = 0.1396
h1 = -0.1800
h2 = 0.0350 
f0 = 1.05 
f1 = 6.39 
f2 = -4.72 
g1 = -0.92 
g2 = 0.57

# trace anomaly: I/T^4
I = np.exp(-h1/t-h2/t**2.)*(h0 + f0*(np.tanh(f1*t+f2)+1.)/(1.+g1*t+g2*t**2.)) 
I = np.asarray(I)
fI_hi = interp1d(Tvec, I, kind='cubic')

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

# pressure: p/T^4
#p = []
#for T in Tvec:
#    p.append(integrate.quad(lambda x: np.exp(-h1/x-h2/x**2.)*(h0 + f0*(np.tanh(f1*x+f2)+1.)/(1.+g1*x+g2*x**2.))/x, 0.0, T/T0)[0])
#p = np.asarray(p)
#fp = interp1d(Tvec, p, kind='cubic')

# energy density: e/T^4
#e = []
#for iT in range(nT):
#    e.append(tr[iT] + 3*p[iT])
#e = np.asarray(e)
#fe = interp1d(Tvec, e, kind='cubic')

# entropy density: s/T^3
#s = []
#for iT in range(nT):
#    T = Tmin + iT*dT;
#    s.append(e[iT]+p[iT])
#s = np.asarray(s)
#fs = interp1d(Tvec, s, kind='cubic')   

# e output mesh: GeV/fm**3
evec = np.arange(1,650000,2)*1e-3
print(evec.size)

# T output mesh
Tinv = []
for ie0, e0 in enumerate(evec):
    Tinv.append(fminbound(lambda x: abs(fe(x)*(x**4.)/(hbarc**3.)-e0),Tmin,Tmax))
Tinv = np.asarray(Tinv)

# plot curves
#plt.plot(Tvec,e,'o')
#plt.plot(Tinv,fe(Tinv),'-')
#plt.plot(Tvec,s,'o')
#plt.plot(Tinv,fs(Tinv),'-')

# output arrays
eEOS = evec
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
#plt.show()
