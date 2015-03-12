#!/usr/bin/env python
import numpy as np
import scipy 
import matplotlib.pyplot as plt
import prettyplotlib as ppl
import fortranformat as ff
from fortranformat import FortranRecordWriter
from prettyplotlib import brewer2mpl
from scipy import interpolate
from scipy.interpolate import interp1d
from scipy.interpolate import UnivariateSpline
from scipy.optimize import fminbound
from scipy import integrate
from scipy.integrate import quad

####################################################################################

# construct low temperature interaction measure 

#T,QS,PN,EN,PEN,SN,CV,SNCV = np.loadtxt("hrg-eos/OUT_2.5.DAT140_2014_bulk",dtype='float',unpack=True,skiprows=2)
#fIHRG = interp1d(T, QS, kind='cubic')

T, I, e, p, Cs = np.loadtxt("hrg-eos/hrg-urqmd-eos.dat",dtype='float',unpack=True,skiprows=1)
fI_lo = interp1d(T, I, kind='cubic')

#T = np.linspace(5,200,100)
#ppl.plot(T/1000,fI_lo(T),'--',color='orange',linewidth=3.0,zorder = 301, label="HRG UrQMD partial reson")
####################################################################################

# construct high temperature interaction measure
 
#T,I,p,e,s,cv,cs = np.loadtxt("hotqcd-eos/EOS-table.txt",dtype='float',unpack=True,skiprows=1)
#T = np.linspace(65,800,1000)
#ppl.plot(T/1000,I,color="red",linewidth=3.0,zorder = 100, label="HotQCD continuum extr.")

# temperature limits [GeV]
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

# pressure over T^4
p = 0.5*(1.+np.tanh(ct*(tb-t0)))*(pid+an/tb+bn/tb**2.+cn/tb**3.+dn/tb**4.)/(1.+ad/tb+bd/tb**2.+cd/tb**3.+dd/tb**4.)
dfp = interp1d(T, interpolate.splev(T,interpolate.splrep(T,p,s=0),der=1), kind='cubic')
fp = interp1d(T, p, kind='cubic')

# trace anomaly over T^4
I = []
for iT in T:
    I.append(dfp(iT)*iT)
fI_hi = interp1d(T*1000, I, kind='cubic')

# plot HotQCD trace anomaly))
#T = np.linspace(130,800,1000)
#ppl.plot(T/1000,fI_hi(T),color="blue",linewidth=3.0,zorder = 99,label="HotQCD parameterization")

#########################################################################################################################

#construct intermediate interaction measure

# define temperature range and gradations
Tmin = 5    # minimum temp
Tmax = 800  # maximum temp
nT = 10000  # partitions
dT = (Tmax-Tmin)/float(nT) # temperature step size
Tvec = np.linspace(Tmin,Tmax,nT)

# import a set of spline coefficients
k1, k4, k2, k3, C5, C6, C7, C8, C9, C10, C11, C12, C13, C14, C15, C16 = np.loadtxt("hotqcd-eos/BootCoeffs.txt",dtype='float',delimiter='\t',unpack=True,usecols=(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16),skiprows=2)

# spline count
N = k1.size

# loop over splines (range fixed to 1 to test first interaction measure)
for ic in range(5):

    # initialize basis functions
    B1 = np.zeros(nT)
    B2 = np.zeros(nT)
    B3 = np.zeros(nT)
    B4 = np.zeros(nT)
    B5 = np.zeros(nT)
    B6 = np.zeros(nT)
    B7 = np.zeros(nT)
    B8 = np.zeros(nT)
    B = (B1,B2,B3,B4,B5,B6,B7,B8)

    # initialize interaction measure
    I   = np.zeros(nT)
    I2  = np.zeros(nT)
    dI  = np.zeros(nT)
    Ilo = np.zeros(nT)
    Imid = np.zeros(nT)
    Ihi = np.zeros(nT)
    dIlo = np.zeros(nT)
    dIhi = np.zeros(nT)
    
    # knot vector
    t = (k1[ic],k1[ic],k1[ic],k2[ic],k3[ic],k4[ic],k4[ic],k4[ic],k4[ic])
    print("knot vector:",t)

    # initialize basis functions
    for i in range(8):
        for iT,T in enumerate(Tvec):
            if t[i] <= T and T < t[i+1]:
                B[i][iT] = 1.0
    
    # iterate recursive algorithm thrice B0 -> B1 -> B2 -> B3
    for k in range(1,4):
        for i in range(5):
            for iT,T in enumerate(Tvec):
                B[i][iT] = (T-t[i])/(t[i+k]-t[i]+1e-15)*B[i][iT] + (t[i+k+1]-T)/(t[i+k+1]-t[i+1]+1e-15)*B[i+1][iT] 

    # construct intermediate interaction measure
    Imid = C6[ic] + C5[ic]*(Tvec-130) + C7[ic]*B[1] + C8[ic]*B[2] + C9[ic]*B[3] + C10[ic]*B[4]
    fI_mid = interp1d(Tvec,Imid)  

    # determine high temperature switching parameters
    Td = np.linspace(350,400,1000)
    dfI_mid = interp1d(Td, interpolate.splev(Td,interpolate.splrep(Td,fI_mid(Td),s=0),der=1))
    dfI_hi = interp1d(Td, interpolate.splev(Td,interpolate.splrep(Td,fI_hi(Td),s=0),der=1))
    deriv = interp1d(Td,abs(dfI_hi(Td)-dfI_mid(Td)))
    Tsw = fminbound(deriv,350,400)
    shift = fI_mid(Tsw) - fI_hi(Tsw)

    # blend low, intermediate and high regions together
    for iT,T in enumerate(Tvec):
        if T < 120:
            I[iT] = fI_lo(T)
        elif T >= 120 and T < 180:
            w1 = lambda x: 5/(1+np.exp((x-176)/2))
            u1 = lambda x: (1+np.tanh((T-150)/w1(T)))/2
            I[iT] = fI_lo(T) + u1(T)*(fI_mid(T)-fI_lo(T))
        elif T >= 180 and T < 500:
            w2 = lambda x: 40*(1/(1+np.exp(-(x-190)/5)) * 1/(1+np.exp((x-490)/5)))
            u2 = lambda x: (1+np.tanh((T-Tsw)/w2(T)))/2
            I[iT] = fI_mid(T) + u2(T)*(fI_hi(T) + shift - fI_mid(T))
        else:
            I[iT] = fI_hi(T) + shift  

    fI = interp1d(Tvec,I)   	            
    
    # calculate p/T**4 from interaction measure
    p = []
    integral = 0.0
    for iT,T in enumerate(Tvec):
        p.append(integral)
        integral += I[iT]/T*dT
    fp = interp1d(Tvec,p)
   
    # calculate e/T**4 
    e = []
    for iT,T in enumerate(Tvec):
        e.append(I[iT] + 3.*p[iT])
    fe = interp1d(Tvec,e)

    # calculate s/T**3
    s = []
    for iT,T in enumerate(Tvec):
        s.append(e[iT] + p[iT])
    fs = interp1d(Tvec,s)

    ppl.plot(Tvec,e,label="energy")
    ppl.plot(Tvec,p,label="pressure")
    ppl.plot(Tvec,s,label="entropy")
    #continue   

    # e output mesh: GeV/fm^3
    emin = 0.1e-2
    emax = 0.310999e3
    #ne = 155500
    ne = 100
    de = (emax-emin)/float(ne)
    evec = np.linspace(emin,emax,ne)

    # T output mesh: GeV
    Tinv = []
    for ie0,e0 in enumerate(evec):
        Tinv.append(fminbound(lambda x: abs(fe(x)*(x/1000)**4./(hbarc**3.)-e0),Tmin,Tmax))
    Tinv = np.asarray(Tinv)

    # output arrays
    eEOS = fe(Tinv)*Tinv**4./(hbarc**3.)
    pEOS = fp(Tinv)*Tinv**4./(hbarc**3.)
    sEOS = fs(Tinv)*Tinv**3./(hbarc**3.)
    TEOS = Tinv

    line = FortranRecordWriter('(4E15.6)')

    # open output file
    with open("eos-tables/HotQCD-EOS-{}.dat".format(ic),'w') as wf:
        for i in range(len(eEOS)):
            wf.write(line.write([eEOS[i],pEOS[i],sEOS[i],TEOS[i]])+"\n")

    # plot that shit!
    ppl.plot(Tvec/1000,p) 

    #dI = interp1d(Tvec, interpolate.splev(Tvec,interpolate.splrep(Tvec,I,s=0),der=2))
    #ppl.plot(Tvec/1000,dI(Tvec))

    # save to text output
    #np.savetxt("realizations/hotqcd-measure/hotqcd-eos-spline_{}.dat".format(ic), np.c_[Tvec,I], fmt='%10.5f')
    

# uncomment to plot basis splines    
#plt.plot(xvec,B[0],label="B1")
#plt.plot(xvec,B[1],label="B2")
#plt.plot(xvec,B[2],label="B3")
#plt.plot(xvec,B[3],label="B4")
#plt.plot(xvec,B[4],label="B5")
#plt.plot(xvec,B[5],label="B6")
#plt.plot(xvec,B[6],label="B7")
#plt.plot(xvec,B[7],label="B8")

#######################################################################



# plot properties
#plt.xlim([0.0,0.8])
#plt.ylim([0.0,5.0])
plt.xlabel("$T$ [GeV]")
plt.ylabel("$(e-3p)/T^4$")
plt.legend(loc='upper right', fontsize=17, frameon=False)
plt.show()
