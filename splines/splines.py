#!/usr/bin/env python
import numpy as np
import scipy 
import matplotlib.pyplot as plt
import prettyplotlib as ppl
from prettyplotlib import brewer2mpl
from scipy import interpolate
from scipy.interpolate import interp1d
from scipy import integrate
from scipy.integrate import quad

####################################################################################

# construct low temperature HRG EOS
T,QS,PN,EN,PEN,SN,CV,SNCV = np.loadtxt("hrg-eos/OUT_2.5.DAT140_2014_bulk",dtype='float',unpack=True,skiprows=2)
fIHRG = interp1d(T, QS, kind='cubic')
N = T.size
ppl.plot(T[0:N/2]*1000.,QS[0:N/2],'r--',color="purple",linewidth=3.0,zorder = 300,label="hadron resonance gas")


####################################################################################

# plot best fit spline from paper table
T,I,p,e,s,cv,cs = np.loadtxt("hotqcd-eos/EOS-table.txt",dtype='float',unpack=True,skiprows=1)
ppl.plot(T,I,color="red",linewidth=3.0,zorder = 100, label="HotQCD continuum extr.")

# plot parameterized best fit spline

# constants
hbarc = 0.19733

# temperature limits in GeV
Tmin = 0.05
Tmax = 0.8
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

# plot
N = T.size
ppl.plot(T[N/3:N]*1000.,fI(T[N/3:N]),color="blue",linewidth=3.0,zorder = 99,label="HotQCD parameterization")

###################################################################################################################################

# define temperature range and gradations
Tmin = 50.0 # lower knot
Tmax = 610.0 # upper knot
nT = 10000   # partitions
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
    print "knot vector:",t

    # initialize basis functions
    for i in range(8):
        for iT in range(nT):
            T = Tmin+iT*dT
            if t[i] <= T and T < t[i+1]:
                B[i][iT] = 1.0


    # iterate recursive algorithm thrice B0 -> B1 -> B2 -> B3
    for k in range(1,4):
        for i in range(5):
            for iT in range(nT):
                T = Tmin+iT*dT
                B[i][iT] = (T-t[i])/(t[i+k]-t[i]+1e-15)*B[i][iT] + (t[i+k+1]-T)/(t[i+k+1]-t[i+1]+1e-15)*B[i+1][iT] 
            
    # curve shifting parameters
    Tlsw = 130.
    Thsw = 400.
    iTlsw = int((Tlsw-Tmin)/dT)
    iThsw = int((Thsw-Tmin)/dT)

    # construct intermediate interaction measure
    Imid = C6[ic] + C5[ic]*(Tvec-130.) + C7[ic]*B[1] + C8[ic]*B[2] + C9[ic]*B[3] + C10[ic]*B[4]

    # construct low temp HRG limit
    for iT in range(nT):
        T = Tmin+iT*dT
        if(T < Tlsw):
            Ilo[iT] = fIHRG(T/1000.)
        else:
            Ilo[iT] = Imid[iT]

    # find where low temp first derivative matches
    iTlb = int((129.-Tmin)/dT)
    iTub = int((131.-Tmin)/dT)
    
    for iT in range(nT-1):
        T = Tmin + iT*dT
        dIlo[iT] = (Ilo[iT+1]-Ilo[iT])/dT
    dImid = (Imid[iTlsw+1]-Imid[iTlsw])/dT
    dI = np.abs(dImid - dIlo)[iTlb:iTub]
   
    iT0 = np.where(dI == dI.min())[0][0] + iTlb
    T0 = Tmin + iT0*dT

    # perform HRG T shift
    Tshift = -(iTlsw-iT0)*dT
    for iT in range(nT):
        T = Tmin+iT*dT
        Tnew = min(max(T-Tshift,50.),180.)
        Ilo[iT] = fIHRG(Tnew/1000.)

    # perform HRG I shift
    Ishift = Imid[iTlsw] - Ilo[iTlsw] 
    for iT in range(nT):
        Ilo[iT] += Ishift

    # blend at cut
    for iT in range(nT):
        T = Tmin+iT*dT
        Imid[iT] = Ilo[iT] + (1.+np.tanh((T-Tlsw)/20.))/2.*(Imid[iT] - Ilo[iT])

    # blend second time into HRG curve
    for iT in range(nT):
        T = Tmin+iT*dT
        if(50. < T and T < 180.):
            Imid[iT] = fIHRG(T/1000.) + (1.+np.tanh((T-100)/40.))/2.*(Imid[iT] - fIHRG(T/1000.))

    # construct high temp limit
    for iT in range(nT):
        T = Tmin+iT*dT
        Ihi[iT] = fI(T/1000.)

    # find where high temp first derivative matches
    iTlb = int((300.-Tmin)/dT)
    iTub = int((500.-Tmin)/dT)

    dImid = (Imid[iThsw+1]-Imid[iThsw])/dT
    for iT in range(nT-1):
        T = Tmin + iT*dT
        dIhi[iT] = (Ihi[iT+1]-Ihi[iT])/dT
    dI = np.abs(dIhi - dImid)[iTlb:iTub]
    iT0 = np.where(dI == dI.min())[0][0] + iTlb
    T0 = Tmin + iT0*dT

    # perform T shift
    Tshift = (iThsw-iT0)*dT
    for iT in range(nT):
        T = Tmin+iT*dT
        Ihi[iT] = fI(max(T-Tshift,50.)/1000.)

    # perform I shift
    Ishift = Ihi[iThsw] - Imid[iThsw]
    for iT in range(nT):
        Ihi[iT] -= Ishift

    for iT in range(nT):
        T = Tmin+iT*dT
        I[iT] = Imid[iT]
        if (T > 400.):
            I[iT] = Ihi[iT]

    # blend at cut
    for iT in range(nT):
        T = Tmin+iT*dT
        I[iT] = Imid[iT] + (1.+np.tanh((T-Thsw)/50.))/2.*(Ihi[iT] - Imid[iT])

    # calc second derivative
    for iT in range(1,nT-1):
        T = Tmin+iT*dT
        I2[iT] = (I[iT+1]-2.*I[iT]+I[iT-1])/(dT*dT)
        

    # plot that shit!
    ppl.plot(Tvec,I)
    

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
plt.xlim([50,600])
plt.ylim([0.0,5.0])
plt.xlabel("T [GeV]")
plt.ylabel("I(T)")
plt.legend(loc='upper right', fontsize=17, frameon=False)
plt.show()
