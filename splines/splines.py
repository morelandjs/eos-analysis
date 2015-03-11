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

# construct low temperature HRG EOS (HotQCD full resonances)
#T,QS,PN,EN,PEN,SN,CV,SNCV = np.loadtxt("hrg-eos/OUT_2.5.DAT140_2014_bulk",dtype='float',unpack=True,skiprows=2)
#fIHRG = interp1d(T, QS, kind='cubic')
#fPHRG = interp1d(T, PN, kind='cubic')
#THRG = np.linspace(30,220,100)
#ppl.plot(THRG/1000,fIHRG(THRG),'--',color="purple",linewidth=3.0,zorder = 300,label="HRG full reson")

# construct low temperature HRG EOS (UrQMD partial resonances)
T, I, e, p, Cs = np.loadtxt("hrg-eos/hrg-urqmd-eos.dat",dtype='float',unpack=True,skiprows=1)
fIHRG = interp1d(T, I, bounds_error=False, fill_value=0.0, kind='cubic')
fPHRG = interp1d(T, p, kind='cubic')
THRG = np.linspace(5,200,100)
ppl.plot(THRG/1000,fIHRG(THRG),'--',color='orange',linewidth=3.0,zorder = 301, label="HRG UrQMD partial reson")

####################################################################################

# plot average HotQCD QGP EOS
T,I,p,e,s,cv,cs = np.loadtxt("hotqcd-eos/EOS-table.txt",dtype='float',unpack=True,skiprows=1)
T /= 1000
#ppl.plot(T,I,color="red",linewidth=3.0,zorder = 100, label="HotQCD continuum extr.")

# constants
hbarc = 0.19733

# temperature limits [GeV]
Tmin = 0.050
Tmax = 0.800
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
fI = interp1d(T*1000, I, bounds_error=False, fill_value=0.0, kind='cubic')

# plot HotQCD trace anomaly
TQGP = np.linspace(130,800,1000)
#ppl.plot(TQGP/1000,fI(TQGP),color="blue",linewidth=3.0,zorder = 99,label="HotQCD parameterization")

###################################################################################################################################

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
    print "knot vector:",t

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

    # construct low/high temp interaction measures
    for iT,T in enumerate(Tvec):
	Ihi[iT] = fI(T)
	if T < 180:
	   Ilo[iT] = fIHRG(T)
	if T >= 180 and T < 200:
	   Ilo[iT] = fIHRG(T) + (1.+np.tanh((T-190)/5))/2.*(fI(T) - fIHRG(T))
	else:
	   Ilo[iT] = fI(T)
    ppl.plot(Tvec/1000,Ihi)

    # blend Ilo into Imid
    Tlsw = 150
    dTlsw = 10

    for iT,T in enumerate(Tvec):
	if T < 600:
    		Imid[iT] = Ilo[iT] + (1.+np.tanh((T-Tlsw)/dTlsw))/2.*(Imid[iT] - Ilo[iT])
	else:
		Imid[iT] = Imid[np.argmax(Tvec >= 600)]
    #ppl.plot(Tvec/1000,Imid)

    # blend Imid into Ihi
    Thsw = 400 
    Tdlo = 350
    Tdhi = 400
    iThsw = int((Thsw-Tmin)/dT)
    iTdlo = int((Tdlo-Tmin)/dT)
    iTdhi = int((Tdhi-Tmin)/dT)

    # calc |dImid/dT - dIhi/dT|
    dImid = np.zeros(nT)
    dIhi = np.zeros(nT)
    for iT in range(nT-1):
    	dImid[iT] = (Imid[iT+1]-Imid[iT])/dT
	dIhi[iT]  = (Ihi[iT+1]-Ihi[iT])/dT
    dslope = np.absolute(dImid-dIhi)[iTdlo:iTdhi]
    	
    # minimize |dImid/dT - dIhi/dT| from dlo to dhi
    iTdmin = np.where(dslope == dslope.min())[0][0] + iTdlo 

    # perform I shift
    Ishift = Ihi[iTdmin] - Imid[iTdmin]
    for iT in range(nT):
        Ihi[iT] -= Ishift

    for iT,T in enumerate(Tvec):
	Thsw = Tmin+iTdmin*dT
        I[iT] = Imid[iT]
        if T > Thsw:
            I[iT] = Ihi[iT]

    # blend at cut
    for iT,T in enumerate(Tvec):
	if T < 600:
        	I[iT] = Imid[iT] + (1.+np.tanh((T-Thsw)/40))/2.*(Ihi[iT] - Imid[iT])

    # calc second derivative
    for iT in range(1,nT-1):
        T = Tmin+iT*dT
        I2[iT] = (Imid[iT+1]-2.*Imid[iT]+Imid[iT-1])/(dT*dT)
    #ppl.plot(Tvec,I2)    

    # plot that shit!
    #ppl.plot(Tvec/1000,I)

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
