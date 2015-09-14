#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
import fortranformat as ff
from fortranformat import FortranRecordWriter
from scipy.interpolate import InterpolatedUnivariateSpline as Spline

import seaborn as sns

# hbarc in [GeV fm]
hbarc = 0.1973269718

# UrQMD (e-3p)/T**4
####################################################################################
urqmd_T, urqmd_I = np.loadtxt("hrg-eos/hrg-urqmd-eos.dat", dtype='float', usecols=(0,1), 
                              unpack=True, skiprows=1)

def urqmd_e3p_T4(T):
    spl = Spline(urqmd_T/1000., urqmd_I)
    return spl(T)

# HotQCD (e-3p)/T**4
#################################################################################
def hotqcd_p_T4(T, Tc=0.154, pid=95*np.pi**2/180, 
             ct=3.8706, an=-8.7704, bn=3.9200, cn=0, dn=0.3419, 
             t0=0.9761, ad=-1.2600, bd=0.8425, cd=0, dd=-0.0475):

    t = T/Tc
    return 0.5 * (1 + np.tanh(ct*(t-t0))) * \
        (pid + an/t + bn/t**2 + cn/t**3 + dn/t**4) / \
        (1   + ad/t + bd/t**2 + cd/t**3 + dd/t**4)

def hotqcd_e3p_T4(T):
    spl = Spline(T, hotqcd_p_T4(T))
    return T*spl(T, nu=1)

# HotQCD (e-3p)/T**4 splines
#########################################################################################################################

# T [GeV] partitions for piecewise function
T = np.linspace(0.005,0.800,1000)
Tlo  = T[T < 0.155]
Tmid = T[(T >= 0.155) & (T < 0.180)]
Thi  = T[T >= 0.180]

# load bootstrap coefficients
k1, k4, k2, k3, C5, C6, C7, C8, C9, C10, C11, C12, C13, C14, C15, C16 = np.loadtxt("hotqcd-eos/BootCoeffs.txt", dtype='float', 
                                                                                   delimiter='\t', unpack=True, usecols=range(1,17), skiprows=2)

# loop over splines (range fixed to 1 to test first interaction measure)
for ic in range(100):

    # initialize basis functions
    B1, B2, B3, B4, B5, B6, B7, B8, hotqcd_e3p_T4_spline = (np.zeros(T.size),)*9
    B = np.array((B1, B2, B3, B4, B5, B6, B7, B8))

    # knot vector
    t = (k1[ic], k1[ic], k1[ic], k2[ic], k3[ic], k4[ic], k4[ic], k4[ic], k4[ic])
    print("knot vector:",t)

    # initialize basis functions
    for i in range(8):
        for iT, T_ in enumerate(T*1000):
            if t[i] <= T_ and T_ < t[i+1]:
                B[i][iT] = 1.0

    # iterate recursive algorithm thrice B0 -> B1 -> B2 -> B3
    for k in range(1, 4):
        for i in range(5):
            for iT, T_ in enumerate(T*1000):
                B[i][iT] = (T_-t[i])/(t[i+k]-t[i]+1e-15)*B[i][iT] + (t[i+k+1]-T_)/(t[i+k+1]-t[i+1]+1e-15)*B[i+1][iT] 
    
    # construct hotqcd spline interaction measure
    hotqcd_e3p_T4_spline = Spline(T, C6[ic] + C5[ic]*(T*1000 - 130) + C7[ic]*B[1] + C8[ic]*B[2] + C9[ic]*B[3] + C10[ic]*B[4])
    
    # evaluate trace anomaly (e-3p)/T**4 on temperature intervals
    e3p_T4_lo  = urqmd_e3p_T4(Tlo)
    z = np.linspace(0, 1, Tmid.size)
    s = 6*z**5 - 15*z**4 + 10*z**3
    e3p_T4_mid = (1-s)*urqmd_e3p_T4(Tmid) + s*hotqcd_e3p_T4_spline(Tmid)
    e3p_T4_hi  = hotqcd_e3p_T4_spline(Thi)
    e3p_T4 = np.concatenate([e3p_T4_lo, e3p_T4_mid, e3p_T4_hi])

    # (e-3p)/T**4, p/T**4, e/T**4 and s/T**3
    p_T4_interp = Spline(T, e3p_T4/T).antiderivative()
    e_T4  = e3p_T4  + 3*p_T4_interp(T)
    s_T3 = e_T4  + p_T4_interp(T)

    # e [GeV/fm**3], p [GeV/fm**3], s [1/fm**3]
    e = e_T4*T**4/hbarc**3
    p = p_T4_interp(T)*T**4/hbarc**3
    s = s_T3*T**3/hbarc**3

    plt.plot(T, e/T**4)
    plt.plot(T, p/T**4)
    plt.plot(T, s/T**3)

    # power law extrapolation
    #for iT, T_ in enumerate(T):
    #    if T_ > 0.4:

    # express as functions of energy density
    p_interp = Spline(e, p)
    s_interp = Spline(e, s) 
    T_interp = Spline(e, T) 

    # e output mesh: GeV/fm**3
#    e_table = np.arange(1,311000,2)*1e-3
    e_table = np.exp(np.linspace(-6.9078,6.476972,1000))
    p_table = p_interp(e_table)
    s_table = s_interp(e_table)
    T_table = T_interp(e_table)

    # print to file in fortran format
#    line = FortranRecordWriter('(4E15.6)')
#    with open('hotqcd.dat','w') as wf:
#        for i in range(len(e_table)):
#            wf.write(line.write([e_table[i], p_table[i], s_table[i], T_table[i]])+"\n")

#plt.plot(T, hotqcd_e3p_T4(T), color='r', lw=4)
#plt.plot(T, urqmd_e3p_T4(T), color='k', lw=4)
plt.xlim(0.050,0.800)
#plt.ylim(0,5)
plt.show()
