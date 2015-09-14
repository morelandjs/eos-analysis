#!/usr/bin/env python


import numpy as np
import fortranformat as ff
from fortranformat import FortranRecordWriter
from scipy.interpolate import InterpolatedUnivariateSpline as Spline

# hbarc in [GeV fm]
hbarc = 0.1973269718

# UrQMD (e-3p)/T**4
################################################################################
urqmd_T, urqmd_I, = np.loadtxt("tables/urqmd/hrg-urqmd-eos.dat", dtype='float', 
                               usecols=(0,1), unpack=True, skiprows=1)
   
def urqmd_e3p_T4(T):
    spl = Spline(urqmd_T/1000., urqmd_I)
    return spl(T) 

# HotQCD p/T**4
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

# Merge the UrQMD trace anomaly with the HotQCD trace anomaly and print the EOS
###############################################################################
def main():

    # T [GeV] partitions for piecewise function
    T = np.linspace(0.005, 0.800, 1000)
    Tlo  = T[T < 0.155]
    Tmid = T[(T >= 0.155) & (T < 0.18)]
    Thi  = T[T >= 0.18]
    
    # evaluate trace anomaly (e-3p)/T**4 on temperature intervals
    e3p_T4_lo  = urqmd_e3p_T4(Tlo)
    z = np.linspace(0, 1, Tmid.size)
    s = 6*z**5 - 15*z**4 + 10*z**3
    e3p_T4_mid = (1-s)*urqmd_e3p_T4(Tmid) + s*hotqcd_e3p_T4(Tmid)
    e3p_T4_hi  = hotqcd_e3p_T4(Thi)
    e3p_T4 = np.concatenate([e3p_T4_lo, e3p_T4_mid, e3p_T4_hi])

    # (e-3p)/T**4, p/T**4, e/T**4 and s/T**3
    p_T4_interp = Spline(T, e3p_T4/T).antiderivative()
    e_T4  = e3p_T4  + 3*p_T4_interp(T)
    s_T3 = e_T4  + p_T4_interp(T)

    # e [GeV/fm**3], p [GeV/fm**3], s [1/fm**3]
    e = e_T4*T**4/hbarc**3
    p = p_T4_interp(T)*T**4/hbarc**3
    s = s_T3*T**3/hbarc**3

    # express as functions of energy density
    p_interp = Spline(e, p)
    s_interp = Spline(e, s) 
    T_interp = Spline(e, T) 

    # e output mesh: GeV/fm**3
    e_table = np.arange(1,311000,2)*1e-3
    #e_table = np.exp(np.linspace(-6.9078,6.476972,1000))
    p_table = p_interp(e_table)
    s_table = s_interp(e_table)
    T_table = T_interp(e_table)

    # print to file in fortran format
    line = FortranRecordWriter('(4E15.6)')
    with open('HotQCD.dat','w') as wf:
        for i in range(len(e_table)):
            wf.write(line.write([e_table[i], p_table[i], s_table[i], T_table[i]])+"\n")

if __name__ == "__main__":
    main()
