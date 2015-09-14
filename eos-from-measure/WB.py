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

# Wuppertal-Budapest (e-3p)/T**4
#################################################################################
def wb_e3p_T4(T, T0=0.2, h0=0.1396, h1=-0.1800, h2=0.0350, 
              f0=1.05, f1=6.39, f2=-4.72, g1=-0.92, g2=0.57):
    t = T/T0
    return np.exp(- h1/t - h2/t**2) * \
        (h0 + f0 * (np.tanh(f1*t + f2) + 1)/(1 + g1*t + g2*t**2))

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
    e3p_T4_mid = (1-s)*urqmd_e3p_T4(Tmid) + s*wb_e3p_T4(Tmid)
    e3p_T4_hi  = wb_e3p_T4(Thi)
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
    with open('WB.dat','w') as wf:
        for i in range(len(e_table)):
            wf.write(line.write([e_table[i], p_table[i], s_table[i], T_table[i]])+"\n")

if __name__ == "__main__":
    main()
