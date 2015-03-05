#!/usr/bin/env python
import numpy as np
import scipy 
import random
import fortranformat as ff
from fortranformat import FortranRecordWriter
import matplotlib.pyplot as plt
from scipy import interpolate
from scipy.interpolate import interp1d
from scipy import integrate
from scipy.integrate import quad

# constants
hbarc = 0.19733
T0 = 130 

# loop over spline files
for ic in range(100):
	
	# load interaction measure (e-3p)/T**4 as function of temperature
	Tvec, I = np.loadtxt("realizations/hotqcd-measure/hotqcd-eos-spline_{}.dat".format(ic),unpack=True)
	fI = interp1d(Tvec, I, bounds_error=True)
	integrand = lambda x: fI(x)/x
	nT = Tvec.size
	Tmin = Tvec[0]
	Tmax = Tvec[nT-1]

	# pressure over T^4
	p = []
	pconst = 0.439122*(1+0.1*random.uniform(-1,1))
	p.append(pconst + quad(integrand, T0, Tmin)[0])
	for iT in range(nT-1):
		Tl = Tvec[iT]
		Th = Tvec[iT+1]
    		p.append(p[iT] + quad(integrand, Tl, Th)[0])
	p = np.asarray(p)
	fp = interp1d(Tvec, p, bounds_error=True)

	# energy density over T^4 
	e = []
	for iT in range(nT):
    		e.append(I[iT] + 3.*p[iT])
	fe = interp1d(Tvec, e, bounds_error=True)

	# entropy density over T^3
	s = []
	for iT in range(nT):
    		s.append(e[iT] + p[iT])
	fs = interp1d(Tvec, s, bounds_error=True)


	# convert e in GeV to temperature in GeV
	def e2T(e):
    		counts = 0
    		error = 1.0
    		Tlb = Tmin
    		Tub = Tmax
    		Tmid = (Tub+Tlb)/2.

    		while error > 1e-12: 
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
	printed = True
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
