#!/usr/bin/env python3

import numpy as np
import itertools
import sys

# Lorentz boost z-direction
def boost_z(p, beta_z):
    E, px, py, pz = p.T
    gamma = 1/np.sqrt(1-beta_z**2)
    E_ = gamma*(E - beta_z*pz)
    px_ = px
    py_ = py
    pz_ = gamma*(pz - beta_z*E)
    return np.array(E_, px_, py_, pz_)


# event generator
def events():
    event = []
    for l in sys.stdin.buffer:
        if len(l) > 200:
            p = l[97:192].replace(b'D', b'E')
            E = p[0:23]
            px = p[24:47]
            py = p[48:71]
            pz = p[72:95]
            ityp = l[218:221]
            chg = l[225:227]
            if int(chg) != 0 and int(ityp) == 101:
                event.append((t, x, y, z, E, px, py, pz)) 
        else:
            if event:
                yield np.array(event, dtype=float)
                event.clear()

    if event: 
        yield np.array(event, dtype=float)
                            
# loop over all hydro events and urqmd oversamples
for particles in events():

    # unpack particle data
    pvec = particles[5:8]
    E, px, py, pz = pvec.T

    # calculate pseudorapidity
    p = np.sqrt(np.square(pvec).sum(axis=1))
    eta = 0.5*np.log((p+pz)/(p-pz))

    # get all p1, p2 momentum pairs from same event
    for particle_1, particle_2 in itertools.combinations(particles,2):

	# unpack pair coordinates
	x1 = np.array(particle_1[0:4])
	x2 = np.array(particle_2[0:4])
	dx = (x2 - x1)[1:4]

	# unpack pair momenta
	p1 = np.array(particle_1[5:9])
	p2 = np.array(particle_2[5:9])

        # boost into longitudinal co-moving frame
        beta_z = (p1[3] + p2[3])/(p1[0] + p2[0])
        p1 = boost_z(p1, beta_z)[1:]
        p2 = boost_z(p2, beta_z)[1:]

        # calculate q and k vectors
        q = p2 - p1
        k = (p1 + p2)/2

        # project out the transverse component of k
        kt = np.sqrt(np.dot(k[0:2].T,k[0:2]))
        
        # project out q_out, q_side, q_long
        qo = np.sqrt(np.dot(q[0:2].T,q[0:2]))
        ql = q[2]
        qs = q - qo - ql

        # add up same event correlations
        C[iqo,iqs,iql,ikt] += 1 + np.cos(np.dot(q.T,dx)) 

