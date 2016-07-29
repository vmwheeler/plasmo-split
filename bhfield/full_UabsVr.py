#!/usr/bin/env python
import numpy as np
from scipy.integrate import simps
import matplotlib.pyplot as plt

data = np.loadtxt("./UabsVr_full.dat",skiprows=3)

#check that this integrates to proper absorption cross section
rads = data[:,0] # in meters
Uabss = data[:,1]

RADCOR = 30.0E-9
RADCOT = 60.0E-9
CC=2.99792458E8 # light speed [m s-1]
EPSVAC=1.0E7/(4.0*np.pi*CC*CC) # eps0[F m-1]
MU=4.0*np.pi*1.0E-7 # assume non-magnetic (MU=MU0=const) [N A-2]

Uabssr2 = []
for rad,Uabs in zip(rads,Uabss):
  Uabssr2.append(Uabs*rad*rad)


UABS_T = simps(Uabssr2,rads)
AP = np.pi*RADCOT*RADCOT
MYQABS = UABS_T/AP

print "MYQABS = " + str(MYQABS)


plt.figure()
plt.plot(rads*1.E9,Uabss)
plt.title('Angle-averaged solar-weighted radial power absorption\n '
  + 'of a Au-CeO_2 core-shell nanoparticle\n ')
plt.xlabel('r [nm]')
plt.ylabel('U_abs [(W m^-3)')
plt.show()

