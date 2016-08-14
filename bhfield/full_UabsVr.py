#!/usr/bin/env python
import numpy as np
from scipy.integrate import simps
import matplotlib.pyplot as plt

absdata = np.loadtxt("./UabsVr.dat",skiprows=3)

#check that this integrates to proper absorption cross section
rads = absdata[:,0] # in meters
Uabss = absdata[:,1]

chosenTindex = 10
emdata = np.loadtxt("./UemVr.dat",skiprows=4)
Uems = emdata[:,chosenTindex]

f = open("./UemVr.dat")
tline = f.readlines()[2]
tline = tline.split()[1:]
temp = float(tline[chosenTindex])


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

cratio = 1
concUabss = cratio*Uabss

plt.figure()
plt.plot(rads*1.E9,concUabss,label="P_abs")
plt.plot(rads*1.E9,Uems,label="P_em")
plt.plot(rads*1.E9,concUabss-Uems,label="P_net")
plt.title('Angle-averaged solar-weighted radial power absorption and emission\n '
  + 'of a ' + str(RADCOR*1.E9) + 'nm-' + str(RADCOT*1.E9) +'nm '
  + 'Au-CeO_2 core-shell nanoparticle\n for '
  + 'solar concentration of ' + str(cratio) + ' and T=' + str(temp) + 'K')
plt.xlabel('r [nm]')
plt.ylabel('U_abs [(W m^-3)')
plt.legend()
plt.show()

