#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import bhcoat_pyed
from scipy.interpolate import interp1d
from scipy.integrate import simps

#NOTE:  use >>f2py -c bhcoat.f -m bhcoat_pyed
# to get the bhcoat_pyed module

h = 6.626e-34
c = 3.0e+8
k = 1.38e-23

def planck(wav, T):
    a = 2.0*h*c**2
    b = h*c/(wav*k*T)
    intensity = a/ ( (wav**5) * (np.exp(b) - 1.0) )
    return intensity/1.0E9


solarDat = np.genfromtxt("./data/ASTMG173.csv",delimiter=",",skip_header = 2)

solI = interp1d(solarDat[:,0]*1E-9,solarDat[:,3],bounds_error=False,fill_value=0)


lammin = 280.E-9
lammax = 30000.E-9
lams = np.linspace(lammin,lammax,2000)

plt.figure()


plt.loglog(lams*1E9,solI(lams),'-k')
#plt.plot(lams*1E9,solI(lams)/np.max(solI(lams)))


#show emission spectrum
temp = 700
#plt.plot(lams*1E9,planck(lams,temp)/np.max(planck(lams,temp)))
plt.loglog(lams*1E9,planck(lams,temp))
plt.loglog(lams*1E9,planck(lams,500))
plt.loglog(lams*1E9,planck(lams,80bla0))
plt.loglog(lams*1E9,planck(lams,1000))
#plt.loglog(lams*1E9,planck(lams,5777)/15000.)
#print planck(lams,800)

plt.xlabel('wavelength [nm]')
plt.ylabel('Intensity [W m^-2 nm^-1]')
plt.xlim(lammin*1E9,lammax*1E9)
plt.ylim(0.01)
plt.show()
