#!/usr/bin/env python
import numpy as np
from scipy.integrate import simps
import matplotlib.pyplot as plt

f = open("./QabsVlam.dat")

data = np.loadtxt("./QabsVlam.dat",skiprows=3)

#check that this integrates to proper absorption cross section
lams = data[:,0] # in micrometers
Qabss = data[:,1]


print lams

integval = simps(Qabss,lams)
print "oh please: " + str(integval)

plt.figure()
plt.plot(lams*1.0E3,Qabss)
plt.xlabel('wavelength [nm]')
plt.ylabel('Q_abs [1]')
plt.show()

