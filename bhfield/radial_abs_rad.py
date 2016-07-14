#!/usr/bin/env python
import numpy as np
from scipy.integrate import tplquad, nquad
from scipy.interpolate import LinearNDInterpolator, Rbf
from math import isnan
#import matplotlib.pyplot as plt
#from mayavi import mlab

data = np.loadtxt("./U_0allf_rad.dat",skiprows=2)

radshell = 0.06
radcore = radshell/2.

#global Uabss
Uabss = data[:,3]

#global pts
pts = data[:,0:3]

ct = 0

interp = Rbf(pts[:,0],pts[:,1],pts[:,2],Uabss)
#interp = LinearNDInterpolator(pts,Uabss)
npts = np.max(Uabss.shape)

print interp(radcore/10.,0,0)

# limits for radius
r1 = 0.
r2 = radcore
r3 = radshell
# limits for theta
t1 = 0.
t2 = np.pi
# limits for phi
p1 = 0.
p2 = 2.*np.pi


def diff_volume(p,t,r):
  return r**2*np.sin(t)
  
def new_diff_Uabs(p,t,r):
  out = interp(r,t,p)*r**2.*np.sin(t)
  #out = r**(1.22)*r**2.*np.sin(t)
  #print out
  #print r, t, p
  if isnan(out):
    print 'oh fuck'
    print r, t, p
    die
  else:
    return out

volume = tplquad(diff_volume, r1, r3, lambda r:   t1, lambda r:   t2,
                                      lambda r,t: p1, lambda r,t: p2)[0]
print 'volume check (' + str(4./3.*np.pi*(r3)**3.) + '): ' + str(volume)


Uabs_T = tplquad(new_diff_Uabs, r1, r3, lambda r:   t1, lambda r:   t2,
                                      lambda r,t: p1, lambda r,t: p2)[0]

eps0 = 8.854187E-12
mu0 = 4.*np.pi*1.E-7
Ap = np.pi*radshell*radshell
Qabs_corrected = Uabs_T/Ap/(1./2.*np.sqrt(eps0/mu0))
print 'my calc divided by sqrt(eps0/mu0) = ' + str(Qabs_corrected*1.E-6)
