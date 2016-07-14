#!/usr/bin/env python
import numpy as np
from scipy.integrate import quad, dblquad, tplquad, nquad
from scipy.interpolate import griddata, LinearNDInterpolator
from math import isnan
import matplotlib.pyplot as plt

data = np.loadtxt("./U_0allf.dat",skiprows=2)

radshell = 0.06  
radcore = radshell/2.

#global Uabss
Uabss = data[:,3]

#global pts
pts = data[:,0:3]

ct = 0
#for r in data[:,0]:
#  if r >= radcore:
#    cut1 = ct - 1
#    cutoff = ct -1
#    print "found: " + str(cut1)
#    print pts[cut1,0]
#    print pts[cut1+1,0]
#    break
#  ct += 1

#print cutoff
#print data[cutoff,0]
#print data[cutoff+1,0]

#pts_core = pts[:cutoff,:]
#print pts_core
#pts_shell = pts[cutoff+1:,:]
#print pts_shell
#Uabss_core = Uabss[:cutoff]
#Uabss_shell = Uabss[cutoff+1:]
#print Uabss_core
#print Uabss_shell


# make an interpolation of the data so it can be called like a function
# for a quadrature routine
#testout  = griddata(pts,Uabss,(0.02,0,0),method='linear')
#print testout
#testout  = griddata(pts_core,Uabss_core,(0.02,0,0),method='linear')
#print testout

##THIS IS NOT EVEN CLOSE!!!:
interp = LinearNDInterpolator(pts,Uabss)
print 'interpolation in phi'
print interp(0.27273E-01,  0.18480E+01,  0.36960E+01)
print interp(0.27273E-01,  0.18480E+01,  0.375E+01)
print interp(0.27273E-01,  0.18480E+01,  0.38808E+01)
print 'interpolation in theta'
print interp(0.27273E-01,  0.18480E+01,  0.36960E+01)
print interp(0.27273E-01,  0.19+01,  0.36960E+01)
print interp(0.27273E-01,  0.20328E+01,  0.36960E+01)
print 'interpolation in r'
print interp(0.27273E-01,  0.18480E+01,  0.36960E+01)
print interp(0.29E-01,  0.18480E+01,  0.36960E+01)
print interp(0.32727E-01,  0.18480E+01,  0.36960E+01)



#print griddata(pts,Uabss,(0.27273E-01,  0.18480E+01,  0.36960E+01),method='linear')
#print griddata(pts,Uabss,(0.27273E-01,  0.18480E+01,  0.377E+01),method='linear')
#print griddata(pts,Uabss,(0.27273E-01,  0.18480E+01,  0.38808E+01),method='linear')


npts = np.max(Uabss.shape)

pi = np.pi

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


radt = np.linspace(0,0.06,100)
#thet = np.linspace(0,np.pi,100)

Uts = interp(radt,np.pi,0)
#Uts = interp(0.05,thet,0.)


print pts[1:5,2]


plt.figure()
#plt.plot(radt,Uts)
#plt.plot(thet,Uts)
plt.plot(pts[1:5,2],Uabss[1:5])
plt.plot(np.linspace(0,0.62832E+01,10),interp(0.015,0.0,np.linspace(0,0.62832E+01,10)))
plt.show()


dphi, dtheta = pi/250.0, pi/250.0
[phi,theta] = np.mgrid[0:pi+dphi*1.5:dphi,0:2*pi+dtheta*1.5:dtheta]
m0 = 4; m1 = 3; m2 = 2; m3 = 3; m4 = 6; m5 = 2; m6 = 6; m7 = 4;
r = np.sin(m0*phi)**m1 + np.cos(m2*phi)**m3 + np.sin(m4*theta)**m5 + np.cos(m6*theta)**m7
x = r*np.sin(phi)*np.cos(theta)
y = r*np.cos(phi)
z = r*np.sin(phi)*np.sin(theta)

# View it.
from mayavi import mlab
s = mlab.mesh(x, y, z)
mlab.show()

moop

def diff_volume(p,t,r):
  return r**2*np.sin(t)

  out = griddata(pts,Uabss,(r,t,p),method='linear')*r**2*np.sin(t)
  print out
  if isnan(out):
    print 'oh fuck'
    die
  else:
    return out
  
def new_diff_Uabs(p,t,r):
  out = interp(r,t,p)*r**2.*np.sin(t)
  #out = r**(1.22)*r**2.*np.sin(t)
  #print out
  if isnan(out):
    print 'oh fuck'
    die
  else:
    return out

volume = tplquad(diff_volume, r1, r2, lambda r:   t1, lambda r:   t2,
                                      lambda r,t: p1, lambda r,t: p2)[0]
print 'volume check (' + str(4./3.*np.pi*(radcore)**3.) + '): ' + str(volume)


Qabs_core = tplquad(new_diff_Uabs, r1, r2, lambda r:   t1, lambda r:   t2,
                                      lambda r,t: p1, lambda r,t: p2,
                                      epsabs=1.49e-04, epsrel=1.49e-04)[0]
print 'Qabs_core: ' + str(Qabs_core)
#Qabs_shell = tplquad(diff_Uabs, r2, r3, lambda r:   t1, lambda r:   t2,
#                                      lambda r,t: p1, lambda r,t: p2)[0]
#print 'Qabs_shell: ' + str(Qabs_shell)
#print 'Qabs (0.65...): ' + str(Qabs_core + Qabs_shell)
#Qabs = tplquad(diff_Uabs, r1, r3, lambda r:   t1, lambda r:   t2,
#                                        lambda r,t: p1, lambda r,t: p2)[0]

#print Qabs
