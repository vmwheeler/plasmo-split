#!/usr/bin/env python
import xalglib
import numpy as np
from scipy.integrate import tplquad, nquad
from math import isnan

plotit = False
splineit = False
integrateit = True

radshell = 0.06
radcore = radshell/2.

if plotit:
  if splineit:
    data_cart = np.loadtxt("./U_0allf_cart.dat",skiprows=2)
    print "making spline... "
    sz = np.max(np.size(data_cart[:,0]))
    print sz
    print data_cart
    print len(data_cart[:,0].tolist())
    splinefit = xalglib.spline3dbuildtrilinearv(data_cart[:,0].tolist(), sz, 
                                                data_cart[:,1].tolist(), sz, 
                                                data_cart[:,2].tolist(), sz, 
                                                data_cart[:,3].tolist(), 1)
    print "...done"
  else:
    print "creating interpolant for cartesian data... "
    data_cart = np.loadtxt("./U_0allf_cart.dat",skiprows=2)
    # make an empty 3 dimensional scalar interpolation model
    model_cart = xalglib.rbfcreate(3, 1)
    xalglib.rbfsetpoints(model_cart, data_cart.tolist())
    xalglib.rbfsetalgomultilayer(model_cart, 0.03, 4, 1.0e-3)
    rep = xalglib.rbfbuildmodel(model_cart)
    print " ...finished."
  
  
  
if integrateit:
  data = np.loadtxt("./U_0allf_rad.dat",skiprows=2)
  print "preview of data sampling... "
  from mayavi import mlab
  xs = data[:,0]*np.sin(data[:,1])*np.cos(data[:,2])
  ys = data[:,0]*np.sin(data[:,1])*np.sin(data[:,2])
  zs = data[:,0]*np.cos(data[:,1])
  mlab.points3d(xs, ys, np.zeros(np.shape(xs)),mode='point')
  mlab.show()
  
  
  # make an empty 3 dimensional scalar interpolation model
  model = xalglib.rbfcreate(3, 1)

  #add dataset (has to be list, not np array)
  xalglib.rbfsetpoints(model, data.tolist())

  

  # Now build the model
  print "creating interpolant for spherical data... "
  #xalglib.rbfsetalgoqnn(model)
  xalglib.rbfsetalgomultilayer(model, 0.03, 1, 1.0e-3)
  rep = xalglib.rbfbuildmodel(model)
  #print(rep.terminationtype) # expected 1

  print " ...finished."

  #v = xalglib.rbfcalc3(model,radcore/10.,0,0)
  #print(v) # expected 2.500

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
  out = xalglib.rbfcalc3(model,r,p,t)*r**2.*np.sin(t)
  #out = r**(1.22)*r**2.*np.sin(t)
 
  if isnan(out):
    print 'oh fuck'
    print r, t, p
    die
  else:
    return out

def integrate_that_shit():
  volume = tplquad(diff_volume, r1, r3, lambda r:   t1, lambda r:   t2,
                                        lambda r,t: p1, lambda r,t: p2)[0]
  print 'volume check (' + str(4./3.*np.pi*(r3)**3.) + '): ' + str(volume)

  print 'now the nasty integral...'
  Uabs_T = tplquad(new_diff_Uabs, r1, r3, lambda r:   t1, lambda r:   t2,
                                        lambda r,t: p1, lambda r,t: p2)[0]

  print '... no problem!'

  eps0 = 8.854187E-12
  mu0 = 4.*np.pi*1.E-7
  Ap = np.pi*radshell*radshell
  Qabs_corrected = Uabs_T/Ap/(1./2.*np.sqrt(eps0/mu0))
  print 'my calc divided by sqrt(eps0/mu0) = ' + str(Qabs_corrected*1.E-6)

def plot_that_shit():
  npts = 100
  spcing = 2*radshell/npts
  xx,yy = np.mgrid[-radshell:radshell:spcing,-radshell:radshell:spcing]
  preds = np.zeros(np.shape(xx))
  print 'evaluating interpolant on a cut of the data... '
  for j in range(len(xx)):
    for i in range(len(yy)):
      if splineit:
        val = xalglib.spline3dcalcv(splinefit, xx[i,j],yy[i,j],0.)
      else:
        val = xalglib.rbfcalc3(model_cart,xx[i,j],yy[i,j],0.)
      print val
      #preds[i,j] = val 
  
  print '... finished'

  print 'plotting...'
  from mayavi import mlab
  
  mlab.surf(xx,yy,preds,warp_scale='auto')
  mlab.show()
  return


if integrateit:
  integrate_that_shit()

if plotit:
  plot_that_shit()
