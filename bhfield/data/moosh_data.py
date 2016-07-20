#!/usr/bin/env python
import numpy as np

#nfile = "nDataPatsalas950.nk"
#kfile = "kDataPatsalas950.nk"
#outfile = "CeO2_patsalas.nk"
nfile = "nAuBabar.nk"
kfile = "kAuBabar.nk"
outfile = "Au_babar.nk"

ndata = np.loadtxt(nfile,delimiter = ',', skiprows=0,)
kdata = np.loadtxt(kfile,delimiter = ',', skiprows=0)

print ndata
mmooo

#fn = open(nfile)
#fk = open(kfile)

#print fn.readlines()

lams = []
ndat = []
kdata = []
for line in fn.readlines():
  stuff = line.split(',')
  print stuff[0]
  print stuff[0].strip()
  
  lams.append(float(stuff[0].strip()))
  ndat.append(float(stuff[1].strip()))
  
for line in fk.readlines():
  stuff = line.split(',')
  kdat.append(float(stuff[1].strip()))
  
#print mooshed



#mooshed = [ndata[:,0],ndata[:,1],kdata[:,1]]


#np.savetxt(outfile,np.transpose(mooshed))
np.savetxt(outfile,np.transpose((lams,ndat,kdat)))
