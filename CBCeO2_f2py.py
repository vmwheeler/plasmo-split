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
    return intensity*1.0E-9

#print bhcoat_pyed.bhcoat.__doc__

nCeO2Dat = np.loadtxt("./data/nDataPatsalas950.txt",delimiter=",")
kCeO2Dat = np.loadtxt("./data/kDataPatsalas950.txt",delimiter=",")
#nAuDat = np.genfromtxt("./data/nAuJohnson.txt",delimiter=",",skip_header = 1)
#kAuDat = np.genfromtxt("./data/kAuJohnson.txt",delimiter=",",skip_header = 1)
#nAuDat = np.genfromtxt("./data/nAuBabar.txt",delimiter=",",skip_header = 1)
#kAuDat = np.genfromtxt("./data/kAuBabar.txt",delimiter=",",skip_header = 1)
nAuDat = np.genfromtxt("./data/nSiC.txt",delimiter=",",skip_header = 1)
kAuDat = np.genfromtxt("./data/kSiC.txt",delimiter=",",skip_header = 1)
solarDat = np.genfromtxt("./data/ASTMG173.csv",delimiter=",",skip_header = 2)

nCeO2 = interp1d(nCeO2Dat[:,0]*1E-9,nCeO2Dat[:,1])
kCeO2 = interp1d(kCeO2Dat[:,0]*1E-9,kCeO2Dat[:,1])
nAu = interp1d(nAuDat[:,0]*1E-6,nAuDat[:,1])
kAu = interp1d(kAuDat[:,0]*1E-6,kAuDat[:,1])
solI = interp1d(solarDat[:,0]*1E-9,solarDat[:,3])

tlams = np.linspace(280E-9,1000E-9,100)


#plt.figure()
#plt.plot(tlams,nCeO2(tlams),tlams,kCeO2(tlams),tlams,nAu(tlams),tlams,kAu(tlams))
#plt.show()

rcore = 80E-9
rrefvac = 1.0

lammin = 280.E-9
lammax = 1200.E-9
nlams = 1200
lams = np.linspace(lammin,lammax,nlams)

rshells = np.linspace(rcore,rcore*2,5)
rrefrefs = np.linspace(1.0,2.35,5)

plt.figure()
plt.plot(lams*1E9,solI(lams)/np.max(solI(lams)),linestyle='--',color='goldenrod')

#show emission spectrum
temp = 1300
#plt.plot(lams*1E9,planck(lams,temp)/np.max(planck(lams,temp)))
plt.plot(lams*1E9,planck(lams,temp)/np.max(solI(lams)),'--r')
#print planck(lams,800)

qabsmax = 0

datout = []
filehead = "wavelength[nm] "
datout.append(np.transpose(lams*1.E9))

for rshell in rshells:
    cabss = []
    qabss = []
    for lam in lams:
        xcore = 2.*np.pi*rcore*rrefvac/lam
        xshell= 2.*np.pi*rshell*rrefvac/lam
        rrefcore = nAu(lam) + kAu(lam)*1j
        rrefshell = nCeO2(lam) + kCeO2(lam)*1j
        [qext, qsca, qback, gsca] = bhcoat_pyed.bhcoat(xcore, xshell, rrefcore, rrefshell)
        cabss.append((qext-qsca)*np.pi*rshell**2.)
        qabss.append(qext-qsca)
  
  #calculate the integral of Qabs over the spectrum
    newmax = np.max(qabss)
    if newmax > qabsmax:
        qabsmax = newmax
    qabsint = simps(qabss,lams)*1.E6
    print("integral over Qabs for rshell = " + str(rshell) + " = " + str(qabsint)) 
    datout.append(qabss)
    filehead += "Qabs--" + str(round((rshell-rcore)*1.E9)) + "nm[1] "
    plt.plot(lams*1E9,qabss, label="t="+str(round((rshell-rcore)*1E9,1))+'nm')


#rshell = max(rshells)
rshell = 45.E-9

cabss = []
qabss = []
for lam in lams:
	xcore = 2.*np.pi*rcore*rrefvac/lam
	xshell= 2.*np.pi*rshell*rrefvac/lam
	rrefcore = nCeO2(lam) + kCeO2(lam)*1j
	rrefshell = nCeO2(lam) + kCeO2(lam)*1j
	[qext, qsca, qback, gsca] = bhcoat_pyed.bhcoat(xcore,xshell,rrefshell,rrefshell)
	cabss.append((qext-qsca)*np.pi*rshell**2.)
	qabss.append(qext-qsca)

datout.append(qabss)
filehead += "Qabs--" + "CeO2 "
qabsint = simps(qabss,lams)*1.E6
print("integral over Qabs for ceria only rad = " + str(rshell) + " = " + str(qabsint))

plt.plot(lams*1E9,qabss,linestyle=':',color='black', label='ceria only')


plt.xlabel('wavelength [nm]')
plt.ylabel('Q_abs [1]')
plt.xlim(lammin*1E9,lammax*1E9)
#plt.ylim(0,1.1*qabsmax)
#plt.title('core radius = '+str(rcore*1E9)+'nm')
plt.legend()
plt.show()

#output data

datout.append(solI(lams))
datout.append(planck(lams,temp))
filehead += "Solar Irradiance[W/m^2/nm] " + "BB intensity@" + str(temp)+"K[W/m^2/nm]"

datout = np.array(datout)
np.savetxt('forLance.dat',np.transpose(datout),delimiter=' ',header=filehead)
