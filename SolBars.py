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

class NCNP:
    def __init__(self,corename,ncore,kcore,shellname,nshell,kshell,rcore,tshell):
        ncore = interp1d(ncore[:,0]*1E-6,ncore[:,1],fill_value='extrapolate')
        kcore = interp1d(kcore[:,0]*1E-6,kcore[:,1],fill_value='extrapolate')
        nshell = interp1d(nshell[:,0]*1E-6,nshell[:,1],fill_value='extrapolate')
        kshell = interp1d(kshell[:,0]*1E-6,kshell[:,1],fill_value='extrapolate')
        self.ncore = ncore
        self.kcore = kcore
        self.nshell = nshell
        self.kshell = kshell
        self.corename = corename
        self.shellname = shellname
        self.rcore = rcore
        self.tshell = tshell
        self.rshell = rcore+tshell
        self.vol = 4./3.*np.pi*(rcore+tshell)**3.
        self.cs = np.pi*(rcore+tshell)**2.

nCeO2 = np.loadtxt("./data/nDataPatsalas950.txt",delimiter=",")
nCeO2[:,0] = nCeO2[:,0]*1E-3
kCeO2 = np.loadtxt("./data/kDataPatsalas950.txt",delimiter=",")
kCeO2[:,0] = kCeO2[:,0]*1E-3

nAu = np.genfromtxt("./data/nAuBabar.txt",delimiter=",",skip_header = 1)
kAu = np.genfromtxt("./data/kAuBabar.txt",delimiter=",",skip_header = 1)
nSiC = np.genfromtxt("./data/nSiC.txt",delimiter=",",skip_header = 1)
kSiC = np.genfromtxt("./data/kSiC.txt",delimiter=",",skip_header = 1)
nB4C = np.genfromtxt("./data/nB4C.txt",delimiter=",",skip_header = 1)
kB4C = np.genfromtxt("./data/kB4C.txt",delimiter=",",skip_header = 1)
nTiN = np.genfromtxt("./data/nTiN.txt",delimiter=",",skip_header = 1)
kTiN = np.genfromtxt("./data/kTiN.txt",delimiter=",",skip_header = 1)
nZrN = np.genfromtxt("./data/nZrN.txt",delimiter=",",skip_header = 1)
kZrN = np.genfromtxt("./data/kZrN.txt",delimiter=",",skip_header = 1)

AuCeO2 = NCNP('Au',nAu,kAu,'CeO2',nCeO2,kCeO2,30E-9,15E-9)
SiCCeO2 = NCNP('SiC',nSiC,kSiC,'CeO2',nCeO2,kCeO2,80E-9,20e-9)
B4CCeO2 = NCNP('B4C',nB4C,kB4C,'CeO2',nCeO2,kCeO2,80E-9,30e-9)
TiNCeO2 = NCNP('TiN',nTiN,kTiN,'CeO2',nCeO2,kCeO2,50E-9,12.5e-9)
ZrNCeO2 = NCNP('ZrN',nZrN,kZrN,'CeO2',nCeO2,kCeO2,35E-9,10E-9)
JustCeO2 = NCNP('CeO2',nCeO2,kCeO2,'CeO2',nCeO2,kCeO2,50E-9,12.5e-9)

candidates = [AuCeO2,SiCCeO2,B4CCeO2,TiNCeO2,ZrNCeO2,JustCeO2]

solarDat = np.genfromtxt("./data/ASTMG173.csv",delimiter=",",skip_header = 2)
solI = interp1d(solarDat[:,0]*1E-9,solarDat[:,3])

tlams = np.linspace(280E-9,1000E-9,100)

which = 4
cand = candidates[which]


plt.figure()
plt.plot(tlams,cand.ncore(tlams),'-r',label='core n')
plt.plot(tlams,cand.kcore(tlams),'--r',label='core k')
plt.plot(tlams,cand.nshell(tlams),'-b',label='shell n')
plt.plot(tlams,cand.kshell(tlams),'--b',label='shell k')
plt.legend()



rrefvac = 1.0

lammin = 280.E-9
lammax = 1200.E-9
nlams = 1200
lams = np.linspace(lammin,lammax,nlams)

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

Ppervols = []
Psolwts = []
CRps = []

for can in candidates:
    cabss = []
    qabss = []
    
    for lam in lams:
        xcore = 2.*np.pi*can.rcore*rrefvac/lam
        xshell= 2.*np.pi*can.rshell*rrefvac/lam
        rrefcore = can.ncore(lam) + can.kcore(lam)*1j
        rrefshell = can.nshell(lam) + can.kshell(lam)*1j
        [qext, qsca, qback, gsca] = bhcoat_pyed.bhcoat(xcore, xshell, rrefcore, rrefshell)
        cabss.append((qext-qsca)*can.cs)
        qabss.append(qext-qsca)
  
    plt.plot(lams*1E9,qabss, label=can.corename + '-' + can.shellname)
    
    #get power absorbed per volume

    Psolwt = simps(cabss*solI(lams),lams)*1.E9
    Phit = simps(can.cs*solI(lams),lams)*1.E9
    Pcheck = simps(solI(lams),lams)*1.e9
    print(Pcheck)
    Psolwts.append(Psolwt)
    Ppervol = Psolwt/can.vol
    Ppervols.append(Ppervol)
    CRp = Psolwt/Pcheck/can.cs
    CRps.append(CRp)

plt.plot(lams*1E9,np.ones(np.shape(lams)),color='darkgray',linestyle='--' , label="max for bulk materials")

#plt.plot(lams*1E9,qabss,linestyle=':',color='black', label='ceria only')


plt.xlabel('wavelength, n')
plt.ylabel('Q_abs, -')
plt.xlim(lammin*1E9,lammax*1E9)
#plt.ylim(0,1.1*qabsmax)
#plt.title('core radius = '+str(rcore*1E9)+'nm')
plt.legend(loc='upper right')

barx = []
for can in candidates:
    barx.append(can.corename)

print(barx)
print(Ppervols)
plt.figure()
plt.bar(barx,Ppervols)
plt.ylabel('P_abs/vol, W/m^3')

'''
plt.figure()
plt.bar(barx,Psolwts)
plt.ylabel('P_{abs,sol}, -')
'''

plt.figure()
plt.bar(barx,CRps)
plt.ylabel('$CR_{particle}$')

plt.show()

