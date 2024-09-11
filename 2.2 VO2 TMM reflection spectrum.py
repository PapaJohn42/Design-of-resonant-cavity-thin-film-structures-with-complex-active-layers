##For use in 'Design of resonant cavity thin...' only.
##Recreating Figure 6, specifically N4 M2 (Solid lines).
##coldVO2
#dataset: 
#nTiO2; Bond 1965: n(o) 0.45–2.4 µm
#nSiO2; Malitson 1965: n 0.21–6.7 µm
#nVO2; Oguntoye et al. 2023: n,k 0.21–2.5 µm; 20 °C
#nAir; Mathar 2007: n 1.3–2.5 µm

#To do: -

import numpy as np
from tmm import tmm
import matplotlib.pyplot as plt

#These  files contains a three columns data: wavelength (in um), Real index and Img index
arr = []
f = open('Refractive Index data/TiO2.data', 'r', encoding='utf-8-sig')	#Bond 1965: n(o) 0.45–2.4 µm
next(f)
for line in f:
	arr.append(map(float, line.split()))
rH = list(zip(*arr))
f.close

arr = []	#reset array
f = open('Refractive Index data/SiO2.data', 'r', encoding='utf-8-sig')	#Malitson 1965: n 0.21–6.7 µm
next(f)
for line in f:
	arr.append(map(float, line.split()))
rL = list(zip(*arr))
f.close

arr = []	#reset array
f = open('Refractive Index data/cVO2.data', 'r', encoding='utf-8-sig')	#Oguntoye et al. 2023: n,k 0.21–2.5 µm; 20 °C
next(f)
for line in f:
	arr.append(map(float, line.split()))
rC = list(zip(*arr))
f.close

rD = rH

arr = []	#reset array
f = open('Refractive Index data/Air.data', 'r', encoding='utf-8-sig')
next(f)
for line in f:
	arr.append(map(float, line.split()))	#Mathar 2007: n 1.3–2.5 µm
ra = list(zip(*arr))
f.close

nlambda = 100
wavelengths = np.linspace(1500, 1600, nlambda)	#[input]
nH = np.interp(wavelengths, np.asarray(rH[0])*1000.0, np.asarray(rH[1]-1j*np.asarray(rH[2])))
nL = np.interp(wavelengths, np.asarray(rL[0])*1000.0, np.asarray(rL[1]-1j*np.asarray(rL[2])))
nC = np.interp(wavelengths, np.asarray(rC[0])*1000.0, np.asarray(rC[1]-1j*np.asarray(rC[2])))
nD = np.interp(wavelengths, np.asarray(rD[0])*1000.0, np.asarray(rD[1]-1j*np.asarray(rD[2])))
ns = np.interp(wavelengths, np.asarray(ra[0])*1000.0, np.asarray(ra[1]-1j*np.asarray(ra[2])))
nf = np.concatenate((nH, nH,nL,nH, nH,nL,nH, nH,nL,nH, nC,nD,nC, nH,nL,nH, nH,nL,nH, nH,nL,nH, nH,nL,nH, nH))
thk_std = 1550.0/4

#thk = np.array([thk_std/2, thk_std,thk_std,thk_std, thk_std,thk_std,thk_std, thk_std,thk_std,thk_std, 0.9,152.1,0.9, thk_std,thk_std,thk_std, thk_std,thk_std,thk_std, thk_std,thk_std,thk_std, thk_std,thk_std,thk_std, thk_std/2 ])	#Thickness of each layer (in nm)
thk = np.array([thk_std/2, thk_std,thk_std,thk_std, thk_std,thk_std,thk_std, thk_std,thk_std,thk_std, 4.38,148.32,4.38, thk_std,thk_std,thk_std, thk_std,thk_std,thk_std, thk_std,thk_std,thk_std, thk_std,thk_std,thk_std, thk_std/2 ])	#Thickness of each layer (in nm)

ts = 500	#Substrate thickness - not relevant unless back=1
na = 1.0	#Incident medium (ie. air) index
q = 0.0     #Incident angle
sp = 'TM'   #Incident polarization (only for off-normal incidence)
back = 0    #back=0 will ignore bakcside of the substrate
[T, R, t, r] = tmm(wavelengths, ns, ts, na, nf, thk, q, sp, back)   #Run TMM

#np.savetxt('tmm_dispersion_complex.data', np.transpose([wavelengths, R]))   #Save to file
#Plot the reflection spectrum
plt.plot(wavelengths, R)
plt.title("Fig.7. coldVO2 Resonant cavity structure with N=4 M=3. \n \
        Designed to resonant at 1550nm")
plt.xlabel('Wavelength (nm)')
plt.ylabel('Reflection')
plt.show()