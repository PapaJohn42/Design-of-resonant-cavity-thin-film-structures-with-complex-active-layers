import numpy as np
from tmm import tmm
import matplotlib.pyplot as plt

#The  files contains a three columns data: wavelength (in um), Real index and Img indexarr = []
arr = []
f = open('Refractive Index data/Si.data', 'r', encoding='utf-8-sig')
next(f)
for line in f:
	arr.append(map(float, line.split()))
r1 = list(zip(*arr))
f.close

arr = []	#reset array
f = open('Refractive Index data/SiO2.data', 'r', encoding='utf-8-sig')
next(f)
for line in f:
	arr.append(map(float, line.split()))
r2 = list(zip(*arr))
f.close

arr = []	#reset array
f = open('Refractive Index data/Cu.data', 'r', encoding='utf-8-sig')
next(f)
for line in f:
	arr.append(map(float, line.split()))
r3 = list(zip(*arr))
f.close

nlambda = 1000
wavelengths = np.linspace(400, 1000, nlambda)
nsi = np.interp(wavelengths, np.asarray(r1[0])*1000.0, np.asarray(r1[1]-1j*np.asarray(r1[2])))
nsio2 = np.interp(wavelengths, np.asarray(r2[0])*1000.0, np.asarray(r2[1]-1j*np.asarray(r2[2])))
ncu = np.interp(wavelengths, np.asarray(r3[0])*1000.0, np.asarray(r3[1]-1j*np.asarray(r3[2])))
ns = nsi
nf = np.concatenate((ncu, nsi, nsio2))
thk = np.array([5.0, 10.0, 240.0])	#Thickness of each layer (in nm)

ts = 500	#Substrate thickness - not relevant unless back=1
na = 1.0	#Incident medium (ie. air) index
q = 0.0     #Incident angle
sp = 'TM'   #Incident polarization (only for off-normal incidence)
back = 0    #back=0 will ignore bakcside of the substrate
[T, R, t, r] = tmm(wavelengths, ns, ts, na, nf, thk, q, sp, back)   #Run TMM

np.savetxt('tmm_dispersion_complex.data', np.transpose([wavelengths, R]))   #Save to file
plt.plot(wavelengths, R)
plt.title("TMM Calculation using three layes of dispersive complex materials\n \
        on a silicon substrate (data from: filmetrics.com)")
plt.xlabel('Wavelength (nm)')
plt.ylabel('Reflection')
plt.show()
