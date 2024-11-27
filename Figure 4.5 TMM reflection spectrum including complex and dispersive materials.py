import numpy as np
import csv
from tmm import tmm
import matplotlib.pyplot as plt

# The  files contains a three columns data: wavelength (in um), Real index and Img indexarr = []
arr = []	# reset array
with open('dataset/Si.csv', 'r', encoding='utf-8-sig') as f:	# Franta et al. 2017: n,k 0.0310–310 µm; 300 K
	reader = csv.reader(f)
	next(reader)	# skip the first line (Header)
	for line in reader:
		arr.append(list(map(float, line)))
r1 = list(zip(*arr))

arr = []	# reset array
with open('dataset/SiO2.csv', 'r', encoding='utf-8-sig') as f:	# Rodríguez-de Marcos et al. 2016: n,k 0.03–1.5 µm
	reader = csv.reader(f)
	next(reader)	# skip the first line (Header)
	for line in reader:
		arr.append(list(map(float, line)))
r2 = list(zip(*arr))

arr = []	# reset array
with open('dataset/Cu.csv', 'r', encoding='utf-8-sig') as f:	# Brimhall et al. 2009: n,k 0.0101–0.0421 µm
	reader = csv.reader(f)
	next(reader)	# skip the first line (Header)
	for line in reader:
		arr.append(list(map(float, line)))
r3 = list(zip(*arr))


npts = 500 # Number of plotted points
wavelengths = np.linspace(400, 1000, npts)
nSi = np.interp(wavelengths, np.asarray(r1[0])*1000.0, np.asarray(r1[1]-1j*np.asarray(r1[2])))
nSiO2 = np.interp(wavelengths, np.asarray(r2[0])*1000.0, np.asarray(r2[1]-1j*np.asarray(r2[2])))
nCu = np.interp(wavelengths, np.asarray(r3[0])*1000.0, np.asarray(r3[1]-1j*np.asarray(r3[2])))

ns = nSi
nf = np.concatenate((nCu, nSi, nSiO2)) # Join a sequence of arrays along an existing axis
thk = np.array([5.0, 10.0, 240.0])	# Thickness of each layer (in nm)

ts = 500	# Substrate thickness - not relevant unless back=1
na = 1.0	# Incident medium (ie. air) index
q = 0.0     # Incident angle
sp = 'TM'   # Incident polarization (only for off-normal incidence)
back = 0    # back=0 will ignore bakcside of the substrate
[T, R, t, r] = tmm(wavelengths, ns, ts, na, nf, thk, q, sp, back)   # Run TMM

#np.savetxt('tmm_dispersion_complex.data', np.transpose([wavelengths, R]))   # Save to file

plt.plot(wavelengths, R)
plt.title("TMM Calculation using three layes of dispersive complex materials\n \
        on a silicon substrate (data from: refractiveindex.info)")
plt.xlabel('Wavelength (nm)')
plt.ylabel('Reflection')
plt.show()
