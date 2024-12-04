""" 
Limiting Optical Diodes Enabled by the Phase Transition of Vanadium Dioxide
Chenghao Wan, et al.
The reference wavelength is 1320 nm
cVO2 = 25°C, hVO2 = 100°C
Au = 10nm, VO2 = 100nm, C-plane Sapphire substrate
*Problem with choosing each layer datasets
"""

import numpy as np
import csv
from tmm import tmm
import matplotlib.pyplot as plt

# The  files contains a three columns data: wavelength (in um), Real index and Img indexarr = []
arr = []	# reset array
with open('dataset/Au.csv', 'r', encoding='utf-8-sig') as f:	# Yakubovsky et al. 2019: 9-nm film; n,k 0.30–3.3 µm
	reader = csv.reader(f)
	next(reader)	# skip the first line (Header)
	for line in reader:
		arr.append(list(map(float, line)))
r1 = list(zip(*arr))

arr = []	# reset array
with open('dataset/cVO2.csv', 'r', encoding='utf-8-sig') as f:	# Beaini et al. 2020: n,k 0.5–25 µm; 25 °C
	reader = csv.reader(f)
	next(reader)	# skip the first line (Header)
	for line in reader:
		arr.append(list(map(float, line)))
r2c = list(zip(*arr))

arr = []	# reset array
with open('dataset/hVO2-100.csv', 'r', encoding='utf-8-sig') as f:	# Beaini et al. 2020: n,k 0.5–25 µm; 100 °C
	reader = csv.reader(f)
	next(reader)	# skip the first line (Header)
	for line in reader:
		arr.append(list(map(float, line)))
r2h = list(zip(*arr))

arr = []	# reset array
with open('dataset/Al2O3.csv', 'r', encoding='utf-8-sig') as f:	# Querry 1985: α-Al2O3 (Sapphire); n,k(o) 0.21–55.6 µm
	reader = csv.reader(f)
	next(reader)	# skip the first line (Header)
	for line in reader:
		arr.append(list(map(float, line)))
r3 = list(zip(*arr))


npts = 2900 # Number of plotted points
wavelengths = np.linspace(500, 3300, npts)
nAu = np.interp(wavelengths, np.asarray(r1[0])*1000.0, np.asarray(r1[1]-1j*np.asarray(r1[2])))
ncVO2 = np.interp(wavelengths, np.asarray(r2c[0])*1000.0, np.asarray(r2c[1]-1j*np.asarray(r2c[2])))
nhVO2 = np.interp(wavelengths, np.asarray(r2h[0])*1000.0, np.asarray(r2h[1]-1j*np.asarray(r2h[2])))
nAl2O3 = np.interp(wavelengths, np.asarray(r3[0])*1000.0, np.asarray(r3[1]-1j*np.asarray(r3[2])))

ns = nAl2O3	# 
nfc = np.concatenate((ncVO2, nAu)) # Join a sequence of arrays along an existing axis
nfh = np.concatenate((nhVO2, nAu)) # Join a sequence of arrays along an existing axis
wl = 1320.0	#reference wavelength
tAu = 10 # Quarterwave thick
tVO2 = 100 # Quarterwave thick
thk = np.array([tVO2, tAu])	# Thickness of each layer (in nm)

ts = 500 * 1000	# Substrate thickness - not relevant unless back=1
na = 1.0	# Incident medium (ie. air) index
q = 0.0     # Incident angle
sp = 'TM'   # Incident polarization (only for off-normal incidence)
back = 0    # back=0 will ignore bakcside of the substrate
[Tc, Rc, tc, rc] = tmm(wavelengths, ns, ts, na, nfc, thk, q, sp, back)   # Run TMM
[Th, Rh, th, rh] = tmm(wavelengths, ns, ts, na, nfh, thk, q, sp, back)   # Run TMM

#np.savetxt('tmm_dispersion_complex.data', np.transpose([wavelengths, R]))   # Save to file


plt.plot(wavelengths, Tc, color="blue", label="VO2 in the insulating phase")
plt.plot(wavelengths, Th, color="red", label="VO2 in the metallic phase")

plt.title("Limiting Optical Diodes Enabled by the Phase Transition of Vanadium Dioxide.\n"
	"Chenghao Wan, et al.\n"
	"Au = 10nm, VO2 = 100nm, Thick c-plane Sapphire", fontsize=11)
plt.xlabel('Wavelength (nm)')
plt.ylabel('Transmission', color="black")
plt.legend(loc='upper right', fontsize=8)

plt.xlim(1300, 1350)
plt.ylim(0.0, 0.35)
plt.show()