""" 
Fig.9a Resonant cavity structure with N = 4 and M = 3 designed for resonance at 1550 nm 
using VO2 in the cold state and TiO2 phasecompensator. 
(a) Absorption spectrum
plot both cVO2& hVO2
"""

import numpy as np
import csv
from tmm import tmm
import matplotlib.pyplot as plt

# The  files contains a three columns data: wavelength (in um), Real index and Img indexarr = []
arr = []	# reset array
with open('dataset/TiO2.csv', 'r', encoding='utf-8-sig') as f:	#TiO2=2.2899; Zhukovsky et al. 2015: Thin film; n 0.211–1.69 µm
	reader = csv.reader(f)
	next(reader)	# skip the first line (Header)
	for line in reader:
		arr.append(list(map(float, line)))
r1 = list(zip(*arr))

arr = []	# reset array
with open('dataset/SiO2.csv', 'r', encoding='utf-8-sig') as f:	# SiO2=1.4657; Lemarchand 2013: n,k 0.25–2.5 µm
	reader = csv.reader(f)
	next(reader)	# skip the first line (Header)
	for line in reader:
		arr.append(list(map(float, line)))
r2 = list(zip(*arr))

arr = []	# reset array
with open('dataset/cVO2.csv', 'r', encoding='utf-8-sig') as f:	# cVO2=3.1462+0.3388j; Beaini et al. 2020: n,k 0.5–25 µm; 25 °C
	reader = csv.reader(f)
	next(reader)	# skip the first line (Header)
	for line in reader:
		arr.append(list(map(float, line)))
r3c = list(zip(*arr))

arr = []	# reset array
with open('dataset/hVO2.csv', 'r', encoding='utf-8-sig') as f:	# hVO2=1.5729+2.7388j; Oguntoye et al. 2023 n,k 0.21–2.5 µm; 80 °C
	reader = csv.reader(f)
	next(reader)	# skip the first line (Header)
	for line in reader:
		arr.append(list(map(float, line)))
r3h = list(zip(*arr))

arr = []	# reset array
with open('dataset/air.csv', 'r', encoding='utf-8-sig') as f:	# na=1.0; Mathar 2007: n 1.3–2.5 µm
	reader = csv.reader(f)
	next(reader)	# skip the first line (Header)
	for line in reader:
		arr.append(list(map(float, line)))
r4 = list(zip(*arr))


npts = 100 # Number of plotted points
wavelengths = np.linspace(1500, 1600, npts)
nTiO2 = np.interp(wavelengths, np.asarray(r1[0])*1000.0, np.asarray(r1[1]-1j*np.asarray(r1[2])))	# np.interpolate(range, wl(um to nm), index n-k)
nSiO2 = np.interp(wavelengths, np.asarray(r2[0])*1000.0, np.asarray(r2[1]-1j*np.asarray(r2[2])))
ncVO2 = np.interp(wavelengths, np.asarray(r3c[0])*1000.0, np.asarray(r3c[1]-1j*np.asarray(r3c[2])))
nhVO2 = np.interp(wavelengths, np.asarray(r3h[0])*1000.0, np.asarray(r3h[1]-1j*np.asarray(r3h[2])))
nAir = np.interp(wavelengths, np.asarray(r4[0])*1000.0, np.asarray(r4[1]-1j*np.asarray(r4[2])))

ns = nAir	# na=1.0; Mathar 2007: n 1.3–2.5 µm
nfc = np.concatenate((nTiO2, nTiO2,nSiO2,nTiO2, nTiO2,nSiO2,nTiO2, nTiO2,nSiO2,nTiO2, nTiO2,nSiO2,nTiO2,\
	 ncVO2,nTiO2,ncVO2,\
	 nTiO2,nSiO2,nTiO2, nTiO2,nSiO2,nTiO2, nTiO2,nSiO2,nTiO2, nTiO2)) # Join a sequence of arrays along an existing axis
nfh = np.concatenate((nTiO2, nTiO2,nSiO2,nTiO2, nTiO2,nSiO2,nTiO2, nTiO2,nSiO2,nTiO2, nTiO2,nSiO2,nTiO2,\
	 nhVO2,nTiO2,nhVO2,\
	 nTiO2,nSiO2,nTiO2, nTiO2,nSiO2,nTiO2, nTiO2,nSiO2,nTiO2, nTiO2)) # Join a sequence of arrays along an existing axis
wl = 1550.0	#reference wavelength
tTiO2 = wl/(4*2.2899) # Quarterwave thick
tSiO2 = wl/(4*1.4657) # Quarterwave thick
thk = np.array([tTiO2/2, tTiO2/2,tSiO2,tTiO2/2, tTiO2/2,tSiO2,tTiO2/2, tTiO2/2,tSiO2,tTiO2/2, tTiO2/2,tSiO2,tTiO2/2,\
	4.92,155.39,4.92,\
	tTiO2/2,tSiO2,tTiO2/2, tTiO2/2,tSiO2,tTiO2/2, tTiO2/2,tSiO2,tTiO2/2, tTiO2/2])	# Thickness of each layer (in nm)

ts = 500	# Substrate thickness - not relevant unless back=1
na = 1.0	# Incident medium (ie. air) index
q = 0.0     # Incident angle
sp = 'TM'   # Incident polarization (only for off-normal incidence)
back = 0    # back=0 will ignore bakcside of the substrate
[Tc, Rc, tc, rc] = tmm(wavelengths, ns, ts, na, nfc, thk, q, sp, back)   # Run TMM
[Th, Rh, th, rh] = tmm(wavelengths, ns, ts, na, nfh, thk, q, sp, back)   # Run TMM

# Calculate Absorption = 1.0 - (Transmission + Reflection)
Ac = np.full(npts, 1.0)
Ah = Ac
Ac = Ac - (Tc + Rc)
Ah = Ah - (Th + Rh)

#np.savetxt('tmm_dispersion_complex.data', np.transpose([wavelengths, R]))   # Save to file



# fig, ax1 = plt.subplots()
# ax2 = ax1.twinx()

plt.plot(wavelengths, Ac, color="purple")
plt.plot(wavelengths, Ah, color="green")
# ax1.plot(wavelengths, Rc, color="purple", label="Reflection")
# ax1.plot(wavelengths, Rh, color="purple", linestyle="dashed")
# ax2.plot(wavelengths, Tc, color="green", label="Transmission")
# ax2.plot(wavelengths, Th, color="green", linestyle="dashed")

plt.title("Resonant cavity structure with N = 4 and M = 3 designed for resonance\n"
	"at 1550 nm using VO2 in the cold state and TiO2 phase compensator.", fontsize=12)
plt.xlabel('Wavelength (nm)')
plt.ylabel('Absorption')
# ax1.set_ylabel('Reflection', color="purple")
# ax2.set_ylabel('Transmission', color="green")

plt.xlim(1500, 1600)
plt.ylim(0, 1.0)
# ax1.set_ylim(0.0, 1.0)	# Reflection
# ax2.set_ylim(0.0, 1.0)	# Transmission
plt.show()