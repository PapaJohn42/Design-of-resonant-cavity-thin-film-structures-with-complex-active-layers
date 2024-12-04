""" 
Limiting Optical Diodes Enabled by the Phase Transition of Vanadium Dioxide
Chenghao Wan, et al.
The reference wavelength is 1320 nm
cVO2 = 25°C, hVO2 = 100°C
nH=nD=TiO2, nL=Al2O3, nC=VO2
*Problem with choosing each layer datasets & substituting Au with TiO2 might cause problem?
"""

import numpy as np
import csv
from tmm import tmm
import matplotlib.pyplot as plt

# The  files contains a three columns data: wavelength (in um), Real index and Img indexarr = []
arr = []	# reset array
with open('dataset/TiO2.csv', 'r', encoding='utf-8-sig') as f:	# TiO2=2.2899; Zhukovsky et al. 2015: Thin film; n 0.211–1.69 µm
	reader = csv.reader(f)
	next(reader)	# skip the first line (Header)
	for line in reader:
		arr.append(list(map(float, line)))
r1 = list(zip(*arr))

arr = []	# reset array
with open('dataset/Al2O3.csv', 'r', encoding='utf-8-sig') as f:	# Al2O3; Querry 1985: α-Al2O3 (Sapphire); n,k(o) 0.21–55.6 µm
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


npts = 2900 # Number of plotted points
wavelengths = np.linspace(500, 1700, npts)
nTiO2 = np.interp(wavelengths, np.asarray(r1[0])*1000.0, np.asarray(r1[1]-1j*np.asarray(r1[2])))	# np.interpolate(range, wl(um to nm), index n-k)
nAl2O3 = np.interp(wavelengths, np.asarray(r2[0])*1000.0, np.asarray(r2[1]-1j*np.asarray(r2[2])))
ncVO2 = np.interp(wavelengths, np.asarray(r3c[0])*1000.0, np.asarray(r3c[1]-1j*np.asarray(r3c[2])))
nhVO2 = np.interp(wavelengths, np.asarray(r3h[0])*1000.0, np.asarray(r3h[1]-1j*np.asarray(r3h[2])))
nAir = np.interp(wavelengths, np.asarray(r4[0])*1000.0, np.asarray(r4[1]-1j*np.asarray(r4[2])))

ns = nAir	# na=1.0; Mathar 2007: n 1.3–2.5 µm
nfc = np.concatenate((nTiO2, nTiO2,nAl2O3,nTiO2, nTiO2,nAl2O3,nTiO2,\
	ncVO2,nTiO2,ncVO2,\
	nTiO2,nAl2O3,nTiO2, nTiO2)) # Join a sequence of arrays along an existing axis
nfh = np.concatenate((nTiO2, nTiO2,nAl2O3,nTiO2, nTiO2,nAl2O3,nTiO2,\
	nhVO2,nTiO2,nhVO2,\
	nTiO2,nAl2O3,nTiO2, nTiO2)) # Join a sequence of arrays along an existing axis

wl = 1320.0	#reference wavelength
tTiO2 = wl/(4*2.2964) # Quarterwave thick
tAl2O3 = wl/(4*1.748) # Quarterwave thick
tVO2 = 36.1		# tcVO2	[input]
tD = 56.9		# TiO2 central layer [input]

thk = np.array([tTiO2/2, tTiO2/2,tAl2O3,tTiO2/2, tTiO2/2,tAl2O3,tTiO2/2,\
	tVO2,tD,tVO2,\
	tTiO2/2,tAl2O3,tTiO2/2, tTiO2/2])	# Thickness of each layer (in nm)

ts = 500*1000	# Substrate thickness - not relevant unless back=1
na = 1.0	# Incident medium (ie. air) index
q = 0.0     # Incident angle
sp = 'TM'   # Incident polarization (only for off-normal incidence)
back = 0    # back=0 will ignore bakcside of the substrate
[Tc, Rc, tc, rc] = tmm(wavelengths, ns, ts, na, nfc, thk, q, sp, back)   # Run TMM
[Th, Rh, th, rh] = tmm(wavelengths, ns, ts, na, nfh, thk, q, sp, back)   # Run TMM

#np.savetxt('tmm_dispersion_complex.data', np.transpose([wavelengths, R]))   # Sav to file



# fig, ax1 = plt.subplots()
# ax2 = ax1.twinx()

plt.plot(wavelengths, Tc, color="blue", label="VO2 in the insulating phase")
plt.plot(wavelengths, Th, color="red", label="VO2 in the metallic phase")

plt.title("Limiting Optical Diodes Enabled by the Phase Transition of Vanadium Dioxide.\n"
	"Using complex resonant cavity. N=4, M=2."
	"Chenghao Wan, et al.\n"
	"TiO2 = 123.0nm, VO2 = 7.6nm, Sapphire = nm", fontsize=11)
plt.xlabel('Wavelength (nm)')
plt.ylabel('Transmission', color="black")
plt.legend(loc='upper right', fontsize=8)

plt.xlim(1300, 1350)
plt.ylim(0.0, 0.5)
plt.show()