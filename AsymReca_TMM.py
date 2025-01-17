""" 
TMM grapher
"""

import numpy as np
import csv
from tmm import tmm
import matplotlib.pyplot as plt

"""
Initial refractive index datasets. Input is .csv files with
[Input] in .csv only: wavelength (in um), Real index and Img index
"""
# Function for reading refractive index datasets
def dataset_reader(path):
	arr = []	# Create an empty array
	with open(path, 'r', encoding='utf-8-sig') as f:
		reader = csv.reader(f)
		next(reader)	# skip the header (1st line)
		for line in reader:
			arr.append(list(map(float, line)))
	return list(zip(*arr))

r1 = dataset_reader("dataset/TiO2.csv")		# [input]
r2 = dataset_reader("dataset/Al2o3.csv")	# [input]	
r3c = dataset_reader("dataset/cVO2.csv")	# [input]
r3h = dataset_reader("dataset/hVO2.csv")	# [input]
r4 = dataset_reader("dataset/air.csv")		# [input]

npts = 2900 # Number of plotted points		# [input]
wavelengths = np.linspace(500, 1700, npts)	# [input]
nH = np.interp(wavelengths, np.asarray(r1[0])*1000.0, np.asarray(r1[1]-1j*np.asarray(r1[2])))		# nTiO2; np.interpolate(range, wl(um -> nm), index n-k)
nL = np.interp(wavelengths, np.asarray(r2[0])*1000.0, np.asarray(r2[1]-1j*np.asarray(r2[2])))		# nAl2O3 (Sapphire)
ncVO2 = np.interp(wavelengths, np.asarray(r3c[0])*1000.0, np.asarray(r3c[1]-1j*np.asarray(r3c[2])))	# ncoldVO2
nhVO2 = np.interp(wavelengths, np.asarray(r3h[0])*1000.0, np.asarray(r3h[1]-1j*np.asarray(r3h[2])))	# nhotVO2
nAir = np.interp(wavelengths, np.asarray(r4[0])*1000.0, np.asarray(r4[1]-1j*np.asarray(r4[2])))		# nAir
ns = nAir


"""
Construct an arrray of resonant cavity structure
"""
# Function for creating Resonant cavity structure in an refractive index array
def AsymReca_index(nH, nL, nf, nD, N, M):
	# [U]^N [C D C] [U]^M
	nU = np.concatenate((nH, nL, nH))	# Unit cell index
	
	# Substrate side
	nPSM = nH
	i=1
	while i <= N:
		nPSM = np.concatenate((nPSM, nU))
		i+=1

	# Central cavity layers
	nPSM = np.concatenate((nPSM, nf, nD, nf))

	# Incident side
	i=1
	while i <= M:
		nPSM = np.concatenate((nPSM, nU))
		i+=1
	nPSM = np.concatenate((nPSM, nH))

	return nPSM

N = 4		# Numbers of trilayer structure (substrate side)	[input]
M = 3		# Numbers of trilayer structures (Incident side)	[input]
nfc = AsymReca_index(nH, nL, ncVO2, nH, N, M)	# Cavity structure in cold VO2 phase
nfh = AsymReca_index(nH, nL, nhVO2, nH, N, M)	# Cavity structure in hot VO2 phase


# Function for creating Resonant cavity structure in an thickness array
def AsymReca_thk(tH, tL, tf, tD, N, M):
	# [U]^N [C D C] [U]^M
	tU = [tH/2, tL, tH/2]
	
	tPSM = [tH/2]
	for N in range(N):
		tPSM += tU

	tPSM += [tf, tD, tf]

	for M in range(M):
		tPSM += tU
	tPSM += [tH/2]

	return tPSM

wl = 1320.0		# Reference wavelength (nm)	[input]
tTiO2 = wl/(4*2.2964) # Quarterwave thick	[input]
tAl2O3 = wl/(4*1.748) # Quarterwave thick	[input]
tVO2 = 7.6		# tVO2						[input]
tD = 123.0		# TiO2 central layer 		[input]
thk = AsymReca_thk(tTiO2, tAl2O3, tVO2, tD, N, M)

ts = 500*1000	# Substrate thickness - not relevant unless back=1
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


"""
Plot a graph(s)
"""

fig, ax1 = plt.subplots()
ax2 = ax1.twinx()

p1, = ax1.plot(wavelengths, Rc, color="purple", label="cold VO2 (resonant)")
p2, = ax1.plot(wavelengths, Rh, color="purple", label="hot VO2 (non-resonant)", linestyle="dashed")

# plt.title("Resonant cavity structure with N = 4 & M = 3 designed for 1550 nm resonance.\n"
# 	"Plot both the Reflection and Transmission spectra of cold & hot VO2 states.", fontsize=11)
ax1.set_xlabel('Wavelength (nm)')
ax1.set_ylabel('Reflection', color="purple")

# choose either T or A
y_graph = "A"	# T: transmission or A: absorption
if y_graph == "T":
	ax2.plot(wavelengths, Tc, color="green")
	ax2.plot(wavelengths, Th, color="green", linestyle="dashed")
	ax2.set_ylabel('Transmission', color="green")
elif y_graph == "A": 
	ax2.plot(wavelengths, Ac, color="red")
	ax2.plot(wavelengths, Ah, color="red", linestyle="dashed")
	ax2.set_ylabel('Absorption', color="red")

ax1.legend(handles=[p1, p2], loc='right', fontsize=8)

plt.xlim(1300, 1350)
ax1.set_ylim(0.0, 1.0)	# Reflection
ax2.set_ylim(0.0, 1.0)	# Transmission
plt.show()