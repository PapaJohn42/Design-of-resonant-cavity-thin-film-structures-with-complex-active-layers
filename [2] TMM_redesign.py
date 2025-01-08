""" 
TMM grapher for "Limiting Optical Diodes Enabled by the Phase Transition of Vanadium Dioxide"
Layout: Substrate --> H/2 U^N [C D C] U^M H/2 --> Air
"""

import numpy as np
import csv
from tmm import tmm
import matplotlib.pyplot as plt


"""
Input
""" 
# Number of unit cells
N = 5		# Numbers of trilayer structure (substrate side)	[input]
M = 2		# Numbers of trilayer structures (Incident side)	[input]

# Layer's thickness
wl = 1320.0		# Reference wavelength (nm)	[input]
tTiO2 = wl/(4*2.2964) # Quarterwave thick	[input]
tAl2O3 = wl/(4*1.4661) # Quarterwave thick	[input]
tVO2 = 87.1		# tVO2						[input]
tD = 48.0		# TiO2 central layer 		[input]

# choose either T or A
y_graph = "T"	# T: transmission or A: absorption	[input]
wl_low = 1300	# lower wavelength limit
wl_high = 1350	# upper wavelength limit


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

# dataset
r1 = dataset_reader("dataset/TiO2.csv")		# Zhukovsky et al. 2015: Thin film; n 0.211–1.69 µm
r2 = dataset_reader("dataset/SiO2.csv")		# SiO2; Lemarchand 2013: n,k 0.25–2.5 µm
r3c = dataset_reader("dataset/VO2_25deg.csv")	# Beaini et al. 2020: n,k 0.5–25 µm; 25 °C
r3h = dataset_reader("dataset/VO2_100deg.csv")	# Beaini et al. 2020: n,k 0.5–25 µm; 100 °C
r4 = dataset_reader("dataset/air.csv")		# Börzsönyi et al. 2008: n 0.4–1.0 µm + Mathar 2007: n 1.3–2.5 µm 

npts = 500 # Number of plotted points		# [input]
wavelengths = np.linspace(wl_low, wl_high, npts)
nTiO2 = np.interp(wavelengths, np.asarray(r1[0])*1000.0, np.asarray(r1[1]-1j*np.asarray(r1[2])))		# nAu; np.interpolate(range, wl(um -> nm), index n-k)
nSiO2 = np.interp(wavelengths, np.asarray(r2[0])*1000.0, np.asarray(r2[1]-1j*np.asarray(r2[2])))	# nAl2O3 (Sapphire)
ncVO2 = np.interp(wavelengths, np.asarray(r3c[0])*1000.0, np.asarray(r3c[1]-1j*np.asarray(r3c[2])))	# ncoldVO2
nhVO2 = np.interp(wavelengths, np.asarray(r3h[0])*1000.0, np.asarray(r3h[1]-1j*np.asarray(r3h[2])))	# nhotVO2
nAir = np.interp(wavelengths, np.asarray(r4[0])*1000.0, np.asarray(r4[1]-1j*np.asarray(r4[2])))		# nAir
ns = nAir
nH = nTiO2
nL = nSiO2


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

	# Central cavity layers [C D C]
	#####nPSM = np.concatenate((nPSM, nf, nD, nf))
	# Central cavity layers [C D]
	nPSM = np.concatenate((nPSM, nD, nf))	# [input]

	# Incident side
	i=1
	while i <= M:
		nPSM = np.concatenate((nPSM, nU))
		i+=1
	nPSM = np.concatenate((nPSM, nH))

	return nPSM

nfc = AsymReca_index(nH, nL, ncVO2, nH, N, M)	# Cavity structure in cold VO2 phase
nfh = AsymReca_index(nH, nL, nhVO2, nH, N, M)	# Cavity structure in hot VO2 phase


# Function for creating Resonant cavity structure in an thickness array
def AsymReca_thk(tH, tL, tf, tD, N, M):
	# [U]^N [C D C] [U]^M
	tU = [tH/2, tL, tH/2]
	
	tPSM = [tH/2]
	for N in range(N):
		tPSM += tU

	#####tPSM += [tf, tD, tf]
	tPSM += [tD, tf]	# [input]

	for M in range(M):
		tPSM += tU
	tPSM += [tH/2]

	return tPSM

thk = AsymReca_thk(tTiO2, tAl2O3, tVO2, tD, N, M)


# TMM
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


""""""
# Find the index of the wavelength closest to 1550 nm
index = np.abs(wavelengths - wl).argmin()

# Extract values for 1550 nm
transmission_cold = Tc[index]
transmission_hot = Th[index]
absorption_cold = Ac[index]
absorption_hot = Ah[index]
reflection_cold = Rc[index]
reflection_hot = Rh[index]

# Print the values
print(f"Data at {wl} nm:")
print(f"Transmission (cold VO2): {transmission_cold:.4f}")
print(f"Transmission (hot VO2): {transmission_hot:.4f}")
print(f"Absorption (cold VO2): {absorption_cold:.4f}")
print(f"Absorption (hot VO2): {absorption_hot:.4f}")
print(f"Reflection (cold VO2): {reflection_cold:.4f}")
print(f"Reflection (hot VO2): {reflection_hot:.4f}")


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

# Plotting Transmission or Absroption
if y_graph == "T":
	ax2.plot(wavelengths, Tc, color="green")
	ax2.plot(wavelengths, Th, color="green", linestyle="dashed")
	ax2.set_ylabel('Transmission', color="green")
elif y_graph == "A": 
	ax2.plot(wavelengths, Ac, color="red")
	ax2.plot(wavelengths, Ah, color="red", linestyle="dashed")
	ax2.set_ylabel('Absorption', color="red")

ax1.legend(handles=[p1, p2], loc='right', fontsize=8)

plt.xlim(wl_low, wl_high)
ax1.set_ylim(0.0, 1.0)	# Reflection
ax2.set_ylim(0.0, 1.0)	# Transmission
plt.show()