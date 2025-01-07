""" 
TMM grapher for "Optical Limiting Sensor Based on Multilayer Optimization of Ag/VO2 Phase Change Material"
Layout: Substrate --> Ag VO2 --> Air
"""

import numpy as np
import csv
from tmm import tmm
import matplotlib.pyplot as plt


"""
Input
""" 
# Layer's thickness
wl = 1550.0		# Reference wavelength (nm)	[input]
tAg = 5			# [input]
tVO2 = 52		# [input]

# choose either T or A
y_graph = "T"	# T: transmission or A: absorption	[input]
wl_low = 1000	# lower wavelength limit
wl_high = 2000	# upper wavelength limit


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
r1 = dataset_reader("dataset/Ag.csv")			# Ciesielski et al. 2017: Ag/SiO2; n,k 0.191–20.9 µm
r2 = dataset_reader("dataset/SiO2.csv")			# SiO2; Lemarchand 2013: n,k 0.25–2.5 µm	
r3c = dataset_reader("dataset/VO2_20deg.csv")	# Oguntoye et al. 2023: n,k 0.21–2.5 µm; 20 °C
r3h = dataset_reader("dataset/VO2_70deg.csv")	# Oguntoye et al. 2023: n,k 0.21–2.5 µm; 70 °C
r4 = dataset_reader("dataset/air.csv")			# Börzsönyi et al. 2008: n 0.4–1.0 µm + Mathar 2007: n 1.3–2.5 µm

npts = 500 # Number of plotted points		# [input]
wavelengths = np.linspace(wl_low, wl_high, npts)	# [input]
nAg = np.interp(wavelengths, np.asarray(r1[0])*1000.0, np.asarray(r1[1]-1j*np.asarray(r1[2])))		# nAg; np.interpolate(range, wl(um -> nm), index n-k)
nSiO2 = np.interp(wavelengths, np.asarray(r2[0])*1000.0, np.asarray(r2[1]-1j*np.asarray(r2[2])))	# nAl2O3 (Sapphire)
ncVO2 = np.interp(wavelengths, np.asarray(r3c[0])*1000.0, np.asarray(r3c[1]-1j*np.asarray(r3c[2])))	# ncoldVO2
nhVO2 = np.interp(wavelengths, np.asarray(r3h[0])*1000.0, np.asarray(r3h[1]-1j*np.asarray(r3h[2])))	# nhotVO2
nAir = np.interp(wavelengths, np.asarray(r4[0])*1000.0, np.asarray(r4[1]-1j*np.asarray(r4[2])))		# nAir
ns = nAir

"""
Construct an arrray of resonant cavity structure
"""
nfc = np.concatenate((nAg, ncVO2))	# Cavity structure in cold VO2 phase
nfh = np.concatenate((nAg, nhVO2))	# Cavity structure in hot VO2 phase


# Function for creating Resonant cavity structure in an thickness array
def AsymReca_thk(tH, tL, tf, tD, N, M):
	# [U]^N [C D C] [U]^M
	tU = [tH/2, tL, tH/2]
	
	tPSM = [tH/2]
	for N in range(N):
		tPSM += tU

	#####tPSM += [tf, tD, tf]
	tPSM += [tD, tf]

	for M in range(M):
		tPSM += tU
	tPSM += [tH/2]

	return tPSM

thk = np.array([tAg, tVO2])


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