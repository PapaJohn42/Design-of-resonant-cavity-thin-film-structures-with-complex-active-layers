""" 
Table 3. Calculated Results for the Resonant Cavity 
Designed for resonant with a VO2 in the Cold State

Result:
N	M	tcVO2(nm)	tTiO2(nm)
4	3	7.6			123.0
4	2	27.6		74.9
4	1	120.2		165.6
3	2	14.5		105.5
3	1	110.8		189.3
2	1	36.1		56.9
"""

import numpy as np
from scipy.optimize import fsolve
from scipy.optimize import minimize
from contour import contour
from unit_cell_contour import unit_cell_contour
import matplotlib.pyplot as plt

"""
Parameters setup
"""
na = 1.0            # The same as air
ns = na            	# na = ns satisfied the resonance condition
nH = 2.2964      	# TiO2; Zhukovsky et al. 2015: Thin film; n 0.211–1.69 µm [input]
nL = 1.748   		# Al2O3; Querry 1985: α-Al2O3 (Sapphire); n,k(o) 0.21–55.6 µm	[input]
#nL = 1.748 -1j*0.019   	# Al2O3; Querry 1985: α-Al2O3 (Sapphire); n,k(o) 0.21–55.6 µm	[input]
nC = 3.1268 -1j*0.3969  # VO2; Beaini et al. 2020: n,k 0.5–25 µm; 25 °C
nD = nH             # TiO2
q = 1.0             # Phase thickness of each films is quarter-wave thick.
N = 4               # Numbers of trilayer structure (substrate side)	[input]
M = 2               # Numbers of trilayer structures (Incident side)	[input]
wl = 1320.0         # Wavelength in nm
npts = 1000			# Numbers of plotted points
q2C = 0.5			# starting estimate [input]
q2D = 0.5			# starting estimate [input]

## The cavity structures are H/2 R3 I R2 H/2; R = (H/2 L H/2); I = (C D)
print("### The cavity structures are [H/2 UN] [nC nD nC] [UM H/2]; U = (H/2 L H/2)")
print('Parameters')
print('nH=', nH)
print('nL=', nL)
print('nC=', nC)
print('nD=', nD)
print('wl=', wl)
print('N=', N)
print('M=', M)
print('===')


"""
Calculate from the substrate side back to the central cavity structure
assuming substrate is air
"""
## The H/2 layer
(cr1a, ci1a, t1a, qr1a, qi1a) = contour(ns, nH, q/2, wl, npts)	# contour(substrate index, film index, film thickness (unit: quarter wave), number of plotted points)
## N = The number unit cells layers
ns = cr1a[-1] + 1j*ci1a[-1]
(cr1b, ci1b, t1b, qr1b, qi1b) = unit_cell_contour(ns, nH, nL, N, wl, 2)


"""
Calculate from the incident side up to the central cavity structure
"""
## The H/2 layer
(cr3a, ci3a, t3a, qr3a, qi3a) = contour(na, -nH, q/2, wl, npts)	# -nH because we calculate the contour in reverse
## M = The number unit cell(s) layers
ns = cr3a[-1] + 1j*ci3a[-1]
(cr3b, ci3b, t3b, qr3b, qi3b) = unit_cell_contour(ns, -nH, -nL, M, wl, 2)    # -nH and -nL are used since we calculate the coutour in reverse


"""
Calculate the thickness of complex film layers and the dielectric layer
"""
# Fucntion for [nC nD nC] structure. Refer to 'two_contour_equations.py' [nC nD]
def three_contour_equations(p, nf1, nf2, ns, na):
	q1, q2 = p				# p are the solved value
	if q1 <= 0 or q2 <= 0:  # Help penalized negative value when sovling for p
		return (1e6, 1e6)
	
	(cr,ci,_,_,_) = contour(ns, nf1, q1, 1, 2)
	ns = cr[-1]+1j*ci[-1]
	(cr,ci,_,_,_) = contour(ns, nf2, q2, 1, 2)
	ns = cr[-1]+1j*ci[-1]
	(cr,ci,_,_,_) = contour(ns, nf1, q1, 1, 2)

	return (cr[-1]-np.real(na), ci[-1]-np.imag(na))

n1b = cr1b[-1] + 1j*ci1b[-1]	# Equivalent to ns
n3b = cr3b[-1] + 1j*ci3b[-1]	# Equivalent to na
q2C, q2D = fsolve(three_contour_equations, (q2C,q2D), (nC, nD, n1b, n3b))

## [nC]
(cr2a, ci2a, t2a, qr2a, qi2a) = contour(n1b, nC, q2C, wl, npts)
ns=cr2a[-1] +1j*ci2a[-1]
# [nC nD]
(cr2b, ci2b, t2b, qr2b, qi2b) = contour(ns, nD, q2D, wl, npts)
ns=cr2b[-1] +1j*ci2b[-1]
# [nC nD nC]
(cr2c, ci2c, t2c, qr2c, qi2c) = contour(ns, nC, q2C, wl, npts)

print('Layers thickness: \
	  \nt2Complex (each)=', round(t2a[-1], 1), 'nm \
	  \nt2Dielectric=', round(t2b[-1], 1), ' nm')


"""
Graph plotter
"""
plt.plot(cr1a, ci1a, color='grey')
plt.plot(cr1b, ci1b, color='purple', linestyle='dashed', marker = 'o')    
plt.plot(cr2a, ci2a, color='green')
plt.plot(cr2b, ci2b, color='yellow')
plt.plot(cr2c, ci2c, color='green')
plt.plot(cr3b, ci3b, color='cyan', linestyle='dashed', marker = 'o')
plt.plot(cr3a, ci3a, color='grey')

plt.xlabel('Real Index')
plt.ylabel('Imaginary Index')
plt.show()