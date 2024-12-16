""" 
Calculate the thickness of each central cavity layers used in asymmetrical resonant cavity structure
H/2 [U]^N [C D C] [U]^M H/2
H = TiO2 = 1.770
L = SiO2 = 1.4657

Result:	[C D C]
N	M	tVO2(nm)	tD (TiO2)(nm)
5	4	2.6			162.6
5	3	9.7			145.3
5	2	36.1		85.1
5	1	159.9		158.2
4	3	6.6			152.6
4	2	29.8		98.9
4	1	157.6		164.6
3	2	18.5		124.5
3	1	151.3		180.3					
2	1	130.8		229.1

Result: [C D]
5	4	5.7			162.1
5	3	27.9		135.9
4	3	17.1		148.3
4	2	142.3		15.7
4	1	188.6		293.7
3	2	128.6		31.3
3	1	182.2		303.1
2	1	166.4		325.0

Result: [D C]	*
5	4	162.6		5.3
5	3	143.4		21.5
4	3	152.6		13.6
4	2	135.8		23.7
4	1	272.0		202.3
3	2	120.3		42.9
3	1	282.1		195.7
2	1	311.0		176.1
"""

import numpy as np
from scipy.optimize import fsolve
from scipy.optimize import minimize
from contour import contour
from unit_cell_contour import unit_cell_contour
from two_contour_equations import three_contour_equations
import matplotlib.pyplot as plt

"""
Parameters setup
"""
na = 1.0            # The same as air
ns = na            	# na = ns satisfied the resonance condition
nH = 2.2899 		# TiO2; Zhukovsky et al. 2015: Thin film; n 0.211–1.69 µm
nL = 1.4657   		# SiO2; Lemarchand 2013: n,k 0.25–2.5 µm	
nC = 2.8593 -1j*0.28114  # Oguntoye et al. 2023: n,k 0.21–2.5 µm; 20 °C
nD = nH             # TiO2
q = 1.0             # Phase thickness of each films is quarter-wave thick.
N = 4               # Numbers of trilayer structure (substrate side)	[input]
M = 2               # Numbers of trilayer structures (Incident side)	[input]
wl = 1550.0         # Wavelength in nm
npts = 1000			# Numbers of plotted points
q2C = 0.2			# starting estimate [input]
q2D = 1			# starting estimate [input]

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
(1) Calculate from the substrate side back to the central cavity structure
assuming substrate is air
"""
## The H/2 layer
(cr1a, ci1a, t1a, qr1a, qi1a) = contour(ns, nH, q/2, wl, npts)	# contour(substrate index, film index, film thickness (unit: quarter wave), number of plotted points)
## N = The number unit cells layers
ns = cr1a[-1] + 1j*ci1a[-1]
(cr1b, ci1b, t1b, qr1b, qi1b) = unit_cell_contour(ns, nH, nL, N, wl, 2)


"""
(2) Calculate from the incident side up to the central cavity structure
"""
## The H/2 layer
(cr3a, ci3a, t3a, qr3a, qi3a) = contour(na, -nH, q/2, wl, npts)	# -nH because we calculate the contour in reverse
## M = The number unit cell(s) layers
ns = cr3a[-1] + 1j*ci3a[-1]
(cr3b, ci3b, t3b, qr3b, qi3b) = unit_cell_contour(ns, -nH, -nL, M, wl, 2)    # -nH and -nL are used since we calculate the coutour in reverse


"""
(3) Calculate the thickness of complex film layers and the dielectric layer
"""
n1b = cr1b[-1] + 1j*ci1b[-1]	# Equivalent to ns
n3b = cr3b[-1] + 1j*ci3b[-1]	# Equivalent to na
def three_contour_equations(p, nf1, nf2, ns, na):
	q1, q2 = p				# p are the solved value
	if q1 <= 0 or q2 <= 0:  # Help penalized negative value when sovling for p
		return (1e6, 1e6)
	
	(cr,ci,_,_,_) = contour(ns, nf1, q1, 1, 2)
	ns = cr[-1]+1j*ci[-1]
	(cr,ci,_,_,_) = contour(ns, nf2, q2, 1, 2)
	ns = cr[-1]+1j*ci[-1]
	# (cr,ci,_,_,_) = contour(ns, nf1, q1, 1, 2)

	return (cr[-1]-np.real(na), ci[-1]-np.imag(na))
q2D, q2C = fsolve(three_contour_equations, (q2D,q2C), (nD, nC, n1b, n3b))

# [nC]
(cr2a, ci2a, t2a, qr2b, qi2b) = contour(n1b, nD, q2D, wl, npts)
ns=cr2a[-1] +1j*ci2a[-1]
# [nC nD]
(cr2b, ci2b, t2b, qr2a, qi2a) = contour(ns, nC, q2C, wl, npts)
ns=cr2b[-1] +1j*ci2b[-1]
# [nC nD nC]
(cr2c, ci2c, t2c, qr2c, qi2c) = contour(ns, nC, q2C, wl, npts)

print('Layers thickness: \
	  \nt2Complex (each)=', round(t2b[-1], 1), 'nm \
	  \nt2Dielectric=', round(t2a[-1], 1), ' nm')


"""
Graph plotter

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
"""