""" 
Table 3. Calculated Results for the Resonant Cavity 
Designed with a VO2 Layer in the Cold State

Result:
N	M	tcVO2(nm)	tTiO2(nm)	tc(paper)	td(paper)
4	3	4.92		155.39		3.51		149.45
4	2	20.01		116.36		15.69		120.48
4	1	135.87		205.21		132.02		192.53
3	2	13.14		133.58		11.05		131.30
3	1	129.10		222.27		126.99		204.12
"""

import numpy as np
from scipy.optimize import fsolve
from scipy.optimize import minimize
from contour import contour
#from two_contour_equations import two_contour_equations
from unit_cell_contour import unit_cell_contour
import matplotlib.pyplot as plt

##Parameters setup [input]
na = 1.0            # The same as air
ns = na            	# na = ns satisfied the resonance condition
nH = 2.2899      	# TiO2=2.2899; Zhukovsky et al. 2015: Thin film; n 0.211–1.69 µm [input]
nL = 1.4657         # SiO2=1.4657; Lemarchand 2013: n,k 0.25–2.5 µm	[input]
nC = 3.1462 -1j*0.3388   # cVO2=3.1462+0.3388j; Beaini et al. 2020: n,k 0.5–25 µm; 25 °C
# nH = 2.45
# nL = 1.44
# nC = 3.02 - 1j*0.325
nD = nH             # TiO2
q = 1.0             # Phase thickness of each films is quarter-wave thick.
N = 3               # Numbers of trilayer structure (substrate side)	[input]
M = 1               # Numbers of trilayer structures (Incident side)	[input]
wl = 1550.0         # Wavelength in nm
npts = 1000			# Numbers of plotted points

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


##Calculating substrate (starting point) to U4, assuming the substrate is air
# substrate/air(start) --> H/2 -> UN
##---From substrate to H/2. '1a'
(cr1a, ci1a, t1a, qr1a, qi1a) = contour(ns, nH, q/2, wl, npts)
##---From H/2 to U3 '1b'
ns = cr1a[-1] + 1j*ci1a[-1]
(cr1b, ci1b, t1b, qr1b, qi1b) = unit_cell_contour(ns, nH, nL, N, wl, 2)
##print('cr1b=', cr1b)
##print('ci1b=', ci1b)
print('> c1b endpoint=', cr1b[-1] + 1j*ci1b[-1])


##Calculating from air (endpoint) back to U2
# <- UM <- H/2 <- air(end)
##---From air to H/2 '3a'
nf3a = -nH          #Use -nH since we want to plot in the reverse direction
(cr3a, ci3a, t3a, qr3a, qi3a) = contour(na, nf3a, q/2, wl, npts)   
##---From H/2 to U2 '3b'
ns = cr3a[-1] + 1j*ci3a[-1]
(cr3b, ci3b, t3b, qr3b, qi3b) = unit_cell_contour(ns, -nH, -nL, M, wl, 2)    #-nH and -nL are used since we calculate the coutour in reverse order (mirror of c1)
##print('cr3b=', cr3b)
##print('ci3b=', ci3b)
print('> c3b endpoint=', cr3b[-1] + 1j*ci3b[-1])


##Fom U4 to Active layers (C D C)
# U4 -> C -> D -> C ->
n1b = cr1b[-1] + 1j*ci1b[-1]
n3b = cr3b[-1] + 1j*ci3b[-1]
q2C = 0.5	#starting estimate [input]
q2D = 1.0	#starting estimate [input]

def three_contour_equations(p, nf1, nf2, ns, na):
	q1,q2 = p
	(cr,ci,_,_,_) = contour(ns, nf1, q1, 1.0, 2)
	ns = cr[-1]+1j*ci[-1]
	(cr,ci,_,_,_) = contour(ns, nf2, q2, 1.0, 2)
	ns = cr[-1]+1j*ci[-1]
	(cr,ci,_,_,_) = contour(ns, nf1, q1, 1.0, 2)
	return (cr[-1]-np.real(na), ci[-1]-np.imag(na))

q2C, q2D = fsolve(three_contour_equations, (q2C,q2D), (nC, nD, n1b, n3b))
print('q2Complex=', q2C, 'pi/2 unit')
print('q2Dielectric=', q2D, 'pi/2 unit')


(cr2a, ci2a, t2a, qr2a, qi2a) = contour(n1b, nC, q2C, wl, npts)
#print('> t2a=', t2a[-1])
#print('> q2a=', qr2a[-1]/(np.pi/2) + 1j*qi2a[-1]/(np.pi/2), 'pi/2 unit')
ns=cr2a[-1] +1j*ci2a[-1]
(cr2b, ci2b, t2b, qr2b, qi2b) = contour(ns, nD, q2D, wl, npts)
#print('> t2b =', t2b[-1])
#print('> q2b=', qr2b[-1]/(np.pi/2) + 1j*qi2b[-1]/(np.pi/2), 'pi/2 unit')
ns=cr2b[-1] +1j*ci2b[-1]
(cr2c, ci2c, t2c, qr2c, qi2c) = contour(ns, nC, q2C, wl, npts)
#print('> t2c =', t2c[-1])
#print('> q2c=', qr2c[-1]/(np.pi/2) + 1j*qi2c[-1]/(np.pi/2), 'pi/2 unit')
print('Layers thickness: \
	  \nt2Complex (each)=', round(t2a[-1], 2), 'nm \
	  \nt2Dielectric=', round(t2b[-1], 2), ' nm')



##Plotting the contour graph
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