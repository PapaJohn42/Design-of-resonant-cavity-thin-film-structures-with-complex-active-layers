##For use in 'Design of resonant cavity thin...' only.
##Recreating Figure 5, specifically N4 M2.

#Unfinished

import numpy as np
from contour import contour
from unit_cell_contour_simplified import unit_cell_contour_simplified

#Parameters setup [input]
na = 1.0            #The same as air
ns = na            	#na = ns satisfied the resonance condition
nH = 2.5            #n1
nL = 1.5            #n2
nC = 2.5 -1j*0.25   #nf1; The central complex layers are (C D C)
nD = 2.5            #nf2; Dielectric layer index
q = 1.0             #Phase thickness of each films is quarter-wave thick.
N = 4               #Numbers of trilayer structure (substrate side)
M = 3               #(Incident side)
wl = 1550.0         #Wavelength in nm
npts = 1000


##Calculate the combined effective refractive index of M layers (air/input side)
#H/2 U4
(cr1a, ci1a, t1a, qr1a, qi1a) = contour(ns, nH, q/2, wl, npts)
ns = cr1a[-1] + 1j*ci1a[-1]
(cr1b, ci1b, t1b, qr1b, qi1b) = unit_cell_contour_simplified(ns, nH, nL, N, wl, 2)
print('> nM=', cr1b[-1] + 1j*ci1b[-1])


##Calculate the combined effective refractive index of N layers (substrate side)
nf3a = -nH          #Use -nH since we want to plot in the reverse direction
(cr3a, ci3a, t3a, qr3a, qi3a) = contour(na, nf3a, q/2, wl, npts) 
ns = cr3a[-1] + 1j*ci3a[-1]
(cr3b, ci3b, t3b, qr3b, qi3b) = unit_cell_contour_simplified(ns, -nH, -nL, M, wl, 2)    #-nH and -nL are used since we calculate the coutour in reverse order (mirror of c1)
print('> nN=', cr3b[-1] + 1j*ci3b[-1])


##Calculate the combined effective refractive index of active layers
n1b = cr1b[-1] + 1j*ci1b[-1]
n3b = cr3b[-1] + 1j*ci3b[-1]
q2C = 0.0406
q2D = 0.9192