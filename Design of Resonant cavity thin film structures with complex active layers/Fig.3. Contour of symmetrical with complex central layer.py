""" 
Fig.3. Contour of the same "[H/2 U4] H [U4 H/2]" cavity structure as in Fig. 2,
except the central layer index is 2.5 - j0.25 
"""

import numpy as np
from contour import contour
from two_contour_equations import two_contour_equations
from unit_cell_contour import unit_cell_contour
from scipy.optimize import fsolve
import matplotlib.pyplot as plt

## Parameters setup
ns = 1.0            # Substrate index. The same as air.
na = ns             # na = ns satisfied the resonance condition
nH = 2.5            # n1
nL = 1.5            # n2
nC = 2.5 -1j*0.25   # nf1; The central I layer is (C D)
nD = nC             # nf2; Dielectric layer index
q = 1.0             # Phase thickness of each films is quarter-wave (pi/2) thick.
N = 4               # Numbers of trilayer structures (Substrate side)
M = N               # Numbers of trilayer structures (Incident side)
wl = 1550.0         # Wavelength in nm
npts = 2000         # Numbers of plotted points

## The cavity structures are H/2 R3 I R2 H/2; R = (H/2 L H/2); I = (C D)
print("### The cavity structures are [H/2 U4] H [U4 H/2]; U = (H/2 L H/2)")


## Calculating substrate (starting point) to U4, assuming the substrate is air
# substrate/air(start) --> H/2 -> U -> U -> U -> U
## From substrate to H/2. '1a'
(cr1a, ci1a, t1a, qr1a, qi1a) = contour(ns, nH, q/2, wl, npts)

## From H/2 to U4 '1b'
ns = cr1a[-1] + 1j*ci1a[-1]
(cr1b, ci1b, t1b, qr1b, qi1b) = unit_cell_contour(ns, nH, nL, N, wl, 2)
# print('cr1b=', cr1b)
# print('ci1b=', ci1b)
print('c1b endpoint =', cr1b[-1] + 1j*ci1b[-1])
# print('-----')


## Fom U4 to Central layers (H)
# R3 -> H/2 -> H/2
n1b = cr1b[-1] + 1j*ci1b[-1]
q2C = q
(cr2a, ci2a, t2a, qr2a, qi2a) = contour(n1b, nC, q2C, wl, npts)
ns = cr2a[-1] + 1j*ci2a[-1]


## Calculating from Central layers (H) to air(end)
# U -> U -> U -> U -> H/2 -> air(end)
## From U4 to H/2 '3b'
(cr3b, ci3b, t3b, qr3b, qi3b) = unit_cell_contour(ns, nH, nL, M, wl, 2)    #-nH and -nL are used since we calculate the coutour in reverse order (mirror of c1)

## From H/2 to air '3a'
ns = cr3b[-1] + 1j*ci3b[-1]
(cr3a, ci3a, t3a, qr3a, qi3a) = contour(ns, nH, q/2, wl, npts)   
# print('cr3b=', cr3b)
# print('ci3b=', ci3b)
print('c3b endpoint =', cr3b[-1] + 1j*ci3b[-1])
print('-----')



## Plotting the contour graph
plt.plot(cr1a, ci1a, color='grey')
plt.plot(cr1b, ci1b, color='purple', linestyle='dashed', marker = 'o')    
plt.plot(cr2a, ci2a, color='green')
plt.plot(cr3b, ci3b, color='cyan', linestyle='dashed', marker = 'o')
plt.plot(cr3a, ci3a, color='grey')

plt.xlabel('Real Index')
plt.ylabel('Imaginary Index')
plt.show()