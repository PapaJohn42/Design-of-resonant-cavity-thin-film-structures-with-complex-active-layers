## Optical Thin Film Design Chapter 13
## Figure 13.9

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
nD = 2.5            # nf2; Dielectric layer index
q = 1.0             # Phase thickness of each films is quarter-wave (pi/2) thick.
N = 3               # Numbers of trilayer structures (Substrate side)
M = 2               # Numbers of trilayer structures (Incident side)
wl = 550.0          # Wavelength in nm
npts = 1000         # Numbers of plotted points

## The cavity structures are H/2 R3 I R2 H/2; R = (H/2 L H/2); I = (C D)
print("### The cavity structures are H/2 R3 I R2 H/2; R = (H/2 L H/2); I = (C D)")

## Calculating substrate (starting point) to R3, assuming the substrate is air
# substrate/air(start) --> H/2 -> R -> R -> R ->
## From substrate to H/2. '1a'
(cr1a, ci1a, t1a, qr1a, qi1a) = contour(ns, nH, q/2, wl, npts)
## From H/2 to R3 '1b'
ns = cr1a[-1] + 1j*ci1a[-1]
(cr1b, ci1b, t1b, qr1b, qi1b) = unit_cell_contour(ns, nH, nL, N, wl, 2)
# print('cr1b=', cr1b)
# print('ci1b=', ci1b)
print('c1b endpoint =', cr1b[-1] + 1j*ci1b[-1])
# print('-----')


## Calculating from air (endpoint) back to R2
# <- R <- R <- H/2 <- air(end)
## From air to H/2 '3a'
nf3a = -nH      #Use -nH since we want to plot in the reverse direction
(cr3a, ci3a, t3a, qr3a, qi3a) = contour(na, nf3a, q/2, wl, npts)   
## From H/2 to R2 '3b'
ns = cr3a[-1] + 1j*ci3a[-1]
(cr3b, ci3b, t3b, qr3b, qi3b) = unit_cell_contour(ns, -nH, -nL, M, wl, 2)    #-nH and -nL are used since we calculate the coutour in reverse order (mirror of c1)
# print('cr3b=', cr3b)
# print('ci3b=', ci3b)
print('c3b endpoint =', cr3b[-1] + 1j*ci3b[-1])
print('-----')


## Fom R3 to Active layers (C D)
# R3 -> C -> D ->
n1b = cr1b[-1] + 1j*ci1b[-1]
n3b = cr3b[-1] + 1j*ci3b[-1]
q2C = 1.0   # starting estimate [input]
q2D = 1.0   # starting estimate [input]

q2C, q2D = fsolve(two_contour_equations, (q2C,q2D), (nC, nD, n1b, n3b))
print('q2C =', q2C, "pi/2")
print('q2D =', q2D, "pi/2")

(cr2a, ci2a, t2a, qr2a, qi2a) = contour(n1b, nC, q2C, wl, npts)
print('> t2Complex =', t2a[-1], "nm")
print('> q2Complex =', qr2a[-1]/(np.pi/2) + 1j*qi2a[-1]/(np.pi/2), 'pi/2 unit')

ns=cr2a[-1] +1j*ci2a[-1]
(cr2b, ci2b, t2b, qr2b, qi2b) = contour(ns, nD, q2D, wl, npts)
print('> t2Dielectric =', t2b[-1], "nm")
print('> q2Dielectric =', qr2b[-1]/(np.pi/2) + 1j*qi2b[-1]/(np.pi/2), 'pi/2 unit')


## Plotting the contour graph
plt.plot(cr1a, ci1a, color='grey')
plt.plot(cr1b, ci1b, color='purple', linestyle='dashed', marker = 'o')    
plt.plot(cr2a, ci2a, color='green')
plt.plot(cr2b, ci2b, color='yellow')
plt.plot(cr3b, ci3b, color='cyan', linestyle='dashed', marker = 'o')
plt.plot(cr3a, ci3a, color='grey')

plt.xlabel('Real Index')
plt.ylabel('Imaginary Index')
plt.show()