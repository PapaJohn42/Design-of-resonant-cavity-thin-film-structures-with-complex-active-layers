##For use in 'Design of resonant cavity thin...' only.
##Based on '15.2.10 Solving for the Metal and Dielectric Thickness...'
import numpy as np
from contour import contour
from unit_cell_contour import unit_cell_contour
import matplotlib.pyplot as plt

##Parameters setup
ns = 1.0    #same as air
nH = 2.5
nL = 1.5
nD = nH     #Central D layer, equal to nH in this case
t = 1.0     #Each films are quarter-wave thick
wl = 1550.0 #Wavelength in nm
npts = 1000

##The cavity structure is H/2 U4 D U4 H/2; U = H/2 L H/2

##From Substrate to H/2, assuming the substrate is air
nf1 = nH
# Since we don't know what physical thickness(t) is used. Assume that qr=quarter-wave/2
#t1 = 0.5*t
#qr1 = (2.0*np.pi/wl)*np.real(nf1)*t1/(np.pi/2.0)    # q=theta=(2pi/wl)*nf*t1
qr1 = t/2   #H/2 = 0.5 of half-wavelength thickness
print('qr1=',qr1)
(cr1, ci1, t1, qr1, qi1) = contour(ns, nf1, qr1, wl, npts)
print('c1 endpoint=', cr1[-1], ci1[-1])
print('-----')

##From H/2 to U^4 (U=H/2 L H/2)
ns = cr1[-1] + 1j*ci1[-1]
(cr2, ci2, t2, qr2, qi2) = unit_cell_contour(ns, nH, nL, 4, wl, npts)
print('cr2=', cr2)
print('ci2=', ci2)
print('c2 endpoint=', cr2[-1], ci2[-1])
print('-----')
####From H/2 to U^4 (U=H/2 L H/2)
##ns = cr1[-1] + 1j*ci1[-1]
##nf2 = +1j*nH    #equivalent index = +-(j*n1); where n1 is the outer layer (nH)
##q2 = np.arccos(complex((-1/2)*((nH/nL)+(nL/nH))))  #based on equation 7
##print('q2=',q2)
##q2U4 = q2 * 4     #Since the unit cells are stack 4 times (U^4)
##print('q2U4=',q2U4)
##q2U4 = q2U4/(np.pi/2)   #Converted to pi/2 unit (Quarter wave)
##(cr2, ci2, t2, qr2, qi2) = contour(ns, nf2, q2U4, wl, npts)
##print('c2U4 endpoint=', cr2[-1], ci2[-1])
##cr2_ep = [cr1[-1], cr2[-1]]
##ci2_ep = [ci1[-1], ci2[-1]]
##print('-----')
##

##From U^4 to H(D) central layer
ns = cr2[-1] + 1j*ci2[-1]
nf3 = nH    #Centarl D layer contain only real part
##nf3 = 2.5 +1j*0.25   #Central D layer cotain both real and imag part. The example use 2.5j*print('nf3=', nf3)
qr3 = t     #Central layer is H
print('qr3=', qr3)
(cr3, ci3, t3, qr3, qi3) = contour(ns, nf3, qr3, wl, npts)
print('c3 endpoint=', cr3[-1], ci3[-1])
print('-----')


##From H(D) to U^4 (U=H/2 L H/2)
ns = cr3[-1] + 1j*ci3[-1]
(cr4, ci4, t4, qr4, qi4) = unit_cell_contour(ns, nH, nL, 4, wl, npts)
print('cr4=', cr4)
print('ci4=', ci4)
print('c4 endpoint=', cr4[-1], ci4[-1])
print('-----')
####From H to U^4
##ns = cr3[-1] + 1j*ci3[-1]
##nf4 = +1j*nH    #equivalent index = +-(j*n1); where n1 is the outer layer (nH)
##q4 = np.arccos(complex((-1/2)*((nH/nL)+(nL/nH))))  #based on equation 7
##print('q4=',q4)
##q4U4 = q4 * 4     #Since the unit cells are stack 4 times (U^4)
##print('q4U4=',q4U4)
##q4U4 = q4U4/(np.pi/2)   #Converted to pi/2 unit (Quarter wave)
##(cr4, ci4, t4, qr4, qi4) = contour(ns, nf4, q4U4, wl, npts)
##print('c4U4 endpoint=', cr4[-1], ci4[-1])
##cr4_ep = [cr3[-1], cr4[-1]]
##ci4_ep = [ci3[-1], ci4[-1]]
##print('-----')

##From U^4 to H/2
ns = cr4[-1] + 1j*ci4[-1]
nf5 = nH
qr5 = t/2
print('qr5=', qr5)
(cr5, ci5, t5, qr5, qi5) = contour(ns, nf5, qr5, wl, npts)
print('c5 endpoint=', cr5[-1], ci5[-1])
print('-----')


##Plotting the contour graph
plt.plot(cr1, ci1, color='grey')
plt.plot(cr2, ci2, color='purple', linestyle='dashed', marker = 'o')    
plt.plot(cr3, ci3, color='green')
plt.plot(cr4, ci4, color='cyan', linestyle='dashed', marker = 'o')
plt.plot(cr5, ci5, color='grey')

plt.xlabel('Real Index')
plt.ylabel('Imaginary Index')
plt.show()
