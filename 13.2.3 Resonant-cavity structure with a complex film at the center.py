##For use in 'Design of resonant cavity thin...' only.
##Recreating Figure 13.9
##The qC and qD answers got are different from the one in pg.195
import numpy as np
from scipy.optimize import fsolve
from contour import contour
#from two_contour_equations import two_contour_equations
from unit_cell_contour import unit_cell_contour
import matplotlib.pyplot as plt

##Parameters setup
ns = 1.0            #The same as air
nH = 2.5            #n1
nL = 1.5            #n2
nC = 2.5 +1j*0.25   #nf1 need to invese the imaginary part; The central I layer is (C D)
nD = 2.5            #nf2; Dielectric layer
t = 1.0             #Thickness of each films is pi/2 thick.
wl = 550.0          #Wavelength in nm
npts = 1000

##The cavity structure is H/2 U3 (C D) U2 H/2; U = H/2 L H/2

##From substrate to U3, assuming the substrate is air
##---From substrate to H/2. '1a'
qr1a = t/2
(cr1a, ci1a, t1a, qr1a, qi1a) = contour(ns, nH, qr1a, wl, npts)
##---From H/2 to U3 '1b'
ns = cr1a[-1] + 1j*ci1a[-1]
(cr1b, ci1b, t1b, qr1b, qi1b) = unit_cell_contour(ns, nH, nL, 3, wl, npts)
print('cr1b=', cr1b)
print('ci1b=', ci1b)
print('> c1b endpoint=', cr1b[-1] + 1j*ci1b[-1])
print('-----')


##From air to U2
##---From air to H/2 '3a'
ns = 1.0
nf3a = -nH      #-nH since c3a[0] is the end of the contour
qr3a = t/2
(cr3a, ci3a, t3a, qr3a, qi3a) = contour(ns, nf3a, qr3a, wl, npts)   
##---From H/2 to U2 '3b'
ns = cr3a[-1] + 1j*ci3a[-1]
(cr3b, ci3b, t3b, qr3b, qi3b) = unit_cell_contour(ns, -nH, -nL, 2, wl, npts)    #-nH and -nL are used since we calculate the coutour in reverse order (mirror of c1)
print('cr3b=', cr3b)
print('ci3b=', ci3b)
print('> c3b endpoint=', cr3b[-1] + 1j*ci3b[-1])
print('-----')


##Fom U3 to Complex film (C D)
n1b = cr1b[-1] + 1j*ci1b[-1]
n3b = cr3b[-1] + 1j*ci3b[-1]
q2C = 0.5   #initial guess
q2D = 0.2  #initial guess
##(cr2,ci2,_,_,_) = contour(n1b, nC, 0.77, wl, npts)
##ns=cr2[-1] +1j*ci2[-1]
##(cr2b,ci2b,_,_,_) = contour(ns, nD, 0.23, wl, npts)
#print('> c2b endpoint=', cr2b[-1] + 1j*ci2b[-1])

def two_contour_equations(p, nf1, nf2, ns, na):
	q1,q2 = p
	(cr,ci,_,qr,qi) = contour(ns, nf1, q1, wl, 2)
	ns = cr[-1]+1j*ci[-1]
	(cr,ci,_,qr,qi) = contour(ns, nf2, q2, wl, 2)
	return (cr[-1]-np.real(na), ci[-1]-np.imag(na))

##n1b=0.093+1j*2.498
##n3b=0.258+1j*2.486
##nC=2.5+1j*0.25
##nD=2.5
q2C, q2D = fsolve(two_contour_equations, (q2C,q2D), (nC, nD, n1b, n3b))
print('q2C=', q2C)
print('q2D=', q2D)
###
(cr2a, ci2a, t2a, qr2a, qi2a) = contour(n1b, nC, q2C, wl, npts)
print('> t2a=', t2a[-1])
print('> q2a=', qr2a[-1]/(np.pi/2) + 1j*qi2a[-1]/(np.pi/2), 'pi/2 unit')
ns=cr2a[-1] +1j*ci2a[-1]
(cr2b, ci2b, t2b, qr2b, qi2b) = contour(ns, nD, q2D, wl, npts)
print('> t2b =', t2b[-1])
print('> q2b=', qr2b[-1]/(np.pi/2) + 1j*qi2b[-1]/(np.pi/2), 'pi/2 unit')


##Plotting the contour graph
plt.plot(cr1a, ci1a, color='grey')
plt.plot(cr1b, ci1b, color='purple', linestyle='dashed', marker = 'o')    
plt.plot(cr2a, ci2a, color='green')
plt.plot(cr2b, ci2b, color='yellow')
plt.plot(cr3b, ci3b, color='cyan', linestyle='dashed', marker = 'o')
plt.plot(cr3a, ci3a, color='grey')

plt.xlabel('Real Index')
plt.ylabel('Imaginary Index')
plt.show()
