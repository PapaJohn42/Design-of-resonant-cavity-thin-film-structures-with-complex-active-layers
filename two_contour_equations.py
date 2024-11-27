import numpy as np
from contour import contour

def two_contour_equations(p, nf1, nf2, ns, na):
	q1,q2 = p
	(cr,ci,_,_,_) = contour(ns, nf1, q1, 1.0, 2)
	ns = cr[-1]+1j*ci[-1]
	(cr,ci,_,_,_) = contour(ns, nf2, q2, 1.0, 2)
	
	return (cr[-1]-np.real(na), ci[-1]-np.imag(na))