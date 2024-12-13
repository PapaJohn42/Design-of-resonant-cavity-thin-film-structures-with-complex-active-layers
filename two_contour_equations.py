import numpy as np
from contour import contour

def two_contour_equations(p, nf1, nf2, ns, na):
	q1, q2 = p				# p are the solved value
	if q1 <= 0 or q2 <= 0:  # Help penalized negative value when sovling for p
		return (1e6, 1e6)
	
	(cr,ci,_,_,_) = contour(ns, nf1, q1, 1.0, 2)
	ns = cr[-1]+1j*ci[-1]
	(cr,ci,_,_,_) = contour(ns, nf2, q2, 1.0, 2)
	
	return (cr[-1]-np.real(na), ci[-1]-np.imag(na))