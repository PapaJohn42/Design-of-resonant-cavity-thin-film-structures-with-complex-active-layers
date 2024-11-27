## Function for calculating equivalent index(n), phase thickness(q)
## and create an endpoint for trilayer unit cells of N size.
## Herpin equivalent principle version, resault might be different from 'unit_cell_contour'
## Can only be use with complex nf

import numpy as np
from contour import contour

# ns = substrate index
# n1 = outer film index
# n2 = central film index
# N = number of trilayer structures
# wl = wavelength (nm)
# npts = number of plotted points

def unit_cell_contour_hep(ns, n1, n2, N, wl, npts):
    # --Create an array of N+1 size. Where cr[0] = ns.
    cr=[np.real(ns)]*(N+1)
    ci=[np.imag(ns)]*(N+1)

    t=[None]*(N+1)
    qr=[None]*(N+1)
    qi=[None]*(N+1)
    nf = +1j * n1   # n = +-jn1
    q2 = np.arccos(complex((-1/2)*((n1/n2)+(n2/n1)))) / (np.pi/2) # in pi/2 unit
    
    i = 1
    while i <= N:
        (cr[i],ci[i],t[i],qr[i],qi[i]) = contour(ns, nf, q2, wl, npts)
        ns = cr[i][-1] + 1j*ci[i][-1]

        cr[i]=cr[i][-1]
        ci[i]=ci[i][-1]
        i += 1
    return cr, ci, t, qr, qi

# c[i] is an array of size N+1
# c[0] is ns
# c[i] is endpoint of each unit cells