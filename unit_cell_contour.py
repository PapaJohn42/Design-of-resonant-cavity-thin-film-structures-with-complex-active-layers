##Function for calculating equivalent index(n), phase thickness(q)
##  and create an endpoint for trilayer unit cells of N size
import numpy as np
from contour import contour

def unit_cell_contour(ns, n1, n2, N, wl, npts):   #n1=outer film index, n2=cetral film index
    #print('--unit cell contour')   #debug
    cr=[np.real(ns)]*(N+1)  #c[0] contains the ns of previous film layer
    #print('cr=',cr)    #debug
    ci=[np.imag(ns)]*(N+1)
    #print('ci=',ci)    #debug
    t=[None]*(N+1)
    qr=[None]*(N+1)
    qi=[None]*(N+1)
    nf = +1j * n1
    q2 = np.arccos(complex((-1/2)*((n1/n2)+(n2/n1)))) / (np.pi/2)
    i = 1
    while i <= N:
        #print('loop ', i)  #debug
        (cr[i],ci[i],t[i],qr[i],qi[i]) = contour(ns, nf, q2, wl, npts)
        ns = cr[i][-1] + 1j*ci[i][-1]
        cr[i]=cr[i][-1]
        ci[i]=ci[i][-1]
        i += 1
    return cr, ci, t, qr, qi
#Both cr, ci are a array with c[0] as ns
#c[i] as endpoint of each unit cells
