##Function for calculating equivalent index(n), phase thickness(q)
##  and create an endpoint for trilayer unit cells of N size
##  In this version, each contours in the trilayer U are calculated.
import numpy as np
from contour import contour

def unit_cell_contour(ns, n1, n2, N, wl, npts):   #n1=outer film index, n2=cetral film index
    #print('--unit cell contour')   #debug
    #--Create an array of N+1 size. Where cr[0] = ns.
    cr=[np.real(ns)]*(N+1)
    #print('cr=',cr)    #debug
    ci=[np.imag(ns)]*(N+1)
    #print('ci=',ci)    #debug
    t=[None]*(N+1)
    qr=[None]*(N+1)
    qi=[None]*(N+1)
    i = 1
    
    while i <= N: #The unit cell structure is H/2, L, H/2
        #print('loop ', i)  #debug
        (cr[i],ci[i],t[i],qr[i],qi[i]) = contour(ns, n1, 0.5, wl, npts)
        ns = cr[i][-1] + 1j*ci[i][-1]
        (cr[i],ci[i],t[i],qr[i],qi[i]) = contour(ns, n2, 1.0, wl, npts)
        ns = cr[i][-1] + 1j*ci[i][-1]
        (cr[i],ci[i],t[i],qr[i],qi[i]) = contour(ns, n1, 0.5, wl, npts)
        ns = cr[i][-1] + 1j*ci[i][-1]
        cr[i]=cr[i][-1]
        ci[i]=ci[i][-1]
        i += 1
    return cr, ci, t, qr, qi

#c[0] is ns
#c[i] is endpoint of each unit cells
