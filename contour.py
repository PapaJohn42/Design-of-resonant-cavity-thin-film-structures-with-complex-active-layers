import numpy as np
#ns = substrate (or equivalent. ie previous film)
#nf = film index
#phase = in units of pi/2 (real part only)
#wl = wavelength in nm
#npts = number of points in arc
def contour(ns, nf, phase, wl, npts):
    #phase = phase*(1.0+1j*np.imag(nf)/np.real(nf))  #np.real(nf))
    #print('phase=', phase)
    phase=complex(phase)    #Use complex() instead. Avoided the problem when phase only contain imag part.
    #print('new_phase=', phase)
    q = np.linspace(0, phase*np.pi/2, npts)
    #print('q=',q)
    #print('ns, nf=',ns,nf)
    nr = nf * ((ns+nf) + (ns-nf)*np.exp(-1j*2*q)) \
         / ((ns+nf) - (ns-nf)*np.exp(-1j*2*q))
    cr = np.real(nr)
    ci = np.imag(nr)
    qr = np.real(q)
    qi = np.imag(q)
    t = np.real(q/(2.0*np.pi/wl)/nf)
    return cr, ci, t, qr, qi
