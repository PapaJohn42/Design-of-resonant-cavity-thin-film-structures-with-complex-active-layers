import numpy as np

### Subroutine for forming the 2x2 matrix
def tmm_matrix(nf, z, k0, q, na, sp):
    nx = na*np.sin(q)
    nz = np.sqrt(nf**2-nx**2)

    if sp == "TE":  # TRANSVERSE ELECTRICE polarization (OBLIQUE INCIDENCE)
        m11 = np.exp(-1j*k0*nz*z)
        m12 = np.exp(1j*k0*nz*z)
        m21 = nz*np.exp(-1j*k0*nz*z)
        m22 = -nz*np.exp(1j*k0*nz*z)
        M = np.array([[m11,m12], [m21,m22]])

    elif sp == "TM":  # TM polarization
        m11 = nz/nf*np.exp(-1j*k0*nz*z)
        m12 = -nz/nf*np.exp(1j*k0*nz*z)
        m21 = nf*np.exp(-1j*k0*nz*z)
        m22 = nf*np.exp(1j*k0*nz*z)
        M = np.array([[m11,m12], [m21,m22]])

    return M


def tmm_n(nf, z1, z2, k0, q, na, sp):
    M1 = tmm_matrix(nf, z1, k0, q, na, sp)
    M2 = tmm_matrix(nf, z2, k0, q, na, sp)
    N = np.matmul(M2, np.linalg.inv(M1))

    return N


def tmm(wl, ns, ts, na, nf, tf, q, sp, back):
    nsi = np.imag(ns)
    tz = np.cumsum(np.insert(np.multiply(tf, -1.0), 0, 0.0))
    nwp = wl.size   # number of wavelength points
    nl = nf.size/nwp        # number of layers
    nl = round(nl)          # convert nl from float to int otherwise it won't be able to use in line 46
    k0 = 2.0*np.pi/wl
    Rr = np.zeros(nwp)
    Rt = np.zeros(nwp)
    Tr = np.zeros(nwp)
    Tt = np.zeros(nwp)
    T = np.zeros(nwp)
    R = np.zeros(nwp)
    t = np.zeros(nwp, dtype = complex)
    r = np.zeros(nwp, dtype = complex)

    for i in np.arange(0, nwp, 1):  # Step through all wavelengths
        Ms = tmm_matrix(ns[i], tz[0], k0[i], q, na, sp)
        N = np.array([[1.0,0.0], [0.0,1.0]])
        for n in np.arange(0, nl, 1):   # Step through all layers, nl must be converted from float to int
            N = np.matmul(tmm_n(nf[n*nwp+i], tz[n], tz[n+1], k0[i], q, na, sp), N)
        Ma = np.linalg.inv(tmm_matrix(na, tz[nl], k0[i], q, na, sp))    # np.linear_algebra.inverse
        M = np.matmul(N, Ms)
        M = np.matmul(Ma, M)
        r[i] = M[1, 0] / M[0, 0]
        R[i] = np.absolute(r[i])**2
        t[i] = 1.0/M[0,0]
        T[i] = np.absolute(t[i])**2
        
        if back == 1:   # If back=1, back interface will be computed. Otherwise ignore.
            Rs = np.absolute((ns[i]-na)/(ns[i]+na))**2
            Ts = np.absolute((1.0-Rs)*ns[i]/na)
            Rr[i] = np.absolute(M[0,1]/M[0,0])**2
            Tr[i] = np.absolute(M[1,1] - M[0,1]*M[1,0]/M[0,0])**2
            Rt[i] = R[i] + Tr[i]*Rs*T[i]*np.exp(4.0*k0[i]*nsi[i]*ts*1000.0)/\
                (1.0-Rr[i]*Rs*np.exp(4.0*k0[i]*nsi[i]*ts*1000.0))
            Tt[i] = Ts*T[i]*np.exp(2.0*k0[i]*nsi[i]*ts*1000.0)/\
                (1.0-Rr[i]*Rs*np.exp(4.0*k0[i]*nsi[i]*ts*1000.0))
        else:
            Rt[i] = R[i]
            Tt[i] = T[i]*np.absolute(ns[i]/na)

    return [Tt, Rt, t, r]