import numpy as np
import math as mth

GM = 3.986005e14         # earth's universal gravitational
                         # parameter m^3/s^2

Omegae_dot = 7.2921151467e-5 # earth rotation rate, rad/s

vel_c = 299792458 # vacuum speed of light, m/s

#==========================================================
def get_binary_eph(navfile) :
    '''
    Eph = get_binary_eph(navfile) 
    
    read in the MATLAB binary file and return a matrix of ephemeris

    input:    
    navfile => brdc0010.12n.nav

    output:
    Eph = [noEphemeris,22[
    '''

    with open(navfile, 'rb') as fid:
        ndata = np.fromfile(fid, np.float64)

    noeph = np.size(ndata) / 22.
    Eph = ndata.reshape(noeph,22)
    return Eph

#==========================================================
def find_eph(Eph,sv,time):
    '''
    %FIND_EPH  Finds the proper column in ephemeris array
    %Kai Borre and C.C. Goad 11-26-96
    %Copyright (c) by Kai Borre
    %$Revision: 1.1 $  $Date: 1998/07/01  $

    find_eph(Eph,sv,time)

    input:
    Eph: Matric returned from get_binary_eph()
    sv = PRN number of satellite
    time = in GPS seconds of week

    output:
    icol: row number of satellite ephemerix data which best matches the time requested
    '''

    icol = 0

    criterion = ( Eph[:,0] == sv )
    isat = np.array(np.where(criterion))

    isat.resize(isat.size,) 

    # Did't find any matches, return
    if np.size(isat) == 0:
        return -1

    dtmin = 86400 #Eph[icol,20] - time

    for t in isat:
        dt = Eph[t,20] - time
        #print "T:", t, "Eph time :",Eph[t,20]," DT:",dt

        if dt < 0.1:
             if np.abs(dt) < np.abs(dtmin):
                #print "ICOL :",icol,"T:",t
                icol = t
                dtmin = dt

    return icol

#==========================================================

def satpos(row,eph,t):
    '''
    function satp = satpos(t,eph);
    %SATPOS Calculation of X,Y,Z coordinates at time t
    %	     for given ephemeris eph

    %Kai Borre 04-09-96
    %Copyright (c) by Kai Borre
    %$Revision: 1.1 $  $Date: 1997/12/06  $
    '''

    #  Units are either seconds, meters, or radians
    #  Assigning the local variables to eph
    svprn       = eph[row,0]
    af2	        = eph[row,1]
    M0	        = eph[row,2]
    roota       = eph[row,3]
    deltan      = eph[row,4]
    ecc	        = eph[row,5]
    omega       = eph[row,6]
    cuc	        = eph[row,7]
    cus	        = eph[row,8]
    crc	        = eph[row,9]
    crs	        = eph[row,10]
    i0	        = eph[row,11]
    idot        = eph[row,12]
    cic	        = eph[row,13]
    cis	        = eph[row,14]
    Omega0      = eph[row,15]
    Omegadot    = eph[row,16]
    toe	        = eph[row,17]
    af0	        = eph[row,18]
    af1	        = eph[row,19]
    toc	        = eph[row,20]

    # Procedure for coordinate calculation
    A = roota*roota
    #tk = check_t(t-toe);
    #MAK
    tk = t-toe
    n0 = np.sqrt(GM/A**3)
    n = n0+deltan
    M = M0+n*tk
    #M = rem(M+2*pi,2*pi);
    #MAK
    E = M
    for i in range(1,10):
        E_old = E
        E = M+ ecc * np.sin(E)
        dE = np.remainder(E-E_old,2*np.pi)
        if abs(dE) < 1.e-12:
            break

    #E = rem(E+2*pi,2*pi);
    #MAK
    v = np.arctan2( np.sqrt(1.-ecc**2)*np.sin(E), np.cos(E)-ecc )
    phi = v + omega
    #phi = rem(phi,2*pi);
    #MAK
    u = phi + cuc*np.cos(2.*phi) + cus*np.sin(2.*phi)
    r = A*(1.-ecc*np.cos(E)) + crc*np.cos(2.*phi)+crs*np.sin(2.*phi)
    i_1 = i0 + idot*tk + cic*np.cos(2.*phi) + cis*np.sin(2.*phi)
    Omega = Omega0 + (Omegadot - Omegae_dot)*tk - Omegae_dot*toe
    #Omega = rem(Omega+2*pi,2*pi)
    #MAK
    x1 = np.cos(u)*r
    y1 = np.sin(u)*r
    #print("R,X1,Y1",r,x1,y1)

    X = x1 * np.cos(Omega) - y1 * np.cos(i_1) * np.sin(Omega)
    Y = x1 * np.sin(Omega) + y1 * np.cos(i_1) * np.cos(Omega)
    Z = y1 * np.sin(i_1)

    satp = np.matrix([X,Y,Z])
    #import pdb
    #pdb.set_trace()
    #return ([X,Y,Z]) 
    return satp

#==========================================================

def togeod(a,finv,X,Y,Z):
    '''
    function [dphi,dlambda,h] = togeod(a,finv,X,Y,Z)
    %TOGEOD   Subroutine to calculate geodetic coordinates
    %         latitude, longitude, height given Cartesian
    %         coordinates X,Y,Z, and reference ellipsoid
    %         values semi-major axis (a) and the inverse
    %         of flattening (finv).

    %  The units of linear parameters X,Y,Z,a must all agree (m,km,mi,ft,..etc)
    %  The output units of angular quantities will be in decimal degrees
    %  (15.5 degrees not 15 deg 30 min).  The output units of h will be the
    %  same as the units of X,Y,Z,a.

    %  Copyright (C) 1987 C. Goad, Columbus, Ohio
    %  Reprinted with permission of author, 1996
    %  Fortran code translated into MATLAB
    %  Kai Borre 03-30-96
    '''
    h = 0
    tolsq = 1.e-10
    maxit = 10
    # compute radians-to-degree factor
    rtd = 180/np.pi
    # compute square of eccentricity
    if finv < 1.e-20:
        esq = 0
    else:
        esq = (2-1/finv)/finv
    oneesq = 1-esq
    # first guess
    # P is distance from spin axix
    P = np.sqrt(X**2+Y**2)
    # direct calculation of longitude
    if P > 1.e-20:
        dlambda = np.arctan2(Y,X)*rtd
    else:
        dlambda = 0

    if (dlambda < 0):
        dlambda = dlambda + 360.
    # r is distance from origin (0,0,0)
    r = np.sqrt(P**2+Z**2)

    if r > 1.e-20:
        sinphi = Z/r
    else:
        sinphi = 0

    dphi = mth.asin(sinphi)
    # initial value of height  =  distance from origin minus
    # approximate distance from origin to surface of ellipsoid
    if r < 1.e-20:
        h = 0
        return
    h = r-a*(1-sinphi*sinphi/finv)
    # iterate
    for i in range(1,maxit):
        sinphi = np.sin(dphi)
        cosphi = np.cos(dphi)
        # compute radius of curvature in prime vertical direction
        N_phi = a/np.sqrt(1-esq*sinphi*sinphi)
        # compute residuals in P and Z
        dP = P - (N_phi + h) * cosphi
        dZ = Z - (N_phi*oneesq + h) * sinphi
        # update height and latitude
        h = h+(sinphi*dZ+cosphi*dP)
        dphi = dphi+(cosphi*dZ-sinphi*dP)/(N_phi + h)
        # test for convergence
        if (dP*dP + dZ*dZ < tolsq):
            break

        # Not Converged--Warn user
        if i == maxit:
            print ' Problem in TOGEOD, did not converge in',i,' iterations'

    dphi = dphi*rtd

    return (dphi,dlambda,h)

def topocent(X,dx):
    '''
    [Az, El, D] = topocent(X,dx)

    TOPOCENT  Transformation of vector dx into topocentric coordinate
          system with origin at X.
          Both parameters are 3 by 1 vectors.
          Output: D    vector length in units like the input
                  Az   azimuth from north positive clockwise, degrees
                  El   elevation angle, degrees

    %Kai Borre 11-24-96
    %Copyright (c) by Kai Borre
    %$Revision: 1.0 $  $Date: 1997/09/26  $
    '''

    #print("TOPOCENT dx",type(dx),dx.shape,dx)
    dtr = np.pi/180.

    [phi,lamb,h] = togeod(6378137.,298.257223563,X[0,0],X[0,1],X[0,2])
    #print("phi ",phi,"lamb",lamb,"h",h)

    cl = np.cos(np.radians(lamb)) 
    sl = np.sin(np.radians(lamb))
    cb = np.cos(np.radians(phi)) 
    sb = np.sin(np.radians(phi))

    F = np.matrix([ [-sl, -sb*cl, cb*cl],
      [cl, -sb*sl, cb*sl],
      [ 0.,    cb,   sb]] )

    local_vector = F.T*dx.T
    #print("Local_vector:",type(local_vector),local_vector.shape,local_vector)

    E = local_vector[0,0]
    N = local_vector[1,0]
    U = local_vector[2,0]

    hor_dis = mth.sqrt(E**2 + N**2)

    if hor_dis < 1.e-20:
        Az = 0.
        El = 90.
    else:
        Az = np.degrees( mth.atan2(E,N) )
        El = np.degrees( mth.atan2(U,hor_dis) )

    if Az < 0:
        Az = Az+360.

    #print("dx",dx.shape)
    D =  mth.sqrt((dx[0,0]**2 + dx[0,1]**2 + dx[0,2]**2))

    return Az,El,D

def earth_rot_corr(traveltime, X_sat):
    '''
    earth_rot_corr  Returns rotated satellite ECEF coordinates
    due to Earth rotation during signal travel time

    %Kai Borre 10-10-96
    %Copyright (c) by Kai Borre
    %$Revision: 1.0 $  $Date: 1997/09/26 $

    '''

    X_sat = np.matrix(X_sat).T
    omegatau = Omegae_dot*traveltime;

    R3 = np.matrix([ [ np.cos(omegatau), np.sin(omegatau), 0],
                    [-np.sin(omegatau), np.cos(omegatau), 0],
                    [ 0., 0., 1.] ])

    X_sat_rot = R3*X_sat

    return X_sat_rot


if __name__ == "__main__":

    ndata = get_binary_eph('test_data/brdc0010.12n.nav')
    # should be 407,22 size 8954
#    print(np.shape(ndata))
#    print(ndata[0,:])

    # Decode 1 set of values for 1 satellite as a test...
    # should return 377
    #icol = find_eph(ndata,1,0)
    #print "ICOL",icol

    # calculate the range
    #satp = satpos(icol,ndata)   
    #print satp 

    import geodetic as geod
    sit1 = {} 
    [a,b,e2] = geod.refell('GRS80')
    sit1['XYZ'] = np.matrix(geod.ell2xyz(-45.0,140.0,0.,a,e2))
    #dx1 = satp - sit1['XYZ'] 

    #(Az,El,Dist) = topocent(sit1['XYZ'],dx1)

    #pos_rot1 = earth_rot_corr( Dist/vel_c, satp )

    # Now try again for sat 6, which is meant to be at elevation 48.17
    icol = find_eph(ndata,10,0)
    satp = satpos(icol,ndata,0)
    dx1 = satp - sit1['XYZ']
    (Az,El,Dist) = topocent(sit1['XYZ'],dx1)

