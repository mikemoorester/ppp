import numpy as np

def togeod(a,finv,X,Y,Z):
    """
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
    """
    h = 0;
    tolsq = 1.e-10;
    maxit = 10;
    # compute radians-to-degree factor
    rtd = 180/pi;
    # compute square of eccentricity
    if finv < 1.e-20:
        esq = 0
    else:
        esq = (2-1/finv)/finv

    oneesq = 1-esq
    # first guess
    # P is distance from spin axix
    P = np.sqrt(X^2+Y^2)
    # direct calculation of longitude
    if P > 1.e-20:
        dlambda = np.atan2(Y,X)*rtd
    else:
        dlambda = 0
    
    if (dlambda < 0):
        dlambda = dlambda + 360
    
    # r is distance from origin (0,0,0)
    r = sqrt(P^2+Z^2)
    if r > 1.e-20:
        sinphi = Z/r
    else:
        sinphi = 0

    dphi = np.asin(sinphi)
    # initial value of height  =  distance from origin minus
    # approximate distance from origin to surface of ellipsoid
    if r < 1.e-20:
        h = 0
        return

    h = r-a*(1-sinphi*sinphi/finv)
    # iterate
    for i = 1:maxit :
        sinphi = sin(dphi);
        cosphi = cos(dphi);
        # compute radius of curvature in prime vertical direction
        N_phi = a/sqrt(1-esq*sinphi*sinphi);
        # compute residuals in P and Z
        dP = P - (N_phi + h) * cosphi;
        dZ = Z - (N_phi*oneesq + h) * sinphi;
        # update height and latitude
        h = h+(sinphi*dZ+cosphi*dP);
        dphi = dphi+(cosphi*dZ-sinphi*dP)/(N_phi + h);
        # test for convergence
        if (dP*dP + dZ*dZ < tolsq):
            break;
        

        # Not Converged--Warn user
        if i == maxit:
            fprintf([' Problem in TOGEOD, did not converge in %2.0f',...
            ' iterations\n'],i)

    dphi = dphi*rtd
    return(dphi,dlambda,h)

def topocent(X,dx) :
    """
    %TOPOCENT  Transformation of vector dx into topocentric coordinate
    %          system with origin at X.
    %          Both parameters are 3 by 1 vectors.
    %          Output: D    vector length in units like the input
    %                  Az   azimuth from north positive clockwise, degrees
    %                  El   elevation angle, degrees

    %Kai Borre 11-24-96
    %Copyright (c) by Kai Borre
    %$Revision: 1.0 $  $Date: 1997/09/26  $
    """
    dtr = pi/180
    [phi,lambda,h] = togeod(6378137,298.257223563,X(1),X(2),X(3))
    cl = cos(lambda*dtr); sl = sin(lambda*dtr)
    cb = cos(phi*dtr); sb = sin(phi*dtr)
    F = [-sl -sb*cl cb*cl;
      cl -sb*sl cb*sl;
       0    cb   sb];

    local_vector = np.transpose(F)*dx

    E = local_vector(1)
    N = local_vector(2)
    U = local_vector(3)

    hor_dis = np.sqrt(E**2 + N**2)
    if hor_dis < 1.e-20:
        Az = 0
        El = 90
    else:
        Az = atan2(E,N)/dtr;
        El = atan2(U,hor_dis)/dtr;
    
    if Az < 0:
        Az = Az+360

    D = np.sqrt(dx(1)**2+dx(2)**2+dx(3)**2)

    return(Az,El,D)
