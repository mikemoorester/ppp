import numpy as np

def refell(type):
    '''
    #def [a,b,e2,finv]=refell(type)
    # ported from:

    # REFELL  Computes reference ellispoid parameters.
    #   Reference: Department of Defence World Geodetic
    #   System 1984, DMA TR 8350.2, 1987.
    #   TOPEX Reference: <http://topex-www.jpl.nasa.gov/aviso/text/general/news/hdbk311.htm#CH3.3>
    # Version: 1 Jun 04
    # Useage:  [a,b,e2,finv]=refell(type)
    # Input:   type - reference ellipsoid type (char)
    #                 CLK66 = Clarke 1866
    #                 GRS67 = Geodetic Reference System 1967
    #                 GRS80 = Geodetic Reference System 1980
    #                 WGS72 = World Geodetic System 1972
    #                 WGS84 = World Geodetic System 1984
    #                 ATS77 = Quasi-earth centred ellipsoid for ATS77
    #                 NAD27 = North American Datum 1927 (=CLK66)
    #                 NAD83 = North American Datum 1927 (=GRS80)
    #                 INTER = International
    #                 KRASS = Krassovsky (USSR)
    #                 MAIRY = Modified Airy (Ireland 1965/1975)
    #                 TOPEX = TOPEX/POSEIDON ellipsoid
    # Output:  a    - major semi-axis (m)
    #          b    - minor semi-axis (m)
    #          e2   - eccentricity squared
    #          finv - inverse of flattening

    # Copyright (c) 2011, Michael R. Craymer
    # All rights reserved.
    # Email: mike@craymer.com
    '''
    type.upper()

    # should return a dictionary 'object'

    if (type=='CLK66' or type=='NAD27'):
        a=6378206.4 
        finv=294.9786982 
    elif type=='GRS67':
        a=6378160.0 
        finv=298.247167427 
    elif (type=='GRS80' or type=='NAD83'):
        a=6378137.0 
        finv=298.257222101 
    elif (type=='WGS72'):
        a=6378135.0 
        finv=298.26 
    elif (type=='WGS84'):
        a=6378137.0 
        finv=298.257223563 
    elif type=='ATS77':
        a=6378135.0 
        finv=298.257 
    elif type=='KRASS':
        a=6378245.0 
        finv=298.3 
    elif type=='INTER':
        a=6378388.0 
        finv=297.0 
    elif type=='MAIRY':
        a=6377340.189 
        finv=299.3249646 
    elif type=='TOPEX':
        a=6378136.3 
        finv=298.257 

    f = 1./finv 
    b = a*(1.-f) 
    e2 = 1.-(1.-f)**2 

    return a,b,e2

def ell2xyz(lat,lon,h,a,e2) :
    '''
    % ELL2XYZ  Converts ellipsoidal coordinates to cartesian.
    %   Vectorized.
    % Version: 2011-02-19
    % Useage:  [x,y,z]=ell2xyz(lat,lon,h,a,e2)
    %          [x,y,z]=ell2xyz(lat,lon,h)
    % Input:   lat - vector of ellipsoidal latitudes (radians)
    %          lon - vector of ellipsoidal E longitudes (radians)
    %          h   - vector of ellipsoidal heights (m)
    %          a   - ref. ellipsoid major semi-axis (m); default GRS80
    %          e2  - ref. ellipsoid eccentricity squared; default GRS80
    % Output:  x \
    %          y  > vectors of cartesian coordinates in CT system (m)
    %          z /

    % Copyright (c) 2011, Michael R. Craymer
    % All rights reserved.
    % Email: mike@craymer.com
    '''
    #if nargin ~= 3 & nargin ~= 5
    #    warning('Incorrect number of input arguments');
    #    return
    #end
    #if nargin == 3
    #    [a,b,e2]=refell('grs80');
    #end

    v=a/np.sqrt(1.-e2*np.sin(np.radians(lat))*np.sin(np.radians(lat)))
    x=(v+h)*np.cos(np.radians(lat))*np.cos(np.radians(lon))
    y=(v+h)*np.cos(np.radians(lat))*np.sin(np.radians(lon))
    z=(v*(1.-e2)+h)*np.sin(np.radians(lat));

    return ([x,y,z])


if __name__ == "__main__":

    (a,b,e2) = refell('GRS80')
    print "GRS80:",a,b,e2

    pos = ell2xyz(-45.,140.,0.,a,e2)
    # should be
    # X: -346675.389, Y: 2903851.443, Z: -4487348.409
    print "POS:",pos

    pos = ell2xyz(-41.,140.,0.,a,e2)
    # X: -3692786.960 , Y: 3098616.176, Z: -4162423.201
    print "POS:",pos

    pos = ell2xyz(-25.,140.,0.,a,e2)
    # X: -4430811.872 Y: 3717892.607, Z: 2679074.463
    print "POS:",pos

    pos = ell2xyz(19.,140.,0.,a,e2)
    # X: -4621383.516, Y: 3877801.204, Z: 2063349.463 
    print "POS:",pos

    pos = ell2xyz(88.,140.,0.,a,e2)
    # X: -171089.653, Y: 143561.265, Z: 6352853.879
    print "POS:",pos

