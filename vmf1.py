#!/usr/bin/env python
from __future__ import division, print_function, absolute_import

import numpy as np
import math as mth

def vmf1map(zenithAngles) :
    """
    vmf1h,vmf1w = vmf1map(e, doy)

    Calculate the hydrostatic and wet mapping values as a function of elevation

    need a default value for ah, aw, usually obtained from a grid file

    """

    # using a typical 'a' coefficient TMH notes
    ah = 1.232e-3
    aw = 0.583e-3
    dlat = -45.
    ht = 1.0
    doy = 001
    vh = []
    vw = []

    bh = 0.0029
    c0h = 0.062
    

    # check to see if it is the Southern Hemisphere
    if (dlat < 0 ) :
        phh = np.pi
        c11h = 0.007
        c10h = 0.002
    else :
        phh = 0.0
        c11h = 0.005
        c10h = 0.001

    #doy = tmp[index,0]
    #ah  = tmp[index,1]
    #aw  = tmp[index,2]
    #dzen = tmp[index,3]

    #     reference day is 28 January
    #     this is taken from Niell (1996) to be consistent
    #doy = dmjd  - 44239.d0 + 1 - 28
    ch = ( c0h + ((np.cos(doy/365.25*2.0*np.pi + phh)+1.0)*c11h/2.0
            + c10h) * (1.0 - np.cos(dlat)) )
 
    #for zd in np.linspace(0,90,181) :
    for zd in zenithAngles :
        zd = mth.radians(zd)
        sine = np.sin(np.pi/2.0 - zd)
        beta = bh/(sine + ch)
        gamma = ah/(sine + beta)
        topcon = ( 1.0 + ah/(1.0 + bh/(1.0 +ch)) )
        vmf1h = topcon/(sine + gamma)

        #
        # assume we need to do a height correction, as we
        # obtained the values from a grid
        #
        a_ht = 2.53e-5
        b_ht = 5.49e-3
        c_ht = 1.14e-3

        hs_km = ht/1000.0
        beta = b_ht/(sine + c_ht)
        gamma = a_ht/(sine + beta)
        topcon = (1.0 + a_ht/(1.0 + b_ht/(1.0+c_ht)))
        ht_corr_coef = 1.0/sine - topcon/(sine + gamma)
        ht_corr      = ht_corr_coef * hs_km
        vmf1h        = vmf1h + ht_corr

        bw = 0.00146
        cw = 0.04391
        beta   = bw/( sine + cw )
        gamma  = aw/( sine + beta)
        topcon = (1.0 + aw/(1.0 + bw/(1.0 + cw)))
        vmf1w   = topcon/(sine+gamma)

        #print data[i,1], data[i,2], data[i,3], data[i,4], vmf1h, vmf1w
        zd = mth.degrees(zd)
        #print tmp[index,0], zd, tmp[index,1], tmp[index,2], tmp[index,3], vmf1h, vmf1w, dzen
        #print( zd, vmf1h, vmf1w)
        vh.append(vmf1h)
        vw.append(vmf1w)

    return vh, vw



if __name__ == "__main__":
    from optparse import OptionParser

    """
    # Expects a concatenated file obtained from grdtab:
    # > grdtab dsymo2.273 2012 001 1 map.grid (or use getMet.pl)
    # > grep ATMOSMAP usymo2.* | awk '{if($2=="MODEL"){}else{print $0}}' > yar2.vmf1a
    #   -which contains the DOY, 'A' hydo, 'A' wet and dzen values
    #
    # This will then calculate the vmf1 mapping function for a range of zenithe values 0:0.5:90
    # every 6 hours
    """
    parser = OptionParser()
    parser.add_option("-f", "--primary", dest="filename",
                                  help="file of met valuse for a station", metavar="FILE")

    parser.add_option("-d", "--doy", dest="doy",type="float",
                   help="Day of Year")

    (option, args) = parser.parse_args()
    #==========================

    #
    # Read in the vmf1a coefficient file
    #
    # for YAR2 :
    dlat = -29.0466
    ht = 241.2965
    #
    # for YAR3 :
    #dlat = -29.0465
    #ht = 242.4516

    if (option.filename) :
        data = np.genfromtxt(option.filename)
        index = 0


    if (option.doy) :
        print( "DOY: ", option.doy)

        criterion = ((data[:,1] >= option.doy ) & (data[:,1] < (option.doy+1.) ))
        inda = np.where(criterion)
        tmp = np.column_stack( ( np.transpose(data[inda,1]), 
                             np.transpose(data[inda,2]),
                             np.transpose(data[inda,3]),
                             np.transpose(data[inda,4])
                         ) )

        dlat = mth.radians(dlat)

    h,w = vmf1map()
    print(h,w)
