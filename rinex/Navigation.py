#!/usr/bin/env python
from __future__ import division, print_function, absolute_import

import numpy as np
import re
import string as s
import datetime as dt
import gpsTime as gpsT

def parseHeader(obs, line):
    #=====================================================
    rinexVersionRGX     = re.compile('RINEX VERSION / TYPE')
    rnxrPgmRGX          = re.compile('PGM / RUN BY / DATE')
    commentRGX          = re.compile('COMMENT')
    endOfHeaderRGX      = re.compile('END OF HEADER')
    #=====================================================

    if endOfHeaderRGX.search(line): 
        #print("EOH",line)
        return 1

def expDtofloat(s):
    return float(s.replace('D','E'))

def parseNavData(nav,line):
    epochRGX            = re.compile(r'\d+\s+\d+\s+\d+\s+\d+\s+\d+\s+\d+\s+')#+[\d.]+')

    line = line.rstrip()

    if epochRGX.search(line):
#       epoch = {}
#       epoch['data'] = []

#       SV = int(line[0:2])
#       YY = int(line[3:5])
#       MM = int(line[6:8])
#       DD = int(line[9:11])
#       hh = int(line[12:14])
#       mm = int(line[15:17])
#       ss = int(line[18:20])
#       ms = int(line[21:22])

#       if YY < 80:
#           YYYY = YY + 2000
#       else:
#           YYYY = YY + 1900

#       epoch['prn'] = SV
#       epoch['time'] = dt.datetime(YYYY,MM,DD,hh,mm,ss,ms)
#       epoch['data'].append( expDtofloat(line[22:41]) )
#       epoch['data'].append( expDtofloat(line[41:60]) )
#       epoch['data'].append( expDtofloat(line[60:79]) )
#       nav['epochs'].append(epoch)
        parseNavEpoch(nav,line)
    else:
        epoch = nav['epochs'][-1]
        l = len(line)
        if l > 21:
            # check that it has the right format, ie D or E for exponential
            # if this is not present assume this field is corrupt and set it to NaN
	    # set the epoch time to 1980, so that this frame or epoch of data
            # will not be selected
            if re.search(r'([DdEe])',line[3:22]) :
                epoch['data'].append(expDtofloat(line[3:22]) )
            else:
                epoch['data'].append(float('NaN'))
    		epoch['time'] = dt.datetime(1980,01,01,00,00,00,00)
        if l > 40:
            if re.search('[DdEe]',line[22:41]):
                epoch['data'].append(expDtofloat(line[22:41]))
            else:
                epoch['data'].append(float('NaN'))
    		epoch['time'] = dt.datetime(1980,01,01,00,00,00,00)
        if l > 59:
            if re.search('[DdEe]',line[41:60]):
                epoch['data'].append(expDtofloat(line[41:60]))
            else:
                epoch['data'].append(float('NaN'))
    		epoch['time'] = dt.datetime(1980,01,01,00,00,00,00)
        if l > 78:
            if re.search('[DdEe]',line[60:79]):
                epoch['data'].append(expDtofloat(line[60:79]))
            else:
                epoch['data'].append(float('NaN'))
    		epoch['time'] = dt.datetime(1980,01,01,00,00,00,00)
        
    return 1

def parseNavEpoch(nav,line):
    """
        parseNavEpoch(nav,line)

        GPS and GLONASS nav files have the same format for the epoch line:

     |PRN / EPOCH / SV CLK| - Satellite number:                      |     I2,    |
     |                    |       Slot number in sat. constellation  |            |
     |                    | - Epoch of ephemerides             (UTC) |            |
     |                    |     - year (2 digits, padded with 0,     |   1X,I2.2, |
     |                    |                if necessary)             |            |
     |                    |     - month,day,hour,minute,             |  4(1X,I2), |
     |                    |     - second                             |    F5.1,   |
     |                    | - SV clock bias (sec)             (-TauN)|   D19.12,  |
     |                    | - SV relative frequency bias    (+GammaN)|   D19.12,  |
     |                    | - message frame time                 (tk)|   D19.12   |
     |                    |   (0 .le. tk .lt. 86400 sec of day UTC)  |            |

    """

    epoch = {}
    epoch['data'] = []

    SV = int(line[0:2])
    YY = int(line[3:5])
    MM = int(line[6:8])
    DD = int(line[9:11])
    hh = int(line[12:14])
    mm = int(line[15:17])
    ss = int(line[18:20])
    ms = int(line[21:22])

    if YY < 80:
        YYYY = YY + 2000
    else:
        YYYY = YY + 1900

    epoch['prn'] = SV
    epoch['time'] = dt.datetime(YYYY,MM,DD,hh,mm,ss,ms)
    epoch['data'].append( expDtofloat(line[22:41]) )
    epoch['data'].append( expDtofloat(line[41:60]) )
    epoch['data'].append( expDtofloat(line[60:79]) )
    nav['epochs'].append(epoch)

    return 1

def getFrame(sat,Ntime,nav):
    '''
    function 
        frame = getFrame(sat,gpssow,nav)
    '''
    ctr = 0
    match = []
    diff = []

    # loop through all of the epochs in the navigation data
    for epoch in nav['epochs'] :
        # look for any recrds that match the requested PRN
        if epoch['prn'] == sat :
            match.append(ctr)
            d = epoch['time'] - Ntime
            #print('Match ',sat,epoch,ctr,d)
            # calculate the difference in time between the satellite epoch time
            # and the time of broadcast in the the nav data
            diff.append(np.abs(d.total_seconds()))
        ctr += 1
    diff = np.array(diff)

    # return the array which has the smallest absolute difference in time 
    if diff.any():
        return(nav['epochs'][match[diff.argmin()]])
    else :
        return -1

def satpos(sat,Ntime,nav):
    '''
    function satp = satpos(t,eph);
    %SATPOS Calculation of X,Y,Z coordinates at time t
    %        for given ephemeris eph

    '''

    GM = 3.986005e14             # earth's universal gravitational m^3/s^2
    Omegae_dot = 7.2921151467e-5 # earth rotation rate, rad/s

    frame = getFrame(sat,Ntime,nav)

    # Check that a frame has been found  
    if frame == -1:
        return
 
    af0         = frame['data'][0]  # SV clock bias (seconds) 
    af1         = frame['data'][1]  # SV clock drift (sec/sec)
    af2         = frame['data'][2]  # clock drift rate (sec/sec2)
    
    iode        = frame['data'][3]  # IODE Issue of data, Ephemeris
    crs         = frame['data'][4]  # (meters)
    deltan      = frame['data'][5]  # (radians/sec)
    M0          = frame['data'][6]  # (radians)
    
    cuc         = frame['data'][7]  # (radians)
    ecc         = frame['data'][8]  # eccentricity
    cus         = frame['data'][9]  # (radians)
    roota       = frame['data'][10]

    toe         = frame['data'][11] # Time of Ephemeris (sec of GPS week)
    cic         = frame['data'][12] # (radians)
    Omega0      = frame['data'][13] # (radians)
    cis         = frame['data'][14] # (radians)

    i0          = frame['data'][15]  # (radians)
    crc         = frame['data'][16]  # (metres)
    omega       = frame['data'][17]  # (radians)
    Omegadot    = frame['data'][18]  # (radians/sec)

    idot        = frame['data'][19]  # (radians/sec)
    codes       = frame['data'][20]  # codes on L2 channel
    weekno      = frame['data'][21]  # gps week # to go with toe, continuous number not mod(1024)
    L2flag      = frame['data'][22]  # L2 P data flag

    svaccuracy  = frame['data'][23]  # sv accuracy
    svhealth    = frame['data'][24]  # sv health
    tgd         = frame['data'][25]  # TGD 
    iodc        = frame['data'][26]  # IODC Issue of Data, Clock

    # The next frame is not always in the navigation message, 
    # might have the z-count, but not the rest of the fields
    tom          = frame['data'][27]  # Transmission time of message from z-count in hand Over Word
    #fitInt      = frame['data'][28]  # Fit Interval (hours)
    #spare1      = frame['data'][29]  # spare
    #spare2      = frame['data'][30]  # spare
    #=============================================================
    # Procedure for coordinate calculation
    A = roota*roota
    w,t = gpsT.dateTime2gpssow(Ntime)
   
    tk = t-toe

    n0 = np.sqrt(GM/A**3)
    n = n0+deltan
    M = M0+n*tk

    E = M
    for i in range(0,10):
        E_old = E
        E = M+ ecc * np.sin(E)
        dE = np.remainder(E-E_old,2*np.pi)
        if abs(dE) < 1.e-12:
            break

    v = np.arctan2( np.sqrt(1.-ecc**2)*np.sin(E), np.cos(E)-ecc )
    phi = v + omega

    u = phi + cuc*np.cos(2.*phi) + cus*np.sin(2.*phi)
    r = A*(1.-ecc*np.cos(E)) + crc*np.cos(2.*phi)+crs*np.sin(2.*phi)
    i_1 = i0 + idot*tk + cic*np.cos(2.*phi) + cis*np.sin(2.*phi)
    Omega = Omega0 + (Omegadot - Omegae_dot)*tk - Omegae_dot*toe

    x1 = np.cos(u)*r
    y1 = np.sin(u)*r

    X = x1 * np.cos(Omega) - y1 * np.cos(i_1) * np.sin(Omega)
    Y = x1 * np.sin(Omega) + y1 * np.cos(i_1) * np.cos(Omega)
    Z = y1 * np.sin(i_1)

    satp = np.matrix([X,Y,Z])

    # Should I be adding in the IERS terrestrial to celestial inertial rotation matrix??
    #import pdb
    #pdb.set_trace()
    return satp    

#==============================================================================

def satpos_iers(sat,Ntime,nav):
    '''

    satp = satpos(t,eph)

    Calculation of X,Y,Z coordinates at time t
    for a given ephemeris obtained from a navigation file
    The ICD earth rotation has been replaced with an iers2003, which include
    polar motion and precession, nutaiton etc. 
    Uses some simplifications

    '''

    GM = 3.986005e14             # earth's universal gravitational m^3/s^2
    Omegae_dot = 7.2921151467e-5 # earth rotation rate, rad/s

    frame = getFrame(sat,Ntime,nav)
   
    af0         = frame['data'][0]  # SV clock bias (seconds) 
    af1         = frame['data'][1]  # SV clock drift (sec/sec)
    af2         = frame['data'][2]  # clock drift rate (sec/sec2)
    
    iode        = frame['data'][3]  # IODE Issue of data, Ephemeris
    crs         = frame['data'][4]  # (meters)
    deltan      = frame['data'][5]  # (radians/sec)
    M0          = frame['data'][6]  # (radians)
    
    cuc         = frame['data'][7]  # (radians)
    ecc         = frame['data'][8]  # eccentricity
    cus         = frame['data'][9]  # (radians)
    roota       = frame['data'][10]

    toe         = frame['data'][11] # Time of Ephemeris (sec of GPS week)
    cic         = frame['data'][12] # (radians)
    Omega0      = frame['data'][13] # (radians)
    cis         = frame['data'][14] # (radians)

    i0          = frame['data'][15]  # (radians)
    crc         = frame['data'][16]  # (metres)
    omega       = frame['data'][17]  # (radians)
    Omegadot    = frame['data'][18]  # (radians/sec)

    idot        = frame['data'][19]  # (radians/sec)
    codes       = frame['data'][20]  # codes on L2 channel
    weekno      = frame['data'][21]  # gps week # to go with toe, continuous number not mod(1024)
    L2flag      = frame['data'][22]  # L2 P data flag

    svaccuracy  = frame['data'][23]  # sv accuracy
    svhealth    = frame['data'][24]  # sv health
    tgd         = frame['data'][25]  # TGD 
    iodc        = frame['data'][26]  # IODC Issue of Data, Clock

    # The next frame is not always in the navigation message, 
    # might have the z-count, but not the rest of the fields
    tom          = frame['data'][27]  # Transmission time of message from z-count in hand Over Word
    #fitInt      = frame['data'][28]  # Fit Interval (hours)
    #spare1      = frame['data'][29]  # spare
    #spare2      = frame['data'][30]  # spare
    #=============================================================
    # Procedure for coordinate calculation
    A = roota*roota
    w,t = gpsT.dateTime2gpssow(Ntime)
   
    tk = t-toe

    n0 = np.sqrt(GM/A**3)
    n = n0+deltan
    M = M0+n*tk

    E = M
    for i in range(0,10):
        E_old = E
        E = M+ ecc * np.sin(E)
        dE = np.remainder(E-E_old,2*np.pi)
        if abs(dE) < 1.e-12:
            break

    v = np.arctan2( np.sqrt(1.-ecc**2)*np.sin(E), np.cos(E)-ecc )
    phi = v + omega

    u = phi + cuc*np.cos(2.*phi) + cus*np.sin(2.*phi)
    r = A*(1.-ecc*np.cos(E)) + crc*np.cos(2.*phi)+crs*np.sin(2.*phi)
    i_1 = i0 + idot*tk + cic*np.cos(2.*phi) + cis*np.sin(2.*phi)
    Omega = Omega0 + (Omegadot - Omegae_dot)*tk - Omegae_dot*toe

    x1 = np.cos(u)*r
    y1 = np.sin(u)*r
    z1 = y1 * np.sin(i_1) # z1 =z according to the ICD

    # now call the IERS terrestrial to celestial inertial rotation matrix
#    RT2C = []
#    RT2C = iers_QRW(mjd)
#    RT2C = np.matrix(RT2C)

#    satp = RT2C.T * np.array([x1 y1 z1])
    #X = x1 * np.cos(Omega) - y1 * np.cos(i_1) * np.sin(Omega)
    #Y = x1 * np.sin(Omega) + y1 * np.cos(i_1) * np.cos(Omega)
    #Z = y1 * np.sin(i_1)

    #satp = np.matrix([X,Y,Z])

    # Should I be adding in the IERS terrestrial to celestial inertial rotation matrix??
    #import pdb
    #pdb.set_trace()
    return satp    



#=========================================================

def parseFile(navfile):
    headerFlag = 0

    nav = {}
    nav['filename'] = navfile #'test_data/brdc21030.12n'
    nav['numObsType'] = 0
    nav['epochs'] = []

    with open(navfile) as f:
        for line in f:
            if headerFlag == 0 :
                rtn = parseHeader(nav,line)
                if rtn == 1:
                    headerFlag = 1
                    #print("Will no longer check for the header")
                    #print(obs)
            else:
                rtn = parseNavData(nav,line)

    return nav
