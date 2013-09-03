#!/usr/bin/env python
import numpy as np
import re
import string as s
import datetime as dt

def parseHeader(obs, line):
    #=====================================================
    rinexVersionRGX     = re.compile('RINEX VERSION / TYPE')
    rnxrPgmRGX          = re.compile('PGM / RUN BY / DATE')
    commentRGX          = re.compile('COMMENT')
    mNameRGX            = re.compile('MARKER NAME')
    mNumberRGX          = re.compile('MARKER NUMBER')
    observerRGX         = re.compile('OBSERVER / AGENCY')
    receiverRGX         = re.compile('REC # / TYPE / VERS')
    antennaRGX          = re.compile('ANT # / TYPE')
    positionRGX         = re.compile('APPROX POSITION XYZ')
    deltaAntRGX         = re.compile('ANTENNA: DELTA H/E/N')
    wavelengthRGX       = re.compile('WAVELENGTH FACT L1/2')
    obsTypeRGX          = re.compile('# / TYPES OF OBSERV')
    intervalRGX         = re.compile('INTERVAL')
    firstObs            = re.compile('TIME OF FIRST OBS')
    endOfHeaderRGX      = re.compile('END OF HEADER')
    #=====================================================

    if endOfHeaderRGX.search(line): 
        #print("EOH",line)
        return 1
    elif obsTypeRGX.search(line):
        #print("obsType",line)
        if obs['numObsType'] == 0 :
            obs['numObsType'] = int(line[0:6]) 
        obs['obsTypes'] =  np.array(s.split(line)[1:-5])
        #obs['obsTypes'] =  tmp.astype(np.int)

def parseSignalFlag(chunk):
    chunk = chunk.strip()
    if chunk != '' :
        return int(chunk)
    else:
        return 0

def parseObservation(chunk):
    chunk = chunk.strip()
    if len(chunk) < 1:
        return(float('NaN'))
    else:
        return(float(chunk))

def parseObsData(obs,line):

    epochRGX            = re.compile(r'\s+\d+\s+\d+\s+\d+\s+\d+\s+\d+\s+[\d.]+\s+')
    epochRGX2           = re.compile(r'^\s{32}\w+')

    line = line.rstrip()

    if epochRGX.search(line):
        epoch = {}
        #print("New Epoch",line)
        YY = int(line[1:3])
        MM = int(line[4:6])
        DD = int(line[7:9])
        hh = int(line[10:12])
        mm = int(line[13:15])
        ss = int(line[16:18])
        ms = int(line[19:26])

        if YY < 80:
            YYYY = YY + 2000
        else:
            YYYY = YY + 1900

        epoch['time'] = dt.datetime(YYYY,MM,DD,hh,mm,ss,ms)

        epoch['flag'] = int(line[27:29]) 
        epoch['numSats'] = int(line[30:32])
        sats = []
        #print("TIME",YY,MM,DD,hh,mm,ss,flag,"NumSats",numSats)
        s = 32
        e = 35
        if epoch['numSats'] > 12 :
            numSats = 12
        else :
            numSats = epoch['numSats']

        for i in range(0,numSats) :            
            sats.append(line[s:e])
            s += 3
            e += 3
        #print("Sats:",sats)
        
        epoch['sats'] = sats
        epoch['satData'] = np.zeros((epoch['numSats'],obs['numObsType']))
        epoch['satDataFlags'] = np.zeros((epoch['numSats'],obs['numObsType']))
        obs['epochs'].append(epoch)

        # reset counter for the number of sat observations parsed for this epoch
        obs['satsParsed'] = 0
        obs['obsParsed']  = 0
        obs['satObsLinesParsed'] = 0
    elif epochRGX2.search(line): 
        #print(line)
        epoch = obs['epochs'][-1]
        numSats = epoch['numSats'] - 12
        s = 32
        e = 35
        for i in range(0,numSats) :            
            epoch['sats'].append(line[s:e])
            s += 3
            e += 3
        #print("Sats:",sats)
    else:
        # check to see if we have gone over the maximum number of satellites expected in this epoch
        epoch = obs['epochs'][-1]
        if obs['satsParsed'] >= epoch['numSats'] :
            print('dodgy data... something wrong with the rinex file, might want to skip the rest')
            l = 0
        else : 
            l = len(line)

        i = obs['satObsLinesParsed']
        j = i * 5

        if l > 13 :
            f1 = parseObservation(line[1:14])
            epoch['satData'][obs['satsParsed']][j] = f1
        if l > 14 : 
            f1_f = parseSignalFlag(line[14:16])
            epoch['satDataFlags'][obs['satsParsed']][j] = f1_f

        if l > 29 :
            f2 = parseObservation(line[17:30])
            epoch['satData'][obs['satsParsed']][j+1] = f2
        if l > 30 :
            f2_f = parseSignalFlag(line[30:32])
            epoch['satDataFlags'][obs['satsParsed']][j+1] = f2_f

        if l > 45 :
            f3 = parseObservation(line[33:46])
            epoch['satData'][obs['satsParsed']][j+2] = f3
        if l > 46 :
            f3_f = parseSignalFlag(line[46:48])
            epoch['satDataFlags'][obs['satsParsed']][j+2] = f3_f

        if l > 61 :
            f4 = parseObservation(line[49:62])
            epoch['satData'][obs['satsParsed']][j+3] = f4
        if l > 62 :
            f4_f = parseSignalFlag(line[62:64])
            epoch['satDataFlags'][obs['satsParsed']][j+3] = f4_f

        if l > 77 :
            f5 = parseObservation(line[65:78])
            epoch['satData'][obs['satsParsed']][j+4] = f5
        #    print("f1",f1,"f2",f2,"f3",f3,"f4",f4,"f5",f5)
        if l > 78 :
            f5_f = parseSignalFlag(line[78:80])
            epoch['satDataFlags'][obs['satsParsed']][j+4] = f5_f

        obs['satObsLinesParsed'] += 1
        #print("length",len(line),"<",line,">",obs['satObsLinesParsed'])

        # check to see if we've finished reading in all of the expected
        # data from this satellite

        if obs['satObsLinesParsed'] > int(obs['numObsType']/5):
            #print("reset",line)
            obs['obsParsed'] = 0
            obs['satsParsed'] += 1
            obs['satObsLinesParsed'] = 0
        
    return 1

def getObs(obsType,obs) :
    #find the row of the obs Type requested
    ctr = 0
    found = 0

    for t in obs['obsTypes']:
        if t == obsType :
            #print("Found:",t)
            found = 1
        elif found == 0:
            ctr +=1

    #data = []
    #for epoch in obs['epochs']:
    #    for sat in epoch['satData']:
    #        #print(epoch['satData'])
    #        #print(np.shape(epoch['satData']))
    #        print(sat[ctr])
    #        #print(epoch['satData'][sat][ctr])
    return ctr
#=========================================================

def parseRinexObsFile(obsfile) :
    headerFlag = 0

    obs = {}
    obs['filename'] = obsfile
    obs['numObsType'] = 0
    obs['epochs'] = []

    with open(obsfile) as f:
        for line in f:
            if headerFlag == 0 :
                rtn = parseHeader(obs,line)
                if rtn == 1:
                    headerFlag = 1
                    #print("Will no longer check for the header")
            else:
                rtn = parseObsData(obs,line)

    return obs

