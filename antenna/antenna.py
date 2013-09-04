#!/usr/bin/env python
import numpy as np
import re
import string as s

#=====================================================
typeRGX             = re.compile('TYPE\s\/\sSERIAL\sNO')
startOfAntennaRGX   = re.compile('START OF ANTENNA')
endOfAntennaRGX     = re.compile('END OF ANTENNA')

daziRGX             = re.compile('DAZI')
dzenRGX             = re.compile('ZEN1\s\/\sZEN2\s\/\sDZEN')
validFromRGX        = re.compile('VALID FROM')
validToRGX          = re.compile('VALID TO')
sinexCodeRGX        = re.compile('SINEX CODE')
pcoRGX              = re.compile('NORTH / EAST / UP')

numFreqRGX          = re.compile('# OF FREQUENCIES')

startOfFrequencyRGX = re.compile('START OF FREQUENCY')
endOfFrequencyRGX   = re.compile('END OF FREQUENCY')

noaziRGX            = re.compile('NOAZI')
typeRGX             = re.compile('TYPE\s\/\sSERIAL\sNO')
#=====================================================

freqFlag = 0

numAntennas = 0
antennas = []

with open('antmod.dat') as f:
    for line in f:
        if typeRGX.search(line):    
            antenna['name'] = line[0:15]+line[16:20]#+serial number
            antenna['type'] = line[0:15]
            antenna['dome'] = line[16:20]

            # might need to put in an exception for satellite antennas
            # if dome is blank set to none
            if antenna['dome'] == '    ':
                antenna['dome'] = 'NONE'
            # define some defaults
            antenna['frequency'] = []
            antenna['data'] = np.zeros(0) 
        elif startOfAntennaRGX.search(line):
            antenna = {}

        elif daziRGX.search(line):
            dazi = np.array(s.split(line))[0]
            antenna['dazi'] = dazi.astype(np.float)

        elif dzenRGX.search(line):
            dzen = np.array(s.split(line))[0:-5]
            antenna['dzen'] = dzen.astype(np.float)

        elif numFreqRGX.search(line):
            numFreq = np.array(s.split(line))[0]
            antenna['numFreq'] = numFreq.astype(np.int)

        elif validFromRGX.search(line):
            # array should have size 6, YYYY MM DD HH MM SS.SSSS
            # maybe less if it has been entered like 2004 01 01, etc...
            validFrom = np.array(s.split(line)[0:-2])
            antenna['validFrom'] = validFrom.astype(np.float)

        elif validToRGX.search(line):
            validTo = np.array(s.split(line)[0:-2])
            antenna['validTo'] = validTo.astype(np.float) 

        elif startOfFrequencyRGX.search(line):
            cFreq = np.array(s.split(line))[0:-3]
            # convert G 1 to G01, etc..
            if cFreq.size == 2:
                tmp = "{}{:02d}".format(cFreq[0],int(cFreq[1]))
                antenna['frequency'].append(tmp)
            else:
                antenna['frequency'].append(cFreq[0])
                
            freqFlag = 1

        elif pcoRGX.search(line):
            pco = np.array(s.split(line))[0:-5]
            name = 'PCO_'+antenna['frequency'][-1]
            antenna[name] = pco.astype(np.float)

            if antenna['data'].size < 1:
                if antenna['dazi'] < 0.0001 :
                    nAZI = 1
                else:
                    nAZI = int(360. / antenna['dazi']) + 1
                nZEN = int((antenna['dzen'][1] - antenna['dzen'][0])/antenna['dzen'][2])+1
                antenna['data'] = np.zeros((antenna['numFreq'],nAZI,nZEN))

        elif noaziRGX.search(line):
            noazi = np.array(s.split(line)[1:])
            antenna['noazi'] = noazi.astype(np.float)

        elif endOfFrequencyRGX.search(line):
            #print("end of frequency")
            freqFlag = 0
            cFreq = ''

        elif endOfAntennaRGX.search(line):
            #print("end of antenna")
            print(antenna)
            antennas.append(antenna)
            numAntennas += 1
            freqFlag = 0
        elif freqFlag :
            tmp = np.array(s.split(line)[1:],dtype=float)
            itr = 0
            f = np.size(antenna['frequency']) - 1
            for v in tmp:
                antenna['data'][f][freqFlag-1][itr] = v 
                itr += 1

            freqFlag += 1     


