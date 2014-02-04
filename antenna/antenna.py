#!/usr/bin/env python
from __future__ import division, print_function, absolute_import

import numpy as np
import re
import string as s

def parseANTEX(atxFile):
    '''
    parseANTEX(antex.atx)
    '''
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

    #print('Opening File<',atxFile,'>')

    with open(atxFile) as f:
        for line in f:
            if typeRGX.search(line):    
                antenna['name'] = line[0:15]+' '+line[16:20]#+serial number
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
                #print(antenna)
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

        return antennas              

#=====================================
def antennaType(antennaType,antennas):
    '''
    antenna = antennaType(antennaType,antennas)
    '''
    for antenna in antennas:
        if antenna['name'] == antennaType:
            return antenna
    print('Coud not find <',antennaType,'>') 
    return -1

#=====================================
if __name__ == "__main__":

    import argparse

    parser = argparse.ArgumentParser(prog='antenna',description='Parse ANTEX files')

    parser.add_argument('-f', '--file', dest='file1', default='./t/antmod.dat')
    #parser.add_argument('-t', '--AntType',dest='AntType', default='ASH701945C_M    NONE')
    parser.add_argument('-t', '--AntType',dest='AntType', default='TRM59800.00     NONE')
    parser.add_argument('-p', '--plot',dest='plot', default=False, action='store_true')
    parser.add_argument('--polar',dest='polar', default=False, action='store_true')
    parser.add_argument('--elevation',dest='elevation', default=False, action='store_true')
    parser.add_argument('--EM',dest='elevationMedian', default=False, action='store_true', help='Plot the Median PCV vs Elevation')

    args = parser.parse_args()

    antennas = parseANTEX(args.file1)
    antenna = antennaType(args.AntType,antennas)

    if args.plot :
        import matplotlib.pyplot as plt
        from matplotlib import cm

        print('Will try to plot the antenna')
        #aData = np.array(antenna['data'][0])
        aData = antenna['data'][0]
        print(aData)

        if args.polar :
            az = np.linspace(0,360,73)
            zz = np.linspace(0,90,19)

            fig = plt.figure(figsize=(3.62, 2.76))

            ax = fig.add_subplot(111,polar=True)
            ax.set_theta_direction(-1)
            ax.set_theta_offset(np.radians(90.))
            ax.set_ylim([0,1])
            ax.set_rgrids((0.00001, np.radians(20)/np.pi*2, np.radians(40)/np.pi*2,np.radians(60)/np.pi*2,np.radians(80)/np.pi*2),labels=('0', '20', '40', '60', '80'),angle=180)

            ma,mz = np.meshgrid(az,zz,indexing='ij')
            ma = ma.reshape(ma.size,)
            mz = mz.reshape(mz.size,)
            polar = ax.scatter(np.radians(ma), np.radians(mz)/np.pi*2., c=aData, s=50, alpha=1., cmap=cm.RdBu,vmin=-15,vmax=15, lw=0)

            cbar = fig.colorbar(polar,shrink=0.75,pad=.10)
            cbar.ax.tick_params(labelsize=8)
            cbar.set_label('PCV (mm)',size=8)

            for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
                    ax.get_xticklabels() + ax.get_yticklabels()):
                item.set_fontsize(8)

            plt.tight_layout()

        if args.elevation :
            
            zz = np.linspace(0,90,19)
            #zz = np.linspace(0,90,181)
            ele = 90. - zz[::-1]

            # Do an elevation only plot
            fig = plt.figure(figsize=(3.62, 2.76))
            ax = fig.add_subplot(111)
            for zen in aData :
                ax.plot(ele,zen[::-1])
            ax.set_xlabel('Elevation Angle (degrees)',fontsize=8)
            ax.set_ylabel('PCV (mm)',fontsize=8)
            #ax.set_ylim([-15, 15])

            for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
                          ax.get_xticklabels() + ax.get_yticklabels()):
                item.set_fontsize(8)

            plt.tight_layout()

        if args.elevationMedian :
            
            zz = np.linspace(0,90,181)
            #zz = np.linspace(0,90,19)
            ele = 90. - zz[::-1]
            # Do an elevation only plot
            fig = plt.figure(figsize=(3.62, 2.76))
            ax = fig.add_subplot(111)
            #med = np.median(aData,axis=0)
            med = np.mean(aData,axis=0)
            ax.plot(ele,med[::-1])
            #for zen in aData :
            #   ax.plot(ele,zen[::-1])
            ax.set_xlabel('Elevation Angle (degrees)',fontsize=8)
            ax.set_ylabel('PCV (mm)',fontsize=8)
            ax.set_xlim([10, 90])
            ax.set_ylim([-10, 5])
            for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
                          ax.get_xticklabels() + ax.get_yticklabels()):
                item.set_fontsize(8)

            plt.tight_layout()

        plt.show()
