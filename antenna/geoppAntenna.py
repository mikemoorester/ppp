#!/usr/bin/env python
from __future__ import division, print_function, absolute_import

import numpy as np
import re
import string as s


def parseGEOPP(geoppFile):
    '''
    parseGEOPP(antmod.are)
    -parse the geopp .are file format
    -parse the geopp .arp file format

    are is an elvation dependent only antenna calibration result, where the
    PCO has been pushed into the PCVs.
    '''

    #=====================================================
    typeRGX            = re.compile('^TYPE=')
    serialRGX          = re.compile('SERIAL NUMBER=')
    calibrationTypeRGX = re.compile('CALIBRATION TYPE=')
    calibrationDateRGX = re.compile('CALIBRATION DATE=')
    numAntennasRGX     = re.compile('NO OF ANTENNAS=')
    numCalibrationRGX  = re.compile('NO OF CALIBRATIONS=')
    gnssTypeRGX        = re.compile('GNSS TYPE=')
    contentRGX         = re.compile('CONTENT TYPE=')
    pcvTypeRGX         = re.compile('PCV TYPE=')
    numFreqRGX         = re.compile('NO OF FREQUENCIES=')
    offsetL1RGX        = re.compile('OFFSETS L1=')
    offsetL2RGX        = re.compile('OFFSETS L2=')
    deleRGX            = re.compile('ELEVATION INCREMENT=')
    daziRGX            = re.compile('AZIMUTH INCREMENT=')
    varL1RGX           = re.compile('VARIATIONS L1=')
    varL2RGX           = re.compile('VARIATIONS L2=')
    #=====================================================

    L1flag = 0
    L2flag = 0
    ObsCtr = 0
    antenna = {}

    with open(geoppFile) as f:
        for line in f:
            line = line.rstrip()
            if typeRGX.search(line) :
                antenna['type'] = line[5:21]
                antenna['dome'] = line[21:25]
                if antenna['dome'] == '    ':
                    antenna['dome'] = 'NONE'
                antenna['name'] = antenna['type'] + antenna['dome']
            elif serialRGX.search(line) and len(line) > 14 :
                antenna['serialnum'] = line[14:]
                antenna['name'] = antenna['name'] + ' ' + antenna['serialnum']
            elif calibrationTypeRGX.search(line) and len(line) > 17 :
                antenna['calType'] = line[17:]
            elif calibrationDateRGX.search(line) and len(line) > 17 :
                antenna['calDate'] = line[17:]
            elif numAntennasRGX.search(line) and len(line) > 15 :
                antenna['numAntennas'] = int(line[15:])
            elif numCalibrationRGX.search(line) and len(line) > 19 :
                antenna['numCalibrations'] = int(line[19:])
            elif gnssTypeRGX.search(line) and len(line) > 10 :
                antenna['gnssType'] = line[10:]
            elif contentRGX.search(line) and len(line) > 13 :
                antenna['content'] = line[13:]
            elif pcvTypeRGX.search(line) and len(line) > 9:
                antenna['pcvType'] = line[9:]
            elif numFreqRGX.search(line) and len(line) > 18 :
                antenna['numFreq'] = int(line[18:])
            elif offsetL1RGX.search(line) and len(line) > 11 :
                antenna['offsetL1'] = line[11:]
            elif offsetL2RGX.search(line) and len(line) > 11 :
                antenna['offsetL2'] = line[11:]
            elif deleRGX.search(line) and len(line) > 20 :
                antenna['dele'] = float(line[20:])
            elif daziRGX.search(line) and len(line) > 18 :
                antenna['dazi'] = float(line[18:])
            elif varL1RGX.search(line) :
                L1flag = 1
                ObsCtr = 0
            elif varL2RGX.search(line) :
                L2flag = 1
                ObsCtr = 0
            elif L1flag == 1:
                tmp = np.array(s.split(line))
                #check that all of the data has been read in before reseting the flag
                if antenna['dazi'] < 0.0001 :
                    antenna['L1PCV'] = tmp.astype(np.float)
                    L1flag = 0
                    ObsCtr = 0
                else :
                    if ObsCtr == 0:
                        rows = int(360./antenna['dazi']) + 1
                        cols = int(90./antenna['dele']) + 1
                        antenna['L1PCV'] = np.zeros((rows,cols))

                    antenna['L1PCV'][ObsCtr,:] = tmp.astype(np.float)
                    ObsCtr += 1

                    if ObsCtr == (int(360./antenna['dazi']) + 1):
                        L1flag = 0
                        ObsCtr = 0
            elif L2flag == 1:
                tmp = np.array(s.split(line))
                if antenna['dazi'] < 0.0001 :
                    antenna['L2PCV'] = tmp.astype(np.float)
                    L2flag = 0
                    ObsCtr = 0
                else :
                    if ObsCtr == 0:
                        rows = int(360./antenna['dazi']) + 1
                        cols = int(90./antenna['dele']) + 1
                        antenna['L2PCV'] = np.zeros((rows,cols))

                    antenna['L2PCV'][ObsCtr,:] = tmp.astype(np.float)
                    ObsCtr += 1

                    if ObsCtr == (int(360./antenna['dazi']) + 1):
                        L2flag = 0
                        ObsCtr = 0
                L2flag = 0

    return antenna

if __name__ == "__main__":

    import matplotlib.pyplot as plt

    from optparse import OptionParser

    parser = OptionParser()
    parser.add_option("-d","--diff", dest="difference",action='store_true',default=False,help="Reference Antenna")
    parser.add_option("-r","--ref", dest="refAntenna", help="Reference Antenna")
    parser.add_option("-p","--print", dest="printFile",action='store_true',default=False,help="Print Difference")

    parser.add_option("-f","--f1", dest="file1", help="Reference Antenna")
    parser.add_option("--f2", dest="file2", help="Reference Antenna")
    parser.add_option("--f3", dest="file3", help="Reference Antenna")
    parser.add_option("--f4", dest="file4", help="Reference Antenna")

    parser.add_option("--l1", dest="legend1", help="Text String for Lengend Label 1")
    parser.add_option("--l2", dest="legend2", help="Text String for Lengend Label 2")
    parser.add_option("--l3", dest="legend3", help="Text String for Lengend Label 3")
    parser.add_option("--l4", dest="legend4", help="Text String for Lengend Label 4")

    (option,args) = parser.parse_args()

    if option.difference:
        refantenna = parseGEOPP(option.refAntenna)
   
        fig = plt.figure(figsize=(3.62, 2.76))
        ax = fig.add_subplot(111) 

        fig2 = plt.figure(figsize=(3.62, 2.76))
        ax2 = fig2.add_subplot(111) 

        if option.file1:
            antenna1 = parseGEOPP(option.file1)
            ax.plot(range(0,91,5),refantenna['L1PCV']-antenna1['L1PCV'])
            ax2.plot(range(0,91,5),refantenna['L2PCV']-antenna1['L2PCV'])
            #if option.printFile:
            #    ctr = 0
            #    f1=open('./fileL1.txt', 'w')
            #    f2=open('./fileL2.txt', 'w')
            #    for e in range(0,91,5):
            #        f1.write("{} {}\n".format(e, (refantenna['L1NoAZI'][ctr] - antenna1['L1NoAZI'][ctr])))
            #        f2.write("{} {}\n".format(e, (refantenna['L2NoAZI'][ctr] - antenna1['L2NoAZI'][ctr])))
            #        ctr += 1
            #    f1.close()
            #    f2.close()
        if option.file2:
            antenna2 = parseARE(option.file2)
            ax.plot(range(0,91,5),refantenna['L1PCV']-antenna2['L1PCV'])
            ax2.plot(range(0,91,5),refantenna['L2PCV']-antenna2['L2PCV'])
        if option.file3:
            antenna3 = parseARE(option.file3)
            ax.plot(range(0,91,5),refantenna['L1PCV']-antenna3['L1PCV'])
            ax2.plot(range(0,91,5),refantenna['L2PCV']-antenna3['L2PCV'])
        if option.file4:
            antenna4 = parseARE(option.file4)
            ax.plot(range(0,91,5),refantenna['L1PCV']-antenna4['L1PCV'])
            ax2.plot(range(0,91,5),refantenna['L2PCV']-antenna4['L2PCV'])

    elif option.printFile:
        antenna1 = parseGEOPP(option.file1)
        if antenna1['dazi'] < 0.001 :
            ctr = 0
            for ele in range(0,91,5) :
                print("{:02d} {}".format(ele,(2.5457*antenna1['L1PCV'][ctr] - 1.5457*antenna1['L2PCV'][ctr])))
                ctr += 1
        else:
            ictr = 0
            jctr = 0
            for az in range(0,361,5):
                jctr = 0
                for ele in range(0,91,5) :
                    print("{:02d} {:02d} {}".format(az,ele,
                            (2.5457*antenna1['L1PCV'][ictr][jctr] - 1.5457*antenna1['L2PCV'][ictr][jctr])))
                    jctr += 1
                ictr +=1
                
    else:
        antenna1 = parseGEOPP(option.file1)
   
        fig = plt.figure(figsize=(3.62, 2.76))
        ax = fig.add_subplot(111) 
        ax.plot(range(0,91,5),antenna1['L1PCV'],'b-d')
        ax.set_xlabel('Elevation Angle (degrees)')
        ax.set_ylabel('PCV (m)')

        if option.file2:
            antenna2 = parseGEOPP(option.file2)
            ax.plot(range(0,91,5),antenna2['L1PCV'],'r-^')
            ax.legend([option.legend1,option.legend2],fontsize=8,loc=0)
        if option.file3:
            antenna3 = parseGEOPP(option.file3)
            ax.plot(range(0,91,5),antenna3['L1PCV'],'k-o')
            ax.legend([option.legend1,option.legend2,option.legend3],fontsize=8,loc=0)
        if option.file4:
            antenna4 = parseGEOPP(option.file4)
            ax.plot(range(0,91,5),antenna4['L1PCV'])
            ax.legend([option.legend1,option.legend2,option.legend3,option.legend4],fontsize=8,loc=0)
    
        for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
                    ax.get_xticklabels() + ax.get_yticklabels()):
            item.set_fontsize(8)
        fig.tight_layout()

        fig2 = plt.figure(figsize=(3.62, 2.76))
        ax2 = fig2.add_subplot(111) 
        ax2.plot(range(0,91,5),antenna1['L2PCV'],'b-d')
        ax2.set_xlabel('Elevation Angle (degrees)')
        ax2.set_ylabel('PCV (mm)')

        if option.file2:
            ax2.plot(range(0,91,5),antenna2['L2PCV'],'r-^')
            ax2.legend([option.legend1,option.legend2],fontsize=8,loc=0)
        if option.file3:
            ax2.plot(range(0,91,5),antenna3['L2PCV'],'k-o')
            ax2.legend([option.legend1,option.legend2,option.legend3],fontsize=8,loc=0)
        if option.file4:
            ax2.plot(range(0,91,5),antenna4['L2PCV'])
            ax2.legend([option.legend1,option.legend2,option.legend3,option.legend4],fontsize=8,loc=0)

        for item in ([ax2.title, ax2.xaxis.label, ax2.yaxis.label] +
                    ax2.get_xticklabels() + ax2.get_yticklabels()):
            item.set_fontsize(8)

        fig2.tight_layout()

        # Now calculate the LC difference
        fig3 = plt.figure(figsize=(3.62, 2.76))
        ax3 = fig3.add_subplot(111) 

        # Only need this bit if I want to output the coords to stdout
        #ctr = 0
        #for ele in range(0,91,5) :
        #    print("{:02d} {}".format(ele,(2.5457*antenna1['L1PCV'][ctr] - 1.5457*antenna1['L2PCV'][ctr])))
        #    ctr += 1

        #ax3.plot(range(0,91,5),(2.5457*antenna1['L1PCV'] - 1.5457*antenna1['L2PCV']),'b-d')
        ax3.plot(range(0,91,5),(2545.7*antenna1['L1PCV'] - 1545.7*antenna1['L2PCV']),'b-d')
        ax3.set_xlabel('Elevation Angle (degrees)')
        ax3.set_ylabel('LC PCV (mm)')
        ax3.set_xlim([0,90])

        if option.file2:
            ax3.plot(range(0,91,5),(2.5457*antenna2['L1PCV'] - 1.5457*antenna2['L2PCV']),'r-^')
            ax3.legend([option.legend1,option.legend2],fontsize=8,loc=0)
        if option.file3:
            ax3.plot(range(0,91,5),(2.5457*antenna3['L1PCV'] - 1.5457*antenna3['L2PCV']),'k-o')
            ax3.legend([option.legend1,option.legend2,option.legend3],fontsize=8,loc=0)
        if option.file4:
            ax3.plot(range(0,91,5),(2.5457*antenna4['L1PCV'] - 1.5457*antenna4['L2NPCV']))
            ax3.legend([option.legend1,option.legend2,option.legend3,option.legend4],fontsize=8, loc=0)

        for item in ([ax3.title, ax3.xaxis.label, ax3.yaxis.label] +
                    ax3.get_xticklabels() + ax3.get_yticklabels()):
            item.set_fontsize(8)

        fig3.tight_layout()

    plt.show()
