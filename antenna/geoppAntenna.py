#!/usr/bin/env python
from __future__ import division, print_function, absolute_import

import numpy as np
import re
import string as s


def parseARE(areFile):
    '''
    parseARE(antmod.are)
    -parse the geopp .are file

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
    antenna = {}

    with open(areFile) as f:
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
            elif varL2RGX.search(line) :
                L2flag = 1
            elif L1flag == 1:
                L1flag = 0
                tmp = np.array(s.split(line))
                antenna['L1NoAZI'] = tmp.astype(np.float)
            elif L2flag == 1:
                L2flag = 0
                tmp = np.array(s.split(line))
                antenna['L2NoAZI'] = tmp.astype(np.float)

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
        refantenna = parseARE(option.refAntenna)
   
        fig = plt.figure(figsize=(3.62, 2.76))
        ax = fig.add_subplot(111) 

        fig2 = plt.figure(figsize=(3.62, 2.76))
        ax2 = fig2.add_subplot(111) 

        if option.file1:
            antenna1 = parseARE(option.file1)
            ax.plot(range(0,91,5),refantenna['L1NoAZI']-antenna1['L1NoAZI'])
            ax2.plot(range(0,91,5),refantenna['L2NoAZI']-antenna1['L2NoAZI'])
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
            ax.plot(range(0,91,5),refantenna['L1NoAZI']-antenna2['L1NoAZI'])
            ax2.plot(range(0,91,5),refantenna['L2NoAZI']-antenna2['L2NoAZI'])
        if option.file3:
            antenna3 = parseARE(option.file3)
            ax.plot(range(0,91,5),refantenna['L1NoAZI']-antenna3['L1NoAZI'])
            ax2.plot(range(0,91,5),refantenna['L2NoAZI']-antenna3['L2NoAZI'])
        if option.file4:
            antenna4 = parseARE(option.file4)
            ax.plot(range(0,91,5),refantenna['L1NoAZI']-antenna4['L1NoAZI'])
            ax2.plot(range(0,91,5),refantenna['L2NoAZI']-antenna4['L2NoAZI'])

    else:
        antenna1 = parseARE(option.file1)
   
        fig = plt.figure(figsize=(3.62, 2.76))
        ax = fig.add_subplot(111) 
        ax.plot(range(0,91,5),antenna1['L1NoAZI'],'b-d')
        ax.set_xlabel('Elevation Angle (degrees)')
        ax.set_ylabel('PCV (m)')

        if option.file2:
            antenna2 = parseARE(option.file2)
            ax.plot(range(0,91,5),antenna2['L1NoAZI'],'r-^')
            ax.legend([option.legend1,option.legend2],fontsize=8,loc=0)
        if option.file3:
            antenna3 = parseARE(option.file3)
            ax.plot(range(0,91,5),antenna3['L1NoAZI'],'k-o')
            ax.legend([option.legend1,option.legend2,option.legend3],fontsize=8,loc=0)
        if option.file4:
            antenna4 = parseARE(option.file4)
            ax.plot(range(0,91,5),antenna4['L1NoAZI'])
            ax.legend([option.legend1,option.legend2,option.legend3,option.legend4],fontsize=8,loc=0)
    
        for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
                    ax.get_xticklabels() + ax.get_yticklabels()):
            item.set_fontsize(8)
        fig.tight_layout()


        fig2 = plt.figure(figsize=(3.62, 2.76))
        ax2 = fig2.add_subplot(111) 
        ax2.plot(range(0,91,5),antenna1['L2NoAZI'],'b-d')
        ax2.set_xlabel('Elevation Angle (degrees)')
        ax2.set_ylabel('PCV (mm)')

        if option.file2:
            ax2.plot(range(0,91,5),antenna2['L2NoAZI'],'r-^')
            ax2.legend([option.legend1,option.legend2],fontsize=8,loc=0)
        if option.file3:
            ax2.plot(range(0,91,5),antenna3['L2NoAZI'],'k-o')
            ax2.legend([option.legend1,option.legend2,option.legend3],fontsize=8,loc=0)
        if option.file4:
            ax2.plot(range(0,91,5),antenna4['L2NoAZI'])
            ax2.legend([option.legend1,option.legend2,option.legend3,option.legend4],fontsize=8,loc=0)

        for item in ([ax2.title, ax2.xaxis.label, ax2.yaxis.label] +
                    ax2.get_xticklabels() + ax2.get_yticklabels()):
            item.set_fontsize(8)

        fig2.tight_layout()

        # Now calculate the LC difference
        fig3 = plt.figure(figsize=(3.62, 2.76))
        ax3 = fig3.add_subplot(111) 
        ctr = 0
        for ele in range(0,91,5) :
            print("{:02d} {}".format(ele,(2.5457*antenna1['L1NoAZI'][ctr] - 1.5457*antenna1['L2NoAZI'][ctr])))
            ax3.scatter(ele,(2.5457*antenna1['L1NoAZI'][ctr] - 1.5457*antenna1['L2NoAZI'][ctr]))
            ctr += 1

        ax3.plot(range(0,91,5),(2.5457*antenna1['L1NoAZI'] - 1.5457*antenna1['L2NoAZI']),'b-d')
        ax3.set_xlabel('Elevation Angle (degrees)')
        ax3.set_ylabel('PCV (mm)')
        ax3.set_xlim([0,90])

        if option.file2:
            ax3.plot(range(0,91,5),(2.5457*antenna2['L1NoAZI'] - 1.5457*antenna2['L2NoAZI']),'r-^')
            ax3.legend([option.legend1,option.legend2],fontsize=8,loc=0)
        if option.file3:
            ax3.plot(range(0,91,5),(2.5457*antenna3['L1NoAZI'] - 1.5457*antenna3['L2NoAZI']),'k-o')
            ax3.legend([option.legend1,option.legend2,option.legend3],fontsize=8,loc=0)
        if option.file4:
            ax3.plot(range(0,91,5),(2.5457*antenna4['L1NoAZI'] - 1.5457*antenna4['L2NoAZI']))
            ax3.legend([option.legend1,option.legend2,option.legend3,option.legend4],fontsize=8, loc=0)

        for item in ([ax3.title, ax3.xaxis.label, ax3.yaxis.label] +
                    ax3.get_xticklabels() + ax3.get_yticklabels()):
            item.set_fontsize(8)

        fig3.tight_layout()

    plt.show()
