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
                print("started the L2 obs")
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
                #L2flag = 0

    return antenna

#==============================================================================

def parseSummaryFile(sumFile) :
    '''

    parseSummaryFile(sumFile) - reads the table output from the program antdpcv

    with the following format:
    # f   el     Min      Max      Mean    AbsMean  Std.dev. AccStd.dev. [m]

    This will then be parsed into a datastructure as follows:

    data['f'][0] => el
    data['f'][1] => Min
    data['f'][2] => Max
    data['f'][3] => Mean
    data['f'][4] => AbsMean
    data['f'][5] => Std.dev.
    data['f'][6] => AccStd.dev

    where f if either L1, L2 or L0 (ionosphere free combination)

    '''

    dataBlock = {}
    dataBlock['L1'] = np.zeros((int(90./5.)+1,7))
    dataBlock['L2'] = np.zeros((int(90./5.)+1,7))
    dataBlock['L0'] = np.zeros((int(90./5.)+1,7))
    ctr = 0

    with open(sumFile) as f:
        for line in f:
            line = line.rstrip()
            data = np.array(s.split(line))
            dataBlock[data[0]][ctr] = data[1:]

            if ctr == int(90./5.):
                ctr = 0
            else:
                ctr += 1

    return dataBlock

#==============================================================================

if __name__ == "__main__":

    import matplotlib.pyplot as plt
    from matplotlib import cm

    from optparse import OptionParser

    parser = OptionParser()
    parser.add_option("-d","--diff", dest="difference",action='store_true',default=False,help="Reference Antenna")
    parser.add_option("-r","--ref", dest="refAntenna", help="Reference Antenna")
    parser.add_option("-p","--print", dest="printFile",action='store_true',default=False,help="Print Difference")
    
    parser.add_option("-q","--quiet", dest="quiet",action='store_true',default=False,help="Don't show the plots")

    parser.add_option("-f","--f1", dest="file1", help="Reference Antenna")
    parser.add_option("--f2", dest="file2", help="Reference Antenna")
    parser.add_option("--f3", dest="file3", help="Reference Antenna")
    parser.add_option("--f4", dest="file4", help="Reference Antenna")

    parser.add_option("--l1", dest="legend1", help="Text String for Lengend Label 1")
    parser.add_option("--l2", dest="legend2", help="Text String for Lengend Label 2")
    parser.add_option("--l3", dest="legend3", help="Text String for Lengend Label 3")
    parser.add_option("--l4", dest="legend4", help="Text String for Lengend Label 4")

    parser.add_option("--dre", dest="dre", help="Difference in Elevation file")
    #parser.add_option("--dre", dest="dre", action='store_true', default=False, help="Difference in Elevation file")
    parser.add_option("--drp", dest="drp", help="Difference in Elevation/Azimuth file")

    parser.add_option("--sum", dest="sum", help="Plot the summary statistics output by antdpcv --sum <file>")

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

    elif option.dre:
        antenna1 = parseGEOPP(option.dre)
        outfile = option.dre
        outfile = re.sub(r'\.dre','',outfile)
        outfile = re.sub(r'\.','_',outfile)

        # draw a horizontal line for the tollerance allowed
        toll = 1

        fig = plt.figure(figsize=(3.62, 2.76))
        ax = fig.add_subplot(111) 
        ax.plot(range(0,91,5),antenna1['L1PCV']*1000,'b-d')
        ax.plot([0,90],[toll,toll],'r--')
        ax.plot([0,90],[0,0],'k-')
        ax.plot([0,90],[-toll,-toll],'r--')
        ax.set_xlabel('Elevation Angle (degrees)')
        ax.set_ylabel('L1 PCV (mm)')
        for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
                    ax.get_xticklabels() + ax.get_yticklabels()):
            item.set_fontsize(8)
        fig.tight_layout()
        plt.savefig(outfile+'_L1.eps')

        fig2 = plt.figure(figsize=(3.62, 2.76))
        ax2 = fig2.add_subplot(111) 
        ax2.plot(range(0,91,5),antenna1['L2PCV']*1000.,'b-d')
        ax2.plot([0,90],[toll,toll],'r--')
        ax2.plot([0,90],[0,0],'k-')
        ax2.plot([0,90],[-toll,-toll],'r--')
        ax2.set_xlabel('Elevation Angle (degrees)')
        ax2.set_ylabel('L2 PCV (mm)')
        for item in ([ax2.title, ax2.xaxis.label, ax2.yaxis.label] +
                    ax2.get_xticklabels() + ax2.get_yticklabels()):
            item.set_fontsize(8)
        fig2.tight_layout()
        plt.savefig(outfile+'_L2.eps')

        # Now calculate the LC difference
        fig3 = plt.figure(figsize=(3.62, 2.76))
        ax3 = fig3.add_subplot(111) 
        ax3.plot(range(0,91,5),(2545.7*antenna1['L1PCV'] - 1545.7*antenna1['L2PCV']),'b-d')
        ax3.set_xlabel('Elevation Angle (degrees)')
        ax3.set_ylabel('LC PCV (mm)')
        ax3.set_xlim([0,90])
        for item in ([ax3.title, ax3.xaxis.label, ax3.yaxis.label] +
                    ax3.get_xticklabels() + ax3.get_yticklabels()):
            item.set_fontsize(8)
        fig3.tight_layout()
        plt.savefig(outfile+'_LC.eps')

        if option.printFile:
            ctr = 0
            f1=open('./fileL1.txt', 'w')
            f2=open('./fileL2.txt', 'w')
            fC=open('./fileLC.txt', 'w')

            for e in range(0,91,5):
                f1.write("{} {}\n".format(e, antenna1['L1PCV'][ctr]) )
                f2.write("{} {}\n".format(e, antenna1['L2PCV'][ctr]) )
                fC.write("{} {}\n".format(e, (2545.7*antenna1['L1PCV'][ctr] - 1545.7*antenna1['L2PCV'][ctr])/1000.))
                ctr = ctr + 1

            f1.close()
            f2.close()
            fC.close()
    elif option.drp:
        antenna1 = parseGEOPP(option.drp)
        outfile = option.drp
        outfile = re.sub(r'\.drp','',outfile)
        outfile = re.sub(r'\.','_',outfile)
        
        az = np.linspace(0,360,73)
        zz = np.linspace(0,90,19)

        cmap = cm.get_cmap('YlOrRd', 3)

        fig = plt.figure(figsize=(3.62, 2.76))
        ax = fig.add_subplot(111,polar=True)
        ax.set_theta_direction(-1)
        ax.set_theta_offset(np.radians(90.))
        ax.set_ylim([0,1])
        ax.set_rgrids((0.00001, np.radians(20)/np.pi*2, np.radians(40)/np.pi*2,np.radians(60)/np.pi*2,np.radians(80)/np.pi*2),labels=('0', '20', '40', '60', '80'),angle=180)

        iCtr = 0
        for a in az:
            jCtr = 0
            for z in zz:
                polar = ax.scatter(np.radians(a), np.radians(90.-z)/np.pi*2., c=np.abs(antenna1['L1PCV'][iCtr,jCtr]*1000), 
                                    s=50, alpha=1., cmap=cmap, vmin=0, vmax=3, lw=0)
                jCtr += 1
            iCtr += 1

        cbar = fig.colorbar(polar,shrink=0.75,pad=.10)
        cbar.ax.tick_params(labelsize=8)
        cbar.set_label('PCV (mm)',size=8)

        for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
            ax.get_xticklabels() + ax.get_yticklabels()):
            item.set_fontsize(8)

        plt.tight_layout()
        plt.savefig(outfile+'_polar_L1.eps')

        #======================================================================
        # L2 plots
        #======================================================================
        fig = plt.figure(figsize=(3.62, 2.76))
        ax = fig.add_subplot(111,polar=True)
        ax.set_theta_direction(-1)
        ax.set_theta_offset(np.radians(90.))
        ax.set_ylim([0,1])
        ax.set_rgrids((0.00001, np.radians(20)/np.pi*2, np.radians(40)/np.pi*2,np.radians(60)/np.pi*2,np.radians(80)/np.pi*2),labels=('0', '20', '40', '60', '80'),angle=180)

        iCtr = 0
        for a in az:
            jCtr = 0
            for z in zz:
                polar = ax.scatter(np.radians(a), np.radians(90.-z)/np.pi*2., c=np.abs(antenna1['L2PCV'][iCtr,jCtr]*1000), 
                                    s=50, alpha=1., cmap=cmap, vmin=0, vmax=3, lw=0)
                jCtr += 1
            iCtr += 1

        cbar = fig.colorbar(polar,shrink=0.75,pad=.10)
        cbar.ax.tick_params(labelsize=8)
        cbar.set_label('PCV (mm)',size=8)

        for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
            ax.get_xticklabels() + ax.get_yticklabels()):
            item.set_fontsize(8)

        plt.tight_layout()
        plt.savefig(outfile+'_polar_L2.eps')

        #======================================================================
        # LC polar plot
        #======================================================================
        fig = plt.figure(figsize=(3.62, 2.76))
        ax = fig.add_subplot(111,polar=True)
        ax.set_theta_direction(-1)
        ax.set_theta_offset(np.radians(90.))
        ax.set_ylim([0,1])
        ax.set_rgrids((0.00001, np.radians(20)/np.pi*2, np.radians(40)/np.pi*2,np.radians(60)/np.pi*2,np.radians(80)/np.pi*2),labels=('0', '20', '40', '60', '80'),angle=180)

        iCtr = 0
        for a in az:
            jCtr = 0
            for z in zz:
                l1 = antenna1['L1PCV'][iCtr,jCtr] * 2545.7
                l2 = antenna1['L2PCV'][iCtr,jCtr] * 1545.7
                lc = l1 - l2
                #print("L1 L2 LC:",l1, l2, lc)
                polar = ax.scatter(np.radians(a), np.radians(90.-z)/np.pi*2., c=np.abs(lc), 
                                    s=50, alpha=1., cmap=cmap, vmin=0, vmax=3, lw=0)
                jCtr += 1
            iCtr += 1

        cbar = fig.colorbar(polar,shrink=0.75,pad=.10)
        cbar.ax.tick_params(labelsize=8)
        cbar.set_label('PCV (mm)',size=8)

        for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
            ax.get_xticklabels() + ax.get_yticklabels()):
            item.set_fontsize(8)

        plt.tight_layout()
        plt.savefig(outfile+'_polar_LC.eps')

        #======================================================================
        # Summary plots
        #======================================================================

        # Plot the mean, minimum and maximum values vs elevation
        # Plot the L1 frequency
        fig1 = plt.figure(figsize=(3.62, 2.76))
        ax = fig1.add_subplot(111)

        # set the tollerance to be 1mm
        toll = 1

        # get the mean for each elevation bin
        antenna1['L1_Ele_Mean'] = np.mean(antenna1['L1PCV'],axis=0)
        antenna1['L1_Ele_Std']  = np.std(antenna1['L1PCV'],axis=0)
        antenna1['L1_Ele_Min']  = np.min(antenna1['L1PCV'],axis=0)
        antenna1['L1_Ele_Max']  = np.max(antenna1['L1PCV'],axis=0)

        antenna1['L2_Ele_Mean'] = np.mean(antenna1['L2PCV'],axis=0)
        antenna1['L2_Ele_Std']  = np.std(antenna1['L2PCV'],axis=0)
        antenna1['L2_Ele_Min']  = np.min(antenna1['L2PCV'],axis=0)
        antenna1['L2_Ele_Max']  = np.max(antenna1['L2PCV'],axis=0)

        ele = range(0,91,5)
        ax.errorbar(ele , antenna1['L1_Ele_Mean']*1000, yerr=antenna1['L1_Ele_Std']*1000., fmt='-o', color='k',ecolor='k',capsize=3)
        ax.plot( ele, antenna1['L1_Ele_Min']*1000,color='k',alpha=0.5 )
        ax.plot( ele, antenna1['L1_Ele_Max']*1000,color='k',alpha=0.5 )

        # plot the tolerances acceptable by the ATWG
        ax.plot([-5,90],[toll,toll],'r--')
        ax.plot([-5,90],[-toll,-toll],'r--')
        ax.set_xlim([-5,90])
        ax.set_xlabel('Elevation Angle (degrees)')
        ax.set_ylabel('Difference in Absolute PCV L1 (mm)')
        for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
            ax.get_xticklabels() + ax.get_yticklabels()):
            item.set_fontsize(8)

        plt.tight_layout()
        plt.savefig(outfile+'_summary_L1.eps')

        #=======================================
        # Plot the L2 frequency
        fig1 = plt.figure(figsize=(3.62, 2.76))
        ax = fig1.add_subplot(111)

        ax.errorbar(ele , antenna1['L2_Ele_Mean']*1000, yerr=antenna1['L2_Ele_Std']*1000., fmt='-o', color='k',ecolor='k',capsize=3)
        ax.plot( ele, antenna1['L2_Ele_Min']*1000,color='k',alpha=0.5 )
        ax.plot( ele, antenna1['L2_Ele_Max']*1000,color='k',alpha=0.5 )

        # plot the tolerances acceptable by the ATWG
        ax.plot([-5,90],[toll,toll],'r--')
        ax.plot([-5,90],[-toll,-toll],'r--')
        ax.set_xlim([-5,90])
        ax.set_xlabel('Elevation Angle (degrees)')
        ax.set_ylabel('Difference in Absolute PCV L2 (mm)')
        for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
            ax.get_xticklabels() + ax.get_yticklabels()):
            item.set_fontsize(8)

        plt.tight_layout()
        plt.savefig(outfile+'_summary_L2.eps')

        #========================================
        # Plot a histogram of the differences
        #========================================
        fig1 = plt.figure(figsize=(3.62, 2.76))
        ax = fig1.add_subplot(111)

        vec = antenna1['L1PCV'].ravel() * 1000
        weights = np.ones_like(vec)/len(vec)
        n, bins, rectangles = ax.hist(vec, bins=11, histtype='bar',weights=weights)

        # Work out how many are above the guidelines for comparisons...
        criterion = ((vec[:] > 1 ) | (vec[:] < -1))        
        ind = np.array(np.where(criterion))
        nonConform = ind.size/vec.size
        titleStr = "Percentage of obs above 1 mm:{:.1f} %".format(nonConform*100)
        ax.set_title(titleStr)

        ax.plot([-1,-1],[0,1],'r--')
        ax.plot([1,1],[0,1],'r--')

        ax.set_xlabel('Difference (mm)')
        ax.set_ylabel('Percentage of Observations')
        ax.set_ylim([0, np.max(n)*1.5])
        for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
            ax.get_xticklabels() + ax.get_yticklabels()):
            item.set_fontsize(8)

        plt.tight_layout()
        plt.savefig(outfile+'_histogram_L1.eps')

        #=======================================
        fig1 = plt.figure(figsize=(3.62, 2.76))
        ax = fig1.add_subplot(111)

        vec = antenna1['L2PCV'].ravel() * 1000
        weights = np.ones_like(vec)/len(vec)
        n, bins, rectangles = ax.hist(vec, bins=11, histtype='bar',weights=weights)

        # Work out how many are above the guidelines for comparisons...
        criterion = ((vec[:] > 1 ) | (vec[:] < -1))        
        ind = np.array(np.where(criterion))
        nonConform = ind.size/vec.size
        titleStr = "Percentage of obs above 1 mm:{:.1f} %".format(nonConform*100)
        ax.set_title(titleStr)

        ax.plot([-1,-1],[0,1],'r--')
        ax.plot([1,1],[0,1],'r--')

        ax.set_xlabel('Difference (mm)')
        ax.set_ylabel('Percentage of Observations')
        ax.set_ylim([0, np.max(n)*1.5])
        for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
            ax.get_xticklabels() + ax.get_yticklabels()):
            item.set_fontsize(8)

        plt.tight_layout()
        plt.savefig(outfile+'_histogram_L2.eps')

        #=======================================
        # Work out the percentage of misfit per elevation angle
        #=======================================
        fig1 = plt.figure(figsize=(3.62, 2.76))
        ax = fig1.add_subplot(111)
        
        # Work out how many are above the guidelines for comparisons...
        jCtr = 0
        for z in zz:
            criterion = ((antenna1['L1PCV'][:,jCtr] > 0.001) | (antenna1['L1PCV'][:,jCtr] < -0.001 ))
            ind = np.array(np.where(criterion))
            percent = (ind.size / az.size) * 100.
            p1 = plt.bar(z, percent, 3.5, color='r')
            p2 = plt.bar(z, 100-percent, 3.5, color='b', bottom=percent) 
            jCtr = jCtr + 1

        ax.legend( (p1[0],p2[0]),('Above 1mm','Below 1mm'),prop={'size':8})
        ax.set_xlabel('Elevation Angle (degrees)')
        ax.set_ylabel('Percentage of Observations')
        ax.set_xlim([0,90])

        for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
            ax.get_xticklabels() + ax.get_yticklabels()):
            item.set_fontsize(8)

        plt.tight_layout()
        plt.savefig(outfile+'_histogram_ele_L1.eps')

        fig1 = plt.figure(figsize=(3.62, 2.76))
        ax = fig1.add_subplot(111)
        
        # Work out how many are above the guidelines for comparisons...
        jCtr = 0
        for z in zz:
            criterion = ((antenna1['L2PCV'][:,jCtr] > 0.001) | (antenna1['L2PCV'][:,jCtr] < -0.001 ))
            ind = np.array(np.where(criterion))
            percent = (ind.size / az.size) * 100.
            p1 = plt.bar(z, percent, 3.5, color='r')
            p2 = plt.bar(z, 100-percent, 3.5, color='b', bottom=percent) 
            jCtr = jCtr + 1

        ax.legend( (p1[0],p2[0]),('Above 1mm','Below 1mm'),prop={'size':8})
        ax.set_xlabel('Elevation Angle (degrees)')
        ax.set_ylabel('Percentage of Observations')
        ax.set_xlim([0,90])

        for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
            ax.get_xticklabels() + ax.get_yticklabels()):
            item.set_fontsize(8)

        plt.tight_layout()
        plt.savefig(outfile+'_histogram_ele_L2.eps')

        #======================================================================
        # need an option to interpolate this onto a 0.1 grid, to be block median'd later in a ppp simulation
        #======================================================================
        if option.printFile:
            ctr = 0
            f1=open('./fileL1.txt', 'w')
            f2=open('./fileL2.txt', 'w')
            fC=open('./fileLC.txt', 'w')
            for a in range(0,361,5):
                for e in range(0,91,5):
                    f1.write("{:.2f} {:.2f} {}\n".format(a, e, antenna1['L1PCV'][ctr]) )
                    f2.write("{:.2f} {:.2f} {}\n".format(a, e, antenna1['L2PCV'][ctr]) )
                    fC.write("{:.2f} {:.2f} {}\n".format(a, e, (2545.7*antenna1['L1PCV'][ctr] - 1545.7*antenna1['L2PCV'][ctr])/1000.))
                    ctr = ctr + 1

            f1.close()
            f2.close()
            fC.close()

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

    elif option.sum:
        outfile = option.sum
        outfile = re.sub(r'\.sum','',outfile)
        outfile = re.sub(r'\.','_',outfile)
        # f   el     Min      Max      Mean    AbsMean  Std.dev. AccStd.dev. [m]
        data = parseSummaryFile(option.sum)

        # set a tollerance of 1 mm
        toll = 1

        # Plot the mean, minimum and maximum values vs elevation
        # Plot the L1 frequency
        fig1 = plt.figure(figsize=(3.62, 2.76))
        ax = fig1.add_subplot(111)

        # 0 => el , 1 -> minumum, 2 -> maximum , 3 -> mean, 4 -> AbsMean, 5 -> Std.Dev, 6-> Acc Std. dev 
        ax.errorbar( data['L1'][:,0], data['L1'][:,3]*1000, yerr=data['L1'][:,5]/2*1000., fmt='-o', color='k',ecolor='k',capsize=3)
        ax.plot( data['L1'][:,0], data['L1'][:,1]*1000,color='k',alpha=0.5 )
        ax.plot( data['L1'][:,0], data['L1'][:,2]*1000,color='k',alpha=0.5 )

        # plot the tolerances acceptable by the ATWG
        ax.plot([-5,90],[toll,toll],'r--')
        ax.plot([-5,90],[-toll,-toll],'r--')
        ax.set_xlim([-5,90])
        ax.set_xlabel('Elevation Angle (degrees)')
        ax.set_ylabel('Difference in Absolute PCV (mm)')
        for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
            ax.get_xticklabels() + ax.get_yticklabels()):
            item.set_fontsize(8)

        plt.tight_layout()
        plt.savefig(outfile+'_summary_L1.eps')

        # Plot the L2 frequency
        fig2 = plt.figure(figsize=(3.62, 2.76))
        ax2 = fig2.add_subplot(111)

        ax2.errorbar( data['L2'][:,0], data['L2'][:,3]*1000, yerr=data['L2'][:,5]/2*1000., fmt='-o', color='k',ecolor='k',capsize=3)
        ax2.plot( data['L2'][:,0], data['L2'][:,1]*1000,color='k',alpha=0.5 )
        ax2.plot( data['L2'][:,0], data['L2'][:,2]*1000,color='k',alpha=0.5 )

        ax2.plot([-5,90],[toll,toll],'r--')
        ax2.plot([-5,90],[-toll,-toll],'r--')
        ax2.set_xlim([-5,90])
        ax2.set_xlabel('Elevation Angle (degrees)')
        ax2.set_ylabel('Difference in Absolute PCV (mm)')
        for item in ([ax2.title, ax2.xaxis.label, ax2.yaxis.label] +
            ax2.get_xticklabels() + ax2.get_yticklabels()):
            item.set_fontsize(8)
        plt.tight_layout()
        plt.savefig(outfile+'_summary_L2.eps')

        # Plot the LC frequency
        fig3 = plt.figure(figsize=(3.62, 2.76))
        ax3 = fig3.add_subplot(111)

        ax3.errorbar( data['L0'][:,0], data['L0'][:,3]*1000, yerr=data['L0'][:,5]/2*1000., fmt='-o', color='k',ecolor='k',capsize=3)
        ax3.plot( data['L0'][:,0], data['L0'][:,1]*1000,color='k',alpha=0.5 )
        ax3.plot( data['L0'][:,0], data['L0'][:,2]*1000,color='k',alpha=0.5 )
        
        ax3.plot([-5,90],[toll,toll],'r--')
        ax3.plot([-5,90],[-toll,-toll],'r--')

        ax3.set_xlim([-5,90])
        ax3.set_xlabel('Elevation Angle (degrees)')
        ax3.set_ylabel('Difference in Absolute PCV (mm)')

        for item in ([ax3.title, ax3.xaxis.label, ax3.yaxis.label] +
            ax3.get_xticklabels() + ax3.get_yticklabels()):
            item.set_fontsize(8)
        plt.tight_layout()
        plt.savefig(outfile+'_summary_L0.eps')

        #==============================
        # Plot the absolute mean and accumulated standard Deviation
        fig4 = plt.figure(figsize=(3.62, 2.76))
        ax4 = fig4.add_subplot(111)
        # 0 => el , 1 -> minumum, 2 -> maximum , 3 -> mean, 4 -> AbsMean, 5 -> Std.Dev, 6-> Acc Std. dev 
        ax4.errorbar( data['L1'][:,0], data['L1'][:,4]*1000, yerr=data['L1'][:,6]*1000., fmt='-o', color='k',ecolor='k',capsize=3)
        ax4.plot([-5,90],[toll,toll],'r--')
        ax4.set_xlim([-5,90])
        ax4.set_xlabel('Elevation Angle (degrees)')
        ax4.set_ylabel('Absolute Difference in Absolute PCV (mm)')
        for item in ([ax4.title, ax4.xaxis.label, ax4.yaxis.label] +
            ax4.get_xticklabels() + ax4.get_yticklabels()):
            item.set_fontsize(8)
        plt.tight_layout()
        plt.savefig(outfile+'_abs_summary_L1.eps')

        fig5 = plt.figure(figsize=(3.62, 2.76))
        ax5 = fig5.add_subplot(111)
        ax5.errorbar( data['L2'][:,0], data['L2'][:,4]*1000, yerr=data['L2'][:,6]*1000., fmt='-o', color='k',ecolor='k',capsize=3)
        ax5.plot([-5,90],[toll,toll],'r--')
        ax5.set_xlim([-5,90])
        ax5.set_xlabel('Elevation Angle (degrees)')
        ax5.set_ylabel('Absolute Difference in Absolute PCV (mm)')
        for item in ([ax5.title, ax5.xaxis.label, ax5.yaxis.label] +
            ax5.get_xticklabels() + ax5.get_yticklabels()):
            item.set_fontsize(8)
        plt.tight_layout()
        plt.savefig(outfile+'_abs_summary_L2.eps')

        fig6 = plt.figure(figsize=(3.62, 2.76))
        ax6 = fig6.add_subplot(111)
        ax6.errorbar( data['L0'][:,0], data['L0'][:,4]*1000, yerr=data['L0'][:,6]*1000., fmt='-o', color='k',ecolor='k',capsize=3)
        ax6.plot([-5,90],[toll,toll],'r--')
        ax6.set_xlim([-5,90])
        ax6.set_xlabel('Elevation Angle (degrees)')
        ax6.set_ylabel('Absolute Difference in Absolute PCV (mm)')
        for item in ([ax6.title, ax6.xaxis.label, ax6.yaxis.label] +
            ax6.get_xticklabels() + ax6.get_yticklabels()):
            item.set_fontsize(8)
        plt.tight_layout()
        plt.savefig(outfile+'_abs_summary_L0.eps')

        # f   el     Min      Max      Mean    AbsMean  Std.dev. AccStd.dev. [m]
        # Save the data block as a LATEX table
        latexFile = outfile + '.tex'
        f1=open(latexFile, 'w') 
        print('\\begin{table}[htbp]',file=f1)
        print('\\centering',file=f1)
        print('\\caption{}',file=f1)
        print('\\label{tab:}',file=f1)
        print('\\begin{footnotesize}',file=f1)
        print('\\begin{tabular}{l c c c c c c c}',file=f1)
        print('\\hline\\noalign{\\smallskip}',file=f1)
        print('f & el  &  Min  &   Max  &   Mean &  AbsMean& Std.dev.&AccStd.dev. [mm]\\\\',file=f1)
        print('\\hline',file=f1)

        for key in ['L1','L2','L0']:
            i = 0

            for el in data[key][:,0]:
                print('{:2s} & {:0.1f} & {:3.2f} & {:3.2f} & {:3.2f} & {:3.2f} & {:3.2f} & {:3.2f} \\\\'.format(key, el, 
                                                                                                           data[key][i,1]*1000.,
                                                                                                           data[key][i,2]*1000.,
                                                                                                           data[key][i,3]*1000.,
                                                                                                           data[key][i,4]*1000.,
                                                                                                           data[key][i,5]*1000.,
                                                                                                           data[key][i,6]*1000.,
                                                                                                           ),file=f1 )
                i += 1

            print('\\hline',file=f1)

        print('\\end{tabular}',file=f1)
        print('\\end{footnotesize}',file=f1)
        print('\\end{table}',file=f1)
        f1.close()
    else:
        antenna1 = parseGEOPP(option.file1)
   
        fig = plt.figure(figsize=(3.62, 2.76))
        ax = fig.add_subplot(111) 
        ax.plot(range(0,91,5),antenna1['L1PCV']*1000,'b-d')
        ax.set_xlabel('Elevation Angle (degrees)')
        ax.set_ylabel('L1 PCV (mm)')

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
        ax2.plot(range(0,91,5),antenna1['L2PCV']*1000.,'b-d')
        ax2.set_xlabel('Elevation Angle (degrees)')
        ax2.set_ylabel('L2 PCV (mm)')

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

    if not option.quiet : 
        plt.show()
