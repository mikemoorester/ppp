#!/usr/bin/env python
from __future__ import division, print_function, absolute_import

import numpy as np
import re

def blockMedian(residualFile, gridSpacing, bType):
    """
        bMedian = blockMedian(residualFile,gridSpacing, bType)

        where,

            residualfile => file path  MOBS_DPH_2012_001_143.all
        
            gridSpacing  => float ie every 0.5 degrees
            
            bType => (int) 0 elevation only
                          1 az/elevation
                          2 azimuth only  

        output:
            bType 0 => bMedian ia a 1d matrix of shape [nZd]
            bType 1 => bMedian is a 2d matrix of shape [nAz,nZd]
            bType 2 => bMedian is a 1d matrix of shape [nAz]

        Example:

            bMedian = blockMedian('./MOBS_DPH_2012_001_143.all',0.5,1)

            bMedian.shape => 721,181
    """

    data = np.genfromtxt(residualFile)

    # Calculate the block median
    zz = np.linspace(0,90,int(90/gridSpacing)+1)
    az = np.linspace(0,360,int(360/gridSpacing)+1)

    # sum of medin diffences used to calculate the rms
    Sdiff = 0.
    SdiffCtr = 0

    if bType == 0:
        bMedian = np.zeros(zz.size)
        bMedianStd = np.zeros(zz.size)
        jCtr = 0
        for j in zz: 
            criterion = (data[:,1] < j + gridSpacing/2.) & (data[:,1] > j - gridSpacing/2. ) 
            ind = np.array(np.where(criterion))
       
            if ind.size > 0 :
                bMedian[jCtr] = np.median(data[ind,2])
                bMedianStd[jCtr] = np.std(data[ind,2])
                Sdiff += bMedianStd[jCtr]**2 
                SdiffCtr += 1
            else :
                bMedian[jCtr] = float('NaN')
            jCtr += 1

        rms = Sdiff * 1./SdiffCtr
        rms = np.sqrt(rms)
        return bMedian,bMedianStd,rms

    elif bType == 1:
        iCtr = 0
        bMedian = np.zeros((az.size,zz.size))
        bMedianStd = np.zeros((az.size,zz.size))

        for i in az:
            jCtr = 0
            criterion = (data[:,0] < i + gridSpacing/2.) & (data[:,0] > i - gridSpacing/2. )
            ind = np.array(np.where(criterion))
    
            if ind.size == 0:
                for j in zz :
                    bMedian[iCtr,jCtr] = float('NaN')
                    jCtr += 1
            else: 
                tmp = data[ind,:]
                tmp = tmp.reshape(tmp.shape[1],tmp.shape[2])

                jCtr = 0
                for j in zz:
                    criterion = (tmp[:,1] < j + 0.25) & (tmp[:,1] > j - 0.25 ) 
                    indZ = np.array(np.where(criterion))
    
                    if indZ.size > 0 :
                        bMedian[iCtr,jCtr] = np.median(tmp[indZ,2])
                        bMedianStd[iCtr,jCtr] = np.std(tmp[indZ,2])
                        if bMedian[iCtr,jCtr] != float('NaN'):
                            Sdiff += np.median(tmp[indZ,2])**2 
                            SdiffCtr += 1
                    else :
                        bMedian[iCtr,jCtr] = float('NaN') 

                    jCtr += 1

                jCtr = 0
                iCtr += 1

        rms = Sdiff * 1./SdiffCtr
        rms = np.sqrt(rms)
        return bMedian, bMedianStd, rms

    return False 

def parseDPH(dphFile) :
    """
    dph = parseDPH(dphFile)
    Read in a GAMIT undifferenced phase residual file.
    Return a DPH structurea

    Will skip any lines in the file which contain a '*' 
    at any column number

    Checks there are no comments in the first column of the file

    """

    asterixRGX = re.compile('\*')

    dph = {}

    obs = {}
    #obs['dphs'] = []
    obs['satsViewed'] = set()
    obs['epochs'] = set()

    debug = 0

    with open(dphFile) as f:
        for line in f:
            dph = {}
            if line[0] != ' ':
                if debug :
                    print('A comment',line)
            elif asterixRGX.search(line):
                if debug :
                    print('Bad observation',line)
            else :
                dph['epoch'] = int(line[1:5])
                dph['l1cyc'] = float(line[6:15])
                dph['l2cyc'] = float(line[16:24])
                dph['p1cyc'] = float(line[25:33])
                dph['p2cyc'] = float(line[34:42])
                dph['lccyc'] = float(line[43:51])
                dph['lgcyc'] = float(line[52:60])
                dph['pccyc'] = float(line[61:69])
                dph['wlcyc'] = float(line[70:78])
                dph['ncyc']  = float(line[79:87])
                dph['lsv']   = int(line[88:91])
                dph['az']    = float(line[94:102])
                dph['el']    = float(line[105:112])
                dph['pf']    = int(line[113:114])
                dph['dataf'] = int(line[115:127])

                if str(line[128:148]).strip() != '' :
                    dph['L1cycles'] = float(line[128:148])
                if str(line[149:169]).strip() != '' :
                    dph['L2cycles'] = float(line[149:169])

                dph['prn']   = int(line[170:172])
                prnSTR = 'prn_'+str(dph['prn'])

                # store the data in lists accessed by the sat prn key
                if dph['prn'] in obs['satsViewed'] :
                    obs[prnSTR].append(dph)
                else:
                    obs[prnSTR] = []
                    obs[prnSTR].append(dph)

                # Keep a record of each satellite viewed at each epoch in a set
                epochStr = str(dph['epoch'])
                if dph['epoch'] in obs['epochs']:
                    obs[epochStr].add(dph['prn'])
                else :
                    obs['epochs'].add(dph['epoch'])
                    obs[epochStr] = set()
                    obs[epochStr].add(dph['prn'])    

                #obs['dphs'].append(dph)

                # keep a record of all the unique satellies which have residuals
                obs['satsViewed'].add(dph['prn'])
                   
    return obs 


#===========================================================================

if __name__ == "__main__":

    from matplotlib import pyplot as plt
    from matplotlib import cm 

    #===================================
    from optparse import OptionParser

    parser = OptionParser()
    parser.add_option("-f", "--filename", dest="filename", help="Result file to plot")
    parser.add_option("-E", "--plotElevation", dest="elevationPlot",action='store_true',default=False,
                        help="Plot Residuals vs Elevation Angle")
    parser.add_option("-P", "--plotPolar", dest="polarPlot",action='store_true',default=False,
                        help="Polar Plot Residuals vs Azimuth & Elevation Angle")
    parser.add_option("--esm","--ESM",dest="esmFilename",help="Example Residual file from which to create an ESM")

    parser.add_option("--dph",dest="dphFilename",help="DPH filename to parse, obtained from GAMIT") 

    (option, args) = parser.parse_args()

    #===================================
    
    if option.dphFilename :
        dphs = parseDPH(option.dphFilename)
    
        fig = plt.figure(figsize=(3.62, 2.76))    
        ax = fig.add_subplot(111)

        for dph in dphs['prn_1']: 
            ax.scatter(dph['epoch'],dph['el'])

        ax.set_xlabel('Elevation Angle (degrees)',fontsize=8)
        ax.set_ylabel('Bias (mm)',fontsize=8)
        ax.set_ylim([0, 90])
        ax.set_xlim([0, 2880])

        for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
                                     ax.get_xticklabels() + ax.get_yticklabels()):
            item.set_fontsize(8)

        plt.tight_layout()
        plt.show() 
    
    # Calculate the block median
    zz = np.linspace(0,90,181)

    if option.elevationMedianPlot :
        med,medStd,medrms = blockMedian(option.filename,0.5,0)
        # flip the array around so that is in order of elevation angle
        ele = 90. - zz[::-1]
        med = med[::-1]

        # Do an elevation only plot
        fig = plt.figure(figsize=(3.62, 2.76))
        ax = fig.add_subplot(111)
        ax.plot(ele,np.median(med))
        ax.set_xlabel('Elevation Angle (degrees)',fontsize=8)
        ax.set_ylabel('Bias (mm)',fontsize=8)
        ax.set_ylim([-7.5, 7.5])
        for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
                                     ax.get_xticklabels() + ax.get_yticklabels()):
            item.set_fontsize(8)

        plt.tight_layout()
    
        #fig.savefig('MOBS_Elevation_Median.eps')


    # Create a polar plot of the residuals
    if option.polarPlot:
        blkMedian,blkMedianStd,rms = blockMedian(option.filename,0.5,1)
        az = np.linspace(0,360,721)
        #fig = plt.figure()
        fig = plt.figure(figsize=(3.62, 2.76))

        ax = fig.add_subplot(111,polar=True)
        ax.set_theta_direction(-1)
        ax.set_theta_offset(np.radians(90.))
        ax.set_ylim([0,1])
        ax.set_rgrids((0.00001, np.radians(20)/np.pi*2, np.radians(40)/np.pi*2,np.radians(60)/np.pi*2,np.radians(80)/np.pi*2), 
                    labels=('0', '20', '40', '60', '80'),angle=180)

        ma,mz = np.meshgrid(az,zz,indexing='ij')
        ma = ma.reshape(ma.size,)
        mz = mz.reshape(mz.size,)
        #polar = ax.scatter(np.radians(ma), np.radians(mz)/np.pi*2., c=blkMedian ,s=1,alpha=1., cmap=cm.RdBu,vmin=-15,vmax=15, lw=0)
        polar = ax.scatter(np.radians(ma), np.radians(mz)/np.pi*2., c=blkMedian ,s=1,alpha=1., cmap=cm.RdBu,vmin=-10,vmax=10, lw=0)
   
        cbar = fig.colorbar(polar,shrink=0.75,pad=.10)
        cbar.ax.tick_params(labelsize=8)
        cbar.set_label('Residuals (mm)',size=8)

        for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
                    ax.get_xticklabels() + ax.get_yticklabels()):
            item.set_fontsize(8)

        plt.tight_layout()
   
        # Print out the ratio if the elvation plot has been selected as well
        if option.elevationPlot:
            ratio = rms/medrms
            print('{} {:.3f} {:.3f} {:.2f}').format(option.filename,medrms,rms,ratio)

    if option.polarPlot | option.elevationPlot :
        plt.show()

    if option.esmFilename :
        esm,esmStd,esmrms = blockMedian(option.esmFilename,0.5,1)