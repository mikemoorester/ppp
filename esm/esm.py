#!/usr/bin/env python
from __future__ import division, print_function, absolute_import

import numpy as np
from scipy import interpolate

import antenna as ant
import residuals as res


#=====================================
if __name__ == "__main__":

    import argparse

    parser = argparse.ArgumentParser(prog='antenna',description='Parse ANTEX files')

    parser.add_argument('-f', '--file', dest='file1', default='./t/antmod.dat')
    parser.add_argument('-g', '--grid', dest='grid', default=5., type=float)
    parser.add_argument('-r', '--resfile', dest='resfile', default='./t/MOBS_DPH_2012_001_143.all')
    parser.add_argument('-t', '--AntType',dest='AntType', default='ASH701945C_M    NONE')

    #parser.add_argument('-p', '--plot',dest='plot', default=False, action='store_true')
    parser.add_argument('--polar',dest='polar', default=False, action='store_true')
    parser.add_argument('--elevation',dest='elevation', default=False, action='store_true')
    
    # Option to create a summary file of all the DPH residuals 
    parser.add_argument('--dph',dest='dphFile', description='Add a DPH file to a Summary file for a station')
    args = parser.parse_args()

    # get the antenna information from an antex file
    antennas = ant.parseANTEX(args.file1)

    # keep the antenna type required
    antenna = ant.antennaType(args.AntType,antennas)
    L1_data = np.array(antenna['data'][0])
    L2_data = np.array(antenna['data'][1])

    # compute the block median from a set of raw residuals
    # 0.5 should be variable that can be input into the program
    med,medStd,medrms = res.blockMedian(args.resfile,0.5,1)

    # add the block median residuals to an interpolate PCV file...
    # args.grid should come from the antenna data based on the grid spacing of the antex file
    x = np.linspace(0,90, int(90./args.grid)+1 )
    y = np.linspace(0,360, int(360./args.grid)+1 )

    L1 = interpolate.interp2d(x, y, L1_data, kind='linear')
    L2 = interpolate.interp2d(x, y, L2_data, kind='linear')

    x = np.linspace(0,90, int(90./0.5)+1 )
    y = np.linspace(0,360, int(360./0.5)+1 )
    esm = np.zeros((y.size,x.size,2))

    i = 0

    for az in y :
        j = 0
        for el in x :
            if med[i,j] > 0.00001 or med[i,j] < -0.00001 :
                esm[i,j,0] = med[i,j] + L1(el,az)[0]
                esm[i,j,1] = med[i,j] + L2(el,az)[0]
            else:
                esm[i,j,0] = L1(el,az)[0]
                esm[i,j,1] = L2(el,az)[0]
            j += 1
        i += 1


    if args.elevation or args.polar :
        import matplotlib.pyplot as plt
        from matplotlib import cm

        if args.polar :
            az = np.linspace(0, 360, int(360./0.5)+1)
            zz = np.linspace(0, 90, int(90./0.5)+1)

            fig = plt.figure(figsize=(3.62, 2.76))

            ax = fig.add_subplot(111,polar=True)
            ax.set_theta_direction(-1)
            ax.set_theta_offset(np.radians(90.))
            ax.set_ylim([0,1])
            ax.set_rgrids((0.00001, np.radians(20)/np.pi*2, np.radians(40)/np.pi*2,np.radians(60)/np.pi*2,np.radians(80)/np.pi*2),labels=('0', '20', '40', '60', '80'),angle=180)

            ma,mz = np.meshgrid(az,zz,indexing='ij')
            ma = ma.reshape(ma.size,)
            mz = mz.reshape(mz.size,)
           
            polar = ax.scatter(np.radians(ma), np.radians(mz)/np.pi*2., c=esm[:,:,0], s=50, alpha=1., cmap=cm.RdBu,vmin=-15,vmax=15, lw=0)

            cbar = fig.colorbar(polar,shrink=0.75,pad=.10)
            cbar.ax.tick_params(labelsize=8)
            cbar.set_label('ESM (mm)',size=8)

            for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
                    ax.get_xticklabels() + ax.get_yticklabels()):
                item.set_fontsize(8)

            plt.tight_layout()

        if args.elevation :
            
            # Do an elevation only plot
            fig = plt.figure(figsize=(3.62, 2.76))
            ax = fig.add_subplot(111)
            ele = np.linspace(0,90, int(90./0.5)+1 )
            for i in range(0,721):
                ax.scatter(90.-ele,esm[i,:,0],s=1,alpha=0.8)

            #ax.plot(ele,esm[:,:,0])
            ax.set_xlabel('Elevation Angle (degrees)',fontsize=8)
            ax.set_ylabel('ESM (mm)',fontsize=8)
            ax.set_xlim([10, 90])

            for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
                          ax.get_xticklabels() + ax.get_yticklabels()):
                item.set_fontsize(8)

            plt.tight_layout()

        plt.show()
