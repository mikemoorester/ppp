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
    parser.add_argument('-p', '--plot',dest='plot', default=False, action='store_true')
    parser.add_argument('--polar',dest='polar', default=False, action='store_true')
    parser.add_argument('--elevation',dest='elevation', default=False, action='store_true')

    args = parser.parse_args()

    # get the antenna information from an antex file
    antennas = ant.parseANTEX(args.file1)

    # keep the antenna type required
    antenna = ant.antennaType(args.AntType,antennas)
    L1_data = np.array(antenna['data'][0])
    L2_data = np.array(antenna['data'][1])

    # compute the block median from a set of raw residuals
    med,medStd,medrms = res.blockMedian(args.resfile,0.5,1)

    # add the block median residuals to an interpolate PCV file...
    x = np.linspace(0,90, int(90./args.grid)+1 )
    y = np.linspace(0,360, int(360./args.grid)+1 )
    print('About to interpolate L1')
    print("x",x.size,x.shape)
    print("y",y.size,y.shape)
    print("L1_data",L1_data.size,L1_data.shape)
    L1 = interpolate.interp2d(y, x, L1_data, kind='linear')
    print('About to interpolate L2')
    L2 = interpolate.interp2d(x, y, L2_data, kind='linear')

    esm = np.zeros((y.size,x.size))

    i = 0
    j = 0

    for az in y :
        j = 0
        for el in x :
            esm[i][j][0] = med[i][j] + L1(az,el)
            esm[i][j][1] = med[i][j] + L2(az,el)  
            j += 1
        i += 1


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

        plt.show()
