#!/usr/bin/env python
from __future__ import division, print_function, absolute_import

import numpy as np
from scipy import interpolate

import antenna as ant
import residuals as res


#=====================================
if __name__ == "__main__":

    import argparse

    parser = argparse.ArgumentParser(prog='consolidateResiduals', description='PARSE DPH files into one summary file for each station')

    #parser.add_argument('-f', '--file', dest='file1', default='./t/antmod.dat')
    #parser.add_argument('-g', '--grid', dest='grid', default=5., type=float)
    #parser.add_argument('-r', '--resfile', dest='resfile', default='./t/MOBS_DPH_2012_001_143.all')
    #parser.add_argument('-t', '--AntType',dest='AntType', default='ASH701945C_M    NONE')

    #parser.add_argument('--polar',dest='polar', default=False, action='store_true')
    #parser.add_argument('--elevation',dest='elevation', default=False, action='store_true')
    
    # Option to create a summary file of all the DPH residuals 
    parser.add_argument('--dph',dest='dphFile')#, description='Add a DPH file to a Summary file for a station')

    parser.add_argument('--year',dest='year')#, description='Year the Phase residuals were observed - will be added to the epoch column to get a correct time stamp')
    parser.add_argument('--doy',dest='doy')#,description='Day-of-year phase residuals were observed - will be added to the epoch column to obtain a correct time stamp')
    args = parser.parse_args()

    # get the antenna information from an antex file
    dphs = res.parseDPH(args.dphFile)
    print('Sats Viewed:',dphs['satsViewed'])
    print('Epochs Recorded:',dphs['epochs'])

    # Iterate over each epoch
    for epoch in dphs['epochs']:
        print('Epoch:',epoch)
        key = str(epoch)
        print(dphs[key])
        for sat in dphs[key]:
            satPRN = 'prn_'+str(sat)
            print(satPRN)
            i = int(epoch)
            print(dphs[satPRN][i]['lccyc'])
            print(dphs['satPRN'][i])
    #print('KEYS:',dphs.keys())
    
    #   for dph in dphs:
#       polar = ax.scatter(np.radians(ma), np.radians(mz)/np.pi*2., c=esm[:,:,0], s=50, alpha=1., cmap=cm.RdBu,vmin=-15,vmax=15, lw=0)
        

#   if args.elevation or args.polar :
#       import matplotlib.pyplot as plt
#       from matplotlib import cm

#       if args.polar :
#           az = np.linspace(0, 360, int(360./0.5)+1)
#           zz = np.linspace(0, 90, int(90./0.5)+1)

#           fig = plt.figure(figsize=(3.62, 2.76))

#           ax = fig.add_subplot(111,polar=True)
#           ax.set_theta_direction(-1)
#           ax.set_theta_offset(np.radians(90.))
#           ax.set_ylim([0,1])
#           ax.set_rgrids((0.00001, np.radians(20)/np.pi*2, np.radians(40)/np.pi*2,np.radians(60)/np.pi*2,np.radians(80)/np.pi*2),labels=('0', '20', '40', '60', '80'),angle=180)

#           ma,mz = np.meshgrid(az,zz,indexing='ij')
#           ma = ma.reshape(ma.size,)
#           mz = mz.reshape(mz.size,)
           
#           polar = ax.scatter(np.radians(ma), np.radians(mz)/np.pi*2., c=esm[:,:,0], s=50, alpha=1., cmap=cm.RdBu,vmin=-15,vmax=15, lw=0)

#           cbar = fig.colorbar(polar,shrink=0.75,pad=.10)
#           cbar.ax.tick_params(labelsize=8)
#           cbar.set_label('ESM (mm)',size=8)

#           for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
#                   ax.get_xticklabels() + ax.get_yticklabels()):
#               item.set_fontsize(8)

#           plt.tight_layout()

#       if args.elevation :
            
#           # Do an elevation only plot
#           fig = plt.figure(figsize=(3.62, 2.76))
#           ax = fig.add_subplot(111)
#           ele = np.linspace(0,90, int(90./0.5)+1 )
#           for i in range(0,721):
#               ax.scatter(90.-ele,esm[i,:,0],s=1,alpha=0.8)

#           #ax.plot(ele,esm[:,:,0])
#           ax.set_xlabel('Elevation Angle (degrees)',fontsize=8)
#           ax.set_ylabel('ESM (mm)',fontsize=8)
#           ax.set_xlim([10, 90])

#           for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
#                         ax.get_xticklabels() + ax.get_yticklabels()):
#               item.set_fontsize(8)

#           plt.tight_layout()

#       plt.show()
