#!/usr/bin/env python
from __future__ import division, print_function, absolute_import

import numpy as np
import datetime as dt
import re

from matplotlib import pyplot as plt
from matplotlib import cm

import geodetic as geod
import gpsTime as gpst

import Navigation as rnxN
import Observation as rnxO

# TODO move the topocent into geodetic
import broadcastNavigation as bcast

import argparse

#================================================================================
parser = argparse.ArgumentParser(prog='plotRinex',description='plot RINEX file')

parser.add_argument('-f', '--file', dest='rnxObsFile', default='./t/yar20010.12o')
parser.add_argument('-n', '--nav', dest='rnxNavFile', default='./t/brdc0010.12n')
parser.add_argument('-g', '--gnav', dest='rnxGlonassNavFile', default='./t/alic0010.13n')

#parser.add_argument('-g', '--grid', dest='grid', default=5., type=float)
#parser.add_argument('--polar',dest='polar', default=False, action='store_true')

args = parser.parse_args()
#================================================================================

nav = rnxN.parseFile(args.rnxNavFile)
obs = rnxO.parseRinexObsFile(args.rnxObsFile)

# TODO: Need to calculate my own position
name   = 'YAR2'
lat    = -29. #0465520472
lon    = 115. #3469787567
h      = 0.
eleAng = 0.

sit1 = { 'name'           : name ,
         'latitude'       : lat,
         'longitude'      : lon,
         'height'         : h,
         'ElCutOff'       : eleAng
         }

[a,b,e2] = geod.refell('GRS80')
sit1['XYZ'] = np.matrix(geod.ell2xyz(sit1['latitude'],sit1['longitude'],0.,a,e2))

# TODO: should be obtained from RINEX obs file
# start time
YYYY = 2012
MM = 01
DD = 01

startymdhms = dt.datetime(YYYY,MM,DD)

# Create an elevation plot for S1
fig = plt.figure()#figsize=(3.62, 2.76))
ax = fig.add_subplot(111)

# create a polar plot of S1
fig2 = plt.figure(figsize=(3.62, 2.76))
ax2 = fig2.add_subplot(111,polar=True)
ax2.set_theta_direction(-1)
ax2.set_theta_offset(np.radians(90.))
ax2.set_ylim([0,1])
ax2.set_rgrids((0.00001, np.radians(20)/np.pi*2, np.radians(40)/np.pi*2,np.radians(60)/np.pi*2,np.radians(80)/np.pi*2),labels=('0', '20', '40', '60', '80'),angle=180)

# get the S1 observations from the RINEX file
s1col = rnxO.getObs('S1',obs)

EL_list = []
AZ_list = []

for epoch in obs['epochs'] :
    dt = epoch['time']
    for i in range(0,epoch['numSats']) :
        sat = epoch['sats'][i]

        # Skip glonass satellites for now
        if re.search('R',sat) :
            break 

        #TODO need to handle GLONASS , QZSS, etc 
        sat = int(sat.replace("G",""))

        #satpos = rnxN.satpos(sat,epoch['time'],nav)
        satpos = rnxN.satpos(sat,dt,nav)

        # calculate the vector dx from the station to the satellite
        dx = satpos - sit1['XYZ']

        # calculate azimuth, elevation to satellite
        # used to compute inter-site satellite visibility only
        Az,El,Dist = bcast.topocent(sit1['XYZ'],dx)

        # get the observations from this satellite at this epoch 
        satData = epoch['satData'][i]

        EL_list.append(El)
        AZ_list.append(Az)
        #print("{:.2f} {:.2f}".format(El,Az))
        if satData[s1col] > 0 :
            ax.scatter(El,satData[s1col])
            polar = ax2.scatter(np.radians(Az), np.radians(90.-El)/np.pi*2., c=satData[s1col] ,s=50,alpha=1., cmap=cm.RdBu,vmin=30,vmax=60, lw=0)
            #print("S1:",satData[s1col])



ax.set_xlabel('Elevation Angle (degrees)',fontsize=8)
ax.set_ylabel('S1',fontsize=8)
ax.set_xlim([0.,90])
#for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
#             ax.get_xticklabels() + ax.get_yticklabels()):
#    item.set_fontsize(8)

#plt.tight_layout()

cbar = fig2.colorbar(polar,shrink=0.75,pad=.10)
#cbar.tick_params(labelsize=8)
#cbar.ax2.tick_params(labelsize=8)
cbar.set_label('S1',size=8)

fig3 = plt.figure(figsize=(3.62, 2.76))
ax3 = fig3.add_subplot(111)
plt.hist(EL_list,bins=180)

# work out the number of observations in each az/el bin
#data = np.zeros(np.size(EL_list),2)
#EL_list = np.array(EL_list)
#AZ_list = np.array(AZ_list)
#data[:,0] = EL_list
#data[:,1] = AZ_list

# create a polar plot of S1
#fig4 = plt.figure(figsize=(3.62, 2.76))
#ax4 = fig4.add_subplot(111,polar=True)
#ax4.set_theta_direction(-1)
#ax4.set_theta_offset(np.radians(90.))
#ax4.set_ylim([0,1])
#ax4.set_rgrids((0.00001, np.radians(20)/np.pi*2, np.radians(40)/np.pi*2,np.radians(60)/np.pi*2,np.radians(80)/np.pi*2),labels=('0', '20', '40', '60', '80'),angle=180)
## Calculate the block median
#zz = np.linspace(0,90,int(90/0.5)+1)
#az = np.linspace(0,360,int(360/0.5)+1)
#
#iCtr = 0
#numObs = np.zeros((az.size,zz.size))
#
#for i in az:
#    jCtr = 0
#    #criterion = (data[:,0] < i + 0.5/2.) & (data[:,0] > i - 0.5/2. )
#    criterion = (AZ_list[:] < i + 0.5/2.) & (AZ_list[:] > i - 0.5/2. )
#    ind = np.array(np.where(criterion))
#    #azimuth = i * 0.5
#
#    if ind.size == 0:
#        for j in zz :
#            polar = ax4.scatter(np.radians(i), np.radians(90.-j)/np.pi*2., c=0 ,s=1,alpha=1., cmap=cm.RdBu,vmin=0,vmax=60, lw=0)
#            numObs[iCtr,jCtr] = 0
#            jCtr += 1
#    else:
#        #tmp = data[ind,:]
#        tmp = EL_list[ind]
#        #tmp = tmp.reshape(tmp.shape[1],tmp.shape[2])
#
#    jCtr = 0
#    for j in zz:
#        #criterion = (tmp[:,1] < j + 0.25) & (tmp[:,1] > j - 0.25 )
#        criterion = (tmp[:] < j + 0.25) & (tmp[:] > j - 0.25 )
#        indZ = np.array(np.where(criterion))
#        #ele = j * 0.5
#
#        if indZ.size > 0 :
#            numObs[iCtr,jCtr] = np.size(tmp[indZ,2])
#            polar = ax4.scatter(np.radians(i), np.radians(90.-j)/np.pi*2., c=np.size(tmp[indZ]) ,s=1,alpha=1., cmap=cm.RdBu,vmin=0,vmax=60, lw=0)
#        else :
#            numObs[iCtr,jCtr] = 0
#            polar = ax4.scatter(np.radians(i), np.radians(90.-j)/np.pi*2., c=0 ,s=1,alpha=1., cmap=cm.RdBu,vmin=0,vmax=60, lw=0)
#
#        jCtr += 1
#
#    jCtr = 0
#    iCtr += 1
#
#cbar = fig4.colorbar(polar,shrink=0.75,pad=.10)
#cbar.tick_params(labelsize=8)
#cbar.ax2.tick_params(labelsize=8)
#cbar.set_label('Number of Obs',size=8)
plt.show()
