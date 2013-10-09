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

# start time
YYYY = 2012
MM = 01
DD = 01

NDAYS = 1

startymdhms = dt.datetime(YYYY,MM,DD)

nav = rnxN.parseFile('t/brdc0010.12n')
obs = rnxO.parseRinexObsFile('t/yar30010.12o')
#obs = rnxO.parseRinexObsFile('t/yar20010.12o')

# Create an elevation plot for S1
fig = plt.figure()#figsize=(3.62, 2.76))
ax = fig.add_subplot(111)

# create a polar plot of S1
#fig2 = plt.figure(figsize=(3.62, 2.76))
#ax2 = fig2.add_subplot(111,polar=True)
#ax2.set_theta_direction(-1)
#ax2.set_theta_offset(np.radians(90.))
#ax2.set_ylim([0,1])
#ax2.set_rgrids((0.00001, np.radians(20)/np.pi*2, np.radians(40)/np.pi*2,np.radians(60)/np.pi*2,np.radians(80)/np.pi*2),labels=('0', '20', '40', '60', '80'),angle=180)

# get the S1 observations from the RINEX file
s1col = rnxO.getObs('S1',obs)

for epoch in obs['epochs'] :
    #print(epoch)
    #epoch['satpos'] = []
    #epoch['Az'] = []
    #epoch['El'] = []
    #epoch['Dist'] = []
    #epoch['dx'] = []
    #ctr = 0
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
        #epoch['satpos'].append(satpos)

        # calculate the vector dx from the station to the satellite
        dx = satpos - sit1['XYZ']
        #epoch['dx'].append(dx)

        # calculate azimuth, elevation to satellite
        # used to compute inter-site satellite visibility only
        Az,El,Dist = bcast.topocent(sit1['XYZ'],dx)
        #print("Az",Az,"El",El,"Dist",Dist)
        #epoch['Az'].append(Az)
        #epoch['El'].append(El)
        #epoch['Dist'].append(Dist)

        # get the observations from this satellite at this epoch 
        satData = epoch['satData'][i]
        print(El, dt, Az, Dist, dx, satData[s1col])

        if satData[s1col] > 0 :
            ax.scatter(El,satData[s1col])
            #polar = ax2.scatter(np.radians(Az), np.radians(90.-El)/np.pi*2., c=satData[s1col] ,s=1,alpha=1., cmap=cm.RdBu,vmin=30,vmax=60, lw=0)
            #print("S1:",satData[s1col])
        #ctr += 1


ax.set_xlabel('Elevation Angle (degrees)',fontsize=8)
ax.set_ylabel('S1',fontsize=8)
ax.set_xlim([0.,90])
#for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
#             ax.get_xticklabels() + ax.get_yticklabels()):
#    item.set_fontsize(8)

#plt.tight_layout()

#cbar = fig2.colorbar(polar,shrink=0.75,pad=.10)
#cbar.tick_params(labelsize=8)
#cbar.ax2.tick_params(labelsize=8)
#cbar.set_label('S1',size=8)

#plt.show()
