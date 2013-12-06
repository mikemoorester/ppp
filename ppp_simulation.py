from __future__ import division, print_function, absolute_import
 
import math as math
import numpy as np

import scipy.io as sio
from scipy import sparse
from scipy import interpolate

import array
import datetime as dt

import matplotlib.pyplot as plt

import gpsTime as gpst
import multipath as mp
import broadcastNavigation as bcast
import geodetic as geod
import vmf1 

import rinex.Navigation as rnxN 
import antenna.residuals as res

#=================================
from optparse import OptionParser

parser = OptionParser()
parser.add_option( "-r", "--randomNoise", dest="randomNoise", 
                    help="Random Noise added to the simulated phase range in mm",
                    type="float", default=1. )

parser.add_option( "-c", "--clockFrequency", dest="clockFrequency", 
                    help="Frequency Clock Term to be estimated in seconds (default every epoch)",
                    type="float", default=120. )

parser.add_option( "-e", "--elevationAngle", dest="elevationAngle", help="Elevation Angle",
                    type="float", default=10.0 )

parser.add_option( "-f", "--filename", dest="filename", help="Output File name for RMS stat" )

parser.add_option( "-p", "--plot", dest="plot", help="Plot a figure",
                    action="store_true", default=False )

parser.add_option( "-n", "--ndays", dest="ndays", help="Number of Days to run Simulation",
                    type="int", default=1 )

parser.add_option( "-s", "--samplingPeriod", dest="sampling", 
                    help="Epoch",
                    type="float", default=120 )

parser.add_option( "-t", "--troposphereFrequency", dest="troposphereFrequency", 
                    help="How often to estimate troposphere frequency in seconds",
                    type="float", default=3600 )

parser.add_option( "-u", "--uniformWeight", dest="uniformWeight", help="Turn Uniform Weighting On",
                    action="store_true", default=False )

parser.add_option( "--gw", "--gamitWeight", dest="gamitWeight", help="Turn GAMIT Weighting Scheme On",
                    action="store_true", default=False )

parser.add_option( "--lat", "--latitude", dest="latitude", help="Latitude of Station",
                    type="float", default=-45.0 )

parser.add_option( "--long", "--longitude", dest="longitude", help="Longitude of Station",
                    type="float",default=140.)

#parser.add_option( "--start", dest="start", help="Start time in YYYMMDD format" )

parser.add_option( "--nt", "--noTropEst", dest="tropEst", help="Turn Off Troposphere Estimation",
                    action="store_false", default=True )

parser.add_option( "--nc", "--noClock", dest="clockEst", help="Turn Off Receiver Clock Estimation",
                    action="store_false", default=True )

parser.add_option( "--nm", "--noMultipath", dest="noMP", help="Don't add any multipath bias",
                    action="store_true", default=False )

parser.add_option( "--YYYY", dest="YYYY", default=2012, help="Start Year for simulation (2012)")
parser.add_option( "--MM", dest="MM", default=01, help="Start Month for simulation (01) range is from 01 to 12")
parser.add_option( "--DD", dest="DD", default=01, help="Start Day for simulation (01), range if from 01 to 31")
parser.add_option( "--DDD", dest="DDD", default=0, help="Start Day for simulation rom the DoY (001), range if from 01 to 366")
#=======================
# Multipath parameters

parser.add_option( "-B", "--bin", dest="BIN", help="grid bin to use for averaging", type="float",default=0.5 )

parser.add_option( "-H", "--height", dest="height", help="Monument Height", type="float",default=1.5 )
parser.add_option( "-S", "--smooth", dest="smooth", help="surface roughness", type="float",default=0.3 )
parser.add_option( "-R", "--refractive", dest="refract", help="Refractive Indices", type="float",default=2 )

parser.add_option( "--MH", dest="herring", action="store_true", default=False, help="Use the herring multipath model")

parser.add_option( "-G", "--gainRate", dest="gainrate", help="gainrate", type="float",default=1.1 )

parser.add_option( "--MRES", dest="residuals", action="store_true", default=False, help="Use the observed MP residuals and calculate a block median")

parser.add_option("--res", dest="residualFile", help="Path to the observed residual file")

parser.add_option("--dph", dest="dphFile", help="Path to the observed DPH residual file")

#
# Add an ESM in
#
parser.add_option("--esm","--ESM", dest="esmFile", help="Create an ESM from a residual file .res or .all to apply" )

#
# Add a radome error in
#
parser.add_option( "--radome", dest="radomeFile", help="Radome Error file to add as a bias in m")
parser.add_option( "--radomeAZ", dest="radomeFileAZ", help="Radome Error filei with azimuth variations to add as a bias in m")
parser.add_option( "--radomeName", dest="radomeName", help="Radome Error file to add as a bias in m")

#
# Print options
#
parser.add_option("--PD", dest='print_daily', default=False, action="store_true", help="Print Daily stats") 
#
# print out the residuals to a file
#
parser.add_option("--PR", dest='print_residuals', default=False, action="store_true", help="Print the Residuals to a file")

(option,args) = parser.parse_args()

#=================================
vel_c = 299792458. # vacuum speed of light, m/s

if option.filename :
    filename = option.filename
else :
    filename = 'ppp_simulation.txt'

#
# Some hardcoded defaults
# to be changed to options
#
name   = 'BASE'
lat    = option.latitude
lon    = option.longitude
h      = 0.0
eleAng = option.elevationAngle

sit1 = { 'name'           : name , 
         'latitude'       : lat, 
         'longitude'      : lon,
         'height'         : h,
         'ElCutOff'       : eleAng
         }

         #'elevationAngle' : eleAng
name   = 'ROVE'

sit2 = { 'name'           : name , 
         'latitude'       : lat, 
         'longitude'      : lon,
         'height'         : h,
         'ElCutOff'       : eleAng
         }

[a,b,e2] = geod.refell('GRS80')
sit1['XYZ'] = np.matrix(geod.ell2xyz(sit1['latitude'],sit1['longitude'],0.,a,e2))
sit2['XYZ'] = np.matrix(geod.ell2xyz(sit2['latitude'],sit2['longitude'],0.,a,e2))

# set the navigation file name
#navfile = { 'name' : 'brdc0010.10n',
#            'binName' : 'brdc0010.10n.nav',
#            'data' : []}
if option.uniformWeight :
    ELE_WEIGHT = 0
    GAMIT_WEIGHT = 0
elif option.gamitWeight :
    ELE_WEIGHT = 0
    GAMIT_WEIGHT = 1
else :
    ELE_WEIGHT = 1
    GAMIT_WEIGHT = 0

GRID = option.BIN

# Check to see if we are applying an ESM to the simulation
if option.esmFile:
    ESM,esmStd,esmrms = res.blockMedian(option.esmFile,GRID,1)
    #ESM = res.blockMedian(option.esmFile,GRID,1)
    NaNs = np.isnan(ESM)
    ESM[NaNs] = 0.
    # convert from mm to m
    ESM = ESM/1000.
    # order the matrix by elevation, rather than zenith angle 
    ESM = ESM[:,::-1]
    esmFlag = 1
    x = np.linspace(0,90, int(90./GRID)+1 )
    y = np.linspace(0,360, int(360./GRID)+1 )
    e = interpolate.interp2d(x, y, ESM, kind='linear')
else:
    esmFlag = 0

# Check to see if we need to import a DPH file
if option.dphFile:
    print("About to parse a DPH file:",option.dphFile)
    dphs = res.parseDPH(option.dphFile)

# Use unhealthy satellites => 1, don't => 0
UNHEALTHY = 0

# number of days to process
NDAYS = option.ndays

# sample rate of simulated observations
# Needs to be divisible by 86400
SAMP = option.sampling #120.

# Solution session length, typically 24 hours
sess=86400

MAXSAT = 32 

# predict the A amd b matrix sizes
MAXOBS = int( sess/SAMP * MAXSAT * NDAYS) #30000

# Set the reolution of the grids
BIN = option.BIN 

# Set start time and end time (based on NDAYS)
YYYY = int(option.YYYY)

if option.DDD:
    MM   = 1 
    DD   = 1
    startDT = dt.datetime(YYYY,MM,DD)
    print("Start Date:",startDT)
    print("Option DDD:",option.DDD)
    startDT = startDT + dt.timedelta(days=(int(option.DDD) -1)) 
    print("After timeDelta startDT:",startDT)
    MM = startDT.strftime("%m")
    DD = startDT.strftime("%d")
else:
    MM   = option.MM
    DD   = option.DD
    startDT = dt.datetime(YYYY,MM,DD)


startymdhms = gpst.cal2jd(YYYY,int(MM),int(DD)+(00/24)+(00/(24*60))+(00.0000/(24*3600))) #user setting
stopymdhms = startymdhms + NDAYS #gpst.cal2jd(YYYY,MM,DD+(((24*NDAYS)-1)/24)+(59/(24*60))+(59.999/(24*3600))) #user setting

amb_count = 0
resFlag = 0

v_El_total = np.array([])
v_full_total = np.array([])
# day iteration counter
ditr = 0

# loop until we've finished all the sessions
# this will ensure we can have multiple days

while startymdhms < stopymdhms :
    # Parameter values
    #
    # users should only need to modify delta and nextest
    #
    # delta - defines how many seconds between estimates, must divede by 86400
    # nest  - defines the numbe rof parameters of this type in this run
    # nextest - defines the next estimation time (JD), setting this value to something
    #           very large or very small ensures the estimate will never happen
    # count - number of paramters of this type estimated to date
    #     

    # Solutions are time-stamped at mid-point of each solution period
    time_deltasec = SAMP
    time_start_jd = startymdhms
    time_stop_jd  = startymdhms + (sess - 1e-3)/86400
    time_deltaday = time_deltasec/86400
    time_epochs   = np.linspace(time_start_jd,time_deltaday,time_stop_jd)

    # coordinate estimates
    north_delta   = 86400 
    north_nest    = round((time_stop_jd - time_start_jd)*(86400/north_delta))
    north_nextest = time_start_jd
    #print("NORTH_NEST :", north_nextest)
    north_count   = 0
    north_time    = []
    A_north = sparse.lil_matrix( (MAXOBS,north_nest) )

    east_delta   = 86400 
    east_nest    = round((time_stop_jd - time_start_jd)*(86400/east_delta))
    east_nextest = time_start_jd
    east_count   = 0
    east_time    = []
    A_east = sparse.lil_matrix( (MAXOBS, east_nest) )

    up_delta   = 86400 
    up_nest    = round((time_stop_jd - time_start_jd)*(86400/up_delta))
    up_nextest = time_start_jd
    up_count   = 0
    up_time    = []
    A_up = sparse.lil_matrix( (MAXOBS, up_nest) )

    # troposphere estimates
    trop_delta   = option.troposphereFrequency 
    trop_nest    = round((time_stop_jd - time_start_jd)*(86400/trop_delta))
    trop_count   = 0 #inilise the column counter for the parameter
    trop_time    = []
    A_trop       = sparse.lil_matrix( (MAXOBS, trop_nest) )
    if option.tropEst:
        trop_nextest = time_start_jd
    else:
        trop_nextest=1e20  # uncomment to not estimate trop

    # ambiguity
    amb_delta    = sess  #user setting
    amb_nest     = round((time_stop_jd-time_start_jd)*(86400/amb_delta))
    amb_nextest  = time_start_jd+amb_delta/86400
    amb_time     = []
    amb_El       = []
    amb_Az       = []
    amb_count    = 0 #initialise the column counter for the parameter
    #amb.count is not used since we need a new ambiguity per sat
    A_amb        = sparse.lil_matrix( (MAXOBS, MAXSAT * int(sess/43200)) )

    # solve for harmonics at one or more frequencies
    harm_delta    = sess    #wouldn't want to change this
    #harm.name=strvcat('m2','s2','n2','t2','k2','r2','k1','o1','p1','q1','s1');
    #harm.w=[28.984104252; 30.0000000; 28.439729568; 29.958933; 28.439729568; 30.041066664; ...
    #       15.041068632; 13.943035584; 14.958931356; 13.398660936; 15.000000]*24;%degrees per hour*24
    #harm.name=strvcat('s2','s1');
    harm_w        = [28.439729568, 30.0000000, 15.041068632, 13.398660936] #degrees per hour*24
    harm_w        = np.transpose(harm_w) * 24.
    #number of constituents to estimate for the above components (actual constituents defined below)
    harm_num      = np.size(harm_w)
    harm_param    = harm_num * 2. #sine and cosine term for each harmonic estimated
    harm_nest     = round((time_stop_jd-time_start_jd)*(86400./harm_delta))*harm_param
    #harm.nextest =time.start.jd  
    harm_nextest  = 10**20   # don't estimate!
    harm_count    = 0.       # initialise the column counter for the parameter
    A_harm        = sparse.lil_matrix( (MAXOBS, harm_nest) )

    sync_delta    = time_deltasec # one clock sync per epoch
    sync_nest     = round((time_stop_jd - time_start_jd)*(86400/sync_delta))
    sync_count    = 0 #initialise the column counter for the parameter
    sync_time     = []
    A_sync        = sparse.lil_matrix( (MAXOBS, sync_nest) )
    if option.clockEst:
        sync_nextest = time_start_jd 
        print('Estimating the clock')
    else:
        sync_nextest = 10**20   #don't estimate!
        print('Not Estimating the clock')

    #=============================================================
    # Define the bias in each coordinate component
    #
    bias_t   = np.linspace(time_start_jd,time_stop_jd,num=86400) 

    #harmonic bias in up component
    biasu_amplitude = 0.000;#metres for all constituents - make a vector for different amps per constit
    biasu_phase     = 0;#degrees for all constituents
    #bias_harmw=[28.984104252; 30.0000000; 28.439729568; 29.958933; 30.082137264; 30.041066664; ...
    #       15.041068632; 13.943035584; 14.958931356; 13.398660936; 15.000000;15.04246526627447]*24;%degrees per hour*24
    #bias_harmw      = [15.04246526627447]*24; #S1
    bias_harmw      =[0]*24 #MJM does not appear to adervsely affect estimation of MP so leave the value in
    # harm.w constituent order is:
    # 'm2','s2','n2','t2','k2','r2','k1','o1','p1','q1','s1','k1-8s(multipath)'
    biasu_angvel    = np.radians(bias_harmw) #angular velocity
    biasu_values    = 0 * bias_t #sum(biasu_amplitude*np.sin(biasu_angvel*np.transpose(bias_t)+np.radians(biasu_phase)),1);

    # harmonic bias in east component
    biase_amplitude = 0.00 #metres
    biase_phase     = 0.0  #degrees
    biase_angvel    = np.radians(bias_harmw[0]) #angular velocity -S1
    biase_values    = biase_amplitude * np.sin( np.radians(biase_angvel*bias_t+biase_phase));

    # harmonic bias in north component
    biasn_amplitude = 0.00 # metres
    biasn_phase     = 0.0  # degrees
    biasn_angvel    = np.radians(bias_harmw[0]) #angular velocity - S1
    biasn_t         = np.linspace(time_start_jd, time_deltaday, time_stop_jd)
    biasn_values    = biasn_amplitude*np.sin(np.radians(biasn_angvel*bias_t+biasn_phase))
    
    #=============================================================
    # Define the multipath bias
    #
    SMOOTH = option.smooth 
    HEIGHT = option.height                  #0.20
    ANTENNATYPE = 'LEIAT504        NONE'
    REFRACTIVE = option.refract 
    GAINRATE = option.gainrate

    if option.herring:
        xi = np.linspace(0,90,(90./GRID)+1)
        tmp = mp.multipath_herring(xi, SMOOTH, HEIGHT, REFRACTIVE,GAINRATE)
        mp_SIGMALC = np.zeros( ((int(360./GRID)+1),int(90./GRID)+1) )
        for i in range(0,360,int((360./GRID))+1):
            mp_SIGMALC[i,:] = tmp[:]
    # 
    # Add an observed residual file as the multipath bias
    elif option.residuals:
        if resFlag == 0:
            mp_SIGMALC = res.blockMedian(option.residualFile,GRID,1)
            NaNs = np.isnan(mp_SIGMALC)
            mp_SIGMALC[NaNs] = 0.
            # convert from mm to m
            mp_SIGMALC = mp_SIGMALC/1000.
            # order the matrix by elevation, rather than zenith angle 
            mp_SIGMALC = mp_SIGMALC[:,::-1]
            resFlag = 1
    #
    # No multipath Option:
    elif option.noMP:
        mp_SIGMALC = np.zeros( ((int(360./GRID)+1),int(90./GRID)+1) )
    else:
        mp_SIGMALC = mp.multipath_moore(SMOOTH,HEIGHT,ANTENNATYPE,REFRACTIVE,GRID,0)

    x = np.linspace(0,90, int(90./GRID)+1 )
    y = np.linspace(0,360, int(360./GRID)+1 )

    f = interpolate.interp2d(x, y, mp_SIGMALC, kind='linear')
    #f = interpolate.RectBivariateSpline(x, y, mp_SIGMALC)
    if option.radomeFile:
        radome = np.genfromtxt(option.radomeFile)
        r = interpolate.interp1d(radome[:,0], radome[:,1], kind='linear')
    elif option.radomeFileAZ:
	radome = np.genfromtxt(option.radomeFileAZ)
    	r = interpolate.interp2d(radome[:,0], radome[:,1], radome[:,2], kind='linear')

    # intialise observation counter
    obs_count = 1
    b = sparse.lil_matrix((MAXOBS,1))

    v_Az = [] # sparse.lil_matrix((MAXOBS,1))
    v_El = []
    v_time = sparse.lil_matrix((MAXOBS,1))

    # work out the navfile name we need to import
    # need to see if there is a routine that will read matlab binary format
    # in python
    #navfile.name = navfileDir + '/brdc' + doy + '012n'
    #navfile.biname = navfile.name + '.nav'
    # convert rinex nav to a binary file
    #rinexe()

    (year, doy) = gpst.jd2doy(startymdhms)
    yy = gpst.yyyy2yy(year)
    #navfile = 'brdc'+str(doy)+'0.'+str(yy)+'n.nav'
    #navfile = 'brdc'+str(doy)+'0.'+str(yy)+'n'
    navfile = 'brdc'+str(doy)+'0.'+ ("%02d" % yy)  +'n'
    #ndata = nav.parseNavData('/Users/michael/code/matlab/king/mikemoore/brdc/'+navfile) #brdc0010.12n.nav')
    #ndata = nav.parseFile('/Users/michael/code/matlab/king/mikemoore/brdc/'+navfile) #brdc0010.12n.nav')
    #ndata = nav.parseFile('/Users/moore/code/matlab/king/mikemoore/brdc/'+navfile) #brdc0010.12n.nav')
    #ndata = bcast.get_binary_eph('/home/547/mjm547/code/ppp/test_data/'+navfile)#brdc0010.12n.nav')
    #ndata = bcast.get_binary_eph('/Users/moore/code/matlab/king/mikemoore/brdc/brdc0010.12n.nav')
    #ndata = bcast.get_binary_eph('/Users/michael/code/matlab/king/mikemoore/brdc/'+navfile)#brdc0010.12n.nav')
    #nav = rnxN.parseFile('/Users/moore/code/matlab/king/mikemoore/brdc/'+navfile)
    #nav = rnxN.parseFile('/Users/moore/code/matlab/king/mikemoore/brdc/'+navfile)
    nav = rnxN.parseFile('/short/dk5/brdc/'+str(YYYY)+'/'+navfile)
    #print('/short/dk5/brdc/'+str(YYYY)+'/'+navfile)

    ep = np.linspace(time_start_jd,time_stop_jd,num=86400./SAMP)
    sat = np.zeros(MAXSAT)

    satD = {}
    for sv in range (1,MAXSAT+1):
        skey = str(sv)
        satD[skey] = {}

    # This naming is stupid...
    # check if epoch_ctr, should really be called obs_ctr
    epoch_ctr = 0
    # E is a real epoch ctr, unlike the variable epoch_ctr!! 
    E = 0

    for epoch_time in ep :
        #E += 1
        epoch_ctr += 1
        # Do some time conversions
        [gpsweek,gpssow,nul] = gpst.jd2gps(epoch_time)
        gpssow = np.round(gpssow)
        epochs = {} 
        epochs['sats'] = [] 

        #
        # Need to update the counters here
        #
        if (epoch_time >= north_nextest):
            print(epoch_time,"epoch_time/86400",(epoch_time/86400),"north_nextest",north_nextest)
            north_count = north_count + 1
            north_time.append(epoch_time+(north_delta/86400.)/2.) #%mid point of batch
            north_nextest=north_nextest+north_delta/86400
  
        #if ((epoch_time+1e-3/86400.) >= east_nextest): <- this results in too many parameters being estimated
        if ( epoch_time >= east_nextest):
            east_count = east_count + 1
            east_time.append(epoch_time+(east_delta/86400.)/2.)
            east_nextest=east_nextest+east_delta/86400
  
        if ( epoch_time >= up_nextest):
            up_count = up_count + 1
            up_time.append(epoch_time+(up_delta/86400.)/2.)
            up_nextest=up_nextest+up_delta/86400
  
        if (epoch_time >= trop_nextest):
            trop_count = trop_count + 1
            trop_time.append(epoch_time+(trop_delta/86400.)/2.)
            trop_nextest=trop_nextest+trop_delta/86400
  
        if ( epoch_time >= harm_nextest):
            harm_count = harm_count + harm_param
            harm_time.append( epoch_time+(harm_delta/86400.)/2.)
            harm_nextest=harm_nextest+harm_delta/86400

        if ( epoch_time >= sync_nextest):
            sync_count = sync_count + 1
            sync_time.append(epoch_time+(sync_delta/86400.)/2.)
            sync_nextest=sync_nextest+sync_delta/86400

        #%% Satellite data
        #sat_dict.append(svn_dict)
        # find out what satellite elevation angles are for the site
        for sv in range(1,MAXSAT+1):  #some of this loop is done more times than needed
            skey = str(sv)
            #satD[skey]['col'] = bcast.find_eph(ndata,sv,gpssow) #FIND_EPH  Finds the proper column in ephemeris array  
            satD[skey]['col'] = 1 

            if (satD[skey]['col'] > -1 ):
                # SATPOS Calculation of X,Y,Z coordinates at time t
                #satD[skey]['pos'] = bcast.satpos(satD[skey]['col'],ndata,gpssow) 
                # dt is a datetime object for the timestamp you want to find the satellites postion
                #[gpsweek,gpssow,nul] = gpst.jd2gps(epoch_time)
                satD[skey]['pos'] = rnxN.satpos(sv,startDT,nav)

		if satD[skey]['pos'] is not None: # == -1 :
		#if satD[skey]['pos'].any(): # == -1 :
		#    print("Can' calc sat pos for sv:",sv)
                #else:
                    dx1 = satD[skey]['pos'] - sit1['XYZ']
                    dx2 = satD[skey]['pos'] - sit2['XYZ']

                    # calculate azimuth, elevation to satellite
                    # used to compute inter-site satellite visibility only
                    (satD[skey]['Az1'],satD[skey]['El1'],satD[skey]['Dist1']) = bcast.topocent(sit1['XYZ'],dx1)
                    (satD[skey]['Az2'],satD[skey]['El2'],satD[skey]['Dist2']) = bcast.topocent(sit2['XYZ'],dx2)

                    # satellite pos needs correction due to earth rotation during travel time
                    # use strang and borre routine e_r_corr
                    satD[skey]['pos_rot1'] = bcast.earth_rot_corr(satD[skey]['Dist1']/vel_c, satD[skey]['pos'])
                    satD[skey]['pos_rot2'] = bcast.earth_rot_corr(satD[skey]['Dist2']/vel_c, satD[skey]['pos'])

                    dx1 = satD[skey]['pos_rot1'].T - sit1['XYZ']
                    dx2 = satD[skey]['pos_rot2'].T - sit2['XYZ']

                    (satD[skey]['Az1'],satD[skey]['El1'],satD[skey]['Dist1']) = bcast.topocent(sit1['XYZ'],dx1)
                    (satD[skey]['Az2'],satD[skey]['El2'],satD[skey]['Dist2']) = bcast.topocent(sit2['XYZ'],dx2)

        amb = {}
        amb['time'] = np.zeros(MAXSAT+1)
        amb['El'] = np.zeros(MAXSAT+1)
        amb['Az'] = np.zeros(MAXSAT+1)
        sv_count = 1

        for sv in range(1,MAXSAT+1):
            skey = str(sv)

            if((epoch_time == time_start_jd) ):
                #initialise the last time this sat was seen
                satD[skey]['lastepoch'] = 0
             
            if(satD[skey]['col'] > -1 and 'El1' in satD[skey] and 'El2' in satD[skey]):
#               %check if the Elevation angle is above the cutoff at both sites AND no obstruction 
#               % - if it is, it is an SD observation
                #print("Made it",epoch_time,skey,gpssow)
                if( (satD[skey]['El2'] > sit2['ElCutOff']) &
                     (satD[skey]['El1'] > sit1['ElCutOff']) ):
                    # check to see if the obs is in the residual file
                    # if this is the option.dphFile
                    if option.dphFile:
                        # DPH is on 30s epochs, need to correlated ths with
                        # the sample time being used in the simulation
                        epochSTR = str(epoch_ctr * int(SAMP/30))
                        #epochSTR = str(epoch_ctr+1)
                        if sv in dphs[epochSTR]:
                            epochs['sats'].append(sv)
                            v_El.append(satD[skey]['El1'])
                            v_Az.append(satD[skey]['Az1'])
                    else:
                        epochs['sats'].append(sv)
                        v_El.append(satD[skey]['El1'])
                        v_Az.append(satD[skey]['Az1'])
                        #v_El[obs_count-1,0] = satD[skey]['El1']

#               %% Two ways to get a new ambiguity parameter
#               % Firstly, a new satellite arc (2 epochs without an obs)
#               % Secondly, a new coordinate solution
                if( ((epoch_time-satD[skey]['lastepoch']) > ( 2. * time_deltasec/86400.))
                     | (epoch_time >= amb_nextest) ):
                    #reference value for Ambiguity columns - caters for missing svs
                    satD[skey]['ref'] = sv_count 

                    # Store the ambiguity informaiton
                    amb['time'][sv_count] = epoch_time
                    amb['El'][sv_count]   = satD[skey]['El1']
                    amb['Az'][sv_count]   = satD[skey]['Az1']

                    # Increment counter for next time
                    sv_count=sv_count+1;

            satD['lastepoch']=epoch_time; #update epoch of last (obs>ElCutOff) for this sat

        # New ambiguities have been assigned so we can wait until the next reset now
        if (epoch_time >= amb_nextest):
            amb_nextest = amb_nextest + amb_delta/86400

        # Loop around the number of satellites at this epoch and form the Matrices
        for i in range(0,len(epochs['sats'])):
            skey = str(epochs['sats'][i])

            # Set up the 'A' (design) Matrix
            # North
            #print("Obs count",obs_count, "North_count:", north_count)
            if(north_count > 0):
                A_north[(obs_count-1, north_count-1)] =-( np.sin(np.radians(90.-satD[skey]['El1'])) * 
                                                           np.cos(np.radians(satD[skey]['Az1'])) )
            # East
            if(east_count > 0):
                A_east[(obs_count-1,east_count-1)]=-( np.sin(np.radians(90.-satD[skey]['El1'])) *
                                                     np.sin(np.radians(satD[skey]['Az1'])) )
            # Up
            if (up_count > 0):
                A_up[(obs_count-1,up_count-1)] =-np.cos(np.radians(90.-satD[skey]['El1']))

            # TZD
            if (trop_count > 0):
                #A_trop[(obs_count-1,trop_count-1)]=1./np.cos(np.radians(90.-satD[skey]['El1']))
                dry, wet = vmf1.vmf1map([90.-satD[skey]['El1']])
                A_trop[(obs_count-1, trop_count-1)] = wet[0] 
                #increase in TZD reduces range

            #  Harmonic terms (local up only at the moment)
            #    if (harm.count > 0)
            #      % cosines all together then sines - cos(wt+phi) - relative to J2000
            #      A.harm(obs_count,(harm.count-harm.param+1):(harm.count-harm.num))=[...
            #        cos(harm.w'*d2r*(epoch.time-2451545.0)+0*d2r).*A.up(obs_count,up.count) ];
            #      A.harm(obs_count,(harm.count-harm.num+1):(harm.count))=[...
            #        sin(harm.w'*d2r*(epoch.time-2451545.0)+0*d2r).*A.up(obs_count,up.count) ];
            #    end

            # Clock sync
            if (sync_count > 0) & option.clockEst:
                A_sync[(obs_count-1,sync_count-1)]=1

            # Ambiguity - 1 per satellite per pass
            #        satD[skey]['ref'] = sv_count 
            #print("Amb:",obs_count,satD[skey]['ref'])
            A_amb[(obs_count-1, satD[skey]['ref']-1)]=1;
            amb_count += 1

            # Set up the 'b' (O-C) Matrix - up bias for this epoch
            #epoch_ctr = epoch_ctr+1
            epochu_bias=np.interp(epoch_time,bias_t,biasu_values)#,'linear')
            epoche_bias=np.interp(epoch_time,bias_t,biase_values)
            epochn_bias=np.interp(epoch_time,bias_t,biasn_values)

            # Santerre (1991) introduces a bias into the fixed station
            # b(obs_count,1)=sin(d2r*(90-sat(epochs.sats(i)).El))*cos(d2r*(90-sat(epochs.sats(i)).El))*...
            #  d2r*(sat(epochs.sats(i)).Az-bias.blAz)*bias.bl*epoch.bias/sat(epochs.sats(i)).Dist;
            # We introduce a bias into the single difference at the free station
            # here the bias is inserted into the height 

            # O-C vector, which is where the biases enter the solution
            # put in e, n bias, using azimuth and length of this bias
            #  see King et al, EPS paper for details
            bias_enaz = (np.pi/2. - np.arctan2(epochn_bias,epoche_bias)) * 180. / np.pi # north is azimuth zero
            bias_enlen = np.sqrt(epochn_bias**2 + epoche_bias**2)

            bias_a = satD[skey]['Dist1'] * np.cos( np.radians(satD[skey]['El1']) )
            bias_b = np.sqrt( bias_a**2 + bias_enlen**2 - 2*bias_a*bias_enlen * 
                        np.cos( np.radians(satD[skey]['Az1'] - bias_enaz + 180.)))
            #%180 because we add it on rather than remove it

            bias_en = np.sqrt((satD[skey]['Dist1'])**2+(bias_b-bias_a)**2 -
                     2*(satD[skey]['Dist1']) * (bias_a-bias_b)*np.cos(np.radians(satD[skey]['El1'])))

            # now have range biased by E, N; now add U bias and subtract "computed" term
            # range bias added - Herring modelled MP term
            #  MJM Add in the AZ/EL depend bias from MP at this point
            bias_mp = f(satD[skey]['El1'], satD[skey]['Az1'])[0]

            # Apply the esm model
            if option.esmFile:
                #bias_mp = bias_mp - e(satD[skey]['El1'], satD[skey]['Az1'])
                bias_mp = bias_mp + e(satD[skey]['El1'], satD[skey]['Az1'])[0]

            if option.radomeFile:
                bias_mp = bias_mp + r(satD[skey]['El1'])
	    elif option.radomeFileAZ:
                bias_mp = bias_mp + r(satD[skey]['El1'],satD[skey]['Az1'])[0]

            b[obs_count-1,0] = ( np.sqrt( bias_en**2 + epochu_bias**2 - 
                                2. * bias_en * epochu_bias * np.cos(np.radians(90.+satD[skey]['El1'])) ) +
                                bias_mp
                                -(satD[skey]['Dist1']))

            # increment obs counter (end of epoch loop)
            obs_count = obs_count +1

        startDT = startDT + dt.timedelta(seconds=option.sampling)

    #======
    # End of the day loop
    startymdhms += 1.

    # Now time to do the math...
    print("coord estimates",north_count,east_count,up_count,"troposphere estimates",trop_count,"OBS count",obs_count)

    # From up the A and Ql matrix without the unused row/columns
    print("Checking Coord matrix")

    # Need to do a more sophisticated check on the numbe rof parameters by checking if it is [], then summing
    nparams = 3 + trop_count + sync_count #+ amb['time'].shape[0]
    print("Number of parameters:",nparams)

    A_full   = sparse.lil_matrix( ((obs_count-1),nparams) )

    A_full[:,0] = A_north[0:obs_count-1,:]
    A_full[:,1] = A_east[0:obs_count-1,:]
    A_full[:,2] = A_up[0:obs_count-1,:]

    v_El = np.array(v_El)
    v_Az = np.array(v_Az)

    #v_time = v_time[0:obs_count-1,:]

    if (trop_count==0):
        A_trop = []
    else:
        A_full[:,3:(trop_count+3)] = A_trop[0:obs_count-1,:]

    if (harm_count==0):
        A_harm = []
    else:
        A_harm = A_harm[0:obs_count-1,:]
         
    if (sync_count==0):
        A_sync = []
    else:
        A_sync = A_sync[0:(obs_count-1),:]
        A_full[:,(trop_count+4):nparams] = A_sync[:,0:(sync_count-1)] 
        #A_full[:,(trop_count+4):nparams] = A_sync[0:(obs_count -1),:] 
        #A_sync = A_sync[0:obs_count-1,:]
                  
    #if (amb['time'].shape[0] == 0):
    #    A_amb = []
    #else:
    #    print("About to sort out the ambiguity matrix")
    #    A_amb = A_amb[0:obs_count-1,:]
    #    A_full[:,(trop_count+4):nparams] = A_amb[:,0:(amb['time'].shape[0] -1)]

    del A_north, A_east, A_up, A_trop, A_amb, A_sync

    b = b[0:obs_count-1,:]

    if(ELE_WEIGHT):
        P=sparse.diags((np.sin(np.radians(v_El)) + 0.000001)**2,0) #;%0.00001 avoids problems when elev = 0
    elif(GAMIT_WEIGHT):
        P=sparse.diags(0.010+0.005**2/(np.sin(np.radians(v_El))**2 + 0.001)**2,0)
    else:
        P=sparse.eye(A_full.shape[0])  

    #print("P shape:",P.shape)

    #%% Add random noise to b matrix
    #%fprintf('Adding 1mm random noise to ranges\n');
    if option.residuals:
        b = b.T
    else :
        if option.randomNoise < 0.0001:
            b = b.T
        else:
            b = b.T + np.random.normal(0,(option.randomNoise/1000.),b.shape[0])
    #print("b shape",b.shape)
    #%% Matrix setup now complete

    #%% Invert Matrix
    Qx = np.linalg.pinv(A_full.todense().T * P.todense() * A_full.todense() )
    dx = Qx *A_full.todense().T * P.todense() * b.T
    
    vfull = A_full.todense() * dx - b.T

    if NDAYS > 1 :
        if ditr > 0:
            v_El_total   = np.concatenate((v_El_total, v_El), axis=1)
            v_full_total = np.concatenate((v_full_total, vfull), axis=0)
        else :
            v_El_total   = v_El
            v_full_total = vfull

	if (option.print_daily) :
	    DITR = str(ditr)
	    fid = open(filename,"a")
            fid.write(str(NDAYS))
            fid.write(' ')
            fid.write(str(HEIGHT))
            fid.write(' ')
            fid.write(DITR)
	    fid.write(' ')
            fid.write(str(dx[0,0]))
            fid.write(' ')
            fid.write(str(dx[1,0]))
            fid.write(' ')
            fid.write(str(dx[2,0]))
            fid.write(' ')
            fid.write('\n')
            fid.close()

        ditr += 1 
        #del vfull, v_El
    # End of the day loop

#print("NDAYS:",NDAYS)

if NDAYS > 1 :
    #print("\nvfull - total",np.shape(v_full_total))
    #print("v_El -total\n\n",np.shape(v_El_total))
    vfull = v_full_total
    v_El = v_El_total
    #print("TOTAL: vfull - total",vfull.shape)
    #print("TOTAL: v_El -total\n\n",v_El.shape)


#%% Now do solution fixing the amb.perc percentage of ambiguities
#%amb.perc=[00 20 40 60 80 100]';
#amb.perc=100
#amb.rand=randperm(size(A.amb,2)).T #%generate random order of the matrices so we get random ambiguity elimination


#loop around each amb.perc
#for i =1:size(amb.perc,1):
#    amb.part=A.amb;
#    amb.part(:,amb.rand(1:floor((amb.perc(i)/100)*size(A.amb,2))))=[];
#    Apart=[A.north A.east A.up A.trop A.harm A.sync amb.part];

#    fprintf('\nNow inverting for the final solution and covariance - %i%% ambiguities fixed\n',amb.perc(i));
#   fprintf('The size of the ambiguity part of the A matrix is now %i\n\n',size(amb.part,2));

  #% Do the computations
#  Qxpart=pinv(full(Apart'*P*Apart));
#    dxpart=Qxpart*Apart'*P*b;
#    #%Npart=Apart'*P*Apart;
#    #%dxpart=Npart\(Apart'*P*b);
#   vpart=Apart*dxpart-b;

  #% Keep a record of this solution
#     perc=sprintf('%02d',amb.perc(i,1));
#    eval(['dxpart_',perc,'=dxpart;']);

#end

# Store the median model into an array
model = {}
m = np.matrix(np.zeros( (int(90/BIN+1),3) ))

# To store the differences from the model and simulated bias
d = np.matrix(np.zeros( (int(90/BIN+1),3) ))

i = 0
sval = 0

for ele in np.linspace(0, 90, int(90/BIN +1)) :
    bind = ((v_El <= ele + BIN ) & (v_El > ele - BIN))
    ind = np.array(np.where(bind)) 
    tmp = np.array(vfull[ind,0])

    res = np.median(tmp) # * -1.
    #sim = f(ele,0.) * -1.
    sim = mp_SIGMALC[0,i] *-1

    m[i,0] = ele
    m[i,1] = res 
    m[i,2] = sim

    d[i,0] = ele
    d[i,1] = sim - res
    d[i,2] = (sim -res)**2
    if np.isnan(d[i,2]):
        sval += 0.
    else :
        sval += (sim -res)**2
    val = np.median(tmp)

    del tmp
    i += 1

#
# print the residuals to a file
#
if option.print_residuals:
    if option.clockEst:
        ppp = 1
    else:
        ppp = 0
    if option.uniformWeight:
        eleW = 0
    else:
        eleW = 1

    h = "{:.2f}".format(option.height)

    if option.tropEst:
        tropName = str(int(option.troposphereFrequency / 3600))
    else:
        tropName = '0'

    if option.radomeName:
        resOut = 'ppp_'+str(ppp)+'_eleW_'+str(eleW)+'_eleA_'+str(int(option.elevationAngle))+'_tropEst_'+tropName+'_ndays_'+str(option.ndays)+'_S_'+str(option.smooth)+'_r_'+str(option.randomNoise)+'_radome_'+option.radomeName+'_h_'+h+'.res'
    else: 
        resOut = 'ppp_'+str(ppp)+'_eleW_'+str(eleW)+'_eleA_'+str(int(option.elevationAngle))+'_tropEst_'+tropName+'_ndays_'+str(option.ndays)+'_S_'+str(option.smooth)+'_r_'+str(option.randomNoise)+'_h_'+h+'.res'

    i = 0
    print("v_El:",np.shape(v_El),"v_Az:",np.shape(v_Az),"vfull:",np.shape(vfull))
    fres = open(resOut,'w')
    
    for r in vfull :
	if i < np.size(v_El) and i < np.size(v_Az) and i < np.size(vfull):
            # need to convert elevation to zenith angle, and the residuals from m to mm
            fres.write('{0:3.4f} {1:3.4f} {2:.2f}\n'.format(v_Az[i],(90.-v_El[i]),(vfull[i,0]*1000.)))
        i += 1
    fres.close()
#
# Calculate the RMS
# mask out the NaN from the calculation
#
criterion = np.logical_not( np.isnan(d[:,2]) )
ind = np.where(criterion)
diff = np.array( d[ind] )

stmp = sum(diff[0,:])
divisor = 1./np.size(diff)
mtmp = sval*divisor
rms = np.sqrt(mtmp)
#print("SUM:",stmp,"Divisor:",divisor,"Mtmp:",mtmp,"RMS:",rms)

fid = open(filename,"a")
fid.write(str(NDAYS))
fid.write(' ')
startDT = startDT - dt.timedelta(days=1)
fid.write(startDT.strftime("%Y %j %H:%M:%S"))
fid.write(' ')
fid.write(str(HEIGHT))
fid.write(' ')

if option.residuals:
    fid.write(str(option.residualFile))

fid.write(' ')
fid.write(str(rms*1000.))
fid.write(' ')
fid.write(str(dx[0,0]))
fid.write(' ')
fid.write(str(dx[1,0]))
fid.write(' ')
fid.write(str(dx[2,0]))
fid.write(' ')

if option.tropEst :
    for i in range(0,trop_count): 
        fid.write(str(dx[i+3,0]))
        fid.write(' ')

fid.write(str(BIN))
fid.write(' ')
fid.write(str(SMOOTH))
fid.write(' ')
fid.write(str(REFRACTIVE))
fid.write(' ')
fid.write(str(lat))
fid.write(' ')
# write out the multipath function called
fid.write('moore')
#fid.write(' ')
#startDT = startDT - dt.timedelta(days=1)
#fid.write(startDT.strftime("%Y %j %H:%M:%S"))
fid.write('\n')
fid.close()

if option.plot:
    #fig1 = plt.figure()
    fig1 = plt.figure(figsize=(3.62, 2.76))
    ax = fig1.add_subplot(111)

    # Plot the residuals
    #ax.plot(v_El, vfull[:,0]*1000.,'b.' )
    # Plot the simulate bias
    ax.plot(m[:,0], m[:,2]*1000.,linewidth=2)
    # Plot the model
    ax.plot(m[:,0], m[:,1]*1000.,linewidth=2)
    # Plot the difference
    ax.plot(m[:,0],(m[:,2]-m[:,1])*1000.,'.')

    ax.set_xlabel('Elevation Angle (degrees)')
    ax.set_ylabel('Bias (mm)')
    ax.set_xlim([eleAng, 90.])
    ax.set_ylim([-20., 20.])
    #txt = 'RMS = (mm)'.format(rms)
    rms = rms*1000.
    ax.text(50, -10, 'RMS = {:.1f} (mm)'.format(rms),fontsize=8)

    ax.legend(['Simulated','Estimated','Difference'],fontsize=8)

    for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
                 ax.get_xticklabels() + ax.get_yticklabels()):
                item.set_fontsize(8)
    plt.tight_layout()

    fig1.savefig('simulation.png')
    plt.show()

