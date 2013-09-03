from __future__ import division, print_function, absolute_import
 
import math as math
import numpy as np

import scipy.io as sio
from scipy import sparse
from scipy import interpolate

import array

import matplotlib.pyplot as plt

import gpsTime as gpst
import multipath as mp
import broadcastNavigation as bcast
import geodetic as geod

#=================================
from optparse import OptionParser

parser = OptionParser()
parser.add_option( "-H", "--height", dest="height", help="Monument Height", 
                   type="float",default=1.5 )

(option,args) = parser.parse_args()
#=================================

vel_c = 299792458 # vacuum speed of light, m/s

#
# Some hardcoded defaults
# to be changed to options
#
name   = 'BASE'
lat    = -45.0
lon    = 140.0
h      = 0.0
eleAng = 10.0

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
#sit1['XYZ'] = np.matrix(sit1['XYZ'])
#sit2['XYZ'] = np.matrix(sit2['XYZ'])
# set the navigation file name
navfile = { 'name' : 'brdc0010.10n',
            'binName' : 'brdc0010.10n.nav',
            'data' : []}

# turn elevation weighting on => 1, off => 0
ELE_WEIGHT = 1
# Use unhealthy satellites => 1, don't => 0
UNHEALTHY = 0

# number of days to process
NDAYS = 1

# sample rate of simulated observations
# Needs to be divisible by 86400
SAMP = 120.

# Solution session length, typically 24 hours
sess=86400

MAXSAT = 32 

# predict the A amd b matrix sizes
MAXOBS = int( sess/SAMP * MAXSAT * NDAYS) #30000

# Set the reolution of the grids
BIN = 0.5

# Set start time and end time (based on NDAYS)
YYYY = 2012
MM = 01
DD = 01

startymdhms = gpst.cal2jd(YYYY,MM,DD+(00/24)+(00/(24*60))+(00.0000/(24*3600))) #user setting
stopymdhms = gpst.cal2jd(YYYY,MM,DD+(((24*NDAYS)-1)/24)+(59/(24*60))+(59.999/(24*3600))) #user setting

amb_count = 0

# loop until we've finished all the sessions
# this will ensure we can have multiple days
while startymdhms < stopymdhms :

    #print np.shape(ndata)

    # class option:
    #nav1 = BroadcastNavigation()
    #print nav1
    #nfile = nav1.navfile()
    #nav1.navfile("brdc1010.12n") # set doesn;t seem to work?
    #print nfile

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
    print("NORTH_NEST :", north_nextest)
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
    trop_delta   = 3600 
    trop_nest    = round((time_stop_jd - time_start_jd)*(86400/trop_delta))
    trop_nextest = time_start_jd
    #trop_nextest=1e20;% uncomment to not estimate trop
    trop_count   = 0 #initialise the column counter for the parameter
    trop_time    = []
    A_trop       = sparse.lil_matrix( (MAXOBS, trop_nest) )

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
    harm_nextest  = 1e20   # don't estimate!
    harm_count    = 0.       # initialise the column counter for the parameter
    A_harm        = sparse.coo_matrix(np.zeros(MAXOBS),np.zeros(harm_nest))

    #sync.delta   = north.delta # only one sync per coordinate solution - as in Santerre
    sync_delta    = time_deltasec # one clock sync per epoch
    sync_nest     = round((time_stop_jd - time_start_jd)*(86400/sync_delta))
    #sync.nextest = time.start.jd  
    sync_nextest  = 1e20   #don't estimate!
    sync_count    = 0 #initialise the column counter for the parameter
    A_sync        = sparse.coo_matrix(np.zeros(MAXOBS),np.zeros(sync_nest))

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
    SMOOTH = 0.3
    HEIGHT = option.height                  #0.20
    ANTENNATYPE = 'LEIAT504        NONE'
    REFRACTIVE = 2
    GRID = 0.5
    mp_SIGMALC = mp.multipath_moore(SMOOTH,HEIGHT,ANTENNATYPE,REFRACTIVE,GRID,0)

    print("\nMP_SIGMALC",type(mp_SIGMALC),np.shape(mp_SIGMALC))

    x = np.linspace(0,90, int(90./GRID)+1 )
    y = np.linspace(0,360, int(360./GRID)+1 )
    print("before interpolate", x.shape, y.shape,x)
    f = interpolate.interp2d(x, y, mp_SIGMALC, kind='linear')
    #f = interpolate.RectBivariateSpline(x, y, mp_SIGMALC)
    print("after interpolate")

    # intialise observation counter
    obs_count = 1
    b = sparse.lil_matrix((MAXOBS,1))

    v_Az = sparse.lil_matrix((MAXOBS,1))
    #v_El = sparse.lil_matrix((MAXOBS,1))
    v_El = []
    v_time = sparse.lil_matrix((MAXOBS,1))

    # work out the navfile name we need to import
    # need to see if there is a routine that will read matlab binary format
    # in python
    #navfile.name = navfileDir + '/brdc' + doy + '012n'
    #navfile.biname = navfile.name + '.nav'
    # convert rinex nav to a binary file
    #rinexe()

    ndata = bcast.get_binary_eph('/Users/moore/code/matlab/king/mikemoore/brdc/brdc0010.12n.nav')
    #ndata = bcast.get_binary_eph('/Users/michael/code/matlab/king/mikemoore/brdc/brdc0010.12n.nav')

    print("Sample rate:",SAMP)
    # begin the epoch loop
    ep = np.linspace(time_start_jd,time_stop_jd,num=86400./SAMP)
    sat = np.zeros(MAXSAT)

    satD = {}
    for sv in range (1,MAXSAT+1):
        skey = str(sv)
        satD[skey] = {}

    epoch_ctr = 0
    for epoch_time in ep :
        # Do some time conversions
        #print(epoch_time)
        [gpsweek,gpssow,nul] = gpst.jd2gps(epoch_time)
        gpssow = np.round(gpssow)
        print("GPSSOW",gpssow)  
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
            satD[skey]['col'] = bcast.find_eph(ndata,sv,gpssow) #FIND_EPH  Finds the proper column in ephemeris array  

            if (satD[skey]['col'] > -1 ):
                # SATPOS Calculation of X,Y,Z coordinates at time t
                satD[skey]['pos'] = bcast.satpos(satD[skey]['col'],ndata,gpssow) 
 
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

        print("2nd epoch_time:",epoch_time)

        for sv in range(1,MAXSAT+1):
            skey = str(sv)

            if((epoch_time == time_start_jd) ):
                #initialise the last time this sat was seen
                satD[skey]['lastepoch'] = 0
             
            if(satD[skey]['col'] > -1):
#               %check if the Elevation angle is above the cutoff at both sites AND no obstruction 
#               % - if it is, it is an SD observation
                #print("Made it",epoch_time,skey,gpssow)
                if( (satD[skey]['El2'] > sit2['ElCutOff']) &
                     (satD[skey]['El1'] > sit1['ElCutOff']) ):
                    epochs['sats'].append(sv)
                    v_El.append(satD[skey]['El1'])
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
                A_trop[(obs_count-1,trop_count-1)]=1./np.cos(np.radians(90.-satD[skey]['El1']))
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
            if (sync_count > 0):
                A_sync[obs_count-1,sync_count-1]=1

            # Ambiguity - 1 per satellite per pass
            #        satD[skey]['ref'] = sv_count 
            #print("Amb:",obs_count,satD[skey]['ref'])
            A_amb[(obs_count-1, satD[skey]['ref']-1)]=1;

            # Set up the 'b' (O-C) Matrix - up bias for this epoch
            #print("epoch_time",epoch_time)
            #print("bias_t",bias_t)
            #print("biasu_values",biasu_values)
            epoch_ctr = epoch_ctr+1
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
            bias_mp = f(satD[skey]['El1'], satD[skey]['Az1'])

            b[obs_count-1,0] = ( np.sqrt( bias_en**2 + epochu_bias**2 - 
                                2. * bias_en * epochu_bias * np.cos(np.radians(90.+satD[skey]['El1'])) ) +
                                bias_mp
                                -(satD[skey]['Dist1']))
                #interp2(mp_ELEV,mp_AZIMUTH,mp_SIGMALC,satD[skey]['El1'], satD[skey]['Az1']) #rangebias - add range bias here

            #v_Az[obs_count-1,0] = satD[skey]['Az1'] #%keep the Azimuth and elevation for plotting the residuals
            #v_El.append(satD[skey]['El1'])
            #v_El[obs_count-1,0] = satD[skey]['El1']
            #v_time[obs_count-1,0] = epoch_time

            # increment obs counter (end of epoch loop)
            obs_count = obs_count +1

    # End of the day loop
    startymdhms += 86400.

# Now time to do the math...
print("coord estimates",north_count,east_count,up_count,"troposphere estimates",trop_count,"OBS count",obs_count)

# From up the A and Ql matrix without the unused row/columns

print("Checking Coord matrix")

# Need to do a more spohisticated check on the numbe rof parameters by checking if it is [], then summing
nparams = 3 + trop_count + amb['time'].shape[0]

A_full   = sparse.lil_matrix( ((obs_count-1),nparams) )

A_full[:,0] = A_north[0:obs_count-1,:]
A_full[:,1] = A_east[0:obs_count-1,:]
A_full[:,2] = A_up[0:obs_count-1,:]

#v_Az = v_Az[0:obs_count-1,:]
#v_El = v_El[0:obs_count-1,:]
v_El = np.array(v_El)
#v_time = v_time[0:obs_count-1,:]

if (harm_count==0):
    A_harm = []
else:
    A_harm = A_harm[0:obs_count-1,:]
         
if (sync_count==0):
    A_sync = []
else:
    A_sync = A_sync[0:obs_count-1,:]
                  
if (trop_count==0):
    A_trop = []
else:
    A_full[:,3:(trop_count+3)] = A_trop[0:obs_count-1,:]

if (amb['time'].shape[0] == 0):
    A_amb = []
else:
    print("About to sort out the ambiguity matrix")
    A_amb = A_amb[0:obs_count-1,:]
    A_full[:,(trop_count+4):nparams] = A_amb[:,0:(amb['time'].shape[0] -1)]

b = b[0:obs_count-1,:]
print("b shape:",np.shape(b))        

#%% Combine A sub-matrices
print("Afull shape:",np.shape(A_full))        

del A_north, A_east, A_up, A_trop, A_amb

#%% Set up the weight matrix (Ql)
print("Setting up weight matrix")
if(ELE_WEIGHT):
    P=sparse.diags((np.sin(np.radians(v_El)) + 0.001)**2,0) #;%0.00001 avoids problems when elev = 0
else:
    P=sparse.eye(A_full.shape[0])  

print("P shape:",P.shape)

#%% Add random noise to b matrix
#%fprintf('Adding 1mm random noise to ranges\n');
b = b.T + np.random.normal(0,0.001,b.shape[0])
print("b shape",b.shape)
#%% Matrix setup now complete

#%% Invert Matrix
print('\nNow inverting for the final solution and covariance\n\n')
Qx = np.linalg.pinv(A_full.todense().T * P.todense() * A_full.todense() )
print("Qx",Qx.shape)
#Qx=pinv(full(Afull'*P*Afull));
print("Now calculating dx")
dx = Qx *A_full.todense().T * P.todense() * b.T
print("dx:",dx.shape)

#%N=Afull'*P*Afull;
#%dx=N\(Afull'*P*b);

print("Calculating the residuals:")
vfull = A_full.todense() * dx - b.T

print("vfull",vfull.shape)
print("\nv_El",v_El.shape)

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

for ele in np.linspace(0,90,int(90/BIN +1)) :
    bind = ((v_El <= ele + BIN ) & (v_El > ele - BIN))
    ind = np.array(np.where(bind)) 
    tmp = np.array(vfull[ind,0])

    res = np.median(tmp) # * -1.
    #sim = f(ele,1.)
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
    #print("ELE:",ele,"Median:",val,"TMP:",tmp.shape)
    del tmp
    i += 1

# mask out the NaN from the calculation
criterion = np.logical_not( np.isnan(d[:,2]) )
ind = np.where(criterion)
diff = np.array( d[ind] )
print("Diff:",np.size(diff),np.shape(diff),np.shape(d))
stmp = sum(diff[0,:])
divisor = 1./np.size(diff)
mtmp = sval*divisor
rms = np.sqrt(mtmp)
rms2 = np.sqrt(sval * divisor )
print("SUM:",stmp,"Divisor:",divisor,"Mtmp:",mtmp,"RMS:",rms,"RMS2:",rms2)

rms = np.sqrt( sum(diff[0,:]) * 1./np.size(diff) )
#rms = np.sqrt( sum(d[:,2]) * 1./(np.size(d[:,2])) )

print("RMS:",rms)
#% calculate the difference from the simualted MP
# try changing vpart with vfull
#for i = 1:size(vpart)
#  diff(i) = b(i) + vpart(i);
#end

#diff(find(isnan(diff))) = 0;

#rms = sqrt( sum((diff.^2)) * 1/(size(diff,2)) )
#mn = mean(diff)
#===

#x = np.linspace(0,90, (90/0.5)+1 )
fig1 = plt.figure()
ax = fig1.add_subplot(111)

ax.plot(v_El, vfull[:,0]*1000.,'b.' )
ax.plot(m[:,0], m[:,1]*1000.,'k.' )
ax.plot(m[:,0], m[:,2]*1000.,'r.')
fig1.savefig('simulation.png')
#plt.show()
