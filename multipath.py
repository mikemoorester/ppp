from __future__ import division, print_function, absolute_import

import numpy as np
from scipy.stats.stats import pearsonr,kendalltau

vel_light = 299792458.0
fL1 = 10.23e6*77.*2. 
fL2 = 10.23e6*60.*2. 
wL1 = vel_light/fL1 
wL2 = vel_light/fL2
lcl1 = 1./(1.-(fL2/fL1)**2)
lcl2 = -(fL2/fL1)/(1.-(fL2/fL1)**2)

def multipath_elosegui(elev,ht,alpha):
    '''

    '''
    elev = np.radians(elev)

    htL1 = ht
    htL2 = ht

    dRL1 = ( wL1/(2.*np.pi) * 
             np.arctan( (alpha * np.sin(4. * np.pi * (htL1/wL1) * np.sin(elev))) / 
                    (1. + alpha * np.cos(4. * np.pi * (htL1/wL1) * np.sin(elev))) )
           )

    dRL2 = ( wL2/(2.*np.pi) * 
             np.arctan( (alpha * np.sin(4. * np.pi * (htL2/wL2) * np.sin(elev)))/
             (1. + alpha * np.cos(4. * np.pi * (htL2/wL2) * np.sin(elev))) )
           )

    LCbias = (2.5457 * dRL1) - (1.5457 * dRL2) 

    return LCbias

def multipath_herring(elev, rough, H, Grate, nP):
    '''% [dsigmaL1,dsigmaL2,dsigmaLc]=multipath_herring(elev, rough, H, Grate, np)
    % Multipath model based on Tom Herring's code
    % uses Fresnel Equations and takes into account antenna gain for direct and
    % reflected signals
    % input H is H above reflector in m
    % input elev is elevation angle of satellite (vectorised)
    % input rough is the surface reflection quality (0<rough<=1)
    % input Grate is how quick the gain goes to zero (modified dipole)
    % input np is the refractive index of the reflecting surface
    %   The dielectric constant is simply the square of the refractive index in a non-magnetic medium 
    %np = sqrt(4); % dieltric constant dry sand (3-5); %concrete is 4; 
    %water is sqrt(88-80.1-55.3-34.5) (0-20-100-200 degC) wikipedia
    '''
    elev = np.radians(elev)
    ht=H  #metres
    surff = rough  # Surface reflection quality (roughness measure
    #Grate = 1.1 ; # Sets how quick the gain goes to zero (modified dipole)

    aL1c = 1.  
    aL2c = 1.

    # Set refractive index of media which determines the refrection
    # coefficients. na is air, np is material below the antenna
    # The dielectric constant is simply the square of the refractive index in a non-magnetic medium 
    na = 1.0 
    #np = sqrt(4); % dieltric constant dry sand (3-5); %concrete is 4; 
    #water is sqrt(88-80.1-55.3-34.5) (0-20-100-200 degC) wikipedia

    gL190 = np.cos(90./Grate*np.pi/180.) 
    gL290 = np.cos(90./Grate*np.pi/180.)

    zd = np.pi/2. - elev
    n = np.floor(zd*10.+1)

    # Additional path length for reflected signal
    dR = 2.*ht*np.sin(elev)

    # Phase differences all reflected path for L1 and L2
    dL1m = np.mod(dR/wL1,1.)*2.*np.pi 
    dL2m = np.mod(dR/wL2,1.)*2.*np.pi

    # Direct gain from antena
    gL1 = np.cos(zd/Grate) 
    gL2 = np.cos(zd/Grate)

    # Reflected gain from antenna
    #aL1 = gL190*aL1c*(1.-np.sin((np.pi/2.-zd))) 
    #aL2 = gL290*aL2c*(1.-np.sin((np.pi/2.-zd)))
    aL1 = gL190*aL1c*(1.-np.sin(elev)) 
    aL2 = gL290*aL2c*(1.-np.sin(elev))

    # E-field perpendicular to plane of incidence
    # this is Fresnel Equation
    Ra = ( (na * np.cos(zd) - np.sqrt(nP**2-(na*np.sin(zd))**2)) / 
           (na * np.cos(zd) + np.sqrt(nP**2-(na*np.sin(zd))**2)) )
    #
    # compute final amplitude of reflected sigals
    aL1 = surff*aL1*Ra 
    aL2 = surff*aL2*Ra

    dsigmaL1 = np.arctan2(gL1+aL1*np.cos(dL1m),aL1*np.sin(dL1m))*wL1/(2*np.pi)
    dsigmaL2 = np.arctan2(gL2+aL2*np.cos(dL2m),aL2*np.sin(dL2m))*wL2/(2*np.pi)

    #dsigmaLc =  -(lcl1*dsigmaL1 + lcl2*dsigmaL2) #units of m
    dsigmaLc =  (lcl1*dsigmaL1 + lcl2*dsigmaL2) #units of m

    return dsigmaLc



def multipath_moore(surff, H, AntennaType, nP, GRID, eleOnly):
    """% [dsigmaL1,dsigmaL2,dsigmaLc]=multipath_moore(surff, H, AntennaType, np)
    % Multipath model based on Tom Herring's code
    % uses Fresnel Equations and takes into account antenna gain for direct and
    % reflected signals
    % input H is H above reflector in m
    % input elev is elevation angle of satellite (vectorised)
    % input surff is the surface reflection quality (0<surff<=1)
    % input Grate is how quick the gain goes to zero (modified dipole)
    % input np is the refractive index of the reflecting surface
    %   The dielectric constant is simply the square of the refractive index in a non-magnetic medium 
    %np = sqrt(4); % dieltric constant dry sand (3-5); %concrete is 4; 
    %water is sqrt(88-80.1-55.3-34.5) (0-20-100-200 degC) wikipedia
    """

    [backL1, backL2]  = interpolateBackLobeGain_504()

    # Set refractive index of media which determines the refrection
    # coefficients. na is air, np is material below the antenna
    # The dielectric constant is simply the square of the refractive index in a non-magnetic medium 
    n1 = 1.0 # refractive coefficient of air
    nP = 2.0 #np.sqrt(4.) # dieletric costant of dry sand (3-5), concrete = 4
    #water is sqrt(88-80.1-55.3-34.5) (0-20-100-200 degC) wikipedia

    #
    # For an AT504 mean = 1.1  
    #
    Grate_L1_504 = 1.1044
    Grate_L2_504 = 1.0931

    #GRID = 0.1
    #GRID = 0.5
    if eleOnly == 0 :
        azs = np.linspace(0,360,(360./GRID) +1 )
        dsigmaLc = np.zeros((360/GRID + 1,90/GRID + 1))
    else :
        azs = [ 0 ] #np.linspace(0,360,(360./GRID) +1 )
        dsigmaLc = np.zeros((1,90/GRID + 1))

    i = 0
    #dsigmaLc = np.zeros((360/GRID + 1,90/GRID + 1))

    for az in azs: 
  
        j = 0
        zds = np.linspace(0,90,(90./GRID)+1 )
        dR = np.zeros(90/GRID + 1)

        for zd in zds : 
            e = 90. - zd
            eRad = e*np.pi/180.
            zRad = zd*np.pi/180.

            # had a bug in previous version N = floor(el*100 + 1) (when incrementing by 0.01 -> produced 
            # a counter that repeated values, so now just use n and perform a fliplr
            n = np.floor( zd * 10 ) # + 1)

            #
            # NB this is using the noazi values
            # Now compute the height of the reflector taking into account the antenna PCV
            #
            htL1 = H # + antL1.neu(1,3)/1000 + pcvL1(i,j)/1000 * sin(eRad);
            htL2 = H # + antL2.neu(1,3)/1000 + pcvL2(i,j)/1000 * sin(eRad);

            Ra = ( n1*np.cos(zRad) - np.sqrt( nP**2 - (n1 * np.sin(zRad))**2  )/
                    n1*np.cos(zRad) + np.sqrt( nP**2 - (n1 * np.sin(zRad))**2  ) )

            # Reflected gain from antenna
            L1_back = surff * backL1[n] * Ra;
            L2_back = surff * backL2[n] * Ra;

            #
            # The direct Gain from the antenna for L1 => gDL1 = cos(z/G)
            # Where Grate is the rate of change of the antenna gain with zenith angle
            #
            gDL1_504 = np.cos( zRad / Grate_L1_504 );
            gDL2_504 = np.cos( zRad / Grate_L2_504 );

            #
            # Reflected Gain rate
            # gR = cos(90/G)(1 - sin(elev))
            #
            dRL1_back = ( wL1/(2.*np.pi) * 
                        np.arctan((L1_back*np.sin(4.*np.pi*(htL1/wL1)*np.sin(eRad))) / 
                            (gDL1_504 + L1_back * np.cos(4*np.pi*(htL1/wL1)*np.sin(eRad))))
                        )
            dRL2_back = (wL2/(2.*np.pi) * 
                        np.arctan((L2_back*np.sin(4*np.pi*(htL2/wL2)*np.sin(eRad))) /
                                (gDL2_504 + L2_back * np.cos(4*np.pi*(htL2/wL2)*np.sin(eRad))))
                        )
            dR[j] = (2.5457 * dRL1_back) - (1.545 * dRL2_back) 

            j = j+1

        # Make the vector a function of elevation, rather than zenith for use in the simulation
        dsigmaLc[i,:] = dR[::-1]

        i = i+1

    return dsigmaLc

#
# scaled from a plot from bedford 2009 for LHCP
# for a LEIAR25
#
def interpolateBackLobeGain_LEIAR25():
    """
    [yi1,yi2] = interpolateBackLobeGain_LEIAR25()
    Interpolate the backlobe gain figures obtained from bedford 2009 for LHCP signals
    for the LEIAR25 antenna.
    """

    G1s = [ -26.0, -26.0, -26.0, -26.0, -26.0, -26.0, -20.0, -19.0, -18.0, -17.0, -16.0, -15.0, -15.35, -15.7, -16.05, -16.4, -16.75, -17.1, -17.5 ]
    G2s = [ -28.5, -28.5, -28.5, -28.5, -28.5, -26.0, -24.0, -22.0, -20.0, -18.0, -17.0, -16.0, -16.25, -16.5, -16.75, -17.0, -17.15, -17.3, -17.5 ]

    for i, item in enumerate(G1s):
        G1s[i] = calcAmplitude(G1s[i],0)
        G2s[i] = calcAmplitude(G2s[i],0)

    x = np.linspace(0,90, (90/5) + 1 )
    xi = np.linspace(0,90, (90/0.1) + 1 )
    yi1 = np.interp(xi,x,G1s)
    yi2 = np.interp(xi,x,G2s)

#
# scaled from a plot from bedford 2009 for LHCP
# for a AT504
#
def interpolateBackLobeGain_504():
    """
    [yi1,yi2] = interpolateBackLobeGain_504()
    Interpolate the backlobe gain figures obtained from bedford 2009 for LHCP signals
    for the LEIAT504 antenna.
    """
    G1s = [ -20.0, -20.5, -21.0, -22.0,   -23.5,   -25.0,   -27.0,   -24.5,   -24.0,   -23.5,   -23.0,   -22.7,   -22.50,   -22.3,   -22.75,   -23.0,   -24.0,    -25.0,   -26.0 ]
    G2s = [ -28.5, -28.5, -28.5, -28.5,   -28.5,   -26.0,   -24.0,   -22.0,   -20.0,   -18.0,   -17.0,   -16.0,   -16.25,   -16.5,   -16.75,   -17.0,   -17.15,   -17.3,   -17.5 ]

    for i, item in enumerate(G1s):
        G1s[i] = calcAmplitude(G1s[i],0)
        G2s[i] = calcAmplitude(G2s[i],0)

    x = np.linspace(0,90, (90/5) + 1 )
    xi = np.linspace(0,90, (90/0.1) + 1 )
    yi1 = np.interp(xi,x,G1s)
    yi2 = np.interp(xi,x,G2s)

    return yi1,yi2

def calcAmplitude(Gain, Gain0):
    """ 
    A = calcAmpiltude(Gain,Gain0)
    Convert the gain from decibels, into amplitude
    """
    A = 10.**(Gain/20.) / 10.**(Gain0/20.)
    return A

#def calcCorrelation(mpArray,weightArray):
#    """
#    corr = calcCorrelation(mpArray)
#    Calculate the correlation of a simulated multipath bias with a 
#    weighting or mapping function
#
#    Where,
#    mpArray is a simulated  multipath signal
#
#    Returns
#    corr the correlation coeeficient
#    """
#    pcc, pval = pearsonr(mpArray,weightArray)
#    #pcc, pval = kendalltau(mpArray,weightArray)
#    return pcc, pval


# profile the code
#import hotshot, hotshot.stats
#prof = hotshot.Profile("multipath.prof")
#prof.runcall("__main__")

#prof.close()

# print the results
#stats = htoshot.stats.load("multipath.prof")
#stats.strip_dirs()
#stats.sort_stats('time','calls')
#stats.print_stats(20)

if __name__ == "__main__":

    #import imp
    #vmf1 = imp.load_source('*','./trop/vmf1.py')
    import vmf1  
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    from matplotlib.widgets import Slider, Button, RadioButtons

    from optparse import OptionParser

    parser = OptionParser()
    parser.add_option("-c","--corr", dest="corr",action="store_true",default=False,
                        help="Calculate the correlation between a multipath signal and a weight function")
 
    parser.add_option("-i","--interactive", dest="inter",action="store_true",default=False,
                        help="Produce an interactive plot of multipath simulation")
 
    parser.add_option("-M","--moore", dest="moore",action="store_true",default=False,
                        help="Produce a plot of simulated multipath using the multipath_moore()")

    parser.add_option("-H","--herring", dest="herring",action="store_true",default=False,
                        help="Produce a plot of simulated multipath using the multipath_herring()")

    parser.add_option("-N","--normalise", dest="normalise",action="store_true",default=False,
                        help="Produce a plot of simulated multipath nromalised with 1/sin(e) ")

    (option,args) = parser.parse_args()

    gainrate = 1.1
    S = 0.3
    grid = 0.5
    nP = 4
    alpha = 0.1
    dry, wet = vmf1.vmf1map([45.2])
    print("DRY",dry,"WET",wet)
    # check to see if we want the interactive multipath plot
    if option.inter :

        ax = plt.subplot(111)
        plt.subplots_adjust(bottom=0.25)
        plt.axis([0,90,-0.035,0.035])
        ax.set_ylabel('Bias (m)')
        ax.set_xlabel('Elevaton Angle (degrees)')
        # Default height of monument
        h = 0.5
        xi = np.linspace(0,90, int(90/grid) + 1 )

        moore_mp = multipath_moore(S, h, 'LEIAT504        NONE', nP,grid,1)     
        moore_h  = multipath_herring(xi, S, h, gainrate, nP)
        moore_e  = multipath_elosegui(xi,h,alpha)

        l, = plt.plot( xi, moore_mp[0,:] ) #, xi, moore_h)  
        m, = plt.plot( xi, moore_h, 'r-' ) #, xi, moore_h)  
        e, = plt.plot( xi, moore_e, 'k-' ) #, xi, moore_h)  

        ax.legend(['enhanced','herring','elosegui'],fontsize=8)
        axcolor = 'lightgoldenrodyellow'
        axheight = plt.axes([0.18, 0.1, 0.65, 0.03], axisbg=axcolor)
        sheight = Slider(axheight, 'Height', 0.001, 2.0, valinit=h)

        def update(val):
            height = sheight.val
            moore_mp = multipath_moore(S, height, 'LEIAT504        NONE', 2,grid,1)     
            moore_h = multipath_herring(xi, S, height, gainrate, 2)
            moore_e  = multipath_elosegui(xi,height,alpha)
            l.set_ydata( moore_mp[0,:] )
            m.set_ydata( moore_h )
            e.set_ydata( moore_e )
            plt.draw()

        sheight.on_changed(update)

        resetax = plt.axes([0.8, 0.025, 0.1, 0.04])
        button = Button(resetax, 'Reset', color=axcolor, hovercolor='0.975')

        def reset(event):
            sheight.reset()

        button.on_clicked(reset)

        plt.show()
    
    xi = np.linspace(0,90, (90/0.5) + 1 )

    if option.moore :

        dsigmaLc_010 = multipath_moore(0.3, 0.17, 'LEIAT504        NONE', 2,0.5,1)
        #dsigmaLc_015 = multipath_moore(0.3, 0.15, 'LEIAT504        NONE', 2,0.5,1)
        #dsigmaLc_050 = multipath_moore(0.3, 0.50, 'LEIAT504        NONE', 2,0.5,1)
        dsigmaLc_150 = multipath_moore(0.3, 1.50, 'LEIAT504        NONE', 2,0.5,1)
        
        fig = plt.figure(figsize=(3.62, 2.76))
        #fig.tight_layout()
        ax = fig.add_subplot(111)

        #ax.plot( xi, dsigmaLc_010[0,:], xi, dsigmaLc_015[0,:], xi, dsigmaLc_050[0,:], xi, dsigmaLc_150[0,:],linewidth=2 )
        ax.plot( xi, dsigmaLc_010[0,:]*1000., xi, dsigmaLc_150[0,:]*1000.,linewidth=2 )
        ax.set_ylim([-25, 25])
        ax.set_xlabel('Elevation Angle (degrees)')
        ax.set_ylabel('Multipath Bias (mm)')
        ax.legend(['h = 0.17 m','h = 1.50 m'],fontsize=8)

        for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
                             ax.get_xticklabels() + ax.get_yticklabels()):
                item.set_fontsize(8)
        plt.tight_layout()
        plt.show()

    if option.herring :
        #multipath_herring(elev, rough, H, Grate, nP):
        herring_010 = multipath_herring(xi, 0.3, 0.10, 1.1, 2)
        herring_015 = multipath_herring(xi, 0.3, 0.17, 1.1, 2)
        herring_050 = multipath_herring(xi, 0.3, 0.5, 1.1, 2)
        herring_150 = multipath_herring(xi, 0.3, 1.5, 1.1, 2)

        fig2 = plt.figure()
        ax2 = fig2.add_subplot(111)
        ax2.plot( xi, herring_010, xi, herring_015, xi, herring_050, xi, herring_150 )
        item = ax2.set_xlabel('Elevation Angle (degrees)')
        ax2.set_ylabel('Multipath Bias (m)')
        ax2.legend(['h = 0.10 m','h = 0.17 m','h = 0.50 m','h = 1.50 m'])

        plt.show()

    if option.normalise :
        dsigmaLc_017 = multipath_moore(0.3, 0.17, 'LEIAT504        NONE', 2,0.5,1)
        dsigmaLc_150 = multipath_moore(0.3, 1.50, 'LEIAT504        NONE', 2,0.5,1)

        # get rid of the divide by 0 error
        small = 0.00000001
        xi[0] = xi[0] + small
        a = 00.0
        b = 3.

        #n = 1./(np.sin(np.radians(xi)))**2
        #n = 1./((b /(np.sin(np.radians(xi)))**2)/1000.)
        #n = 1./ ( ( a**2 + b**2 / ((np.sin(np.radians(xi)))**2)) /1000.)   
        n = 1./ ( ( a**2 + b**2 / ((np.sin(np.radians(xi)))**2)) /1000.)   
        #n = 1./ (a**2 + b**2 /(np.sin(np.radians(xi)))**2)   
        print(n.shape)
        fig3 = plt.figure()
        ax3 = fig3.add_subplot(111)
        ax3.plot(xi, dsigmaLc_017[0,:],'b-',linewidth=2)
        ax3.plot(xi, (dsigmaLc_017[0,:]/n),'b-.',linewidth=2)
        #ax3.plot(xi, dsigmaLc_150[0,:],'r-',linewidth=2)
        #ax3.plot(xi, (dsigmaLc_150[0,:]/n),'r-.',linewidth=2)
        #ax3.plot(xi, (n*dsigmaLc_150[0,:]),'r-.',linewidth=2)
        plt.show()

    if option.corr :
        GRID = 0.5 
        h = 0.0

        x = []
        y = []
        yy = []

        wetW = []
        wetWW = []
        dryD = []
        dryDD = []
        tropP = []

        xi = np.linspace(0,90, (90/GRID) + 1 )
        
        # get the vmf1 mapping function for some random default values
        dryZD, wetZD = vmf1.vmf1map(xi)
        # make the functions in order of elevation, rather than zenith
        dry = dryZD[::-1]
        wet = wetZD[::-1]

        for hs in range(0,200) :
            h += 0.01
            
            herring = multipath_herring(xi, 0.3, h, 1.1, 2)
            mp = multipath_moore(0.3, h, 'LEIAT504        NONE', 2,0.5,1)
            mpArray = np.array(mp[0,:])
            weight = []
            trop = []

            for i in range(0, int(90./GRID)+1 ) :
                weight.append( np.sin(np.radians(i*GRID))**2 )
                trop.append( np.sin(np.radians(i*GRID)) )
                i += 1

            weight = np.array(weight)
            #print(mpArray.shape,weight.shape)
            #p, z = calcCorrelation(herring,weight)

            p, z = pearsonr(herring,weight)
            pp, zz = pearsonr(mpArray,weight)
            y.append(p)
            yy.append(pp)

            #d, dz = pearsonr(herring,dry)
            #dd, ddz = pearsonr(mpArray,dry)
            #dryD.append(d)
            #dryDD.append(dd)
            
            w, wz = pearsonr(herring,wet)
            ww, wwz = pearsonr(mpArray,wet)
            wetW.append(w)
            wetWW.append(ww)

            t, tz = pearsonr(mpArray,trop)
            tropP.append(t)

            x.append(h)


        # Plot the correlation with the elevation mapping function
        fig4 = plt.figure(figsize=(3.62, 2.76))
        ax4 = fig4.add_subplot(111)
        #ax4.plot(x,y,x,yy)
        ax4.plot(x,yy)
        #ax4.legend(['herring','enhanced'])
        ax4.set_xlabel('Monument Height (m)')
        ax4.set_ylabel('Correlation Co-efficient')
        ax4.set_ylim([0, 1])
        ax4.set_xlim([0, 2])
        for item in ([ax4.title, ax4.xaxis.label, ax4.yaxis.label] +
                             ax4.get_xticklabels() + ax4.get_yticklabels()):
                item.set_fontsize(8)
        plt.tight_layout()
        #ax4.set_title('Correlation with weighting function ')

        # Plot the Correlation with troposphere
        fig5 = plt.figure(figsize=(3.62, 2.76))
        ax5 = fig5.add_subplot(111)
        #ax4.plot(x,y,x,yy)
        ax5.plot(x,tropP)
        #ax4.legend(['herring','enhanced'])
        ax5.set_xlabel('Monument Height (m)')
        ax5.set_ylabel('Correlation Co-efficient')
        ax5.set_ylim([0, 1])
        ax5.set_xlim([0, 2])
        for item in ([ax5.title, ax5.xaxis.label, ax5.yaxis.label] +
                             ax5.get_xticklabels() + ax5.get_yticklabels()):
                item.set_fontsize(8)
        plt.tight_layout()
        #ax4.set_title('Correlation with weighting function ')


        # Plot the correlation with the dry mapping function
        #fig5 = plt.figure()
        #ax5 = fig5.add_subplot(111)
        #ax5.plot(x,dryD,x,dryDD)
        #ax5.set_ylim([-1, 1])
        #ax5.legend(['herring','enhanced'])
        #ax5.set_title('Correlation with VMF1 Dry mapping function ')

        # plot the correlation with the wet mapping function
        fig6 = plt.figure(figsize=(3.62, 2.76))
        #fig6 = plt.figure()
        ax6 = fig6.add_subplot(111)
        #ax6.plot(x,wetW,x,wetWW)
        ax6.plot(x,wetWW)
        ax6.set_xlabel('Monument Height (m)')
        ax6.set_ylabel('Correlation Co-efficient')
        ax6.set_ylim([-0.5, 0.5])
        ax6.set_xlim([0, 2])
        for item in ([ax6.title, ax6.xaxis.label, ax6.yaxis.label] +
                             ax6.get_xticklabels() + ax6.get_yticklabels()):
                item.set_fontsize(8)
        plt.tight_layout()
        #ax6.legend(['herring','enhanced'])
        #ax6.set_title('Correlation with VMF1 Wet mapping function ')

        fig7 = plt.figure(figsize=(3.62, 2.76))
        ax7 = fig7.add_subplot(111)
        #ax7.plot(x,yy,x,tropP,x,wetWW)
        ax7.plot(x,yy,x,tropP)

        ax7.set_xlabel('Monument Height (m)')
        ax7.set_ylabel('Correlation Co-efficient')
        #ax7.set_ylim([-0.5, 1])
        ax7.set_ylim([0, 1])
        ax7.set_xlim([0, 2])
        #ax7.legend(['Elevation Weighting','Simplified Trop','VMF1 (wet)'],fontsize=8)
        ax7.legend(['Elevation Weighting','Simplified Trop'],fontsize=8)
        for item in ([ax7.title, ax7.xaxis.label, ax7.yaxis.label] +
                      ax7.get_xticklabels() + ax7.get_yticklabels()):
            item.set_fontsize(8)
        plt.tight_layout()

        #fig8 = plt.figure()
        #ax8 = fig8.add_subplot(111)
        #ax8.plot(xi,wet)
        #ax8.set_title('Wet mapping function ')

        plt.show()
