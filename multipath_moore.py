import numpy as np

def multipath_moore(surff, H, AntennaType, np):
    """% [dsigmaL1,dsigmaL2,dsigmaLc]=multipath_herring(elev, rough, H, Grate, np)
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

    vel_light = 299792458.0
    fL1 = 10.23e6*77.*2. 
    fL2 = 10.23e6*60.*2. 
    wL1 = vel_light/fL1 
    wL2 = vel_light/fL2
    lcl1 = 1./(1.-(fL2/fL1)**2)
    lcl2 = -(fL2/fL1)/(1.-(fL2/fL1)**2)

    # Set refractive index of media which determines the refrection
    # coefficients. na is air, np is material below the antenna
    # The dielectric constant is simply the square of the refractive index in a non-magnetic medium 
    n1 = 1.0 # refractive coefficient of air
    np = np.sqrt(4) # dieletric costant of dry sand (3-5), concrete = 4
    #water is sqrt(88-80.1-55.3-34.5) (0-20-100-200 degC) wikipedia

    #
    # For an AT504 mean = 1.1  
    #
    Grate_L1_504 = 1.1044
    Grate_L2_504 = 1.0931

    GRID = 0.1
    #GRID = 0.5
    azs = np.linspace(0,360,(360./GRID) +1 )
    i = 1

    for az in azs: 
  
        j = 1
        zds = np.linspace(0,90,(90./GRID)+1 )

        for zd in zds : 
            e = 90. - zd
            eRad = e*np.pi/180.
            zRad = zd*np.pi/180.

            # had a bug in previous version N = floor(el*100 + 1) (when incrementing by 0.01 -> produced 
            # a counter that repeated values, so now just use n and perform a fliplr
            n = np.floor(zd*10+1)

            #
            # NB this is using the noazi values
            # Now compute the height of the reflector taking into account the antenna PCV
            #
            htL1 = H # + antL1.neu(1,3)/1000 + pcvL1(i,j)/1000 * sin(eRad);
            htL2 = H # + antL2.neu(1,3)/1000 + pcvL2(i,j)/1000 * sin(eRad);

            Ra = ( n1*np.cos(zRad) - np.sqrt( np**2 - (n1 * np.sin(zRad))**2  )/
                    n1*np.cos(zRad) + np.sqrt( np**2 - (n1 * np.sin(zRad))**2  ) )

            # Reflected gain from antenna
            L1_back = surff * backL1(n) * Ra;
            L2_back = surff * backL2(n) * Ra;

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
                        atan((L1_back*np.sin(4.*np.pi*(htL1/wL1)*np.sin(eRad))) / 
                            (gDL1_504 + L1_back * np.cos(4*np.pi*(htL1/wL1)*np.sin(eRad))))
                        )
            dRL2_back = ( wL2/(2.*np.pi) * 
                        np.atan((L2_back*np.sin(4*np.pi*(htL2/wL2)*np.sin(eRad))) /
                                (gDL2_504 + L2_back * np.cos(4*np.pi*(htL2/wL2)*np.sin(eRad))))
                        )
            #dR(j) = (lcl1 * dRL1_back) + (lcl2 * dRL2_back) ;
            dR(j) = (2.5457 * dRL1_back) - (1.545 * dRL2_back) 

            j = j+1

        # Make the vector a function of elevation, rather than zenith for use in the simulation
        dsigmaLc(i,:) = fliplr(dR) # ./1000;
        #dsigmaLc(i,:) = dR;

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

    G1s = [ -26.0 -26.0 -26.0 -26.0   -26.0   -26.0   -20.0   -19.0   -18.0   -17.0   -16.0   -15.0   -15.35   -15.7   -16.05   -16.4   -16.75   -17.1   -17.5 ]
    G2s = [ -28.5 -28.5 -28.5 -28.5   -28.5   -26.0   -24.0   -22.0   -20.0   -18.0   -17.0   -16.0   -16.25   -16.5   -16.75   -17.0   -17.15   -17.3   -17.5 ]

    for i=1:1:19
        G1s(i) = calcAmplitude(G1s(i),0)
        G2s(i) = calcAmplitude(G2s(i),0)

    yi1 = interp1(0:5:90,G1s,0:0.1:90,'linear')
    yi2 = interp1(0:5:90,G2s,0:0.1:90,'linear')

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
    G1s = [ -20.0 -20.5 -21.0 -22.0   -23.5   -25.0   -27.0   -24.5   -24.0   -23.5   -23.0   -22.7   -22.50   -22.3   -22.75   -23.0   -24.0    -25.0   -26.0 ]
    G2s = [ -28.5 -28.5 -28.5 -28.5   -28.5   -26.0   -24.0   -22.0   -20.0   -18.0   -17.0   -16.0   -16.25   -16.5   -16.75   -17.0   -17.15   -17.3   -17.5 ]

    for i=1:1:19 :
        G1s(i) = calcAmplitude(G1s(i),0)
        G2s(i) = calcAmplitude(G2s(i),0)

    yi1 = interp1(0:5:90,G1s,0:0.1:90,'linear')
    yi2 = interp1(0:5:90,G2s,0:0.1:90,'linear')

  return yi1,yi2

def function calcAmplitude(Gain, Gain0):
    """ 
    A = calcAmpiltude(Gain,Gain0)
    Convert the gain from decibels, into amplitude
    """
    A = 10^(Gain/20) / 10^(Gain0/20)
    return A

