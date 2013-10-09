import numpy as np

def cal2jd(yr,mn,dy) :
    """
    CAL2JD  Converts calendar date to Julian date using algorithm
    from "Practical Ephemeris Calculations" by Oliver Montenbruck
    (Springer-Verlag, 1989). Uses astronomical year for B.C. dates
    (2 BC = -1 yr). 

    Input:
        yr  : (int) YYYY
        mn  : (int) MM 01 to 12
        day : (int) DD 01 to 31

    Output:
        jd : julian date (float)
    """

#    if mn < 1 | mn > 12
#        warning('Invalid input month');
#        return
#    end

#    if dy < 1
#        if (mn == 2 & dy > 29) | (any(mn == [3 5 9 11]) & dy > 30) | (dy > 31)
#        warning('Invalid input day');
#        return
#        end
#    end

    if mn > 2:
        y = yr
        m = mn
    else:
        y = yr - 1
        m = mn + 12

    date1=4.5+31.*(10.+12.*1582.)   # Last day of Julian calendar (1582.10.04 Noon)
    date2=15.5+31.*(10.+12.*1582.)  # First day of Gregorian calendar (1582.10.15 Noon)
    date=dy+31.*(mn+12.*yr)

    if date <= date1:
        b = -2
    elif date >= date2 :
        b = np.fix(y/400.) - np.fix(y/100.)
    else:
        #warning('Dates between October 5 & 15, 1582 do not exist');
        return

    if y > 0:
        jd = np.fix(365.25*y) + np.fix(30.6001*(m+1)) + b + 1720996.5 + dy
    else:
        jd = np.fix(365.25*y-0.75) + np.fix(30.6001*(m+1)) + b + 1720996.5 + dy

    return jd

def jd2gps(jd):
    """
    % JD2GPS  Converts Julian date to GPS week number (since
    %   1980.01.06) and seconds of week. 
    % Usage:   [gpsweek,sow,rollover]=jd2gps(jd)
    % Input:   jd       - Julian date
    % Output:  gpsweek  - GPS week number
    %          sow      - seconds of week since 0 hr, Sun.
    %          rollover - number of GPS week rollovers (modulus 1024)

    % Copyright (c) 2011, Michael R. Craymer
    % All rights reserved.
    % Email: mike@craymer.com
    """
#    if jd < 0
#        warning('Julian date must be greater than or equal to zero');
#        return;
#    end

    jdgps = cal2jd(1980,1,6);    # beginning of GPS week numbering
    nweek = int(np.fix((jd-jdgps)/7.))
    sow = (jd - (jdgps+nweek*7)) * 3600*24
    rollover = np.fix(nweek/1024)  # rollover every 1024 weeks
    #gpsweek = mod(nweek,1024);
    gpsweek = int(nweek)

# rollover is being returned as an array?
# should just be an int

    return(gpsweek,sow,rollover)


def jd2cal(jd):
    """
    % JD2CAL  Converts Julian date to calendar date using algorithm
    %   from "Practical Ephemeris Calculations" by Oliver Montenbruck
    %   (Springer-Verlag, 1989). Must use astronomical year for B.C.
    %   dates (2 BC = -1 yr). Non-vectorized version. See also CAL2JD,
    %   DOY2JD, GPS2JD, JD2DOW, JD2DOY, JD2GPS, JD2YR, YR2JD.
    % Version: 24 Apr 99
    % Usage:   [yr, mn, dy]=jd2cal(jd)
    % Input:   jd - Julian date
    % Output:  yr - year of calendar date
    %          mn - month of calendar date
    %          dy - day of calendar date (including decimal)

    % Copyright (c) 2011, Michael R. Craymer
    % All rights reserved.
    % Email: mike@craymer.com
    """

    #if nargin ~= 1
    #  warning('Incorrect number of arguments');
    #    return;
    #end
    #if jd < 0
    #      warning('Julian date must be greater than or equal to zero');
    #     return;
    #end

    a = np.fix(jd+0.5)

    if a < 2299161. :
        c = a + 1524.
    else:
        b = np.fix( (a-1867216.25) / 36524.25 )
        c = a + b - np.fix(b/4.) + 1525.

    d = np.fix( (c-122.1)/365.25 )
    e = np.fix(365.25*d)
    f = np.fix( (c-e) / 30.6001 )
    dy = c - e - np.fix(30.6001*f) + np.remainder((jd+0.5),a)
    mn = f - 1. - 12. * np.fix(f/14.)
    yr = d - 4715. - np.fix( (7.+mn)/10. )

    return (yr,mn,dy)


def jd2doy(jd):
    """
    % JD2DOY  Converts Julian date to year and day of year.
    % . Non-vectorized version. See also CAL2JD, DOY2JD,
    %   GPS2JD, JD2CAL, JD2DOW, JD2GPS, JD2YR, YR2JD.
    % Version: 24 Apr 99
    % Usage:   [doy,yr]=jd2doy(jd)
    % Input:   jd  - Julian date
    % Output:  doy - day of year
    %          yr  - year

    % Copyright (c) 2011, Michael R. Craymer
    % All rights reserved.
    % Email: mike@craymer.com
    """
    #if nargin ~= 1
    #  warning('Incorrect number of arguments');
    #    return;
    #end
    #if jd < 0
    #      warning('Julian date must be greater than or equal to zero');
    #        return;
    #    end

    [yr,mn,dy] = jd2cal(jd)
    doy = jd - cal2jd(yr,1,0)

    # MM ensure the doy is 0 padded
    doy = "%03d" % doy
    return yr, doy

def yyyy2yy(year):
    """
      yy = yyyy2yy(YYYY)

      return the yy form of YYYY
        yy - last two digits of YYYY
           - returned as an int

      very messy hack
    """
    yy = int( str(int(year))[-2] + str(int(year))[-1] )
    return(yy)

def dateTime2gpssow(dt):
    """
    Usage: week,sow = dateTime2gpssow(dateTime)

    Input:
        dt          is a python datetime object

    Output:
        week        gps week (int)
        sow         seconds into gpsweek since 0 hr, Sunday (float)


    """
    jd = cal2jd(dt.year,dt.month,dt.day)
    week, sow, rollover = jd2gps(jd)

    return week, sow
#=========================

if __name__ == "__main__":

    startymdhms = cal2jd(2012,01,01+(00/24)+(00/(24*60))+(00.0000/(24*3600)))
    (year,doy) = jd2doy(startymdhms)
    print(year,doy)
    # now obatain the yy version of YYYY
    yy = yyyy2yy(year)
    print(yy)
    yy = yyyy2yy(2012.7)
    print(yy)

#print(jd)
#2455927.5

#(gpsweek,sow,rollover) = gt.jd2gps(jd)
#print(gpsweek,sow,rollover)
#(1669, 0.0, array(1.0))
#=> should be (1669, 0.0, 1.0)

