import numpy as np

from matplotlib import pyplot as plt
from matplotlib.dates import YearLocator, MonthLocator, DateFormatter

import datetime as dt

from optparse import OptionParser

# match a float 
def mfloat(a,b):
    if abs(a-b)<0.00000001:
        return 1
    else:
        return 0



parser = OptionParser()

parser.add_option("-f", "--filename", dest="filename",
                                  help="Result file to plot", metavar="FILE")

parser.add_option("-n","--north",dest="north", action="store_true", default=False, help="Plot the north component") 
parser.add_option("-e","--east",dest="east", action="store_true", default=False, help="Plot the east component") 
parser.add_option("-u","--up",dest="up", action="store_true", default=False, help="Plot the up component") 
parser.add_option("-t","--trop",dest="trop",type="int",default = 0,help="Plot n troposphere estimates")

parser.add_option("--MH",dest="MH", default = 0.2, type="string", help="Monument height to plot")


parser.add_option("--l1", "--legend1", dest="legend1", default="Data 1",
                                  help="Label for data being plotted to appear in the legend")

parser.add_option("--lat", dest="lat", type='string', help="Latitude of station to Plot")
parser.add_option("--yyyy",dest="YYYY",help="YYYY Time stamp to plot (One day only)")
parser.add_option("--ddd",dest="DDD", type=int, help="DDD Day-of-Year Time stamp to plot, make sure to use the YYYY option as well")
parser.add_option("-o", "--outfile", dest="outfile",
                                  help="output name for plot", metavar="FILE")

parser.add_option("-d", "--diff", dest="diff", action="store_true", default=False, help = "Calculate a difference")

parser.add_option("-q", "--quiet", dest="verbose", action="store_false", default=True,
                                    help = "Don't display plots")
(option, args) = parser.parse_args()

data = np.genfromtxt(option.filename)
ind = np.argsort(data[:,1])

# Check to see if we re only plotting on epoch
# then initialse a default value if not set
if option.YYYY :
    if option.DDD:
        DDD = option.DDD
    else:
        DDD = 001


# Plot the first file as a function of height
fig1 = plt.figure(figsize=(3.62, 2.76))
ax = fig1.add_subplot(111)
ax.set_xlabel('Time (years)')
ax.set_ylabel('RMS  (mm)')

# Plot a time series from a particular latitude
if option.lat :
    print("Option lat:",option.lat)
    if option.YYYY :
        print("Latitude vs feature",int(option.YYYY),DDD)
        # look for data that occurs on the same day as requested
        #criterion = ((data[:,1] == int(option.YYYY) ) & ((data[:,2]) == DDD))
        criterion = ((data[:,1] == int(option.YYYY) ) & ((data[:,2] -DDD) < 1.) )
        ind = np.array(np.where(criterion))

        ax.plot(data[ind,-2].T, data[ind,8].T*1000.)
        ax.set_xlabel('Latitude')
        ax.set_ylabel('Height Bias (mm)')
    else: 
        criterion = (data[:,-2] == float(option.lat))
        ind = np.array(np.where(criterion))

        FEATURE = 8
        ax.set_ylabel('Height Bias (mm)')

        time = [] 
        ctr = 0

        for i in data[:,1] :
            time.append( dt.datetime( int(data[ctr,1]), 1, 1 ) )
            time[ctr] = time[ctr] + dt.timedelta(days=(int(data[ctr,2]) -1))
            ctr += 1

        # format the ticks
        time = np.array(time)
        ax.plot_date(time[ind,], data[ind,FEATURE]*1000.) 
        #ax.set_ylim([1,4])
        #ax.suptitle("Latitude:",option.lat)

        years    = YearLocator()   # every year
        months   = MonthLocator()  # every month
        yearsFmt = DateFormatter('%Y')

        ax.xaxis.set_major_locator(years)
        ax.xaxis.set_major_formatter(yearsFmt)
        ax.xaxis.set_minor_locator(months)
        ax.autoscale_view()

        ax.fmt_xdata = DateFormatter('%Y-%m-%d')
        ax.grid(True)

        fig1.autofmt_xdate()

else:
    #
    # year = data[ind,1]
    # DoY = data[ind,2]
    # HH:MM:SS data[ind,3]
    # h = data[ind,4]       height of monument
    # rms = data[ind,5]
    # north_bias = data[ind,6]
    # east_bias  = data[ind,7]
    # up_bias    = data[ind,8]
    #
    criterion = (data[:,4] == float(option.MH))

    ind = np.array(np.where(criterion))
    time = [] #np.zeros(data.shape[0])
    ctr = 0
    for i in data[:,1] :
        time.append( dt.datetime( int(data[ctr,1]), 1, 1 ) )
        time[ctr] = time[ctr] + dt.timedelta(days=(int(data[ctr,2]) -1))
        ctr += 1

    # set the default for the rms of residuals?
    FEATURE = 5
    if option.up:
        FEATURE = 8
        ax.set_ylabel('Height Bias(mm)')
        ax.set_ylim([0., 3.5])
 
    # format the ticks
    time = np.array(time)
    ax.plot_date(time[ind,], data[ind,FEATURE]*1000.) 

    years    = YearLocator()   # every year
    months   = MonthLocator()  # every month
    yearsFmt = DateFormatter('%Y')

    ax.xaxis.set_major_locator(years)
    ax.xaxis.set_major_formatter(yearsFmt)
    ax.xaxis.set_minor_locator(months)
    ax.autoscale_view()

    ax.fmt_xdata = DateFormatter('%Y-%m-%d')
    ax.grid(True)

    fig1.autofmt_xdate()


for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
                         ax.get_xticklabels() + ax.get_yticklabels()):
    item.set_fontsize(8)

plt.tight_layout()

if option.outfile:
    fig1.savefig(option.outfile)

if option.verbose:
    plt.show()
