import numpy as np

from matplotlib import pyplot as plt
from matplotlib.dates import YearLocator, MonthLocator, DateFormatter

import datetime as dt

from optparse import OptionParser

parser = OptionParser()

parser.add_option("-f", "--filename", dest="filename",
                                  help="Result file to plot", metavar="FILE")

parser.add_option("--f1", "--filename1", dest="filename",
                                  help="Result file to plot", metavar="FILE")

parser.add_option("-n","--north",dest="north", action="store_true", default=False, help="Plot the north component") 
parser.add_option("-e","--east",dest="east", action="store_true", default=False, help="Plot the east component") 
parser.add_option("-u","--up",dest="up", action="store_true", default=False, help="Plot the up component") 
parser.add_option("-t","--trop",dest="trop",type="int",default = 0,help="Plot n troposphere estimates")

parser.add_option("--MH",dest="MH", default = 0.2, type=float, help="Monument height to plot")


parser.add_option("--l1", "--legend1", dest="legend1", default="Data 1",
                                  help="Label for data being plotted to appear in the legend")

parser.add_option("--f2", "--filename2", dest="filename2",
                                  help="Result file to plot", metavar="FILE")

parser.add_option("--l2", "--legend2", dest="legend2", default="Data 2",
                                  help="Label for data being plotted to appear in the legend")

parser.add_option("--f3", "--filename3", dest="filename3",
                                  help="Result file to plot", metavar="FILE")

parser.add_option("--l3", "--legend3", dest="legend3", default="Data 3",
                                  help="Label for data being plotted to appear in the legend")

parser.add_option("--lat", dest="lat", action='store_true',default=False)

parser.add_option("-o", "--outfile", dest="outfile",
                                  help="output name for plot", metavar="FILE")

parser.add_option("-d", "--diff", dest="diff", action="store_true", default=False, help = "Calculate a difference")

parser.add_option("-q", "--quiet", dest="verbose", action="store_false", default=True,
                                    help = "Don't display plots")
(option, args) = parser.parse_args()

data = np.genfromtxt(option.filename)
ind = np.argsort(data[:,1])

# Plot the first file as a function of height
fig1 = plt.figure(figsize=(3.62, 2.76))
ax = fig1.add_subplot(111)
ax.set_xlabel('Monument Height (m)')
ax.set_ylabel('RMS  (mm)')

# If we have two files then we want to do a comparison of results
if option.lat :
    data2 = np.genfromtxt(option.filename2)

    ind  = np.argsort(data[:,-2])
    ind2 = np.argsort(data2[:,-2])

    ax.plot(data[ind,-2], data[ind,1], data2[ind2,-2], data2[ind2,1])#,linewidth=2)
    ax.legend([option.legend1,option.legend2],fontsize=8)
    ax.set_xlabel('Latitude (degrees)')
    ax.set_xlim([-90, 90])
    ax.set_ylabel('RMS  (mm)')
elif option.filename3 :
    data2 = np.genfromtxt(option.filename2)
    data3 = np.genfromtxt(option.filename3)
    ind2  = np.argsort(data2[:,1])
    ind3  = np.argsort(data3[:,1])

    ax.plot(data[ind,1],   data[ind,2],'b',linewidth=2)
    ax.plot(data2[ind2,1], data2[ind2,2],'r',linewidth=2)
    ax.plot(data3[ind3,1], data3[ind3,2],'k',linewidth=2)

    ax.legend([option.legend1,option.legend2,option.legend3],fontsize=8)
elif option.filename2 :
    data2 = np.genfromtxt(option.filename2)
    ind = np.argsort(data[:,1])
    ind2 = np.argsort(data2[:,1])
    if option.diff:
        ax.plot(data[ind,1], (np.abs(data[ind,2])-np.abs(data2[ind2,2])),linewidth=2)
    else:
        ax.plot(data[ind,1], data[ind,2], data2[ind2,1], data2[ind2,2],linewidth=2)
        ax.legend([option.legend1,option.legend2],fontsize=8)
elif option.north :
    ax.plot(data[ind,1], data[ind,3]*1000.,linewidth=2)
    ax.set_ylabel('North Bias (mm)')
elif option.east :
    ax.plot(data[ind,1], data[ind,4]*1000.,linewidth=2)
    ax.set_ylabel('East Bias (mm)')
elif option.trop > 0:
    #idx = 4 + option.trop
    for idx in range(6,option.trop+5):
        ax.plot(data[ind,1], data[ind,idx]*100,'b.',alpha=0.7)
    ax.set_ylabel('Troposphere Bias')
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
    #criterion = (data[:,4] == 0.20)
    criterion = (data[:,4] == option.MH)
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
        ax.set_ylabel('Height Bias(m)')
 
    # format the ticks
    time = np.array(time)
    ax.plot_date(time[ind,], data[ind,FEATURE]) 

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
