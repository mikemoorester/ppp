import numpy as np

from matplotlib import pyplot as plt
from optparse import OptionParser


parser = OptionParser()
parser.add_option("-f", "--filename", dest="filename",
                                  help="Result file to plot", metavar="FILE")

parser.add_option("--f1", "--filename1", dest="filename",
                                  help="Result file to plot", metavar="FILE")

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

parser.add_option("-o", "--outfile", dest="outfile",
                                  help="output name for plot", metavar="FILE")

parser.add_option("-q", "--quiet", dest="verbose", action="store_false", default=True,
                                    help = "Don't display plots")
(option, args) = parser.parse_args()

data = np.genfromtxt(option.filename)
ind  = np.argsort(data[:,-2])

# Plot the first file as a function of height
fig1 = plt.figure(figsize=(3.62, 2.76))
#fig1 = plt.figure()
ax = fig1.add_subplot(111)
ax.set_xlabel('Monument Height (m)')
ax.set_ylabel('RMS  (mm)')

# If we have two files then we want to do a comparison of results
if option.filename3 :
    data2 = np.genfromtxt(option.filename2)
    data3 = np.genfromtxt(option.filename3)
    ind2 = np.argsort(data2[:,-2])
    ind3  = np.argsort(data3[:,-2])

    ax.plot(data[ind,-2], data[ind,1], data2[ind2,-2], data2[ind2,1], data3[ind3,-2], data3[ind3,1])#,linewidth=2)
    ax.legend([option.legend1,option.legend2,option.legend3],fontsize=8)

elif option.filename2 :
    data2 = np.genfromtxt(option.filename2)
    ind2 = np.argsort(data2[:,-2])

    ax.plot(data[ind,-2], data[ind,1], data2[ind2,-2], data2[ind2,1],linewidth=2)
    ax.legend([option.legend1,option.legend2],fontsize=8)
else:
    ax.plot(data[ind,-2], data[ind,1],linewidth=2)

# vertical line from 143 to 
ax.plot([-60, -60], [0, 6], 'k--', lw=1)
ax.plot([60, 60], [0, 6], 'k--', lw=1)

#plt.xticks([-90,-60,-30,0,30,60,90],['-90','-60','-30','0','30','60','90'])
plt.xticks([-90,-60,-30,0,30,60,90])

ax.set_xlabel('Latitude (degrees)')
ax.set_xlim([-90, 90])
ax.set_ylabel('RMS  (mm)')

for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
                         ax.get_xticklabels() + ax.get_yticklabels()):
                    item.set_fontsize(8)

plt.tight_layout()

if option.outfile :
    fig1.savefig(option.outfile)

if option.verbose:
    plt.show()
