import numpy as np
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 15})

ts = 480
replicates = 20
cell_diameter = 3
Names = ['Time','avgI1rd','avgI2rd','avgDrd','Beta','Beff']
mean = False

def plot_Radius(variable, color = None, mean = False):
    avgrd = np.zeros((ts, replicates))
    for r in range(1, replicates + 1):
        f = np.genfromtxt('Data/PlaqueAssay_%i.txt' % (r), skip_header=1, delimiter=',', names=Names, max_rows=ts)
        avgrd[:, r - 1] = f['avg%srd' % variable]
    if not mean:
        plt.plot(f['Time'], np.percentile(avgrd, 50, axis = 1) / cell_diameter, color=color, linewidth=4.0, label=variable)
        lower_bound = np.percentile(avgrd, 5, axis = 1) / cell_diameter
        upper_bound = np.percentile(avgrd, 95, axis=1) / cell_diameter
    else:
        plt.plot(f['Time'], np.mean(avgrd, axis = 1) / cell_diameter, color=color, linewidth=4.0, label=variable)
        lower_bound = (np.mean(avgrd, axis = 1) - np.std(avgrd, axis = 1)) / cell_diameter
        upper_bound = (np.mean(avgrd, axis = 1) + np.std(avgrd, axis = 1)) / cell_diameter
    plt.fill_between(f['Time'], lower_bound, upper_bound, facecolor=color, interpolate=True, alpha=0.2)
plot_Radius('I1',color= 'orange', mean = mean)
plot_Radius('I2',color= 'red', mean = mean)
plot_Radius('D',color= 'purple', mean = mean)
plt.grid()
plt.title('Plaque Radius')
plt.xlabel('Time (hrs)')
plt.ylabel('Cell Diameters')
plt.ylim([0,100])
plt.xlim([0,80])
plt.legend(loc=2)
plt.gca().set_aspect(1.0/plt.gca().get_data_ratio()*1.0)
if not mean:
    plt.savefig('Fig.median.Radius.pdf', transparent = True)
else:
    plt.savefig('Fig.mean.Radius.pdf', transparent=True)
plt.clf()