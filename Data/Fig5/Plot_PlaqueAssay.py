import numpy as np
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 15})

names = ['Time','avgI1rd','avgI2rd','avgDrd']
ts = 480
replicates = 20
parameter = 'k31'
p = [1.0,5.0,10.0, 50.0, 100.0]
# p = [1.0,3.98, 5.0, 6.31, 7.94, 10.0, 12.59, 15.85, 19.95, 25.12, 39.81, 50.0, 63.10, 79.43, 100.0, 125.89]

bpI1 = np.zeros(len(p))
bpI1Min = np.zeros(len(p))
bpI1Max = np.zeros(len(p))
cell_diameter = 3.0

mean = True
plot_Individual = False

def plot_Radius(T,data,name, color = None, mean= False):
    if not mean:
        plt.plot(T, np.percentile(data, 50, axis = 1) / cell_diameter, color=color, linewidth=4.0, label=name)
        lower_bound = np.percentile(data, 5, axis = 1) / cell_diameter
        upper_bound = np.percentile(data, 95, axis = 1) / cell_diameter
    else:
        plt.plot(T, np.mean(data, axis = 1) / cell_diameter, color=color, linewidth=4.0, label=name)
        lower_bound = (np.mean(data, axis = 1) - np.std(data, axis = 1)) / cell_diameter
        upper_bound = (np.mean(data, axis = 1) + np.std(data, axis = 1)) / cell_diameter
    plt.fill_between(T, lower_bound, upper_bound,facecolor=color, interpolate=True, alpha=0.2)

for i in range(len(p)):
    avgI1rd = np.zeros((ts, replicates))
    avgI2rd = np.zeros((ts, replicates))
    avgDrd = np.zeros((ts, replicates))
    bpI1r = np.zeros(replicates)
    for r in range(1,replicates+1):
        f = np.genfromtxt('Data/PlaqueAssay_%s_%.2f_%i.txt' % (parameter,p[i],r), skip_header=1, delimiter=',', names=names,max_rows=ts)
        avgI1rd[:, r - 1] = f['avgI1rd']
        avgI2rd[:, r - 1] = f['avgI2rd']
        avgDrd[:, r - 1] = f['avgDrd']
        [mI1, y0I1] = np.polyfit(f['Time'][-10:],f['avgI1rd'][-10:], 1)
        bpI1r[ r - 1] = mI1
        if plot_Individual:
            plt.plot(f['Time'], f['avgI1rd'] / cell_diameter, color='orange', linewidth=4.0, label='I1')
            plt.plot(f['Time'], f['avgI2rd'] / cell_diameter, color='red', linewidth=4.0, label='I2')
            plt.plot(f['Time'], f['avgDrd'] / cell_diameter, color='purple', linewidth=4.0, label='D')
            plt.grid()
            plt.ylabel('Cell Diameters')
            plt.xlabel('Time (hrs)')
            plt.ylim([0, 100])
            plt.xlim([0, 80])
            plt.gca().set_aspect(1.0 / plt.gca().get_data_ratio() * 1.0)
            plt.savefig('Fig.%s_%.2f_%i.pdf' % (parameter, p[i], r), transparent=True)
            plt.clf()
    plot_Radius(f['Time'], avgI1rd, 'I1', color='orange', mean=mean)
    plot_Radius(f['Time'], avgI2rd, 'I2', color='red', mean=mean)
    plot_Radius(f['Time'], avgDrd, 'D', color='purple', mean=mean)
    plt.grid()
    plt.xlabel('Time (hrs)')
    if i == 0:
        plt.ylabel('Cell Diameters')
    plt.ylim([0, 100])
    plt.xlim([0, 80])
    # plt.legend(loc=2)
    if i > 0 :
        plt.gca().set_yticklabels([])
    plt.gca().set_aspect(1.0 / plt.gca().get_data_ratio() * 1.0)
    if not mean:
        plt.savefig('Fig.median.%s_%.2f.pdf' % (parameter, p[i]), transparent = True)
    else:
        plt.savefig('Fig.mean.%s_%.2f.pdf' % (parameter, p[i]), transparent = True)
    plt.clf()

