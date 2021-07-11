import numpy as np
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 15})

names = ['Time','avgI1rd','avgI2rd','avgDrd','Beta','Beff']
ts = 480
replicates = 20
parameter = 'k11'
p = [0.0 , 0.25, 0.5, 0.75, 1.0]
bpI1 = np.zeros(len(p))
bpI1Min = np.zeros(len(p))
bpI1Max = np.zeros(len(p))
cell_diameter = 3.0
for i in range(len(p)):
    avgI1rd = np.zeros((ts, replicates))
    avgI2rd = np.zeros((ts, replicates))
    avgDrd = np.zeros((ts, replicates))
    bpI1r = np.zeros(replicates)
    for r in range(1,replicates+1):
        f = np.genfromtxt('PlaqueAssay_%s_%.2f_%i.txt' % (parameter,p[i],r), skip_header=1, delimiter=',', names=names,max_rows=ts)
        avgI1rd[:, r - 1] = f['avgI1rd']
        avgI2rd[:, r - 1] = f['avgI2rd']
        avgDrd[:, r - 1] = f['avgDrd']
        [mI1, y0I1] = np.polyfit(f['Time'][len(f['Time'])//2-10:len(f['Time'])//2+10],f['avgI1rd'][len(f['avgI1rd'])//2-10:len(f['avgI1rd'])//2+10], 1)
        bpI1r[ r - 1] = mI1
    plt.plot(f['Time'], np.mean(avgI1rd, 1) / cell_diameter, color='orange', linewidth=4.0, label='I1')
    plt.fill_between(f['Time'], np.amin(avgI1rd, 1) / cell_diameter, np.amax(avgI1rd, 1) / cell_diameter,
                     facecolor='orange', interpolate=True, alpha=0.2)
    plt.plot(f['Time'], np.mean(avgI2rd, 1) / cell_diameter, color='red', linewidth=4.0, label='I2')
    plt.fill_between(f['Time'], np.amin(avgI2rd, 1) / cell_diameter, np.amax(avgI2rd, 1) / cell_diameter,
                     facecolor='red', interpolate=True, alpha=0.2)
    plt.plot(f['Time'], np.mean(avgDrd, 1) / cell_diameter, color='purple', linewidth=4.0, label='D')
    plt.fill_between(f['Time'], np.amin(avgDrd, 1) / cell_diameter, np.amax(avgDrd, 1) / cell_diameter,
                     facecolor='purple', interpolate=True, alpha=0.2)
    plt.grid()
    plt.xlabel('Time (hrs)')
    if i == 0:
        plt.ylabel('Cell Diameters')
    plt.ylim([0, 100])
    plt.xlim([0, 80])
    # plt.legend(loc=2)
    if i> 0 :
        plt.gca().set_yticklabels([])
    plt.gca().set_aspect(1.0 / plt.gca().get_data_ratio() * 1.0)
    plt.savefig('Fig.%s_%.2f.pdf' % (parameter, p[i]), transparent = True)
    plt.clf()
    bpI1[i] = np.mean(bpI1r)
    bpI1Min[i] = np.amin(bpI1r)
    bpI1Max[i] = np.amax(bpI1r)

plt.plot(p,bpI1,color='black',linewidth = 4.0, label = 'Plaque Growth Rate')
plt.fill_between(p, bpI1Min,bpI1Max, facecolor='black', interpolate=True, alpha=0.2)
plt.title('Plaque Growth Rate')
plt.ylabel('Plaque Growth Rate (cell diameter/hrs)')
plt.xlabel('Parameter Multiplier')
plt.gca().set_aspect(1.0 / plt.gca().get_data_ratio() * 1.0)
plt.gca().set_xticks(p)
plt.legend(loc=2)
plt.grid()
plt.savefig('Fig.Growth_Rate_%s.pdf' % parameter, transparent=True)