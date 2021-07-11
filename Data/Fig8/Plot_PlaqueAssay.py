import numpy as np
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 15})

names = ['Time','avgI1rd','avgI2rd','avgDrd','Beta','Beff']
ts = 468
replicates = 20
parameter = 'dcs'
p_I = [5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0]
p_V = [1.00, 2.00, 3.00, 4.00]
cell_diameter = 3.0

mean = False

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


for i in range(len(p_I)):
    for j in range(len(p_V)):
        avgI1rd = np.zeros((ts, replicates))
        avgI2rd = np.zeros((ts, replicates))
        avgDrd = np.zeros((ts, replicates))
        # plt.title('Plaque Growth\n%s *= %.2f,\n%s *= %.2f' % ('IFN dc', p_I[i],'Virus dc', p_V[j]))
        for r in range(1,replicates+1):
            f = np.genfromtxt('Data/PlaqueAssay_%s_%.2f_%.2f_%i.txt' % (parameter,p_I[i],p_V[j],r), skip_header=1, delimiter=',', names=names,max_rows=ts)
            avgI1rd[:, r - 1] = f['avgI1rd']
            avgI2rd[:, r - 1] = f['avgI2rd']
            avgDrd[:, r - 1] = f['avgDrd']
        plot_Radius(f['Time'], avgI1rd, 'I1', color='orange', mean=mean)
        plot_Radius(f['Time'], avgI2rd, 'I2', color='red', mean=mean)
        plot_Radius(f['Time'], avgDrd, 'D', color='purple', mean=mean)
        plt.grid()
        plt.ylim([0, 80])
        plt.xlim([0, 80])
        # plt.legend(loc=2)
        if j == 0:
            plt.ylabel('Cell Diameters')
        if j > 0:
            plt.gca().set_yticklabels([])
        if i == 0:
            plt.xlabel('Time (hrs)')
        if i > 0:
            plt.gca().set_xticklabels([])
        plt.gca().set_aspect(1.0 / plt.gca().get_data_ratio() * 1.0)
        if not mean:
            plt.savefig('Fig.median.%s_%.2f_%.2f.pdf' % (parameter, p_I[i],p_V[j]), transparent = True)
        else:
            plt.savefig('Fig.mean.%s_%.2f_%.2f.pdf' % (parameter, p_I[i],p_V[j]), transparent = True)
        plt.clf()
