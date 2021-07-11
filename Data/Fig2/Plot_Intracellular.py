import numpy as np
import matplotlib.pyplot as plt

plt.rcParams.update({'font.size': 15})
plt.rcParams['svg.fonttype'] = 'none'

CC3Dnames = ['Time', 'CC3DV', 'CC3DH', 'CC3DP', 'CC3DIFNe_Scalar', 'CC3DIFNe_Field', 'CC3DSTATP', 'CC3DIRF7',
             'CC3DIRF7P', 'CC3DIFN']
ODEnames = ['Time', 'ODEV', 'ODEH', 'ODEP', 'ODEIFNe', 'ODESTATP', 'ODEIRF7', 'ODEIRF7P', 'ODEIFN']

CC3Dts = 180
replicates = 20

mean = True

def Plot_Intracellular(variable, color=None, title=None, xlabel='Time (hrs)', ylabel=None , mean = False):
    ODEf = np.genfromtxt('Data/JordanOriginalODE_1.txt', skip_header=1, delimiter=',', names=ODEnames, max_rows=CC3Dts)
    V = np.zeros((CC3Dts, replicates))
    for r in range(1, replicates + 1):
        f = np.genfromtxt('Data/JordanOriginalCC3D_%i.txt' % (r), skip_header=1, delimiter=',', names=CC3Dnames,
                          max_rows=CC3Dts)
        T = f['Time']
        V[:, r - 1] = f['CC3D%s' % variable]
    if variable == 'IFNe_Field':
        plt.plot(ODEf['Time'][::2], ODEf['ODEIFNe'][::2], '.', color=color, markersize=10, label='ODE')
    else:
        plt.plot(ODEf['Time'][::2], ODEf['ODE%s' % variable][::2], '.', color=color, markersize=10, label='ODE')
    if not mean:
        plt.plot(f['Time'],np.percentile(V, 50, axis =1),color=color,linewidth = 3.0)
        lower_bound = np.percentile(V, 5, axis = 1)
        upper_bound = np.percentile(V, 95, axis = 1)
    else:
        plt.plot(T, np.mean(V, axis=1), color=color, linewidth=3.0)
        lower_bound = np.mean(V, axis=1) - np.std(V, axis=1)
        upper_bound = np.mean(V, axis=1) + np.std(V, axis=1)
    plt.fill_between(T, lower_bound, upper_bound, facecolor=color, interpolate=True,
                     alpha=0.2)
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(r'%s' % ylabel)
    plt.legend(loc=1)
    plt.gca().set_aspect(1.0 / plt.gca().get_data_ratio() * 1.0)
    plt.grid()
    if not mean:
        plt.savefig('Fig.median.%s.pdf' % variable, transparent=True)
    else:
        plt.savefig('Fig.mean.%s.pdf' % variable, transparent=True)
    plt.clf()

Plot_Intracellular('V', color='#6600FF', title='Virus Level',ylabel='Virus Level', mean = mean)
Plot_Intracellular('H', color='#339900', title='Cell Viability', ylabel='Viability/cell', mean = mean)
Plot_Intracellular('P', color='red', title='Fraction of Live Cells', ylabel='Fraction of Live Cells', mean = mean)
Plot_Intracellular('IFNe_Field', color='#993300', title='Extracellular IFN', ylabel='[IFNe] $\mu$M/cell', mean = mean)
Plot_Intracellular('STATP', color='#0033FF', title='STATP', ylabel='[STATP] $\mu$M/cell', mean = mean)
Plot_Intracellular('IRF7', color='#FF6600', title='IRF7', ylabel='[IRF7] $\mu$M/cell', mean = mean)
Plot_Intracellular('IRF7P', color='#666666', title='IRF7P', ylabel='[IRF7P] $\mu$M/cell', mean = mean)
Plot_Intracellular('IFN', color='#CC0033', title='Intracellular IFN', ylabel='[IFN] $\mu$M/cell', mean = mean)
