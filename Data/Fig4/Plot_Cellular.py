import numpy as np
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 15})

CellNames = ['Time','U','I1','I2','D','Ve','IFNe']
ts = 480
replicates = 20
mean = False

def plot_CellType(type , mean = False, color = None):
    T = np.zeros((ts, replicates))
    for r in range(1, replicates + 1):
        f = np.genfromtxt('Data/FullModelCellular_%i.txt' % (r), skip_header=1, delimiter=',', names=CellNames,max_rows = ts)
        T[:, r - 1] = f[type]
    if not mean:
        plt.plot(f['Time'], np.percentile(T, 50, axis=1), color=color, linewidth=4.0)
        lower_bound = np.percentile(T, 5, axis=1)
        upper_bound = np.percentile(T, 95, axis=1)
    if mean:
        plt.plot(f['Time'], np.mean(T, axis=1), color=color, linewidth=4.0)
        lower_bound = np.mean(T, axis=1) - np.std(T, axis=1)
        upper_bound = np.mean(T, axis=1) + np.std(T, axis=1)
    plt.fill_between(f['Time'], lower_bound, upper_bound, facecolor=color,interpolate=True, alpha=0.2)
    plt.title('Cells Types')
    plt.xlabel('Time (hrs)')
    plt.ylabel('Fraction of Cells by Type')

plot_CellType('U', color = 'blue', mean = mean)
plot_CellType('I1', color = 'orange', mean = mean)
plot_CellType('I2', color = 'red', mean = mean)
plot_CellType('D', color = 'purple', mean = mean)
plt.grid()
plt.gca().set_aspect(1.0 / plt.gca().get_data_ratio() * 1.0)
if mean:
    plt.savefig('Fig.mean.Cells.pdf',transparent=True)
else:
    plt.savefig('Fig.median.Cells.pdf', transparent=True)
plt.clf()

def plot_Cellular(variable, mean = False, color = None, title = None, xlabel = 'Time (hrs)', ylabel = None):
    V = np.zeros((ts, replicates))
    for r in range(1, replicates + 1):
        f = np.genfromtxt('Data/FullModelCellular_%i.txt' % (r), skip_header=1, delimiter=',', names=CellNames,
                          max_rows=ts)
        V[:, r - 1] = f[variable]
    epsilon = 1E-20
    if not mean:
        plt.plot(f['Time'],np.percentile(np.log10(V+epsilon), 50, axis =1),color=color,linewidth = 3.0,label='Ve')
        lower_bound = np.percentile(np.log10(V+epsilon), 5, axis = 1)
        upper_bound =  np.percentile(np.log10(V+epsilon), 95, axis = 1)
    else:
        plt.plot(f['Time'], np.mean(np.log10(V + epsilon), axis=1), color=color, linewidth=3.0, label='Ve')
        lower_bound = np.mean(np.log10(V + epsilon), axis=1) - np.std(np.log10(V + epsilon), axis=1)
        upper_bound = np.mean(np.log10(V + epsilon), axis=1) + np.std(np.log10(V + epsilon), axis=1)
    plt.fill_between(f['Time'], lower_bound, upper_bound, facecolor=color, interpolate=True, alpha=0.2)
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(r'%s' % ylabel)
    plt.gca().set_aspect(1.0 / plt.gca().get_data_ratio() * 1.0)
    plt.grid()
    if not mean:
        plt.savefig('Fig.median.%s.pdf' % variable, transparent=True)
    else:
        plt.savefig('Fig.mean.%s.pdf' % variable, transparent=True)
    plt.clf()

plot_Cellular('Ve',color = '#666699', title = 'Viral Load', ylabel = 'log$_{10}$ PFU/ml', mean = mean)
plot_Cellular('IFNe',color = '#993300', title = 'Extracellular IFN', ylabel = '[IFNe] $\mu$M', mean = mean)