import numpy as np
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 15})

names = ['Time','U','I1','I2','D','Ve','IFNe']
t = 480
replicate = 20
parameter = 'beta'
# p = [0.01,0.1,1.0,10.0,100.0]
p = [0.01, 0.08, 0.1, 0.13, 0.16, 0.20, 0.25, 0.32, 0.40, 0.50, 0.63, 0.79, 1.0, 1.26, 10.0, 100.0]
dt = 1.666667e-01

VeAUC = np.zeros(len(p))
VeAUCMin = np.zeros(len(p))
VeAUCMax = np.zeros(len(p))
VeMax = np.zeros(len(p))
VeMaxMin = np.zeros(len(p))
VeMaxMax = np.zeros(len(p))
IFNeAUC = np.zeros(len(p))
IFNeAUCMin = np.zeros(len(p))
IFNeAUCMax = np.zeros(len(p))
IFNeMax = np.zeros(len(p))
IFNeMaxMin = np.zeros(len(p))
IFNeMaxMax = np.zeros(len(p))
for i in range(len(p)):
    TVe = np.zeros((replicate))
    TVeMax = np.zeros((replicate))
    TIFNe = np.zeros((replicate))
    TIFNeMax = np.zeros((replicate))
    for r in range(1,replicate+1):
        f = np.genfromtxt('Data/FullModelCellular_%s_%.2f_%i.txt' % (parameter,p[i],r), skip_header=1, delimiter=',', names=names,max_rows=t)
        TVe[r-1] = np.trapz(f['Ve'], dx = dt)
        TVeMax[r - 1] = np.amax(f['Ve'])
        TIFNe[r - 1] = np.trapz(f['IFNe'],dx = dt)
        TIFNeMax[r - 1] = np.amax(f['IFNe'])
    VeAUC[i] = np.median(TVe)
    VeAUCMin[i] = np.percentile(TVe,10)
    VeAUCMax[i] = np.percentile(TVe,90)
    IFNeAUC[i] = np.median(TIFNe)
    IFNeAUCMin[i] = np.percentile(TIFNe,10)
    IFNeAUCMax[i] = np.percentile(TIFNe,90)

plt.plot(p,VeAUC,color='#666699',linewidth = 4.0, label = 'Virus AUC')
plt.fill_between(p, VeAUCMin,VeAUCMax, facecolor='#666699', interpolate=True, alpha=0.2)
plt.title('Virus AUC')
plt.ylabel(r'PFU hr/ml')
plt.xlabel('Parameter Multiplier')
plt.xscale('log')
plt.yscale('log')
plt.grid()
plt.gca().set_aspect(1.0 / plt.gca().get_data_ratio() * 1.0)
plt.gca().set_xticks([0.01,0.1,1.0,10.0,100.0])
plt.gca().set_yticks([2E9, 3E9,4E9,6E9,1E10])
# plt.legend(loc=2)
plt.savefig('Fig.Virus_AUC_%s.pdf' % parameter,transparent=True)
plt.clf()

plt.plot(p,IFNeAUC,color='#993300',linewidth = 4.0, label = 'IFN AUC')
plt.fill_between(p, IFNeAUCMin,IFNeAUCMax, facecolor='#993300', interpolate=True, alpha=0.2)
plt.title('IFN AUC')
plt.ylabel(r'$\mu$M hr')
plt.xlabel('Parameter Multiplier')
plt.xscale('log')
plt.yscale('log')
plt.grid()
plt.gca().set_aspect(1.0 / plt.gca().get_data_ratio() * 1.0)
plt.gca().set_xticks([0.01,0.1,1.0,10.0,100.0])
# plt.legend(loc=2)
plt.savefig('Fig.IFNe_AUC_%s.pdf' % parameter,transparent=True)
plt.clf()