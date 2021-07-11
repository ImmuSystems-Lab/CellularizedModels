import matplotlib.pyplot as plt
import numpy as np

vDC = np.arange(1,5)
Bifurcation = [7.5,12.5,22.5,32.5]

sampling_x = [1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,4,4,4,4,4,4,4,4]
sampling_y = [5,10,15,20,25,30,35,40,5,10,15,20,25,30,35,40,5,10,15,20,25,30,35,40,5,10,15,20,25,30,35,40,]
plt.plot(vDC,Bifurcation,color = 'black')
plt.fill_between(vDC, 0, Bifurcation, facecolor='#666699', interpolate=True)
plt.fill_between(vDC, Bifurcation, 40, facecolor='#993300', interpolate=True)
plt.scatter(sampling_x,sampling_y,marker='x',color = 'black')
plt.ylim(1,40)
plt.xlim(1,4)
plt.xlabel('Virus Diffusion Coefficient')
plt.ylabel('IFN Diffusion Coefficient')
plt.title('Phase Portrait \n Virus vs IFN Diffusion Coefficient')
plt.savefig('PhasePortrait.svg')
plt.show()