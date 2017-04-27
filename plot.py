import numpy as np
import matplotlib.pyplot as plt

radios=np.genfromtxt('rad_dat.txt')
densidad=np.genfromtxt('dens_dat.txt')

rads =  np.unique(radios)
dens=np.zeros(len(rads))
i=0
for r in rads:
    a= np.where(radios==rads[i])
    dens[i]=np.mean(densidad[a])
    i+=1

fig=plt.figure()
plt.plot(rads,dens)
plt.xlabel(r'Radio al centro de la explosion ($r$)')
plt.ylabel(r'Densidad ($\rho$)')
plt.title('Densidad radial de la explosion')
plt.savefig('Densidad.pdf', fomat='pdf')
plt.close()
