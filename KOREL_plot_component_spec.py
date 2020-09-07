import numpy as np
import matplotlib.pyplot as plt

plt.figure(figsize=(7,5), dpi=300)
# plt.figure(dpi=300)
ax = plt.gca()

spe_file = 'KOREL/6678/korel.spe'

wav, spe1, spe2 = np.loadtxt(spe_file, unpack=True)
#wav, spe1 = np.loadtxt(spe_file, unpack=True)

ax.plot(wav, spe1, color='black', lw='0.5')
ax.plot(wav, spe2-0.025, color='black', lw='0.5')

ax.set_xlabel('$\lambda$ ($\AA$)')
ax.set_ylabel('normalized flux')

# ax.set_xlim(6658, 6698)
# ax.set_ylim(0.93, 1.02)


plt.savefig(spe_file+'.png', format='png')
#plt.savefig('KOREL/4713/korel_4713_spe.png', format='png')
#plt.savefig('KOREL/6678/korel_6678_spe.png', format='png')

#plt.show()
plt.close()
