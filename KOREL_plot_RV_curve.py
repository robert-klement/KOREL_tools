import numpy as np
import matplotlib.pyplot as plt
from PyAstronomy import pyasl
import radvel.kepler

plt.figure(figsize=(7,5), dpi=300)
ax = plt.gca()

# ROBERT ORBIT
P = 53.762 # days
t0 = 51006.25
ecc = 0.076
small_omega = 147.5
K1 = 48.17
K2 = 52.44
gamma = 0

t = np.linspace(0, P, 100)
orbel1 = [P, 0, ecc, np.radians(small_omega), K1]
orbel2 = [P, 0, ecc, np.radians(small_omega), K2]

rv1 = radvel.kepler.rv_drive(t, orbel1)
rv2 = radvel.kepler.rv_drive(t, orbel2)

phase = t/P

ax.plot(phase, rv1, color='red')
ax.plot(phase, -rv2, color='red', alpha=0.5)


# rv_file = 'KOREL/HeI4009_HeI4026/korel.rv'
rv_file = 'KOREL/6678/korel.rv_6678_2048_merged'

#N, t, vr1, oc1 = np.loadtxt(rv_file, unpack=True, skiprows=1)
N, t, vr1, oc1, vr2, oc2 = np.loadtxt(rv_file, unpack=True, skiprows=1)
# N, t, vr1, oc1, vr2, oc2, vr3, oc3 = np.loadtxt(rv_file, unpack=True, skiprows=1)

# with open('RVs_korel_C.dat', 'a') as f:
#     for i in range(len(N)):
#         f.write('t, ')

phases = pyasl.foldAt(t, P, T0=t0, getEpoch=False)

ax.axhline(y=0, ls=':', color='black')
ax.scatter(phases, vr1, marker='.', color='black')
# ax.scatter(phases, vr1-oc1, marker='.', color='red')
ax.scatter(phases, vr2, marker='.', color='grey')
# ax.scatter(phases, vr2-oc2, marker='.', color='red', alpha=0.5)

# plt.scatter(phases, vr3, marker='+', color='blue')
# plt.scatter(phases, vr3-oc3, marker='+', color='green')

ax.set_xlim(-0.1, 1.1)
ax.set_ylim(-60, 60)
ax.set_xlabel('phase')
ax.set_ylabel('RV (km/s)')

plt.savefig(rv_file+'.png', format='png')
#plt.savefig('KOREL/4713/korel_4713_rv.png', format='png')
#plt.savefig('KOREL/6678/korel_6678_rv.png', format='png')


#plt.show()
plt.close()
