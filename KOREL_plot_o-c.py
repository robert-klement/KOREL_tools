import numpy as np
from matplotlib import pyplot as plt
from PyAstronomy import pyasl
import sys
sys.path.append('/home/robert/Documents/python/python_tools')
from normalize_tools import normalize_line_linear

mjds = []
residual_green = 0
residual_red = 0
counter = 0
name = '6678'

# SPECTRA WITH RVs calculated from orbital solution
N, mjd_rv, vr1, oc1, vr2, oc2 = np.loadtxt('KOREL/'+name+'/korel.rv', unpack=True, skiprows=1)
vr1_C = vr1-oc1
vr2_C = vr2-oc2

spe_file = 'KOREL/'+name+'/korel.spe'
wav_spe, flux_spe1, flux_spe2 = np.loadtxt(spe_file, unpack=True)

#### plot input
with open('KOREL/'+name+'/korel.dat', 'r') as fdat:

    with open('KOREL/'+name+'/korel.o-c', 'r') as foc:

        fdat_contents = fdat.readline()

        nbins = int(fdat_contents[-5:-1])

        while len(fdat_contents) > 0:
            if int(fdat_contents[1]) >1:
                mjds.append(fdat_contents[1:12])
            fdat_contents = fdat.readline()

        fdat.seek(0)

        for m in range(len(mjds)):

            plt.close()

            wav = np.empty(nbins)
            flux = np.empty(nbins)

            fdat_contents = fdat.readline()

            wav0 = float(fdat_contents[13:22])
            rv_step = float(fdat_contents[24:29])
            weight = float(fdat_contents[31:37])

            for i in range(nbins):

                wav[i] = (wav0 + (i*rv_step/3e5)*wav0)

                if ((i+1)%10) == 0:
                    fdat.seek(fdat.tell()+1)
                flux[i] = fdat.read(8)

            fdat.seek(fdat.tell()+1)

            flux_oc = np.empty(nbins)

            foc_contents = foc.readline()

            for i in range(nbins):

                if ((i+1)%10) == 0:
                    foc.seek(foc.tell()+1)
                flux_oc[i] = foc.read(8)

            foc.seek(foc.tell()+1)

            # if wav0 == 6653:
            # if wav0 == 4688:

            # cont_wavs = (4700, 4720)
            cont_wavs = (6665, 6688)

            nflux_spe1, nwav_spe1 = pyasl.dopplerShift(wav_spe, flux_spe1, vr1_C[m], edgeHandling="firstlast")
            nflux_spe2, nwav_spe2 = pyasl.dopplerShift(wav_spe, flux_spe2, vr2_C[m], edgeHandling="firstlast")
            cont1 = normalize_line_linear(wav_spe, nflux_spe1, cont_wavs[0], cont_wavs[1])
            cont2 = normalize_line_linear(wav_spe, nflux_spe2, cont_wavs[0], cont_wavs[1])

            nflux_spe1_1 = nflux_spe1 / cont1(wav_spe)
            nflux_spe2_1 = nflux_spe2 / cont2(wav_spe)

            nflux_spe1, nwav_spe1 = pyasl.dopplerShift(wav_spe, flux_spe1, vr1[m], edgeHandling="firstlast")
            nflux_spe2, nwav_spe2 = pyasl.dopplerShift(wav_spe, flux_spe2, vr2[m], edgeHandling="firstlast")
            cont1 = normalize_line_linear(wav_spe, nflux_spe1, cont_wavs[0], cont_wavs[1])
            cont2 = normalize_line_linear(wav_spe, nflux_spe2, cont_wavs[0], cont_wavs[1])

            nflux_spe1_2 = nflux_spe1 / cont1(wav_spe)
            nflux_spe2_2 = nflux_spe2 / cont2(wav_spe)


            counter += 1

            plt.figure(figsize=(10,5), dpi=300)
            ax = plt.gca()
            #fig, (ax,ax2) = plt.subplots(nrows=1,ncols=2, figsize=(10,5), dpi=300)

            ax.plot(wav, flux, color='blue', lw='1.0', label='O', alpha=0.6)

            ax.plot(wav_spe, nflux_spe1_1-0.08, color='green', lw='1.07', alpha=0.5)
            ax.plot(wav_spe, nflux_spe2_1-0.08, color='green', lw='1.0', alpha=0.5)
            ax.plot(wav_spe, nflux_spe1_1 + nflux_spe2_1 -1, color='green', lw='1.0', alpha=0.8, label='C from orbital solution')
            ax.plot(wav_spe, nflux_spe1_2-0.08, color='red', lw='1.0', alpha=0.5)
            ax.plot(wav_spe, nflux_spe2_2-0.08, color='red', lw='1.0', alpha=0.5)
            ax.plot(wav_spe, nflux_spe1_2 + nflux_spe2_2 -1, color='red', lw='1.0', alpha=0.8, label='C from RVs')

            ax.text(wav[0]+1, 1.085, 'RV1 from orb. solution = {:.2f} km/s'.format(vr1[m]-oc1[m]), fontsize=8)
            ax.text(wav[0]+1, 1.07, 'RV2 from orb. solution = {:.2f} km/s'.format(vr2[m]-oc2[m]), fontsize=8)
            ax.text(wav[0]+1, 1.055, 'RV1 from korel RVs = {:.2f} km/s'.format(vr1[m]), fontsize=8)
            ax.text(wav[0]+1, 1.04, 'RV2 from korel RVs = {:.2f} km/s'.format(vr2[m]), fontsize=8)
            ax.text(wav[-1]-10, 1.07, '(O-C)1 = {:.2f} km/s'.format(oc1[m]), fontsize=8)
            ax.text(wav[-1]-10, 1.05, '(O-C)2 = {:.2f} km/s'.format(oc2[m]), fontsize=8)

            ax.scatter(wav, flux - (nflux_spe1_1 + nflux_spe2_1 -1) -- 0.85, color='green', marker='.', s=0.2)
            ax.scatter(wav, flux - (nflux_spe1_2 + nflux_spe2_2 -1) -- 0.77, color='red', marker='.', s=0.2)

            ax.set_title('mjds[m]')
            ax.text(wav0+1, 0.7, 'weight = {}'.format(weight), fontsize=8)
            ax.set_ylim(0.65, 1.1)

            ax.axhline(y=1.0, lw=0.4, ls=':', color='black')
            ax.axhline(y=0.92, lw=0.4, ls=':', color='black')
            ax.axhline(y=0.85, lw=0.4, ls=':', color='black')
            ax.axhline(y=0.77, lw=0.4, ls=':', color='black')
            ax.legend(loc='lower right', fontsize=8)

            plt.savefig('KOREL/{}/O-C_plots/{}.png'.format(name, mjds[m]), format='png')
            # plt.show()
            plt.close()


            residual_green += flux - (nflux_spe1_1 + nflux_spe2_1 -1)
            residual_red += flux - (nflux_spe1_2 + nflux_spe2_2 -1)

            if mjds[m] == mjds[-1]:

                plt.figure(figsize=(5,5), dpi=300)
                ax = plt.gca()

                print(counter)

                for i in range(len(residual_green)):

                    residual_green[i] = residual_green[i] / counter
                    residual_red[i] = residual_red[i] / counter


                ax.plot(wav, residual_green + 1.0, color='green', lw='0.5', label='residual')
                ax.plot(wav, residual_red + 0.95, color='red', lw='0.5', label='residual')
                ax.axhline(y=1.0, lw=0.4, ls=':', color='black')
                ax.axhline(y=0.95, lw=0.4, ls=':', color='black')

                ax.set_ylim(0.9, 1.05)

                # plt.legend(loc='lower right')
                plt.savefig('KOREL/'+name+'/O-C_plots/RESIDUAL_'+name+'.png', format='png')
                # plt.show()
                plt.close()

