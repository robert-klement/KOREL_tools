import numpy as np
from glob import glob
import os
import sys
import matplotlib.pyplot as plt
import astropy.io.fits as fits
import astropy.units as unit
import astropy.constants as const
sys.path.append('/home/robert/Documents/python/python_tools')
from spec_tools import *
from normalize_tools import *

normalize = True

region = '6678'
nbins = 2048

fits_files = sorted(glob('spectra/ondrejov_COUDE700/ha/*.fit'))
files = sorted(glob('spectra/ondrejov_COUDE700/ha/*.dat'))

print('number of files: ', len(files))

datfile = 'KOREL/' + region + '/korel.dat_ondrejov_COUDE700_'+ str(nbins)

if os.path.isfile(datfile):
    os.remove(datfile)

SNRs = []
MJDs = []

for i, fits_file in enumerate(fits_files):

    wav, flux, MJD, dateobs, datereduc, fitsfile = loadfits(fits_file)
    MJDs.append(MJD)
    #print(MJD, wav[0], wav[-1], wav[0]/(wav[1]-wav[0]))

    SNRs.append(get_spec_snr(wav, flux, (6658, 6698)))

print('maximum SNR: ', np.amax(SNRs))

for i, file in enumerate(files):

    wav, flux = np.loadtxt(file, unpack=True)

    flux = flux[(wav>6653) & (wav<6703)]
    wav = wav[(wav>6653) & (wav<6703)]

    print(SNRs[i])

    if normalize:

        cont = normalize_line_polynom(wav, flux, 6668, 6688, polynom=1)
        flux = flux / cont(wav)

    plt.plot(wav, flux, alpha=0.5, color='black', lw=0.5)
    plt.text(6655, 1.2, 'SNR = {:.2f}'.format(SNRs[i]))
    plt.axhline(y=1, lw=0.5, color='grey', ls=':')
    # ax.set_xlim(3995, 4045)
    # ax.set_xlim(4450, 4498)
    plt.ylim(0.5, 1.3)
    plt.savefig(file+'.png', format='png', dpi=300)

    plt.close()

    # plt.text(wav[-1]-10, 0.97, 'SNR = {:.2f}'.format(SNRs[i]))
    weight = (SNRs[i]/854.4600241890643)**2
    # plt.text(wav[-1]-10, 0.93, 'weight = {:.2f}'.format(weight))

    print('datapoints in region: ', len(wav))

    wav_new = np.linspace(6653., 6702.5, num=nbins)

    x = 3e5 * np.log(wav_new / wav_new[0])
    delta_x_min = x[-1] - x[-2]
    print('smallest RV step: ', delta_x_min)
    rv_step = delta_x_min

    equidistant_log_scale = np.empty(nbins)
    for j in range(len(equidistant_log_scale)):
        equidistant_log_scale[j] = x[0] + j * rv_step

    wav_equi = np.exp(equidistant_log_scale / 3e5) * wav_new[0]
    flux_equi = np.interp(wav_equi, wav, flux)

    # for i in range(1, len(wav_equi)):
    #     print(3e5 * (np.log(wav_equi[i]) - np.log(wav_equi[i-1])))

    wav, flux = wav_equi, flux_equi

    flux[0], flux[-1] = 1.0, 1.0

    with open(datfile, 'a') as koreldat:
        if i != 0:
            koreldat.write('\n')
        # koreldat.write('{:12.5f}{:10.4f}{:7.3f}  {:5.3f}     {}\n'.format(MJD, wav[0], rv_step, weight, nbins))
        koreldat.write('{:12.5f}{:10.4f}{:7.3f}  {:5.3f}     {}\n'.format(
            MJDs[i], wav[0], rv_step, weight, nbins))
        for j in range(len(flux)):
            if j != 0:
                if j % 10 == 0:
                    koreldat.write('\n')
            koreldat.write(' {:7.5f}'.format(flux[j]))

print(np.amax(SNRs))
