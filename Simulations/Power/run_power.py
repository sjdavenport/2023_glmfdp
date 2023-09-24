# %%  Scripts to run the power calculations
import numpy as np
import pyperm as pr


def run_power_sims(run, nbootstraps=1000, simtype=1, niters=5000, FWHM=4):
    # simtype = -1 (simes/ARI), 1 (boot/boostrap stepdown)
    # import sys
    # sys.path.insert(0, '/storage/store2/work/sdavenpo/drago/Results/Power_Results/')
    # from run_power import run_power_sims
    # run_power_sims(2, nbootstraps = 1, simtype = -2, niters = 10, FWHM = 4)

    C = np.array([[1, -1, 0], [0, 1, -1]])
    saveloc = '/storage/store2/work/sdavenpo/drago/Results/Power_Results/'

    # Set the dimension of the image
    runmod2 = run % 2
    if runmod2 == 0:
        dimside = 50
    else:
        dimside = 100

    dim = (dimside, dimside)

    # Set the proportion of true null hypotheses
    runmod3 = run % 3
    if runmod3 == 0:
        pi0 = 0.5
    elif runmod3 == 1:
        pi0 = 0.8
    else:
        pi0 = 0.9

    # Calculate and save the power
    runmod5 = run % 5
    nsubj = (runmod5+1)*10
    power, power_sd = pr.bootpower(
        dim, nsubj, C, FWHM, 0, nbootstraps, niters, pi0, simtype)
    if simtype > -1:
        np.savez(saveloc + 'nsubj_' + str(nsubj) + '_pi0_' + str(int(100*pi0)) + '_nbootstraps_' + str(nbootstraps) + '_niters_' +
                 str(niters) + '_dim_' + str(dimside) + str(dimside) + '_simtype_' + str(simtype), power=power, power_sd=power_sd)
    else:
        np.savez(saveloc + 'nsubj_' + str(nsubj) + '_pi0_' + str(int(100*pi0)) + '_niters_' + str(niters) +
                 '_dim_' + str(dimside) + str(dimside) + '_simtype_' + str(simtype), power=power, power_sd=power_sd)

    nsubj = 50 + (runmod5+1)*10
    power, power_sd = pr.bootpower(
        dim, nsubj, C, FWHM, 0, nbootstraps, niters, pi0, simtype)
    if simtype > -1:
        np.savez(saveloc + 'nsubj_' + str(nsubj) + '_pi0_' + str(int(100*pi0)) + '_nbootstraps_' + str(nbootstraps) + '_niters_' +
                 str(niters) + '_dim_' + str(dimside) + str(dimside) + '_simtype_' + str(simtype), power=power, power_sd=power_sd)
    else:
        np.savez(saveloc + 'nsubj_' + str(nsubj) + '_pi0_' + str(int(100*pi0)) + '_niters_' + str(niters) +
                 '_dim_' + str(dimside) + str(dimside) + '_simtype_' + str(simtype), power=power, power_sd=power_sd)


def run_power_boot_sims(nbootstraps, niters=5000, FWHM=4):
    C = np.array([[1, -1, 0], [0, 1, -1]])
    saveloc = '/storage/store2/work/sdavenpo/drago/Results/Power_Results/Variable_bootstrap_sims/'

    dimside = 50
    dim = (dimside, dimside)
    pi0 = 0.8
    nsubj = 50

    # Calculate and save the power
    power, power_sd = pr.bootpower(
        dim, nsubj, C, FWHM, 0, nbootstraps, niters, pi0, 1)
    np.savez(saveloc + 'nsubj_' + str(nsubj) + '_pi0_' + str(int(100*pi0)) + '_nbootstraps_' +
             str(nbootstraps) + '_niters_' + str(niters) +
             '_dim_' + str(dimside) + str(dimside),
             power=power, power_sd=power_sd)


def run_power_sims_unrunsims(run, nbootstraps=100, niters=5000, FWHM=4):
    C = np.array([[1, -1, 0], [0, 1, -1]])
    saveloc = '/storage/store2/work/sdavenpo/drago/Results/Power_Results/'

    # Set the dimension of the image
    if run < 6:
        dimside = 50
        nsubj_vec = np.array((20, 70, 50, 100, 30, 80))
        pi0_vec = np.array((0.5, 0.5, 0.8, 0.8, 0.9, 0.9))
        nsubj = nsubj_vec[run]
        pi0 = pi0_vec[run]
    else:
        run = run - 6
        dimside = 100
        nsubj_vec = np.array((40, 90, 100, 20, 70, 10, 60))
        pi0_vec = np.array((0.5, 0.5, 0.5, 0.8, 0.8, 0.9, 0.9))
        nsubj = nsubj_vec[run]
        pi0 = pi0_vec[run]

    dim = (dimside, dimside)

    # Calculate and save the power
    simtype = 1  # Only missed the bootstraps

    power = pr.bootpower(dim, nsubj, C, FWHM, 0,
                         nbootstraps, niters, pi0, simtype)
    np.savez(saveloc + 'nsubj_' + str(nsubj) + '_pi0_' + str(int(100*pi0)) + '_nbootstraps_' + str(nbootstraps) +
             '_niters_' + str(niters) + '_dim_' + str(dimside) + str(dimside) + '_simtype_' + str(simtype), power=power)

# %%
#import numpy as np
#import pyrft as pr
#contrast_matrix = np.array([[1,-1,0],[0,1,-1]]); nbootstraps = 100; run = 0; nsubj_vec = np.array((40, 90, 100, 20, 70, 10, 60)); pi0_vec = np.array((0.5, 0.5, 0.5, 0.8, 0.8, 0.9, 0.9)); nsubj = nsubj_vec[run]; pi0 = pi0_vec[run]; dim = (50,50); simtype = 1; niters = 1; fwhm = 4;
#power = pr.bootpower(dim, nsubj, contrast_matrix, fwhm, 0, nbootstraps, niters, pi0, simtype)
