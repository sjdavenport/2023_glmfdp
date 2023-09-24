"""
Functions to calculate the empirical JER
"""
# Import statements
import numpy as np
import pyperm as pr


def run_nsubj_dim(run, B=100, simtype=1, niters=5000):
    # Set the location to save the results
    saveloc = '/storage/store2/work/sdavenpo/drago/Results/FPR_results/EJER_nsubj_dim/'

    # Set the fwhm
    fwhm = 2*(run % 5)

    # Initialize the contrast matrix
    if run % 2 == 0:
        # Single contrast
        C = np.array([[1, -1, 0]])
    else:
        # 2 contrasts to test
        C = np.array([[1, -1, 0], [0, 1, -1]])
        # Set the save file name

    nsubj_vec = np.arange(10, 101, 10)
    dim_sides = np.arange(20, 201, 10)

    saveloc += 'nsubjdim2D_' + str(C.shape[0]) + 'contrasts_fwhm_' + str(fwhm)
    saveloc += '_B_' + str(B) + '_niters_' + str(niters)

    # Initialize matrices to store the estimated FPRs
    store_JER = np.zeros((len(nsubj_vec), len(dim_sides)))
    store_FWER = np.zeros((len(nsubj_vec), len(dim_sides)))

    # Save the coverage
    for J in np.arange(len(dim_sides)):
        print('J:', J)
        Dim = (dim_sides[J], dim_sides[J])

        # In the single voxel case use this dimension
        if Dim == (1, 1):
            Dim = 1

        pi0 = 1
        for I in np.arange(len(nsubj_vec)):
            print('I:', I)
            FWER_FPR, JER_FPR = pr.bootfpr(
                Dim, nsubj_vec[I], C, fwhm, 0, B, niters, pi0, simtype)
            store_JER[I, J] = JER_FPR
            store_FWER[I, J] = FWER_FPR
            np.savez(saveloc + '.npz', JER_FPR=store_JER, FWER_FPR=store_FWER)


def run_nsubj_dim_mean(run, B=100, simtype=1, niters=5000):
    # Set the location to save the results
    saveloc = '/storage/store2/work/sdavenpo/drago/Results/FPR_results/EJER_nsubj_dim/'

    # Set the fwhm
    fwhm = 2*(run % 5)

    if run % 2 == 0:
        pi0 = 0.8
    else:
        pi0 = 0.9

    C = np.array([[1, -1, 0], [0, 1, -1]])

    nsubj_vec = np.arange(10, 101, 10)
    dim_sides = np.array((50, 100, 150))

    saveloc += 'meanpropdim2D' + '_fwhm_' + \
        str(fwhm) + '_pi0_' + str(int(100*pi0))
    saveloc += '_B_' + str(B) + '_niters_' + \
        str(niters) + '_simtype_' + str(simtype)

    # Initialize matrices to store the estimated FPRs
    store_JER = np.zeros((len(nsubj_vec), len(dim_sides)))
    store_FWER = np.zeros((len(nsubj_vec), len(dim_sides)))

    # Save the coverage
    for J in np.arange(len(dim_sides)):
        print('J:', J)
        Dim = (dim_sides[J], dim_sides[J])

        # In the single voxel case use this dimension
        if Dim == (1, 1):
            Dim = 1

        for I in np.arange(len(nsubj_vec)):
            print('I:', I)
            FWER_FPR, JER_FPR = pr.bootfpr(
                Dim, nsubj_vec[I], C, fwhm, 0, B, niters, pi0, simtype)
            store_JER[I, J] = JER_FPR
            store_FWER[I, J] = FWER_FPR
            np.savez(saveloc + '.npz', JER_FPR=store_JER, FWER_FPR=store_FWER)


def run_nsubj_dim_update(run, B=100, simtype=1, niters=5000):
    # Set the location to save the results
    saveloc = '/storage/store2/work/sdavenpo/drago/Results/FPR_results/EJER_nsubj_dim/'

    # Set the fwhm
    # fwhm_vec = np.arange(0,9,4)
    runmod3 = run % 3
    # fwhm = fwhm_vec[runmod3]
    if runmod3 == 2:
        fwhm = 4
    else:
        return 0

    # Initialize the contrast matrix
    C = np.array([[1, -1, 0], [0, 1, -1]])

    # Set the number of subjects
    nsubj = 10*(run % 11 + 1)
    if nsubj == 110:
        return 0  # Skip 110 only included to make things prime!

    pi0_vec = np.array((1, 0.9, 0.8, 0.5))
    pi0 = pi0_vec[run % 4]

    dim_sides = np.array((25, 50, 100))

    saveloc += 'nsubjdim2D_fwhm_' + \
        str(fwhm) + '_pi0_' + str(int(100*pi0)) + '_nsubj_' + str(nsubj)
    saveloc += '_B_' + str(B) + '_niters_' + \
        str(niters) + '_simtype_' + str(simtype)

    # Initialize matrices to store the estimated FPRs
    store_JER = np.zeros((1, len(dim_sides)))
    store_FWER = np.zeros((1, len(dim_sides)))
    store_JER_SD = np.zeros((1, len(dim_sides)))
    store_FWER_SD = np.zeros((1, len(dim_sides)))

    # Save the coverage
    for J in np.arange(len(dim_sides)):
        print('J:', J)
        Dim = (dim_sides[J], dim_sides[J])

        # In the single voxel case use this dimension
        if Dim == (1, 1):
            Dim = 1

        FWER_FPR, JER_FPR, FWER_FPR_SD, JER_FPR_SD = pr.bootfpr(
            Dim, nsubj, C, fwhm, 0, B, niters, pi0, simtype, do_sd=1)
        store_JER[0, J] = JER_FPR
        store_FWER[0, J] = FWER_FPR
        store_JER_SD[0, J] = JER_FPR_SD
        store_FWER_SD[0, J] = FWER_FPR_SD
        np.savez(saveloc + '.npz', JER_FPR=store_JER, FWER_FPR=store_FWER,
                 JER_FPR_SD=store_JER_SD, FWER_FPR_SD=store_FWER_SD)


def run_nsubj_dim_review(run, B=100, simtype=1, niters=5000):
    # run should be 1 to 60;

    # Set the location to save the results
    saveloc = 'C:\\Users\\12SDa\\davenpor\\davenpor\\Toolboxes\\Papers\\lmfdp\\Simulations\\JER_control\\EJER_runs_data\\Bootstrap\\'

    # Set the fwhm
    fwhm_vec = np.arange(0, 9, 4)
    runmod3 = run % 3
    fwhm = fwhm_vec[runmod3]

    # Initialize the contrast matrix
    C = np.array([[1, -1, 0], [0, 1, -1]])

    # Set the number of subjects
    nsubj = 100 + 10*(run % 5 + 1)

    pi0_vec = np.array((1, 0.9, 0.8, 0.5))
    pi0 = pi0_vec[run % 4]

    # dim_sides = np.array((25, 50, 100))
    dim_sides = np.array((50))

    saveloc += 'nsubjdim2D_fwhm_' + \
        str(fwhm) + '_pi0_' + str(int(100*pi0)) + '_nsubj_' + str(nsubj)
    saveloc += '_B_' + str(B) + '_niters_' + \
        str(niters) + '_simtype_' + str(simtype)

    # Initialize matrices to store the estimated FPRs
    store_JER = 0
    store_FWER = 0
    store_JER_SD = 0
    store_FWER_SD = 0

    # Save the coverage
    Dim = (50, 50)

    # In the single voxel case use this dimension
    if Dim == (1, 1):
        Dim = 1

    FWER_FPR, JER_FPR, FWER_FPR_SD, JER_FPR_SD = pr.bootfpr(
        Dim, nsubj, C, fwhm, 0, B, niters, pi0, simtype, do_sd=1)
    store_JER = JER_FPR
    store_FWER = FWER_FPR
    store_JER_SD = JER_FPR_SD
    store_FWER_SD = FWER_FPR_SD
    np.savez(saveloc + '.npz', JER_FPR=store_JER, FWER_FPR=store_FWER,
             JER_FPR_SD=store_JER_SD, FWER_FPR_SD=store_FWER_SD)
