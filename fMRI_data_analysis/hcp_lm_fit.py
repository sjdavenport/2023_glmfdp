"""
A file that fits linear models to the HCP data using the covariates,
see the readme file for information on how to acquire this data which is
freely available online.
"""
from hcp_inference import lm_setup, obtain_hcp_data
import numpy as np
import pyperm as pr
from nilearn import plotting
from nilearn.image import load_img
import matplotlib.pyplot as plt

import sys
sys.path.insert(0, '/storage/store2/work/sdavenpo/drago/LinearModelFitting/')

df_unrelated, df_ur_unrelated, imgs = obtain_hcp_data()

# Select the variables to use in the regression
restricted_vars2use = ['Age_in_Yrs', 'Height', 'Weight',
                       'BMI', 'BPSystolic', 'BPDiastolic', 'Handedness']
unrestricted_vars2use = ['Gender', 'PMAT24_A_CR']

# Select the variables for which to test contrasts
contrast_vars = ['Gender', 'PMAT24_A_CR']
pair = contrast_vars[0] + '_' + contrast_vars[1]

# Calculate the number of contrast variables
n_contrasts = len(contrast_vars)

# Obtain the images
nsubj = len(imgs)
imgs = imgs[0:nsubj]

# Obtain the design and contrast matrices
design, contrast_matrix = lm_setup(
    restricted_vars2use, unrestricted_vars2use, contrast_vars, df_unrelated, df_ur_unrelated, nsubj)

nan_entries = np.argwhere(np.isnan(design))
nan_subjects = np.unique(nan_entries[:, 0])
nonnan_subjects = np.setdiff1d(np.arange(nsubj), nan_subjects)

# Select the images from the subjects that have all the required data
imgs = np.ndarray.tolist(np.asarray(imgs)[nonnan_subjects])

# Restrict the design to these subjects
design = design[nonnan_subjects, :]

mask = '/storage/store2/work/sdavenpo/drago/Images/mask_GM_forFunc.nii'

n_bootstraps = 100
for simtype in np.array((-1, 1)):
    tdp_bounds, tdp_bounds_sd, lambda_quant, lambda_quant_sd, masker = pr.cluster_tdp_brain(
        imgs, design, contrast_matrix, mask, n_bootstraps=n_bootstraps, min_cluster_size=1, simtype=simtype)

    if simtype == 1:
        method = 'boot'
    else:
        method = 'parametric'

    saveloc = '/storage/store2/work/sdavenpo/drago/LinearModelFitting/tdp_bounds_' + \
        pair + '_' + method

    if method == 'boot':
        saveloc = saveloc + '_B_' + str(n_bootstraps)

    np.savez(saveloc + '.npz', tdp_bounds=tdp_bounds, tdp_bounds_sd=tdp_bounds_sd,
             lambda_quant=lambda_quant, lambda_quant_sd=lambda_quant_sd)

    # Obtain the boolean mask
    mask_bool = load_img(mask).get_fdata()
    mask_bool = mask_bool > 0

    # Plot the data
    image_dir = '/storage/store2/work/sdavenpo/drago/Images/'

    for L in np.arange(n_contrasts):
        tdp_bounds_nifti = masker.inverse_transform(tdp_bounds[mask_bool, L])
        plotting.plot_stat_map(
            tdp_bounds_nifti,
            display_mode='z', vmax=1, colorbar=True,
            title='Lower bounds on TDP', cut_coords=30)
        if method == 'boot':
            plt.savefig(image_dir + contrast_vars[L] + '_TDP_' +
                        pair + '_' + method + '_B_' + str(n_bootstraps) + '.pdf')
        else:
            plt.savefig(
                image_dir + contrast_vars[L] + '_TDP_' + pair + '_' + method + '.pdf')

#for L in np.arange(n_contrasts): tdp_bounds_nifti = masker.inverse_transform(tdp_bounds[mask, L]); plotting.plot_stat_map(tdp_bounds_nifti,display_mode='z', vmax=1, colorbar=True,title='Lower bounds on TDP', cut_coords = 30); plt.savefig(image_dir + contrast_vars[L] + '_TDP'+method+'.pdf')
