"""
Creating the plots in the paper corresponding for HCP analyses
"""
import numpy as np
from nilearn import plotting
from nilearn.input_data import NiftiMasker
from nilearn.image import load_img
import matplotlib.pyplot as plt
import pyperm as pr

lm_type = 'iq_sex'
# 'iq_hand'
for method in ['parametric', 'boot']:

    if lm_type == 'bmi':
        # Select the variables for which to test contrasts
        contrast_vars = ['Height', 'Weight', 'BMI']
    elif lm_type == 'iq_hand':
        contrast_vars = ['Handedness', 'PMAT24_A_CR']
    elif lm_type == 'iq_age':
        contrast_vars = ['Age_in_Yrs', 'PMAT24_A_CR']
    elif lm_type == 'iq_sex':
        contrast_vars = ['Gender', 'PMAT24_A_CR']
    else:
        raise Exception("wrong lm_type input used")

    # Obtain the masker
    # the fwhm that was used to generate things (shouldn't affect things here I don't think)
    fwhm = 4
    mask = 'C:\\Users\\12SDa\\global\\Intern\\drago\\Images\\mask_GM_forFunc.nii'
    masker = NiftiMasker(smoothing_fwhm=fwhm, mask_img=mask, memory='./').fit()

    # Obtain the boolean mask
    mask = load_img(mask).get_fdata()
    mask = mask > 0

    # Plot the data
    image_dir = 'C:\\Users\\12SDa\\global\\Intern\\drago\\Images\\'

    if lm_type == 'bmi':
        saveloc = 'C:\\Users\\12SDa\\global\\Intern\\drago\\LinearModelFitting\\bmi_tdp_bounds_' + method
    else:
        saveloc = 'C:\\Users\\12SDa\\global\\Intern\\drago\\LinearModelFitting\\lm_run_store\\tdp_bounds' \
            + '_' + contrast_vars[0] + '_' + contrast_vars[1]

    if method == 'boot':
        saveloc += '_boot_B_100'
    else:
        saveloc += '_parametric'

    load_results = np.load(saveloc + '.npz')

    if method == 'boot':
        tdp_bounds_boot = load_results['tdp_bounds']
        tdp_bounds_boot_sd = load_results['tdp_bounds_sd']
    else:
        tdp_bounds_simes = load_results['tdp_bounds']
        tdp_bounds_ari = load_results['tdp_bounds_sd']

    #print(np.unique(tdp_bounds[mask, index]))

# %%
saveloc = 'C:\\Users\\12SDa\\global\\Intern\\drago\\Figures\\HCP\\'

for index in np.array((0, 1)):
    tdp_bounds_nifti = masker.inverse_transform(tdp_bounds_ari[mask, index])
    plotting.plot_stat_map(
        tdp_bounds_nifti,
        display_mode='z', vmax=1, colorbar=True,
        title='ARI TDP bounds', cut_coords=[-27, 48, 69])

    plt.savefig(saveloc + 'TDP_ari_' + lm_type + '_' + str(index) + '.pdf')

    tdp_bounds_nifti = masker.inverse_transform(tdp_bounds_simes[mask, index])
    plotting.plot_stat_map(
        tdp_bounds_nifti,
        display_mode='z', vmax=1, colorbar=True,
        title='Simes TDP bounds', cut_coords=[-27, 48, 69])

    plt.savefig(saveloc + 'TDP_simes_' + lm_type + '_' + str(index) + '.pdf')

    tdp_bounds_nifti = masker.inverse_transform(tdp_bounds_boot[mask, index])
    plotting.plot_stat_map(
        tdp_bounds_nifti,
        display_mode='z', vmax=1, colorbar=True,
        title='Bootstrap TDP bounds', cut_coords=[-27, 48, 69])

    plt.savefig(saveloc + 'TDP_boot_' + lm_type + '_' + str(index) + '.pdf')

# %%
cluster_image, c_sizes = pr.find_clusters(tdp_bounds_boot[..., 1], 0.01)
_, c_sizes_sex = pr.find_clusters(tdp_bounds_boot[..., 0], 0.01)
c_sizes_including_sex = np.append(c_sizes, c_sizes_sex[0])
# Get the sort indices in decreasing order
c_size_index = np.argsort(-c_sizes_including_sex)
c_sizes_including_sex[c_size_index]

# %% Find cluster indices
nclusters = np.max(cluster_image)
boot_cluster_tdps = np.zeros((1, nclusters+1))[0]
simes_cluster_tdps = np.zeros((1, nclusters+1))[0]
ari_cluster_tdps = np.zeros((1, nclusters+1))[0]
for i in np.arange(1, nclusters+1):
    sample_cluster_voxel_idx = np.where(np.ravel(cluster_image) == i)[0][0]
    ravelled_tdp_boot = np.ravel(tdp_bounds_boot[..., 1])
    boot_cluster_tdps[i-1] = ravelled_tdp_boot[sample_cluster_voxel_idx]
    ravelled_tdp_ari = np.ravel(tdp_bounds_ari[..., 1])
    ari_cluster_tdps[i-1] = ravelled_tdp_ari[sample_cluster_voxel_idx]
    ravelled_tdp_simes = np.ravel(tdp_bounds_simes[..., 1])
    simes_cluster_tdps[i-1] = ravelled_tdp_simes[sample_cluster_voxel_idx]

boot_cluster_tdps[nclusters] = np.unique(tdp_bounds_boot[..., 0])[1]
ari_cluster_tdps[nclusters] = np.unique(tdp_bounds_ari[..., 0])[1]
simes_cluster_tdps[nclusters] = np.unique(tdp_bounds_simes[..., 0])[1]

boot_cluster_tps = boot_cluster_tdps*c_sizes_including_sex
ari_cluster_tps = ari_cluster_tdps*c_sizes_including_sex
simes_cluster_tps = simes_cluster_tdps*c_sizes_including_sex

# %% Compare simes, ARI and bootstrap bounds
saveloc = 'C:\\Users\\12SDa\\global\\TomsMiniProject\\Latex\\MyPapers\\FDP_control_via_the_bootstrap\\Article\\Figures\\HCP\\'

for greyscale in np.arange(2):
    if greyscale == 1:
        paracolor = 'grey'
        bootcolor = 'black'
        biomaddon = '_biometrika_' + str(greyscale)
    else:
        paracolor = 'red'
        bootcolor = 'blue'
        biomaddon = ''

    font = {'size': 16}
    plt.rc('font', **font)
    nclusters = len(np.unique(tdp_bounds_boot)) - 1

    plt.plot(np.arange(nclusters),
             ari_cluster_tdps[c_size_index], '*', color=paracolor, label='ARI')
    plt.plot(np.arange(nclusters),
             simes_cluster_tdps[c_size_index], '+', color=paracolor, label='Simes')
    plt.plot(np.arange(nclusters),
             boot_cluster_tdps[c_size_index], '*', color=bootcolor, label='bootstrap')

    plt.xlabel('Suprathreshold Clusters', labelpad=14.0)
    plt.ylabel('TDP')
    plt.title('Lower bounds on the TDP per cluster')
    plt.xticks(np.arange(1, 14))
    plt.tick_params(axis='x', labelbottom=False)
    plt.legend()
    plt.subplots_adjust(left=0.15, right=0.95, bottom=0.15, top=0.9)
    plt.savefig(saveloc + 'TDP_cluster_plot_' + lm_type + biomaddon + '.pdf')
    plt.clf()

    plt.plot(np.arange(nclusters),
             ari_cluster_tps[c_size_index], '*', color=paracolor, label='ARI')
    plt.plot(np.arange(nclusters),
             simes_cluster_tps[c_size_index], '+', color=paracolor, label='Simes')
    plt.plot(np.arange(nclusters),
             boot_cluster_tps[c_size_index], '*', color=bootcolor, label='bootstrap')

    plt.xlabel('Suprathreshold Clusters', labelpad=14.0)
    plt.ylabel('True Positives (TP)')
    plt.title('Lower bounds on the TP per cluster')
    plt.xticks(np.arange(14))
    plt.tick_params(axis='x', labelbottom=False)
    plt.legend()
    plt.subplots_adjust(left=0.15, right=0.95, bottom=0.15, top=0.9)
    plt.savefig(saveloc + 'TP_cluster_plot_' + lm_type + biomaddon + '.pdf')
    plt.clf()

# %%
tdp_bounds_nifti = masker.inverse_transform(tdp_bounds_simes[mask, 0])
plotting.plot_stat_map(
    tdp_bounds_nifti,
    display_mode='z', vmax=1, colorbar=True,
    title='Cluster for the sex contrast', cut_coords=[24])

plt.savefig(saveloc + 'sex_cluster.pdf')

# %% Identify same and different clusters
bootsgreaterthansimes = 1 * \
    (tdp_bounds_boot[..., index] > 0) - (tdp_bounds_simes[..., index] > 0)
boot_extra_clusters = tdp_bounds_boot[np.array(
    bootsgreaterthansimes, dtype='bool'), index]
print(np.unique(boot_extra_clusters))

# %% Find and compare clusters
cluster_image, _ = pr.find_clusters(tdp_bounds_simes[..., index] > 0, 0.5)

simes_tdps = np.zeros(len(np.unique(cluster_image))-1)
boot_tdps = np.zeros(len(np.unique(cluster_image))-1)

for i in np.arange(len(np.unique(cluster_image))):
    if i > 0:
        x, y, z = np.where(cluster_image == i)

        simes_tdps[i-1] = tdp_bounds_simes[x[0], y[0], z[0], index]
        boot_tdps[i-1] = tdp_bounds_boot[x[0], y[0], z[0], index]

print(simes_tdps)
print(boot_tdps)
