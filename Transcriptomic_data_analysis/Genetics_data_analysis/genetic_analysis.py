"""
Perform the analysis on the transcriptomic dataset
"""
import numpy as np
import pandas as pd
import sanssouci as sa
import pyperm as pr
import matplotlib.pyplot as plt

Xpath = "Data/GSE42057_X.csv"
Xdata = pd.read_csv(Xpath)

Ypath = "Data/GSE42057_Y.csv"
Ydata = pd.read_csv(Ypath)

# %%
finalgold = Xdata.finalgold.values
knownsubjects = finalgold > -2 # Just obtains a vector of TRUEs
knownsubjects[48] = False  # As contains a NaN value!
knownsubject_idx = np.arange(len(Xdata))[knownsubjects]
Xdata_known = Xdata[knownsubjects]

nsubj = len(Xdata_known)

# Obtain the gene expression data
snp_data = Ydata.values[:, knownsubjects]
snp_data = pr.make_field(snp_data)

X = np.ones((135,8))
X[:,1:4] = Xdata_known.values[:, 1:4]
X[:,4] = Xdata_known.values[:, 5]
X[:,5:] = Xdata_known.values[:, 8:]

#np.hstack((Xdata_known.values[:, 1:3], Xdata_known.values[:, 5],
#            Xdata_known.values[:, 8:]))

X = np.array(X, dtype='float')

np.linalg.inv(X.T @ X)

C = np.zeros((1, X.shape[1]))
nparams = len(C)
C[0, 4] = 1
print(C)

# %% Calculate the p-values
design = X
lat_data = snp_data
#tstat_field, residuals = pr.contrast_tstats(snp_data, X, C, check_error=1)

tstat_field, residuals, Cbeta_field = pr.contrast_tstats_noerrorchecking(snp_data, X, C)

df = snp_data.fibersize - nparams
orig_pvalues = pr.tstat2pval( tstat_field, df )

# Save data as csv for exporting to R
np.savetxt("./Cbeta.csv", Cbeta_field.field, delimiter=",")
np.savetxt("./orig_pvalues.csv", orig_pvalues.field, delimiter=",")

#Volcano plot (the R version is nicer!)
plt.plot(Cbeta_field.field, -np.log10(orig_pvalues.field), '*')

# %% Plot a histogram of the pvalues
plt.hist(np.ravel(orig_pvalues.field))
alpha = 0.05
ARImestimate = pr.compute_hommel_value(np.ravel(orig_pvalues.field), alpha)

# %% Run the bootstrap algorithm
minp_perm, orig_pvalues_boot, pivotal_stats, bootstore = pr.boot_contrasts(snp_data, X, C, display_progress=1, store_boots = 1)
saveloc = './'

np.savez(saveloc + 'store_genetics_results.npz', minp_perm = minp_perm, orig_pvalues = orig_pvalues_boot.field, pivotal_stats = pivotal_stats, bootstore = bootstore)

# %%
saveloc = './'
stored_data = np.load(saveloc + 'store_genetics_results.npz')

pivotal_stats = stored_data['pivotal_stats']
pvals = stored_data['orig_pvalues']
boot_store = stored_data['bootstore']

alpha = 0.1
lambda_quant = np.quantile(pivotal_stats, alpha)

# Calculate the number of voxel-contrasts
m = np.prod(pvals.shape)

# Obtain the template threshold
thr_boot = sa.linear_template(lambda_quant, m, m)
lambda_quant_sd, _ = pr.step_down( boot_store, alpha = alpha, do_fwer = 0)
thr_boot_sd = sa.linear_template(lambda_quant_sd, m, m)

ARImestimate = pr.compute_hommel_value(np.ravel(pvals), alpha)
thr_ARI = sa.linear_template(alpha, ARImestimate, ARImestimate)
thr_simes = sa.linear_template(alpha, m, m)

# Calculate the bound on the number of false positives within the chosen set
pvalue_set = pvals[pvals < 0.001]
npvals = len(pvalue_set)

max_FP_bound = sa.max_fp(pvalue_set, thr_boot)
print(max_FP_bound)
print((npvals - max_FP_bound)/npvals)

max_FP_bound = sa.max_fp(pvalue_set, thr_ARI)
print(max_FP_bound)
print((npvals - max_FP_bound)/npvals)

# %% Compute and plot the lower bound on the number of false positives
max_FP = sa.curve_max_fp(np.sort(pvalue_set), thr_boot)
one2npvals = np.arange(1, npvals + 1)
plt.plot(max_FP)

# %%
orig_pvalues_sorted = np.sort(np.ravel(orig_pvalues))
pvalue_set = orig_pvalues_sorted[:100]
max_FP = sa.curve_max_fp(pvalue_set, thr_boot)

# %% Run BH on the pvalues and compute the bounds on the TDP and FDP
rejection_ind, n_rejections, rejection_locs = pr.fdr_bh( np.ravel(pvals), alpha = 0.05)
print(n_rejections)
print(pvals[rejection_locs])
max_FP_bound_boot = sa.max_fp(pvals[rejection_locs], thr_boot)
min_TP_bound_boot = len(pvals[rejection_locs]) - max_FP_bound_boot
print(min_TP_bound_boot)

max_FP_bound_boot_sd = sa.max_fp(pvals[rejection_locs], thr_boot)
min_TP_bound_boot_sd = len(pvals[rejection_locs]) - max_FP_bound_boot_sd 
print(min_TP_bound_boot_sd)

max_FP_bound_ARI = sa.max_fp(pvals[rejection_locs], thr_ARI)
min_TP_bound_ARI = len(pvals[rejection_locs]) - max_FP_bound_ARI
print(min_TP_bound_ARI)

max_FP_bound_simes = sa.max_fp(pvals[rejection_locs], thr_simes)
min_TP_bound_simes = len(pvals[rejection_locs]) - max_FP_bound_simes
print(min_TP_bound_simes)

# %% Including BMI as a contrast of interest as well!
C = np.zeros((2, X.shape[1]))
C[0, 4] = 1
C[1, 3] = 1
print(C)

# %% Calculate the p-values
design = X
lat_data = snp_data
tstat_field, residuals = pr.contrast_tstats(snp_data, X, C, check_error=1)

df = snp_data.fibersize - C.shape[1]
orig_pvalues = pr.tstat2pval( tstat_field, df )

# %% Run BH on the pvalues
rejection_ind, n_rejections, rejection_locs = pr.fdr_bh( np.ravel(orig_pvalues.field), alpha = 0.05)
print(n_rejections)

# %% Run the bootstrap algorithm
minp_perm, orig_pvalues, pivotal_stats, bootstore = pr.boot_contrasts(snp_data, X, C, display_progress=1)
saveloc = './'

np.savez(saveloc + 'store_genetics_results_fevandbmi.npz', minp_perm = minp_perm, orig_pvalues = orig_pvalues.field, pivotal_stats = pivotal_stats, bootstore = bootstore)

# %%
saveloc = './'
stored_data = np.load(saveloc + 'store_genetics_results_fevandbmi.npz')

pivotal_stats = stored_data['pivotal_stats']
orig_pvalues = stored_data['orig_pvalues']

alpha = 0.1
lambda_quant = np.quantile(pivotal_stats, alpha)

# Calculate the number of voxel-contrasts
m = np.prod(orig_pvalues.shape)

# Obtain the template threshold
thr_boot = sa.linear_template(lambda_quant, m, m)
thr_ARI = sa.linear_template(alpha, m, m)

# Calculate the bound on the number of false positives within the chosen set
pvalue_set = orig_pvalues[orig_pvalues < 0.001]
npvals = len(pvalue_set)

max_FP_bound = sa.max_fp(pvalue_set, thr_boot)
print(max_FP_bound)
print((npvals - max_FP_bound)/npvals)

max_FP_bound = sa.max_fp(pvalue_set, thr_ARI)
print(max_FP_bound)
print((npvals - max_FP_bound)/npvals)

# Calculate the bound on the number of false positives within the chosen set
pvalue_set = orig_pvalues[orig_pvalues < 0.01]
max_FP_bound = sa.max_fp(pvalue_set, thr_boot)
print(max_FP_bound)

max_FP_bound = sa.max_fp(pvalue_set, thr_ARI)
print(max_FP_bound)
