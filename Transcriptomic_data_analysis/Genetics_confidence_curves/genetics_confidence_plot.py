"""
Calculate and plot lower bounds on the true postives and TDP for all sets of  
the form {p_(i), i < k} where p_(i) is the ith smallest pvalue and 1 <= k <= m
where m is the number of genes.
"""
# %% Import necessary packages
import numpy as np
import sanssouci as sa
import pyperm as pr

import matplotlib.pyplot as plt

# %%
resultsloc = '../Genetics_data_analysis/'
stored_data = np.load(resultsloc + 'store_genetics_results.npz')

pivotal_stats = stored_data['pivotal_stats']
orig_pvalues = stored_data['orig_pvalues']
boot_store = stored_data['bootstore']

# %%
# Set the desired alpha level
alpha = 0.1

# Calculate the lambda quantile
lambda_quant = np.quantile(pivotal_stats, alpha)

# Calculate the number of voxel-contrasts (only one contrast in this case)
m = len(orig_pvalues)

# Obtain the bootstrap template threshold
thr_boot = sa.linear_template(lambda_quant, m, m)
lambda_quant_sd, _ = pr.step_down( boot_store, alpha = alpha, do_fwer = 0)
thr_boot_sd = sa.linear_template(lambda_quant_sd, m, m)

ARImestimate = pr.compute_hommel_value(np.ravel(orig_pvalues), alpha)
thr_ARI = sa.linear_template(alpha, ARImestimate, ARImestimate)
thr_simes = sa.linear_template(alpha, m, m)

# Plot the FDP and TDP curves
saveloc = 'C:\\Users\\12SDa\\davenpor\\davenpor\\Toolboxes\\Papers\\lmfdp\\'+ \
    'Transcriptomic_data_analysis\\Genetics_confidence_curves\\fdp_curves'

font = {'size' : 15}
plt.rc('font', **font)

for greyscale in np.array((0,1)):
    if greyscale == 1:
        biomaddon = '_greyscale_' + str(greyscale)
        saveloc += biomaddon
    for number2plot in np.array((1000,2300,3000,3500,5000,m)):
        if greyscale == 1:
            paracolor = 'grey'
            bootcolor = 'black'
        else:
            paracolor = 'red'
            bootcolor = 'blue'
                
        if number2plot == 3000:
            dolegend = 1
        else:
            dolegend = 0
            
        pr.fdp_plot(orig_pvalues, [thr_boot, thr_boot_sd, thr_simes, thr_ARI], 
                        ['bootstrap', 'bootstrap step down', 'Simes', 'ARI'], 
                        [bootcolor, bootcolor, paracolor, paracolor], 
                        ['dashed', 'solid', 'dashed', 'solid'],
                        number2plot = number2plot, saveloc = saveloc,
                        dolegend = dolegend, vertline = 1745)

# %%
orig_pvalues_sorted = np.sort(np.ravel(orig_pvalues))
print(orig_pvalues_sorted[2300])