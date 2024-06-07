#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 31 12:09:09 2024

@author: sdavenport
"""
savedir = '/Users/sdavenport/Documents/Code/Python/Papers/2023_glmfdp/Simulations/Review_results/'
# pltsaveloc = '/Users/sdavenport/Documents/Code/Python/Papers/2023_glmfdp/Simulations/JER_control/Review_plots/'
pltsaveloc = '/Users/sdavenport/Documents/MyPapers/0papers/bootspatlm/Article/Figures/ReviewFigures/'

file = 'nsubjdim2D_fwhm_0_pi0_80_nsubj_70_B_100_niters_1000_simtype_-1_df_10_dim_25.npz'

data = np.load(savedir + file)

data['JER_FPR_SD']

# %%
greyscale = 0

if greyscale == 1:
    paracolor = 'grey'
    bootcolor = 'black'
    biomaddon = '_greyscale_' + str(greyscale)
else:
    paracolor = 'red'
    bootcolor = 'blue'
    biomaddon = ''

# set the font size
font = {'size': 18}
plt.rc('font', **font)

nsubj_vec = np.arange(20, 101, 10)
dim_sides = np.asarray((25))
# lendimsides = len(dim_sides)
lendimsides = 1

# Calculate the error bars (using the normal approximation to the binomial)
interval = pr.bernstd(0.1, 1000, 0.95)[0]
ones_vec = np.ones(len(nsubj_vec))

FWHM_vec = np.arange(0, 9, 4)
# alpha_vec = (0.4, 0.5, 0.6, 0.75, 0.9)
alpha_vec = (0.4, 0.6, 0.9)

lw = 3

for df in np.array((3, 10)):
    for pi0 in np.array((0.5, 0.8, 0.9, 1)):
        print(pi0)
        pi0invec = int(100*pi0)
        for j in np.arange(lendimsides):
            dim = dim_sides
            # max_eJER = 0
            for i in np.arange(len(FWHM_vec)):
                # Plot the line alpha = 0.1
                plt.plot(nsubj_vec, 0.1*ones_vec, 'k-')

                # Plot the error bars
                plt.plot(nsubj_vec, interval[0]*ones_vec, 'k--')
                plt.plot(nsubj_vec, interval[1]*ones_vec, 'k--')

                # Set the axis labels
                plt.xlabel('Number of Subjects')
                plt.ylabel('Empirical JER')

                # Initialize vectors to store the jer results
                observed_jer_simes = np.zeros((1, len(nsubj_vec)))
                observed_jer_ari = np.zeros((1, len(nsubj_vec)))
                observed_jer_boot = np.zeros((1, len(nsubj_vec)))
                observed_jer_boot_sd = np.zeros((1, len(nsubj_vec)))

                # Initialize vectors to store the fwer results
                observed_fwer_simes = np.zeros((1, len(nsubj_vec)))
                observed_fwer_ari = np.zeros((1, len(nsubj_vec)))
                observed_fwer_boot = np.zeros((1, len(nsubj_vec)))
                observed_fwer_boot_sd = np.zeros((1, len(nsubj_vec)))

                # Load the results for each subject
                for k in np.arange(len(nsubj_vec)):
                    # Load the parametric results
                    parametric_results = np.load(savedir + 'nsubjdim2D_fwhm_' + str(FWHM_vec[i]) + '_pi0_' + str(
                        pi0invec) + '_nsubj_' + str(nsubj_vec[k]) + '_B_100_niters_1000_simtype_-1_df_' + str(df) + '_dim_' + str(dim_sides) + '.npz')

                    # Load the bootstrap results
                    bootstrap_results = np.load(savedir + 'nsubjdim2D_fwhm_' + str(FWHM_vec[i]) + '_pi0_' + str(
                        pi0invec) + '_nsubj_' + str(nsubj_vec[k]) + '_B_100_niters_1000_simtype_1_df_' + str(df) + '_dim_' + str(dim_sides) + '.npz')

                    # Calculates the Simes JER
                    jer_fpr_simes = parametric_results['JER_FPR']
                    observed_jer_simes[0, k] = jer_fpr_simes

                    # Calculates the ARI JER
                    jer_fpr_ari = parametric_results['JER_FPR_SD']
                    observed_jer_ari[0, k] = jer_fpr_ari

                    # Calculates the bootstrap JER
                    jer_fpr_boot = bootstrap_results['JER_FPR']
                    observed_jer_boot[0, k] = jer_fpr_boot

                    # Calculates the bootstrap stepdown JER
                    jer_fpr_boot_sd = bootstrap_results['JER_FPR_SD']
                    observed_jer_boot_sd[0, k] = jer_fpr_boot_sd

                    # Calculates the Simes FWER
                    fwer_fpr_simes = parametric_results['FWER_FPR']
                    observed_fwer_simes[0, k] = fwer_fpr_simes

                    # Calculates the ARI FWER
                    fwer_fpr_ari = parametric_results['FWER_FPR_SD']
                    observed_fwer_ari[0, k] = fwer_fpr_ari

                    # Calculates the bootstrap FWER
                    fwer_fpr_boot = bootstrap_results['FWER_FPR']
                    observed_fwer_boot[0, k] = fwer_fpr_boot

                    # Calculates the bootstrap stepdown FWER
                    fwer_fpr_boot_sd = bootstrap_results['FWER_FPR_SD']
                    observed_fwer_boot_sd[0, k] = fwer_fpr_boot_sd

                # max_eJER = np.max((max_eJER, np.max(observed_jer_simes)))
                # max_eJER = np.max((max_eJER, np.max(observed_jer_ari)))
                # Plot the error bars

                plt.plot(nsubj_vec, observed_jer_simes[0], label='Simes',
                         color=paracolor, linestyle='dashed', linewidth=lw)
                plt.plot(
                    nsubj_vec, observed_jer_ari[0], label='ARI', color=paracolor, linewidth=lw)
                plt.plot(nsubj_vec, observed_jer_boot[0], label='Bootstrap',
                         color=bootcolor, linestyle='dashed', linewidth=lw)
                plt.plot(
                    nsubj_vec, observed_jer_boot_sd[0], label='Bootstrap Step Down', color=bootcolor, linewidth=lw)

                plt.plot(nsubj_vec, interval[0]*ones_vec, 'k--')
                plt.plot(nsubj_vec, interval[1]*ones_vec, 'k--')
                plt.yticks((0, 0.05, 0.1, 0.12))
                # plt.ylim(0, np.max((0.15, max_eJER + 0.01)))
                # if (pi0 == 1.0) and (FWHM_vec[i] == 0):
                if (pi0 == 1.0) and (FWHM_vec[i] == 0):
                    plt.legend(loc="lower right")

                plt.title(
                    "FWHM = " + str(FWHM_vec[i]) + " and $\pi_0$ = " + str(pi0))
                plt.subplots_adjust(left=0.17, right=0.95,
                                    bottom=0.17, top=0.9)
                # saveloc = 'C:\\Users\\12SDa\\global\\Intern\\drago\\Figures\\FinalEJerPlots\\'

                # saveloc = 'C:\\Users\\12SDa\\global\\TomsMiniProject\\Latex\\MyPapers\\FDP_control_via_the_bootstrap\\Current_Draft\\Figures\\EJerFinalPlots\\'
                plt.savefig(pltsaveloc + 'JER_dim_' + str(dim) + '_fwhm_' +
                            str(FWHM_vec[i]) + '_pi0_' + str(pi0invec) + str(biomaddon) + '_df_' + str(df) + '.pdf')
                plt.close()

                plt.plot(nsubj_vec, 0.1*ones_vec, 'k-')

                # Plot the error bars
                # plt.plot(nsubj_vec, interval[0]*ones_vec, 'k--')
                # plt.plot(nsubj_vec, interval[1]*ones_vec, 'k--')

                # Set the axis labels
                plt.xlabel('Number of Subjects')
                plt.ylabel('Empirical FWER')

                plt.plot(nsubj_vec, observed_fwer_simes[0], label='Simes',
                         color=paracolor, linestyle='dashed', linewidth=lw)
                plt.plot(
                    nsubj_vec, observed_fwer_ari[0], label='Holm', color=paracolor, linewidth=lw)
                plt.plot(nsubj_vec, observed_fwer_boot[0], label='Bootstrap',
                         color=bootcolor, linestyle='dashed', linewidth=lw)
                plt.plot(
                    nsubj_vec, observed_fwer_boot_sd[0], label='Bootstrap Step Down', color=bootcolor, linewidth=lw)
                plt.plot(nsubj_vec, interval[0]*ones_vec, 'k--')
                plt.plot(nsubj_vec, interval[1]*ones_vec, 'k--')

                if pi0 == 1.0:
                    plt.yticks((0, 0.05, 0.1, 0.15))
                elif pi0 == 0.9:
                    plt.yticks((0, 0.05, 0.1, 0.14))
                else:
                    plt.yticks((0, 0.05, 0.1, 0.12))

                # plt.ylim(0, np.max((0.15, max_eJER + 0.01)))
                if (pi0 == 1.0) and (FWHM_vec[i] == 0):
                    plt.legend(loc="lower right")

                plt.title(
                    "FWHM = " + str(FWHM_vec[i]) + " and $\pi_0$ = " + str(pi0))
                # saveloc = 'C:\\Users\\12SDa\\global\\Intern\\drago\\Figures\\FinalEJerPlots\\'
                plt.subplots_adjust(left=0.17, right=0.95,
                                    bottom=0.17, top=0.9)
                plt.savefig(pltsaveloc + 'FWER_dim_' + str(dim) + '_fwhm_' +
                            str(FWHM_vec[i]) + '_pi0_' + str(pi0invec) + str(biomaddon) + '_df_' + str(df) + '.pdf')
                plt.close()
