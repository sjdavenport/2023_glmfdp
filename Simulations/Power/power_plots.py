"""
Code to plot the power simulations figures from the paper
"""
import numpy as np
import matplotlib.pyplot as plt

font = {'size': 18}

plt.rc('font', **font)

para_results_dir = '.\\Power_results\\Parametric_Power\\'
boot_result_dir = '.\\Power\\Power_results\\Bootstrap_Power\\'

saveloc = '.\\Power\\Power_plots\\'

# %% Test load
nsubj = 100
pi0 = 0.5
dim = 50
boot = np.load(boot_result_dir + 'nsubj_' + str(nsubj) + '_pi0_' + str(int(100*pi0)) +
               '_nbootstraps_1000_niters_5000_dim_' + str(dim) + str(dim) + '_simtype_1.npz')
print(boot['power'])
print(boot['power_sd'])

ARI = np.load(para_results_dir + 'nsubj_' + str(nsubj) + '_pi0_' + str(int(100*pi0)
                                                                       ) + '_niters_5000_dim_' + str(dim) + str(dim) + '_simtype_-1.npz')
ARIpower = ARI['power']

simes = np.load(para_results_dir + 'nsubj_' + str(nsubj) + '_pi0_' +
                str(int(100*pi0)) + '_niters_5000_dim_' + str(dim) + str(dim) + '_simtype_-2.npz')
simespower = simes['power']

print(simespower)
print(ARIpower)
# %%
nsubj_vec = np.arange(10, 101, 10)

greyscale = 0

if greyscale == 1:
    paracolor = 'grey'
    bootcolor = 'black'
    biomaddon = '_greyscale_' + str(greyscale)
else:
    paracolor = 'red'
    bootcolor = 'blue'
    biomaddon = ''

for powertype in np.arange(5):
    print(powertype)
    for pi0 in np.array((0.5, 0.8, 0.9)):
        for dim in np.array((50,)):
            ARI_power_store = np.zeros(10)
            simes_power_store = np.zeros(10)
            boot_power_store = np.zeros(10)
            boot_power_store_sd = np.zeros(10)
            for I in np.arange(10):
                nsubj = 10*(I+1)
                try:
                    ARI = np.load(para_results_dir + 'nsubj_' + str(nsubj) + '_pi0_' + str(
                        int(100*pi0)) + '_niters_5000_dim_' + str(dim) + str(dim) + '_simtype_-1.npz')
                    ARIpower = ARI['power']
                    ARI_power_store[I] = ARIpower[powertype]
                except:
                    a = 1

                try:
                    simes = np.load(para_results_dir + 'nsubj_' + str(nsubj) + '_pi0_' + str(
                        int(100*pi0)) + '_niters_5000_dim_' + str(dim) + str(dim) + '_simtype_-2.npz')
                    simespower = simes['power']
                    simes_power_store[I] = simespower[powertype]
                except:
                    a = 1

                try:
                    boot = np.load(boot_result_dir + 'nsubj_' + str(nsubj) + '_pi0_' + str(int(100*pi0)) +
                                   '_nbootstraps_1000_niters_5000_dim_' + str(dim) + str(dim) + '_simtype_1.npz')
                    bootpower = boot['power']
                    bootpower_sd = boot['power_sd']
                    if nsubj < 60:
                        boot_power_store[I] = bootpower[powertype]
                        boot_power_store_sd[I] = bootpower_sd[powertype]
                    else:
                        boot_power_store[I] = bootpower[0, powertype]
                        boot_power_store_sd[I] = bootpower[1, powertype]
                except:
                    a = 1
                    print(a)

                plt.plot(nsubj_vec, simes_power_store, label='Simes',
                         color=paracolor, linestyle='dashed')
                plt.plot(nsubj_vec, ARI_power_store,
                         label='ARI', color=paracolor)
                plt.plot(nsubj_vec, boot_power_store, label='Bootstrap',
                         color=bootcolor, linestyle='dashed')
                plt.plot(nsubj_vec, boot_power_store_sd,
                         label='Bootstrap Step Down', color=bootcolor)

                plt.xlabel('Number of Subjects')
                plt.ylabel('Power')
                plt.title('$\pi_0$ = ' + str(pi0))
                if pi0 == 0.9 and powertype > 0:
                    plt.legend(loc='lower right')
                # plt.tight_layout()
                plt.subplots_adjust(left=0.15, right=0.95,
                                    bottom=0.17, top=0.9)
                plt.savefig(saveloc + 'dim_' + str(dim) + str(dim) + '_pi0_' +
                            str(int(100*pi0)) + '_powertype_' + str(powertype) + str(biomaddon))
                plt.close()

# %% Test pyplot.adjust
plt.plot(nsubj_vec, simes_power_store, label='Simes',
         color=paracolor, linestyle='dashed')
plt.plot(nsubj_vec, ARI_power_store, label='ARI', color=paracolor)
plt.plot(nsubj_vec, boot_power_store, label='Bootstrap',
         color=bootcolor, linestyle='dashed')
plt.plot(nsubj_vec, boot_power_store_sd,
         label='Bootstrap Step Down', color=bootcolor)

plt.xlabel('Number of Subjects')
plt.ylabel('Power')
plt.title('$\pi_0$ = ' + str(pi0))
plt.subplots_adjust(left=0.3, right=0.9, bottom=0.3, top=0.9)

# %% Testing
result_dir = 'C:\\Users\\12SDa\\global\\Intern\\drago\\Results\\Power_Results\\'
dim = 50
pi0 = 0.5
nsubj = 90
boot = np.load(result_dir + 'nsubj_' + str(nsubj) + '_pi0_' + str(int(100*pi0)) +
               '_nbootstraps_1000_niters_5000_dim_' + str(dim) + str(dim) + '_simtype_1.npz')
bootpower = boot['power']
print(bootpower)
