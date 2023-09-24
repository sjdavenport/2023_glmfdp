"""
Functions used for analysing the HCP data
"""

import os
import pandas as pd
import numpy as np

def obtain_hcp_data():
    df = pd.read_csv('/storage/store/data/HCP900/RESTRICTED_bthirion_5_26_2021_14_53_29.csv')
    df_ur = pd.read_csv('/storage/store/data/HCP900/unrestricted_hcp_s900.csv')

    # Obtain the ids of the subjects for which the image data is available
    subjects = os.listdir('/storage/store/data/HCP900/glm/')
    subject_image_exist_idx = np.zeros(len(subjects), dtype = 'bool')
    for I in np.arange(len(subjects)):
        if os.path.exists('/storage/store/data/HCP900/glm/' + subjects[I] + '/WM/level2/z_maps/z_2BK-0BK.nii.gz'):
            subject_image_exist_idx[I] = 1

    subjects = np.array(subjects, dtype = 'int')
    subjects_with_images = subjects[subject_image_exist_idx]

    # Restrict the data frame to this subset
    df_ur = df_ur[df_ur.Subject.isin(subjects_with_images)].reset_index()
    df = df[df.Subject.isin(subjects_with_images)].reset_index()

    # Restrict to the unrelated subjects
    npstringidx = np.array(['NotTwin'])
    df_unrelated = df[df.ZygositySR.isin(npstringidx)]

    df_ur_unrelated = df_ur[df_ur.Subject.isin(df_unrelated.Subject)].reset_index()
    df_unrelated = df[df.Subject.isin(df_ur_unrelated.Subject)].reset_index()

    nsubj = len(df_unrelated)
    imgs = [None] * len(df_unrelated)
    for I in np.arange(nsubj):
        imgs[I] = '/storage/store/data/HCP900/glm/' + str(df_unrelated.Subject[I]) + '/WM/level2/z_maps/z_2BK-0BK.nii.gz'
        
    return df_unrelated, df_ur_unrelated, imgs

def lm_setup(restricted_vars2use, unrestricted_vars2use, contrast_vars, df_unrelated, df_ur_unrelated, nsubj):
    all_vars2use = restricted_vars2use + unrestricted_vars2use
    
    # Compute the number of parameters
    p_restricted = len(restricted_vars2use)
    p_unrestricted = len(unrestricted_vars2use)
    n_params = p_restricted + p_unrestricted

    # Initialize the contrast matrix
    n_contrasts = len(contrast_vars)
    contrast_matrix = np.zeros((n_contrasts, n_params))
    for L in np.arange(n_contrasts):
        contrast_matrix[L, all_vars2use.index(contrast_vars[L])] = 1

    design = np.zeros([nsubj, n_params])

    # Replace M/F in the sex register as 0/1
    df_ur_unrelated.Gender =  df_ur_unrelated.Gender.replace(['M', 'F'], [0,1])

    # Obtain the design matrix
    for I in np.arange(p_restricted):
        design[:, I] = df_unrelated[restricted_vars2use[I]].values[:nsubj]
    
    for I in np.arange(p_unrestricted):
        design[:, I + p_restricted] = df_ur_unrelated[unrestricted_vars2use[I]].values[:nsubj]
    
    return design, contrast_matrix
    
