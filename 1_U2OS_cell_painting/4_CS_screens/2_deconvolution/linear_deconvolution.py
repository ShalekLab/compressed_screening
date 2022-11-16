import numpy as np
import pandas as pd
from sklearn.manifold import TSNE
import umap
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
import os
import sys
import random
import pickle
import multiprocess as mp
from itertools import repeat
from sklearn.utils import shuffle
from sklearn import linear_model
import copy

### Funtions to operate on compressed_screening object for linear deconvolution

def LM_decon_shuffle(CompressedScreen):    
    """
    Function randomly shuffles features matrix row indicies
    
    Input:       CompressedScreen object, which has been subsetted to a single screen
                 / condition for deconvolution
    
    Output:     a new CompressedScreen object called 'shuffled', where shuffled.features is randomly shuffled
    
    """
    shuffled=copy.deepcopy(CompressedScreen)
    shuffled.features = shuffle(shuffled.features).reset_index(drop=True)
    shuffled.features['index'] = shuffled.metadata.index    
    shuffled.features.set_index('index', drop=True, inplace=True)
        
    return shuffled
        
def LM_decon_reg(CompressedScreen, method='ridge', top_alpha=None, top_l1_ratio=None): 
    """
    Function runs linear regression using specified method to deconvolute compressed screen
    
    Input:      CompressedScreen object, which has been subsetted to a single screen
                / condition for deconvolution
                 
                method, options:   'ols' -
                                   'ridge' -
                                   'lasso' -
                                   'elasticnet' -
                 
                top_alpha, optional if method == 'ridge', 'lasso', 'elasticnet', 
                uses as regression parameter
                 
                top_l1_ratio, optional if method == 'elasticnet', 
                uses as regression parameter
    
    Output:     CompressedScreen object, where linear deconvolution model coefficients have been 
                saved in .decon_result['method']
                
                top_alpha, if method == 'ridge', 'lasso', 'elasticnet', returns CV-id'd or user-input value
                
                top_l1_ratio, if method == 'elasticnet', returns CV-id'd or user-input value
    
    """
    # OLS (ordinary least squares) regression
    if method == 'ols':
        reg = linear_model.LinearRegression()
        reg.fit(CompressedScreen.design_matrix, CompressedScreen.features)
                
        model_coef = pd.DataFrame(reg.coef_.T, columns = CompressedScreen.features.columns.values, 
                                  index = CompressedScreen.design_matrix.columns.values)
        
        CompressedScreen.save_deconvolution(model_coef, method)
        
        return CompressedScreen
    
    # Ridge regression
    if method == 'ridge':
        if top_alpha != None:
            reg = linear_model.Ridge(alpha=top_alpha)
            reg.fit(CompressedScreen.design_matrix, CompressedScreen.features)
            
        else:
            reg = linear_model.RidgeCV(alphas=np.logspace(-6, 6, 100))
            reg.fit(CompressedScreen.design_matrix, CompressedScreen.features)
            top_alpha = reg.alpha_
            
        model_coef = pd.DataFrame(reg.coef_.T, columns = CompressedScreen.features.columns.values, 
                                  index = CompressedScreen.design_matrix.columns.values)
        
        CompressedScreen.save_deconvolution(model_coef, method)
        
        return CompressedScreen, top_alpha
        
    # Lasso regression
    elif method == 'lasso':
        if top_alpha != None:
            reg = linear_model.MultiTaskLasso(alpha=top_alpha)
            reg.fit(CompressedScreen.design_matrix, CompressedScreen.features)
            
        else:
            reg = linear_model.MultiTaskLassoCV(n_jobs=-1)
            reg.fit(CompressedScreen.design_matrix, CompressedScreen.features)
            top_alpha = reg.alpha_
            
        model_coef = pd.DataFrame(reg.coef_.T, columns = CompressedScreen.features.columns.values, 
                                  index = CompressedScreen.design_matrix.columns.values)
        
        CompressedScreen.save_deconvolution(model_coef, method)
        
        return CompressedScreen, top_alpha
    
    # ElasticNet regression
    elif method == 'elasticnet':
        if top_alpha != None:
            reg = linear_model.MultiTaskElasticNet(l1_ratio=top_l1_ratio, alpha=top_alpha)
            reg.fit(CompressedScreen.design_matrix, CompressedScreen.features)
            
        else:
            reg = linear_model.MultiTaskElasticNetCV(l1_ratio=[0.01, 0.05, .1, 0.3, .5, .7, .9, .95, .99, 1], 
                                                     n_jobs=-1)
            reg.fit(CompressedScreen.design_matrix, CompressedScreen.features)
            top_alpha = reg.alpha_
            top_l1_ratio = reg.l1_ratio_
            
        model_coef = pd.DataFrame(reg.coef_.T, columns = CompressedScreen.features.columns.values, 
                                  index = CompressedScreen.design_matrix.columns.values)
        
        CompressedScreen.save_deconvolution(model_coef, method)
        
        return CompressedScreen, top_alpha, top_l1_ratio
    
    else:
        print("Improper method chosen, try again")
        return None

#### Function to run LM decon shuffle and return significant p-values
###### WORK IN PROGRESS
def LM_shuff_test(CompressedScreen, screen_name, num_shuff=10000, p_cutoff=0.05, method='ridge'):
    ### NOTES ABOUT WHAT THIS DOES    

    i = 0
    
    if method == 'ols':
        # Calculate non-shuffled model_coef
        CompressedScreen_non_shuff = LM_decon_reg(CompressedScreen, method=method)
        model_coef=CompressedScreen_non_shuff.get_deconvolution(method)
        pval_table = np.zeros(shape=model_coef.shape)
        
        # shuffle pool assignments to compute p-val
        print('starting shuffle test')
        while i<num_shuff:

            shuff_CompressedScreen = LM_decon_shuffle(CompressedScreen_non_shuff)
            shuff_CompressedScreen= LM_decon_reg(shuff_CompressedScreen, method=method)
            shuff_model_coef=shuff_CompressedScreen.get_deconvolution(method)
            pval_table = pval_table + (abs(model_coef.values) < abs(shuff_model_coef.values))*1
            i+=1
                
    elif method == 'elasticnet':
        # Calculate non-shuffled model_coef
        CompressedScreen_non_shuff, top_alpha, top_l1_ratio = LM_decon_reg(CompressedScreen, method=method)
        model_coef=CompressedScreen_non_shuff.get_deconvolution(method)
        pval_table = np.zeros(shape=model_coef.shape)
        
        # shuffle pool assignments to compute p-val
        print('starting shuffle test')
        while i<num_shuff:
        
            shuff_CompressedScreen = LM_decon_shuffle(CompressedScreen_non_shuff)
            shuff_CompressedScreen, _, _ = LM_decon_reg(shuff_CompressedScreen,
                                                  top_alpha=top_alpha,
                                                  top_l1_ratio=top_l1_ratio,
                                                  method=method)
            shuff_model_coef=shuff_CompressedScreen.get_deconvolution(method)
            pval_table = pval_table + (abs(model_coef.values) < abs(shuff_model_coef.values))*1
            i+=1
            
    else:
        # Calculate non-shuffled model_coef
        CompressedScreen_non_shuff, top_alpha = LM_decon_reg(CompressedScreen, method=method)
        model_coef=CompressedScreen_non_shuff.get_deconvolution(method)
        pval_table = np.zeros(shape=model_coef.shape)
        
        # shuffle pool assignments to compute p-val
        print('starting shuffle test')
        while i<num_shuff:
        
            shuff_CompressedScreen = LM_decon_shuffle(CompressedScreen_non_shuff)
            shuff_CompressedScreen, _ = LM_decon_reg(shuff_CompressedScreen,
                                               top_alpha=top_alpha,
                                               method=method)
            shuff_model_coef=shuff_CompressedScreen.get_deconvolution(method)
            pval_table = pval_table + (abs(model_coef.values) < abs(shuff_model_coef.values))*1
            i+=1
        
    pval_table = pval_table/num_shuff
    print('done with shuffle')
    
    # plot p val distribution
    plt.hist(pval_table.flatten(),bins='auto',color='grey')
    plt.axvline(p_cutoff)
    
    if not os.path.isdir('plots/'):
        os.mkdir('plots')
    plt.savefig('plots/'+method+ '_%s_p_val_dist.pdf'%screen_name)
    # if PCA_trans == False:
    
    #     plt.savefig('plots/' +pop+'_'+method+ '_run2_p_val_dist.pdf')
        
    # else:
        
    #     plt.savefig('plots/' +pop+'_'+method+ 'PCA_run2_p_val_dist.pdf')
        
    plt.close()
    
    # drop coefficients for p>0.05
    permute_model_coef = (pval_table < p_cutoff)*model_coef
    pval_table = pd.DataFrame(pval_table, 
                              columns=list(model_coef.columns.values),
                              index=list(model_coef.index.values))
    
    print('done!')
    return permute_model_coef, model_coef, pval_table

### add on Conner's functions to evaluate deconvolution


### Helper functions for metrics evaluating the performance of deconvolution (aside from vanilla euclidean distance)
def get_drug_consistency(compressed_screen, ground_truth):
    """
    Drug consistency here is define for each drug as 
    the percentage of features where the 'directionality' (positive, negative, none) of 
    drug effect is consistent between compressed and singular drug screens.

    Output:     a numpy array with the same length as the number of perturbations
    NOTE: the two input matrices/dataframes should have the same indices (perturbation) and 
    the same columns (features)
    """
    def check_single_effect(a, b):
        if a*b>0 or (a==0 and b==0):
            return True
        else: return False
    vectorize_check=np.vectorize(check_single_effect)
    all_consistency=vectorize_check(compressed_screen, ground_truth)
    drug_consistency=np.mean(all_consistency, axis=1)
    return drug_consistency

def get_vector_norm_diff(compressed_screen, ground_truth, norm_ord=1):
    """
    compute the difference between L1 norms of each drug's coeff vector
    NOTE: the two input matrices/dataframes should have the same indices (perturbation) and 
    the same columns (features)

    Output:     a numpy array with the same length as the number of perturbations
    """
    ground_truth_vector_norm=np.linalg.norm(ground_truth.values,ord=norm_ord, axis=1)
    cs_vector_norm=np.linalg.norm(compressed_screen.values,ord=1, axis=1)
    vector_norm_diff=cs_vector_norm-ground_truth_vector_norm
    return vector_norm_diff

def get_featureNormalized_distance(compressed_screen,ground_truth):
    """ 
    Normalizes both input matrix within each feature, such that the max absolute value
    of effect for each feature is 1.  Then, just outputs the Euclidean distance for 
    each drugs given the transformed vectors. 
    NOTE: the two input matrices/dataframes should have the same indices (perturbation) and 
    the same columns (features)

    Output:     a numpy array with the same length as the number of perturbations
    """
    compressed_screen_copy=compressed_screen/compressed_screen.abs().max()
    compressed_screen_copy.fillna(0, inplace=True)
    ground_truth_copy=ground_truth/ground_truth.abs().max()
    ground_truth_copy.fillna(0, inplace=True)
    diff=compressed_screen_copy.values-ground_truth_copy.values
    modified_euc_dist=np.linalg.norm(diff, ord=2, axis=1)
    return modified_euc_dist

def strong_drug_recovery(compressed_screen,ground_truth, n_top=10):
    """ 
    Evaluates how much of strong drugs are 'correctly' (same direction) identified in compressed setting
    compared to ground truth.

    NOTE: the two input matrices/dataframes should have the same indices (perturbation) and 
    the same columns (features)

    Output:     a numpy array with the same length as the number of features
    """
    def preserve_effect_in_feature(feature, n_top=10):
        """
        helper for comparing how many of the top n strongest drugs are preserved
        in the compressed setting. Inputs are two columns in the
        corresponding drug x feature matrix
        """
        effect_compressed, effect_gt=compressed_screen[feature], ground_truth[feature]
        top_n_drugs=effect_gt.abs().rank(axis=0)[:n_top].index
        signs=list(effect_gt[top_n_drugs]>=0)
        top_n_drugs_comp=effect_compressed.abs().rank(axis=0)[:n_top].index
        signs_comp=list(effect_compressed[top_n_drugs_comp]>=0)
        count=sum([1 for i in range(n_top) if (top_n_drugs[i] in top_n_drugs_comp \
                                     and signs_comp[list(top_n_drugs_comp).index(top_n_drugs[i])]==signs[i])])
        return count/float(n_top)
    effect_in_feature=np.vectorize(preserve_effect_in_feature)
    strong_drugs_preserved=effect_in_feature(ground_truth.columns, n_top)
