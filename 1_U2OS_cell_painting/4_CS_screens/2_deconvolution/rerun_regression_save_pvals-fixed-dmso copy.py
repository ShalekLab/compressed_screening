import scanpy as sc
import os, glob
import sys
from math import floor
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import umap
import random
import pickle
import multiprocess as mp
from itertools import repeat
from sklearn.decomposition import PCA
from sklearn.utils import shuffle
from sklearn import linear_model
from ingest_screening_data import *
from CompressedScreen import *
from linear_deconvolution import *

LM = ['Daunorubicin','Doxorubicin','Fluvastatin','Mycophenolic Acid',
          'Riluzole','Ticlopidine','Tropicamide','Vinblastine Sulfate']
controls = ['DMSO']

def main():
    #loop over preprocessing
    input_dir = "../../../cell_painting_data_lock/4_CS_OBJECTS_median_aggregated/"
    metadata_files= glob.glob(input_dir+"*metadata.csv")
    # runs=['CS_run1','CS_run2','CS_run3']
    runs=['CS_run2','CS_run3']
    for run in runs:
        # for layer in ['PCH','raw']:
        for layer in ['raw']:
            print(run+"_"+layer)
            # Find the right metadata file for this run & layer
            metadata_file = [x for x in metadata_files if run+"_"+layer in x][0]
            print(metadata_file)
            features_file = metadata_file.split("metadata.csv")[0] + "features.csv"
            design_matrix_file = metadata_file.split("metadata.csv")[0] + "design_matrix.csv"

            metadata = pd.read_csv(metadata_file,index_col=0)
            features = pd.read_csv(features_file,index_col=0)
            design_matrix = pd.read_csv(design_matrix_file,index_col=0)
            # design matrices should be all fixed now
            

            outdir=os.path.join('~/Dropbox (MIT)/lets_do_drugs/CK/compressed_cell_painting_screen_analysis/regression_deconvolution_December2021/decon_out',"_".join(metadata_file.split("/")[-1].split("metadata.csv")[0].split("_")[2:])+layer+'_decon_out/')

            if run!='CS_run3':
        #         optimization_schemes = ['random','weighted_cosine','phenograph_all','mahalanobis']
                optimization_schemes = ['random']
                pooled_meta = metadata.loc[np.isin(metadata['Metadata_perturbation'],optimization_schemes)]
            else:
                optimization_schemes='random'
                pooled_meta = metadata.loc[metadata['Metadata_perturbation'].str.contains(optimization_schemes)]
            pooled_meta.Metadata_perturbation = pooled_meta.Metadata_perturbation.astype('category')
            pooled_meta.loc[:,'Metadata_perturbation']=pooled_meta.Metadata_perturbation.cat.set_categories(pooled_meta.Metadata_perturbation.unique()).values
            scheme_data = pooled_meta.groupby(['Metadata_compression','Metadata_replicates','Metadata_perturbation']).size().reset_index()
            for i in range(len(scheme_data)):
                optimization = scheme_data['Metadata_perturbation'][i]
                compression = scheme_data['Metadata_compression'][i]
                replicates = scheme_data['Metadata_replicates'][i]
                indices=metadata.index[(metadata.Metadata_perturbation==optimization)&(metadata.Metadata_compression==compression)&(metadata.Metadata_replicates==replicates)]
                cs=CompressedScreen(metadata.loc[indices,:], design_matrix.loc[indices,:], features.loc[indices,:])
                scheme_name=str(compression)+"x_"+str(replicates)+"r_"+optimization
                print(scheme_name)
                permute_model_coef, model_coef, pval_table=LM_shuff_test(cs, run+'_'+scheme_name, num_shuff=1000, p_cutoff=0.05, method='elasticnet')
                permute_model_coef.to_csv(os.path.join(outdir,'%s_%s_permute_model_coef.csv'%(run, scheme_name)))
                model_coef.to_csv(os.path.join(outdir,'%s_%s_model_coef.csv'%(run, scheme_name)) )
                pval_table.to_csv(os.path.join(outdir,'%s_%s_pval_table.csv'%(run, scheme_name)))
                
    runs=['CS_run3']
    for run in runs:
        for layer in ['PCH']:
            print(run+"_"+layer)
            # Find the right metadata file for this run & layer
            metadata_file = [x for x in metadata_files if run+"_"+layer in x][0]
            print(metadata_file)
            features_file = metadata_file.split("metadata.csv")[0] + "features.csv"
            design_matrix_file = metadata_file.split("metadata.csv")[0] + "design_matrix.csv"

            metadata = pd.read_csv(metadata_file,index_col=0)
            features = pd.read_csv(features_file,index_col=0)
            design_matrix = pd.read_csv(design_matrix_file,index_col=0)
            # design matrices should be all fixed now
            

            outdir=os.path.join('~/Dropbox (MIT)/lets_do_drugs/CK/compressed_cell_painting_screen_analysis/regression_deconvolution_December2021/decon_out',"_".join(metadata_file.split("/")[-1].split("metadata.csv")[0].split("_")[2:])+layer+'_decon_out/')

            if run!='CS_run3':
        #         optimization_schemes = ['random','weighted_cosine','phenograph_all','mahalanobis']
                optimization_schemes = ['random']
                pooled_meta = metadata.loc[np.isin(metadata['Metadata_perturbation'],optimization_schemes)]
            else:
                optimization_schemes='random'
                pooled_meta = metadata.loc[metadata['Metadata_perturbation'].str.contains(optimization_schemes)]
            pooled_meta.Metadata_perturbation = pooled_meta.Metadata_perturbation.astype('category')
            pooled_meta.loc[:,'Metadata_perturbation']=pooled_meta.Metadata_perturbation.cat.set_categories(pooled_meta.Metadata_perturbation.unique()).values
            scheme_data = pooled_meta.groupby(['Metadata_compression','Metadata_replicates','Metadata_perturbation']).size().reset_index()
            for i in range(len(scheme_data)):
                optimization = scheme_data['Metadata_perturbation'][i]
                compression = scheme_data['Metadata_compression'][i]
                replicates = scheme_data['Metadata_replicates'][i]
                indices=metadata.index[(metadata.Metadata_perturbation==optimization)&(metadata.Metadata_compression==compression)&(metadata.Metadata_replicates==replicates)]
                cs=CompressedScreen(metadata.loc[indices,:], design_matrix.loc[indices,:], features.loc[indices,:])
                scheme_name=str(compression)+"x_"+str(replicates)+"r_"+optimization
                print(scheme_name)
                permute_model_coef, model_coef, pval_table=LM_shuff_test(cs, run+'_'+scheme_name, num_shuff=1000, p_cutoff=0.05, method='elasticnet')
                permute_model_coef.to_csv(os.path.join(outdir,'%s_%s_permute_model_coef.csv'%(run, scheme_name)))
                model_coef.to_csv(os.path.join(outdir,'%s_%s_model_coef.csv'%(run, scheme_name)) )
                pval_table.to_csv(os.path.join(outdir,'%s_%s_pval_table.csv'%(run, scheme_name)))


if __name__=='__main__':
    main()
    