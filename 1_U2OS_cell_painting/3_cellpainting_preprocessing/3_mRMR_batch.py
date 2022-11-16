### BEM 09252021 Run pymRMR on filtered datasets

import os
import pandas as pd
import numpy as np
from pymrmr import mRMR
import matplotlib.pyplot as plt
import seaborn as sns

def corr_heatmap(data, features, plotname):

    feature_all_values = {} # Dict with keys as features and values as ALL well values
    for f in features:
        feature_all_values[f] = data[f].values

    heatmap_data = []
    for f in feature_all_values.keys():
        feature1 = feature_all_values[f]

        tmp = []
        for f2 in feature_all_values.keys():
            feature2 = feature_all_values[f2]
            cor = np.corrcoef(feature1, feature2)[0, 1]
            tmp.append(cor)
        heatmap_data.append(tmp)

    sns.clustermap(data=heatmap_data, cmap="RdBu")
    plt.savefig(plotname +'_feature_CorrelationHeatMap.jpg')
    plt.close()

def run_mrmr(path_to_data, data_out_name, plot_out_name):
    
    data = pd.read_csv(path_to_data, low_memory = False)

    features = [col for col in \
                    data.columns if 'Metadata' not in col and \
                    'Number_of_Cells' not in col]

    ### Run mRMR to select features after signal / noise filtering
    # reformat datatable to only features + perturbations in GT/landmark
    perturbs = list(set(data['Metadata_perturbation'].values) - {'fulldeckrandom1','fulldeckrandom2',
                                                                 'mahalanobis','phenograph_all','random',
                                                                 'random1','random2','weighted_cosine'})

    data_mRMR = data.loc[data.Metadata_perturbation.isin(perturbs),
                         (['Metadata_perturbation']+list(features))]

    i = 0
    for p in perturbs:
        data_mRMR.at[data_mRMR.Metadata_perturbation == p,'Metadata_perturbation'] = i
        i += 1

    feats_mRMR = mRMR(data_mRMR, 'MIQ', 150)

    metadata_all = [col for col in data.columns if 'Metadata' in col]

    data_mrmr = data[metadata_all+['Number_of_Cells']+list(feats_mRMR)]

    data_mrmr.to_csv(data_out_name, index=False, compression='gzip')

    corr_heatmap(data_mrmr, feats_mRMR, plot_out_name)


#### Main Script starts here
if __name__ == '__main__':
    
    dataname = input("Enter data name: ")
    
    path_to_data = '09242021_QC_both_'+dataname+'_feature_table.gz'
    data_out_name = '09252021_QC_mrmr_'+dataname+'_feature_table.gz'
    plot_out_name = '09252021_mrmr_'+dataname

    run_mrmr(path_to_data, data_out_name, plot_out_name)
