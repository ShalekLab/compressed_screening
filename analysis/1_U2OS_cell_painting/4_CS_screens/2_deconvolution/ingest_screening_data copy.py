import os
import sys
import numpy as np
import pandas as pd
import random
import pickle
from math import floor
import warnings

####################################################################################################
####################################  # Readme / How to use ########################################
####################################################################################################

# Script for ingesting cell painting

# For each cell painting run, writes a
    # Metadata pandas table as a csv
    # Features pandas table as a csv
    # Design matrix pandas table as a csv

# Currently set up to use 'Metadata_run','Metadata_Plate', and 'Metadata_Well' to generate unique ids for
    # "perbutations" (where a perturbation can be DMSO, a single drug, or a pool of drugs)

# Wrote the code in two ways in the main function below
    
    # Case 1: Ingesting a single cell painting run
        # Specify 
            # 1) (pandas table) The pandas table with metadata and features in it & 
            # 2) (string)the path to the well_drug_dicts 
            # 3) (list of strings) The landmarks used for the screen 
            # 4) (list of strings) The vehicle control / controls 

    # Case 2: Ingesting several cell painting runs at a time from a single big hunk o' data
        # Specify 
            # 1) (pandas table) the pandas table with metadata and features in it &
            # 2) (dictionary) a dictionary of paths to all of the well drug dictionary paths (as strings) for all of the cell painting runs
            # 3) (dictionary) a dictionary of the list of strings of landmarks used for each of the cell painting runs
            # 4) (dictionary) a dictionary of the list of strings of vehicle controls used for each of the cell paitning runs

####################################################################################################
####################################  # Helper functions ###########################################
####################################################################################################


def load_metadata_features_cell_painting(metadata_features,landmarks,controls):
    """
    Loads in the metadata and features from a cell painiting screen

    Inputs: 
        metadata_features - pandas table with metadata and features in it
        landmarks - list of the landmakr perturbations
        controls  - list of the control perturbation
    """


    # NOTE: plate 10 has 2 wells with improper # of drugs, all plates have 8X DMSO level

    # split into features and metadata    
    features_name = [col for col in \
                metadata_features.columns if 'Metadata' not in col and \
                'Number_of_Cells' not in col]

    metadata_name = [col for col in \
                metadata_features.columns if 'Metadata' in col or \
                'Number_of_Cells' in col]

    #STAN hi-cor features + metadata
    features = metadata_features.loc[:,features_name].reset_index(drop=True)
    metadata = metadata_features.loc[:,metadata_name].reset_index(drop=True)

    # add compression scheme to metadata
    metadata['Metadata_comp_method'] = metadata['Metadata_perturbation']

    for perturb in landmarks+controls:
        metadata.loc[metadata.Metadata_comp_method==perturb,'Metadata_comp_method']='control'

    metadata['Metadata_comp_method'] = metadata['Metadata_compression'].astype(str)+'x_'+\
                                       metadata['Metadata_replicates'].astype(str)+'r_'+\
                                       metadata['Metadata_comp_method']

    return metadata, features

def load_sample_perturbation_dictionaries(metadata,landmarks,sample_perturbation_dict_path):
    """
    Loads in the sample perturbation dictionaries for all plates in a given screen
    
    Inputs: metadata, pd.Dataframe
            landmarks, list of strings for the landmark perturbations
    
    """
    

    sample_perturbation_dicts = {}

    for p in set(metadata['Metadata_Plate']):
        # if '_' in p: #if in the format of CS_run1_plate4
        #     plate=p.split('_')[-1]
        # else:
        #     plate = p
        plate = p
        with open(sample_perturbation_dict_path+plate+'.pickle','rb') as fp:
            sample_perturbation_dicts[p] = pickle.load(fp)
            
    # fix well names in sample_perturbation_dicts
    new_sample_perturbation_dicts = {}

    for p in set(metadata['Metadata_Plate']):
        plate_sample_perturbation_dicts = {}
        for w in list(sample_perturbation_dicts[p].keys()):
        
            string = w
            newString = string[0]+string[1:] if len(string) == 3 else string[0]+'0'+string[1]
        
            plate_sample_perturbation_dicts[newString] = sample_perturbation_dicts[p][w]
        
        new_sample_perturbation_dicts[p] = plate_sample_perturbation_dicts
        
    sample_perturbation_dicts = new_sample_perturbation_dicts         

    # denote landmark (LM) compounds in sample_perturbation_dicts
    new_sample_perturbation_dicts = {}
    for p in set(metadata['Metadata_Plate']):
        wells = set(metadata.loc[metadata.Metadata_Plate == p, 'Metadata_Well'])

        plate_sample_perturbation_dicts = {}
        for w in wells:
            perturb = metadata.loc[(metadata.Metadata_Plate==p)&(metadata.Metadata_Well==w),
                                  'Metadata_perturbation'].values

#             if perturb in LM:
#                 replace = 'LM_'+perturb
#                 plate_sample_perturbation_dicts[w] = replace.tolist()

            plate_sample_perturbation_dicts[w] = sample_perturbation_dicts[p][w]

        new_sample_perturbation_dicts[p] = plate_sample_perturbation_dicts

    sample_perturbation_dicts = new_sample_perturbation_dicts
    return(sample_perturbation_dicts) 

def make_unique_sample_ids(metadata,column_names_for_id):
    """
    Inputs
    metadata: pd.Dataframe, each row is a sample and each column is a metadata column
    column_names_for_id: list of strings specifying the metadata column names to use for the unique ID
    sep: strring to use when separating column names in the unique ID
    
    Output: numpy array of the unique ids
    """
    if len(column_names_for_id) > 1:
        ids = metadata[column_names_for_id[0]].values.astype(str).astype(object) + "_" + metadata[column_names_for_id[1]].values.astype(str).astype(object)
        for i in range(2,len(column_names_for_id)):
            ids = ids + "_" + metadata[column_names_for_id[i]].values.astype(str).astype(object)
    else:
        ids = metadata[column_names_for_id[0]].values.astype(str).astype(object)
    return(ids)


####################################################################################################
####################################    Data ingesting   ###########################################
####################################################################################################
# Function for ingesting one screen at a time when 
def ingest_cell_painting(run_name,landmarks, controls,metadata,features,metadata_columns_for_unique_ids,
                        sample_perturbation_dict_path):
    """
    Main function for ingesting the cell painting data

    
    Inputs
        - run_name: Unique string indicating the compressed screening run (i.e. Nov2020)
        - landmarks: list of strings of the landmark perturbations
        - controls: list of strings of the control perturbations
        - metadata_features_path: path to the file with the metadata and features data
        - metaadata_columns_for_unique_ids: list of metadata columns to concatenate together when making the unique ids
        - sample_perturbation_dict_path: path to the dir containing the pickled dictionaries for the well_drug_dicts for the plates
     
     Outputs:
         metadata
         features
         design_matrix
    """
    
    # Loading in features and metadata
    metadata['Metadata_run'] = run_name
    
    # Generating unique sample ids and attaching them to features and metadata
    unique_ids = make_unique_sample_ids(metadata,metadata_columns_for_unique_ids)
    metadata.index = unique_ids
    features.index = unique_ids
    
    # Loading in the sample perturbation dictionaries
    sample_perturbation_dicts = load_sample_perturbation_dictionaries(metadata=metadata,
                                                                    landmarks=landmarks,
                                                                     sample_perturbation_dict_path=sample_perturbation_dict_path)
    
    # getting the unique drugs out off the sample perturbation dict
    perturbations = []
    for plate in sample_perturbation_dicts:
        for well in sample_perturbation_dicts[plate]:
            new_perturbations = sample_perturbation_dicts[plate][well]
            for new_perturbation in new_perturbations:
                if new_perturbation not in perturbations:
                    perturbations.append(new_perturbation)

    # initializing the design matrix
    design_matrix = pd.DataFrame(data = np.zeros((len(unique_ids),len(perturbations))),
                                index = unique_ids,
                                columns=perturbations)
    # assigning the drug names in the design matrix
    for plate in sample_perturbation_dicts:
        for well in sample_perturbation_dicts[plate]:
            id_to_use = run_name +"_"+plate+"_" + well
            design_matrix.loc[id_to_use,sample_perturbation_dicts[plate][well]] = 1
            
    return metadata, features, design_matrix
    


    

def main():

    # Example code use cases:

    ################################################################################################
    # Case 1: Ingesting the data straight out of the compressed screens
    # Useful sort of setup to use when ingesting one screen at a time
    # Could be useful for importing any future cell painting data
    ################################################################################################

    ##########################
    # Example: Compressed screen run 1
    ##########################

    # Specifying the run name, landmarks and controls
 #    run = "CS_run1"
	# LM = ['Daunorubicin','Doxorubicin','Fluvastatin','Mycophenolic Acid',
 #      'Riluzole','Ticlopidine','Tropicamide','Vinblastine Sulfate']
 #    controls =['DMSO']
	

 #    # Specifying the file path for the metadata and features pandas table
 #    data= pd.read_csv('/Users/connerkummerlowe/Dropbox (MIT)/lets_do_drugs//BEM/CellPainting_compressed_screen/3_Norm_FeatureSelection/2_Feature_Selection_re-run1+2/Output/04122021_variFeat_compressed_screen_run1_STAN_ALLcor_feature_table.gz')
 #    metadata, features = load_metadata_features_cell_painting(data,LM,controls)

 #    # Specifying the well drug dictionary path
 #    well_drug_dict_dir = '/Users/connerkummerlowe/Dropbox (MIT)/lets_do_drugs//BEM/CellPainting_compressed_screen/4_Analysis//0_RAW_data/plate_dicts_run1/'
	

 #    metadata, features, design_matrix = ingest_cell_painting(run_name=run,landmarks=LM,controls=controls,metadata=metadata,features=features,
 #                                                        metadata_columns_for_unique_ids=['Metadata_run','Metadata_Plate','Metadata_Well'],
 #                                                        sample_perturbation_dict_path=well_drug_dict_dir)
	

    # Saving the resulting files
    # metadata.to_csv("")
    # features.to_csv("")
    # design_matrix.to_csv("")


    ################################################################################################
    # Case 2: Running on the data all cell painting data, going to go with this from here on
    ################################################################################################

    # Loading in the dataset with all of the cell painting runs in it
    data_preprocessing_method = "denoised_data" #use this to specify the preprocessing that's been used
    data = pd.read_csv("~/Dropbox (MIT)/lets_do_drugs/cell_painting_data_lock/2_median_aggregated/08122021_QC_both_all_data_feature_table.gz")


    LM = ['Daunorubicin','Doxorubicin','Fluvastatin','Mycophenolic Acid',
            'Riluzole','Ticlopidine','Tropicamide','Vinblastine Sulfate']
    controls = ['DMSO']

    # Splitting data into metadata nd features
    all_metadata, all_features = load_metadata_features_cell_painting(data,LM,controls)

    # Renaming the plate names to be more consistent with my previous setup
    for plate in np.unique(all_metadata['Metadata_Plate'].values):
        for run in np.unique(all_metadata['Metadata_run'].values):
            if run in plate:
                all_metadata['Metadata_Plate'] = all_metadata['Metadata_Plate'].replace(plate,plate.replace(run+"_",""))

    # Setting the file paths for all of the well drug dict dictionaries
    CS_run1_well_drug_dict_dir = '/Users/connerkummerlowe/Dropbox (MIT)/lets_do_drugs//BEM/CellPainting_compressed_screen/4_Analysis//0_RAW_data/plate_dicts_run1/'
    CS_run2_well_drug_dict_dir = '/Users/connerkummerlowe//Github/lets_do_drugs/cell_painting_screen_design/Picklist design (Second screen)/well_drug_dicts/'
    CS_run3_well_drug_dict_dir = '/Users/connerkummerlowe/Github/lets_do_drugs/Picklist design (Third Screen)/well_drug_dicts_without_ions/'

    GT_run1_batch1_well_drug_dict_dir = '/Users/connerkummerlowe/GitHub/lets_do_drugs/Picklist design (Ground Truth run1 batch1)/well_drug_dicts/'
    GT_run1_batch2_well_drug_dict_dir = '/Users/connerkummerlowe/GitHub/lets_do_drugs/Picklist design (Ground Truth run1 batch 2)/well_drug_dicts/'
    GT_run2_well_drug_dict_dir = '/Users/connerkummerlowe/GitHub/lets_do_drugs/Picklist design (Ground truth run2)/well_drug_dicts/'

    
    well_drug_dict_dir_map = {"CS_run1":CS_run1_well_drug_dict_dir,
                              "CS_run2":CS_run2_well_drug_dict_dir,
                              "CS_run3":CS_run3_well_drug_dict_dir,
                              "GT_run1_batch1":GT_run1_batch1_well_drug_dict_dir,
                              "GT_run1_batch2":GT_run1_batch2_well_drug_dict_dir,
                              "GT_run2":GT_run2_well_drug_dict_dir}

    controls_map = {"CS_run1":['DMSO'],
                              "CS_run2":['DMSO'],
                              "CS_run3":['DMSO'],
                              "GT_run1_batch1":['DMSO'],
                              "GT_run1_batch2":['DMSO'],
                              "GT_run2":['DMSO']}

    landmarks_map = {"CS_run1":LM,
                              "CS_run2":LM,
                              "CS_run3":LM,
                              "GT_run1_batch1":[],
                              "GT_run1_batch2":[],
                              "GT_run2":[]}


    # Loop through each of the cell painting runs and ingest the data into a 
        # metadata csv
        # features csv
        # design matrix csv
    for run in np.unique(all_metadata['Metadata_run'].values):
        well_drug_dict_dir = well_drug_dict_dir_map[run]
        landmarks = landmarks_map[run]
        controls = controls_map[run]
        metadata = all_metadata.loc[all_metadata.Metadata_run==run,:]
        features = all_features.loc[all_metadata.Metadata_run==run,:]
        metadata, features, design_matrix = ingest_cell_painting(run_name=run,landmarks=landmarks,controls=controls,metadata=metadata,features=features,
                                                            metadata_columns_for_unique_ids=['Metadata_run','Metadata_Plate','Metadata_Well'],
                                                            sample_perturbation_dict_path=well_drug_dict_dir)

        metadata.to_csv(run+"_" + data_preprocessing_method+"_metadata.csv")
        features.to_csv(run+"_" + data_preprocessing_method+"_features.csv")
        design_matrix.to_csv(run+"_" + data_preprocessing_method+"_design_matrix.csv")  
     


if __name__ == '__main__':
	main()























