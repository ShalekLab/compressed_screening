### This script will reformat data prior to initializing as compressed_screening object

## req inputs: features, metadata, sample_perturbation_dictionary, metadata_cols_to_create_unique_ID
## optional inputs:

## what this should do: make sure all samples + perturbations have unique IDs
## reformats into:


# Todo
	# format the code for the cell painitng stuff
		# Make a function for writing unique IDs
		# load in metadata and features
		# Assign unique ids to metadata and features
		# Load in well drug dicts and make the final design matrix
	# Also get the code going for all of the ligand formats
	# Include checks
		# Should throw errors for: non-unique sample_ID / perturbation_ID


#   feratures
#     sample_ID x features cols
 
#   metadata
#     sample_ID x metadata cols


#   sample_perturbation_design_matrix
#       sample_ID x perturbation ID

# Note for ben, he is going to need to specific the landmark and control ids when doing the deconvolution


def load_metadata_features_cell_painting(path_data,landmarks,controls):
	"""
	Loads in the metadata and features from a cell painiting screen

	Inputs: 
		path_data - the path to data.frame csv with metadata and features in it
		landmarks - list of the landmakr perturbations
		controls  - list of the control perturbation
	"""

    #### load in STAN features
    data= pd.read_csv(path_data)
 
    # NOTE: plate 10 has 2 wells with improper # of drugs, all plates have 8X DMSO level

    # split into features and metadata    
    features_name = [col for col in \
                data.columns if 'Metadata' not in col and \
                'Number_of_Cells' not in col]

    metadata_name = [col for col in \
                data.columns if 'Metadata' in col or \
                'Number_of_Cells' in col]

    #STAN hi-cor features + metadata
    features = data.loc[:,features_name].reset_index(drop=True)
    metadata = data.loc[:,metadata_name].reset_index(drop=True)

    # add compression scheme to metadata
    metadata['Metadata_comp_method'] = metadata['Metadata_perturbation']

    for perturb in landmarks+controls:
        metadata.loc[metadata.Metadata_comp_method==perturb,'Metadata_comp_method']='control'
        
    metadata['Metadata_comp_method'] = metadata['Metadata_compression'].astype(str)+'x_'+\
                                       metadata['Metadata_replicates'].astype(str)+'r_'+\
                                       metadata['Metadata_comp_method']

    return metadata, features


def main():

	LM = ['Daunorubicin','Doxorubicin','Fluvastatin','Mycophenolic Acid','Riluzole','Ticlopidine','Tropicamide','Vinblastine Sulfate']
	controls = ['DMSO']
	path = '/Users/connerkummerlowe/Dropbox (MIT)/lets_do_drugs//BEM/CellPainting_compressed_screen/3_Norm_FeatureSelection/2_Feature_Selection_re-run1+2/Output/04122021_variFeat_compressed_screen_run1_STAN_ALLcor_feature_table.gz'
	metadata, features = load_metadata_features_cell_painting(path,LM,controls)



if __name__ == '__main__':
	main()


