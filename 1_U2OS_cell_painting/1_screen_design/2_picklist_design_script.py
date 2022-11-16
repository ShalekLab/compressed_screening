#!/usr/bin/env python3

# Comp Screening - Example code for designing picklist for Echo acoustic liquid handler
## For designing compressed screens
#### _CK Nov. 15, 2022_

### This code is currently highly bespoke to our drug library and data formats.

### What this does: makes filesfor running a compressed screens with an Echo acoustic liquid handler
### With an Echo, need to tell it to flies drops from source plates to destination plates
### - Inputs
###		- Drug library plate layout csv
###		- Npy files of the random pool layour from make_random_pools.py
### - Outputs for each desintation plate:
###		- picklist: csv input to Echo, tells it what drops come from  which source plates and go to which wells in destination plate
###		- Metadata: csv that records for each well in the desintation plate what the dose, time, comrpession, replicates, etc are
###		- Well_drug_dict: dictionary describing wells in destation plate, keys are wells, values are lists of the perturbations in that well
###		- Pool_key: intermediate file used for building above, not used in downstream analysis


import sys
import numpy as np
import pandas as pd
import pickle
import math

##############################################################################
##############################################################################
######################## Loading functions ###################################
##############################################################################
##############################################################################

# Returns a numpy array of all well numbers in a 384 well plate
	# Does not have leading zero bc there is no leading zero in the echo
	# source map
def generate_wells_384_picklist():
	well_rows = np.hstack((np.repeat("A",24),
		  np.repeat("B",24),
		  np.repeat("C",24),
		  np.repeat("D",24),
		  np.repeat("E",24),
		  np.repeat("F",24),
		  np.repeat("G",24),
		  np.repeat("H",24),
		  np.repeat("I",24),
		  np.repeat("J",24),
		  np.repeat("K",24),
		  np.repeat("L",24),
		  np.repeat("M",24),
		  np.repeat("N",24),
		  np.repeat("O",24),
		  np.repeat("P",24)))
	well_cols = np.tile(np.arange(24),16)
	well_cols += 1
	well_cols = [str(int(number)) for number in well_cols]
	wells_384 = np.core.defchararray.add(well_rows, well_cols)
	return(wells_384)


# Returns a numpy array of all well numbers in a 384 well plate
# Wells have numeric + two digit long number (i.e. )
def generate_wells_384_leading_zero_not_for_picklist():
	well_rows = np.hstack((np.repeat("A",24),
		  np.repeat("B",24),
		  np.repeat("C",24),
		  np.repeat("D",24),
		  np.repeat("E",24),
		  np.repeat("F",24),
		  np.repeat("G",24),
		  np.repeat("H",24),
		  np.repeat("I",24),
		  np.repeat("J",24),
		  np.repeat("K",24),
		  np.repeat("L",24),
		  np.repeat("M",24),
		  np.repeat("N",24),
		  np.repeat("O",24),
		  np.repeat("P",24)))
	well_cols = np.tile(np.arange(24),16)
	well_cols += 1
	well_cols = [str(int(number)).zfill(2) for number in well_cols]
	wells_384 = np.core.defchararray.add(well_rows, well_cols)
	return(wells_384)

def load_plate_map_all_drugs(id_for_plate,plate_map_path,
	fda_deck_metadata="BML-2843 FDA-Approved Library.KI-HTS.384well_singlept.csv"):
	# Loadin the the echo plates, add convert to table formate
	plate = pd.read_csv(plate_map_path,index_col=0)
	wells = []
	plate_drugs = []
	plate_id = [id_for_plate]* plate.shape[0]*plate.shape[1]
	for i in range(plate.shape[0]):
		for j in range(plate.shape[1]):
			wells.append(plate.index.values[i] + plate.columns.values[j])
			plate_drugs.append(plate.iloc[i,j])

	source_plate_map = pd.DataFrame({"plate_id":plate_id,
							 "well":wells,
							 "drug":plate_drugs})

	source_plate_map = source_plate_map.loc[source_plate_map['drug']!='EMPTY']

	fda_metadata = pd.read_csv(fda_deck_metadata)
	compound_names = []
	for drug in source_plate_map['drug']:
		if drug =='DMSO':
			compound_names.append("DMSO")
		elif not drug in fda_metadata['KI-ID'].values:
			compound_names.append(drug)
		else:
			compound_names.append(fda_metadata['Name'].loc[fda_metadata['KI-ID']==drug].values[0])
	source_plate_map['Name']= compound_names

	return(source_plate_map)

##############################################################################
##############################################################################
######################## Writing functions ###################################
##############################################################################
##############################################################################


def write_pool_key_and_dict(plate_number, source_plate_map, well_layout,
	compression_schemes,drugs_in_order,positive_controls,num_DMSO_controls,
	pools_path,key_path,pool_dict_path):

	pool_dict = {} # used to acces the pools when writing picklist
	pool_source = []
	pool_index = []

	# Add in the compression schemes on the plate
	for compression_scheme in compression_schemes:
		print(compression_scheme)
		scheme_type = compression_schemes[compression_scheme]['scheme_type']
		compression_factor = compression_schemes[compression_scheme]['compression_factor']
		num_replicates = compression_schemes[compression_scheme]['num_replicates']

		num_drugs_per_pool = (compression_factor/2)*num_replicates
		num_DSMO_flies = num_drugs_per_pool * 2
		num_pos_control_dmso_flies = num_DSMO_flies - 2

		if scheme_type == "random":
			pools = np.load(pools_path+scheme_type+"/"+str(compression_factor)+"x_"+str(num_replicates)+"r_"+scheme_type+".npy")
		elif scheme_type == "random1":
			pools = np.load(pools_path+scheme_type+"/"+str(compression_factor)+"x_"+str(num_replicates)+"r_"+scheme_type+".npy")
		elif scheme_type == "random2":
			pools = np.load(pools_path+scheme_type+"/"+str(compression_factor)+"x_"+str(num_replicates)+"r_"+scheme_type+".npy")
		elif scheme_type == "fulldeckrandom1":
			pools = np.load(pools_path+scheme_type+"/"+str(compression_factor)+"x_"+str(num_replicates)+"r_"+scheme_type+".npy")
		elif scheme_type == "fulldeckrandom2":
			pools = np.load(pools_path+scheme_type+"/"+str(compression_factor)+"x_"+str(num_replicates)+"r_"+scheme_type+".npy")

		else:# pickle different with the optimized pools
			pools = np.load(pools_path+scheme_type+"/"+str(compression_factor)+"x_"+str(num_replicates)+"r_"+scheme_type+".npy",allow_pickle=True)
			energies = np.array([y for x,y in pools])
			pools = [x for x,y in pools]
			pools = pools[np.where(energies==min(energies))[0][0]] # add in pool with lowest energy


		pool_dict[compression_scheme] = pools
		# add the source of the pools to the pool source list
		pool_source = pool_source + [compression_scheme]*pools.shape[0]
		# add in the index of the pools in the pool matrix
		pool_index = pool_index + np.arange(pools.shape[0]).tolist()

	# Add in the positive controls	
	for positive_control in positive_controls:
		for i in range(5): # adding in five replicates of each positive control
			pool_source.append(positive_control)
			pool_index.append(np.nan)

	# Add in the DMSO wells
	for i in range(num_DMSO_controls):
		pool_source.append("DMSO")
		pool_index.append(np.nan)


	# adding in the destination plate and well to the pool key
	destination_plate = [plate_number]*len(pool_source)
	destination_wells = np.random.choice(well_layout,len(pool_source),replace=False)

	pool_key = pd.DataFrame({"pool_source":pool_source,
						"pool_index":pool_index,
						 "destination_plate":destination_plate,
						 "destination_well":destination_wells
						})
	pool_key.to_csv(key_path+plate_number+".csv")

	with open(pool_dict_path+plate_number+'.pickle', 'wb') as handle:
		pickle.dump(pool_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)


	# lets also save the specific final pool used

	return pool_key, pool_dict,num_DSMO_flies,num_pos_control_dmso_flies


def pool_key_dict_to_picklist(pool_key,pool_dict,num_DSMO_flies,num_pos_control_dmso_flies,plate_number,source_plate_map, well_layout,drugs_in_order,positive_controls,
	num_DMSO_controls,picklist_path):

	# define lists that will become the columns of the picklist
	picklist_source_plate = []
	picklist_source_well = []
	picklist_destination_plate = []
	picklist_destination_well = []
	picklist_transfer_volume = []
	picklist_compound_name = []

	# Loop through the rows of the pool key
	for pool_key_index in range(pool_key.shape[0]):

		# What type of pool is this? Scheme or positive control or DMSO
		pool_source = pool_key.iloc[pool_key_index]['pool_source']
	
		# If that row is a positive control
		if pool_source in positive_controls:
	
			# Add 5 uL of positive control
			picklist_source_plate.append(source_plate_map['plate_id'][source_plate_map['Name']==pool_source].values[0])
			picklist_source_well.append(source_plate_map['well'][source_plate_map['Name']==pool_source].values[0])
			picklist_destination_plate.append(pool_key.iloc[pool_key_index]['destination_plate'])
			picklist_destination_well.append(pool_key.iloc[pool_key_index]['destination_well'])
			picklist_transfer_volume.append(5.0)
			picklist_compound_name.append(pool_source)

			# Choose a random source plate to add DMSO from
			random_plate_id = np.random.choice(np.unique(source_plate_map['plate_id']))
			# Choose random DMSO well in source plate
			temp_plate_map = source_plate_map.loc[source_plate_map['plate_id']==random_plate_id,:]
			dmso_source_well = np.random.choice(temp_plate_map['well'][temp_plate_map['Name']=="DMSO"])

			# add 2.5 uL * num_pos_control_dmso_flies flies DMSO
			picklist_source_plate.append(random_plate_id)
			picklist_source_well.append(dmso_source_well)
			picklist_destination_plate.append(pool_key.iloc[pool_key_index]['destination_plate'])
			picklist_destination_well.append(pool_key.iloc[pool_key_index]['destination_well'])
			picklist_transfer_volume.append(2.5*num_pos_control_dmso_flies)
			picklist_compound_name.append("DMSO")


		elif pool_source=="DMSO":
			# Choose a random source plate to add DMSO from
			random_plate_id = np.random.choice(np.unique(source_plate_map['plate_id']))
			# Choose random DMSO well in source plate
			temp_plate_map = source_plate_map.loc[source_plate_map['plate_id']==random_plate_id,:]
			dmso_source_well = np.random.choice(temp_plate_map['well'][temp_plate_map['Name']=="DMSO"])
			# add 2.5 nL * num_dmso_flies 
			picklist_source_plate.append(random_plate_id )
			picklist_source_well.append(dmso_source_well)
			picklist_destination_plate.append(pool_key.iloc[pool_key_index]['destination_plate'])
			picklist_destination_well.append(pool_key.iloc[pool_key_index]['destination_well'])
			picklist_transfer_volume.append(2.5*num_DSMO_flies)
			picklist_compound_name.append("DMSO")
		# if we got a pool on our hands
		else:
			# get the pool matrix for the scheme at that row
			pool_matrix = pool_dict[pool_source]

			# Find the drugs in that pool
			matrix_row = pool_matrix[int(pool_key.iloc[pool_key_index]['pool_index']),:]
			non_zero_drug_indices = np.nonzero(matrix_row)[0]
			print(drugs_in_order)
			non_zero_drugs = drugs_in_order[non_zero_drug_indices]
			# Add those drugs to the picklistt
			for non_zero_drug in non_zero_drugs:
				picklist_source_plate.append(source_plate_map['plate_id'][source_plate_map['Name']==non_zero_drug].values[0])
				picklist_source_well.append(source_plate_map['well'][source_plate_map['Name']==non_zero_drug].values[0])
				picklist_destination_plate.append(pool_key.iloc[pool_key_index]['destination_plate'])
				picklist_destination_well.append(pool_key.iloc[pool_key_index]['destination_well'])
				picklist_transfer_volume.append(5.0)
				picklist_compound_name.append(non_zero_drug)
	picklist = pd.DataFrame({"source_plate":picklist_source_plate,
							"source_well":picklist_source_well,
							 "destination_plate":picklist_destination_plate,
							"destination_well":picklist_destination_well,
							"transfer_volume":picklist_transfer_volume,
							"compound_name":picklist_compound_name})
	picklist = picklist.sort_values(by=['source_plate'])
	picklist.to_csv(picklist_path+plate_number+".csv")
	return(picklist)
	
def write_well_drug_dict(pool_key,picklist,plate_number,well_drug_path):

    final_dict = {}
    for i in range(picklist.shape[0]):
        dest_well = picklist['destination_well'].iloc[i]
        if dest_well not in final_dict:
            final_dict[dest_well] = [picklist['compound_name'].iloc[i]]
        else:
            final_dict[dest_well].append(picklist['compound_name'].iloc[i])
    with open(well_drug_path+plate_number+'.pickle', 'wb') as handle:
        pickle.dump(final_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)
    return(final_dict)

def write_metadata(pool_key,plate_number,positive_controls,metadata_path):

	# Metadata_Plate = []
	Metadata_Well = []
	Metadata_perturbation = []
	Metadata_compression = []
	Metadata_replicates = []
	Metadata_conc_uM = []
	Metadata_time_hr = []

	for i in range(pool_key.shape[0]):
		# Metadata_Plate.append(pool_key['destination_plate'].iloc[i])

		well = pool_key['destination_well'].iloc[i]

		# add leading zero single digit numbers
		if len(well)==2:
			well = well[0] +"0"+well[1]

		Metadata_Well.append(well)
		

		pool_source = pool_key['pool_source'].iloc[i]
		if not pool_source in positive_controls and not pool_source == "DMSO":
			perturbation = pool_source.split("_")[0]
			Metadata_perturbation.append(perturbation)
			scheme = pool_source.split("_")[1]
			if len(pool_source.split("_")) < 2:
				sys.exit('Need to reformat scheme name')
			if len(scheme)==4:
				compression = scheme[0]
				Metadata_compression.append(int(compression))
				replicates = scheme[2]
				Metadata_replicates.append(int(replicates))
			elif len(scheme) == 5:
				compression = scheme[0:2]
				Metadata_compression.append(int(compression))
				replicates = scheme[3]
				Metadata_replicates.append(int(replicates))
			else:
				sys.exit('Need to reformat scheme name')
		else:
			Metadata_perturbation.append(pool_source)
			Metadata_compression.append(compression)
			Metadata_replicates.append(replicates)

		Metadata_conc_uM.append(1.0)
		Metadata_time_hr.append(24)

	metadata = pd.DataFrame({"Metadata_Well":Metadata_Well,
							"Metadata_perturbation":Metadata_perturbation,
							"Metadata_compression":Metadata_compression,
							"Metadata_replicates":Metadata_replicates,
							"Metadata_conc_uM":Metadata_conc_uM,
							"Metadata_time_hr":Metadata_time_hr})
	metadata.to_csv(metadata_path+plate_number+".csv")
	return(metadata)


# Tester function that makes sure that everything is in order
def check_plate(plate_number, pool_key,pool_dict,picklist,well_drug_dict,metadata,
	compression_schemes,num_DMSO_controls,positive_controls,drugs_in_order):

	print("Testing "+plate_number)

	# Calculate the total number of needed pools
	num_pools = 0
	num_schemes = 0
	for compression_scheme in compression_schemes:
		num_schemes +=1
		scheme_type = compression_schemes[compression_scheme]['scheme_type']
		compression_factor = compression_schemes[compression_scheme]['compression_factor']
		num_replicates = compression_schemes[compression_scheme]['num_replicates']
		num_pools += len(drugs_in_order)/(compression_factor/6)

		# also calculate the useful number of drugs per pool
		num_drugs_per_pool = [(float(compression_factor)/6.0)*float(num_replicates)]
		num_per_pool_if_constant_float = len(drugs_in_order) / (compression_factor/2.0)
		print("is this even?")
		print(num_per_pool_if_constant_float)
		if num_per_pool_if_constant_float != math.ceil(num_per_pool_if_constant_float):
			num_drugs_per_pool.append(num_drugs_per_pool[0] - 1)





	num_pools += num_DMSO_controls
	num_pools += 5*len(positive_controls)
	num_pools = int(math.ceil(num_pools))

	print("total number of compression schemes on plate: "+str(num_schemes))


	# checking that everything has the correct number of pools
	print("Number of expected pools: "+str(num_pools))
	if num_pools != metadata.shape[0] != pool_key.shape[0] != len(well_drug_dict) !=len(np.unique(picklist['destination_well'])):
		print("Pool numbers are not equal between files")
		return(False)


	# calculate the correct number of drugs per pool (picklist and well drug dict)
	lengths = []
	for well in well_drug_dict:
	    lengths.append(len(well_drug_dict[well]))
	num_in_dict = max(set(lengths), key=lengths.count)

	lengths = picklist['destination_well'].value_counts().tolist()
	num_in_picklist = max(set(lengths), key=lengths.count)

	print("Number of expected drugs per pool: ")
	print(num_drugs_per_pool)
	if len(num_drugs_per_pool) == 1:
		if int(num_drugs_per_pool[0]) != num_in_dict != num_in_picklist:
			print("num drugs per pool is messed up")
			return(False)
	else:
		if int(num_drugs_per_pool[0]) != num_in_dict != num_in_picklist:
			if int(num_drugs_per_pool[1]) != num_in_dict != num_in_picklist:
				print("num drugs per pool is messed up")
				return(False)

	# check to make sure that there are the correct number of replicates
	drug_counts_in_picklist = picklist['compound_name'].value_counts().tolist()
	# print(drug_counts_in_picklist)
	num_reps_in_picklist = max(set(drug_counts_in_picklist), key=drug_counts_in_picklist.count)
	# print(num_reps_in_picklist)

	print("Number of expected replicates: "+str(num_replicates))
	if int(num_replicates) != num_reps_in_picklist / num_schemes:
		print("replicate numbers messed up")
		return(False)

	# calculate the needed volumes and flies for things
	if len(num_drugs_per_pool)==1:
		total_volume = [num_drugs_per_pool[0] * 5.0]
		print("Expected total volume: "+str(total_volume[0]))
	else:
		total_volume = [num_drugs_per_pool[0] * 5.0,num_drugs_per_pool[1] * 5.0]

	volumes = picklist.groupby(['destination_well']).sum()['transfer_volume']
	print(volumes.value_counts())
	# if not np.all(volumes == total_volume):
	if not np.sum(np.isin(volumes,total_volume)) == len(volumes):
		print("total volumes are wrong")
		return(False)

	total_flies = num_schemes*(len(drugs_in_order)*num_replicates) + 2*(len(positive_controls)*5) + num_DMSO_controls
	print("Expected total number of echo flies: "+str(total_flies))
	if total_flies != picklist.shape[0]:
		print("total fly number wrong")
		return(False)

	# checking to make sure no repeated wells in pool key
	if len(np.unique(pool_key['destination_well']))!=len(pool_key['destination_well']):
		print("repeated wells in pool key")
		print(len(np.unique(pool_key['destination_well'])))
		print(len(pool_key['destination_well']))
		return(False)
	
	# Checking that all landmarks exist
	# each landmark should have 5 + num_replicates instances
	for positive_control in positive_controls:
		num_for_positive = 5+ num_replicates*num_schemes
		actual_num = picklist.loc[picklist['compound_name']==positive_control].shape[0]
		if actual_num != num_for_positive:
			print("Wrong number of positive controls")
			return(False)

	# Checking that DMSO exists
	num_for_dmso = num_DMSO_controls + 5*len(positive_controls)
	actual_dmso = picklist.loc[picklist['compound_name']=="DMSO"].shape[0]
	if num_for_dmso != actual_dmso:
		print("wrong number of DMSO")
		return(False)

	# Doing the 3 well spot check 

	# choose 3 random wells (not positive or landmark), and make sure that they
	# line up with the well drug dictionary

	for compression_scheme in compression_schemes:

		scheme_type = compression_schemes[compression_scheme]['scheme_type']
		compression_factor = compression_schemes[compression_scheme]['compression_factor']
		num_replicates = compression_schemes[compression_scheme]['num_replicates']
		pool_scheme = pool_key.loc[pool_key['pool_source']==compression_scheme,:]
		test_wells = np.random.choice(pool_scheme['destination_well'].values,3,replace=False)

		for well in test_wells:
			picklist_well = picklist.loc[picklist['destination_well']==well,:]
			picklist_drugs = picklist_well['compound_name'].values

			well_drugs = well_drug_dict[well]
			well_drugs = np.array(well_drugs)

			if not np.array_equal(well_drugs,picklist_drugs):
				print("Picklist and well drugs don't line up")
				return(False)


	return(True)

# Inputs
	# Plate number
	# source plate map (echo plate for this instance), in nice dataframe format
	# well layout (384 wells in this instance)
	# compression scheme dictionary of dictionaries 
		# E.g.
		# random2x3r: {scheme:"random", compression_factor:2, num_replicates:3}
	# Drugs in order
	# Positive controls (list of strings of positive control names)
	# Number of DMSO controls (integer, number of DMSO only control wells)

def write_plate(plate_number, source_plate_map, well_layout,
	compression_schemes,drugs_in_order,positive_controls,num_DMSO_controls,
	pools_path,key_path,pool_dict_path,picklist_path,well_drug_path,metadata_path):

	pool_key, pool_dict,num_DSMO_flies,num_pos_control_dmso_flies = write_pool_key_and_dict(plate_number=plate_number, source_plate_map=source_plate_map, 
		well_layout=well_layout,compression_schemes=compression_schemes,
		drugs_in_order=drugs_in_order,positive_controls=positive_controls,
		num_DMSO_controls=num_DMSO_controls,pools_path=pools_path,key_path=key_path,pool_dict_path=pool_dict_path)

	picklist = pool_key_dict_to_picklist(pool_key=pool_key,pool_dict=pool_dict,plate_number=plate_number,
		num_DSMO_flies=num_DSMO_flies,num_pos_control_dmso_flies=num_pos_control_dmso_flies,
		source_plate_map=source_plate_map, well_layout=well_layout,drugs_in_order=drugs_in_order,
		positive_controls=positive_controls,num_DMSO_controls=num_DMSO_controls,picklist_path=picklist_path)

	well_drug_dict = write_well_drug_dict(pool_key=pool_key,picklist=picklist,plate_number=plate_number,
		well_drug_path=well_drug_path)

	#check that the wells in the well drug dictionary line up with the pool keys

	metadata = write_metadata(pool_key=pool_key,plate_number=plate_number,positive_controls=positive_controls,
		metadata_path=metadata_path)

	plate_is_good = check_plate(plate_number=plate_number,pool_key=pool_key,pool_dict=pool_dict,
		picklist=picklist,well_drug_dict=well_drug_dict,metadata=metadata,
		compression_schemes=compression_schemes,num_DMSO_controls=num_DMSO_controls,
		positive_controls=positive_controls,drugs_in_order=drugs_in_order)

	if not plate_is_good:
		print("^^^^Bummer^^^^")
		sys.exit(1)


	print("wrote "+plate_number)


##############################################################################
##############################################################################
######################## Main ################################################
##############################################################################
##############################################################################


def main():
	wells_384 = generate_wells_384_picklist()

	# Plate maps for the full FDA deck
	fda_plate_map = load_plate_map_all_drugs("KIECHO000655",plate_map_path="../echo_plate_maps/plate1_KIECHO000655.csv")
	
	# Reads in a csv of drugs in a specified order
	# Order is arbitraay but needs to the same in all screens
	drugs_in_order = pd.read_csv("../Optimization/drugs_in_order.txt",sep="\t",header=None).values.flatten()
	
	positive_controls = ['Daunorubicin', 'Doxorubicin', 'Fluvastatin', 
		'Mycophenolic Acid','Riluzole', 'Ticlopidine', 'Tropicamide', 
		'Vinblastine Sulfate']

	# Example for two schemes:
		# Add in more entries to dictionary for more schemes

	plate_schemes = {"plate1":{"random1_2x6r":{"scheme_type":"random1",
											"compression_factor":6,
											"num_replicates":3}},
					"plate2":{"random2_2x6r":{"scheme_type":"random2",
											"compression_factor":6,
											"num_replicates":3}}}
	for plate in plate_schemes_redo_and_higher:
		write_plate(plate_number=plate,source_plate_map=fda_plate_map,
			well_layout=wells_384, compression_schemes=plate_schemes_redo_and_higher[plate],drugs_in_order=drugs_in_order,
			positive_controls=positive_controls,num_DMSO_controls=28,pools_path="../Optimization/pools/",
			key_path="pool_keys/",pool_dict_path='pool_dicts/',picklist_path="picklists/",
			well_drug_path="well_drug_dicts/",metadata_path="metadata/")



if __name__ == '__main__':
	main()