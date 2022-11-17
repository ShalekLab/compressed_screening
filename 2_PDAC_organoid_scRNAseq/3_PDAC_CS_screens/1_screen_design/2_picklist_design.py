#!/usr/bin/env python3

# CK 2021
# Designing the PDAC compressed screens, makes filesfor running a compressed screens with an Echo acoustic liquid handler.With an Echo, need to tell it to flies drops from source plates to destination plates
# - Inputs
#     - Drug library plate layout csv
#     - Npy files of the random pool layour from make_random_pools.py
# - Outputs for each desintation plate:
#     - picklist: csv input to Echo, tells it what drops come from  which source plates and go to which wells in destination plate
#     - Metadata: csv that records for each well in the desintation plate what the dose, time, comrpession, replicates, etc are
#     - Well_drug_dict: dictionary describing wells in destation plate, keys are wells, values are lists of the perturbations in that well
#     - Pool_key: intermediate file used for building above, not used in downstream analysis

import sys
import numpy as np
import pandas as pd
import pickle
import math
import string


##############################################################################
##############################################################################
######################## Loading functions ###################################
##############################################################################
##############################################################################


def get_source_plate_wells_for_drugs(source_plate_layout):
    perturbations = []
    wells=[]
    for i in range(source_plate_layout.shape[0]):
        for j in range(source_plate_layout.shape[1]):
            if source_plate_layout.iloc[i,j] !="EMPTY":
                perturbations.append(source_plate_layout.iloc[i,j])
                wells.append(source_plate_layout.index.values[i] + source_plate_layout.columns.values[j])
    return(pd.DataFrame({"perturbation":perturbations,"well":wells}))

def make_source_plate_map(fly_design,screen,source_plate_id,source_plate_wells):
    plate_ids = []
    drugs = []
    names = []
    amount_used = []
    num_flies = []
    wells = []
    for i in range(fly_design.shape[0]):
        plate_ids.append(source_plate_id)
        drugs.append(fly_design['ID'].iloc[i])
        names.append(fly_design['Ligand'].iloc[i])
        num_flies.append(fly_design["# of ligand flys"].iloc[i])
        amount_used.append(0)
        wells.append(source_plate_wells['well'].loc[source_plate_wells['perturbation']==fly_design['Ligand'].iloc[i]].values[0])
    source_plate_map = pd.DataFrame({'plate_id':plate_ids,
                                    'well':wells,
                                    'ID':drugs,
                                    'Name':names,
                                    'num_flies':num_flies,
                                    'amount_used':amount_used})
    source_plate_map.to_csv("source_plate_maps/source_plate_map_"+screen+".csv")

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

# generate_wells_384_no_outer_wells_picklist
def generate_well_384_no_outer_wells_picklist():
    well_rows = np.hstack((np.repeat("B",22),
          np.repeat("C",22),
          np.repeat("D",22),
          np.repeat("E",22),
          np.repeat("F",22),
          np.repeat("G",22),
          np.repeat("H",22),
          np.repeat("I",22),
          np.repeat("J",22),
          np.repeat("K",22),
          np.repeat("L",22),
          np.repeat("M",22),
          np.repeat("N",22),
          np.repeat("O",22)))
    well_cols = np.tile(np.arange(1,23),14)
    well_cols += 1
    well_cols = [str(int(number)) for number in well_cols]
    wells_384_no_outer = np.core.defchararray.add(well_rows, well_cols)
    return(wells_384_no_outer)


# Creates lists defining the quadrant positions
def create_quadrant_positions():
    g1 = []
    g2 = []
    g3 = []
    g4 = []
    for i in [0,2,4,6,8,10,12,14]:
        for j in [0,2,4,6,8,10,12,14,16,18,20,22]:
            g1.append(i*24 + j)
    for i in [0,
              2,4,6,8,10,12,14]:
        for j in [1,3,5,7,9,11,13,15,17,19,21,23]:
            g2.append(i*24 + j)
    for i in [1,3,5,7,9,11,13,15]:
        for j in [0,2,4,6,8,10,12,14,16,18,20,22]:
            g3.append(i*24 + j)
    for i in [1,3,5,7,9,11,13,15]:
        for j in [1,3,5,7,9,11,13,15,17,19,21,23]:
            g4.append(i*24 + j)
    g1 = np.array(g1)
    g2 = np.array(g2)
    g3 = np.array(g3)
    g4 = np.array(g4)
    
    return g1, g2, g3, g4

# generate wells 384 quadrant 1
def generate_wells_384_quadrant1_picklist():
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
    quadrant1,quadrant2,quadrant3,quadrant4 = create_quadrant_positions()
    return(wells_384[quadrant1])

    # generate wells 384 quadrant 2
def generate_wells_384_quadrant2_picklist():
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
    quadrant1,quadrant2,quadrant3,quadrant4 = create_quadrant_positions()
    return(wells_384[quadrant2])


def generate_wells_384_quadrant_picklist_no_outer_well_b2_is_q4(quadrant):
    wells_384 = generate_wells_384_picklist()
    quadrant1,quadrant2,quadrant3,quadrant4 = create_quadrant_positions()
    
    if quadrant == "quadrant1":
        quad = wells_384[quadrant1]
    if quadrant == "quadrant2":
        quad = wells_384[quadrant2]
    if quadrant == "quadrant3":
        quad = wells_384[quadrant3]
    if quadrant == "quadrant4":
        quad = wells_384[quadrant4]
    
    final_wells = []
    for well in quad:
        if well[0]!="A" and well[0]!="P":
            if len(well) ==2:
                if well[1]!="1":
                    final_wells.append(well)
            else:
                if well[1:3]!="24":
                    final_wells.append(well)
    return(np.array(final_wells))

# Creates lists defining the quadrant positions
def create_quadrant_positions_no_outer_wells_b2_is_q1():
    g1 = []
    g2 = []
    g3 = []
    g4 = []
    for i in [0,2,4,6,8,10,12]:
        for j in [0,2,4,6,8,10,12,14,16,18,20]:
            g1.append(i*22 + j)
    for i in [0,2,4,6,8,10,12]:
        for j in [1,3,5,7,9,11,13,15,17,19,21]:
            g2.append(i*22 + j)
    for i in [1,3,5,7,9,11,13]:
        for j in [0,2,4,6,8,10,12,14,16,18,20]:
            g3.append(i*22 + j)
    for i in [1,3,5,7,9,11,13]:
        for j in [1,3,5,7,9,11,13,15,17,19,21]:
            g4.append(i*22 + j)
    g1 = np.array(g1)
    g2 = np.array(g2)
    g3 = np.array(g3)
    g4 = np.array(g4)
    
    return g1, g2, g3, g4

def generate_wells_384_quadrant_picklist_no_outer_well_b2_is_q1(quadrant):
    wells_308 = generate_well_384_no_outer_wells_picklist()
    quadrant1,quadrant2,quadrant3,quadrant4 = create_quadrant_positions_no_outer_wells_b2_is_q1()
    
    if quadrant == "quadrant1":
        return(wells_308[quadrant1])
    if quadrant == "quadrant2":
        return(wells_308[quadrant2])
    if quadrant == "quadrant3":
        return(wells_308[quadrant3])
    if quadrant == "quadrant4":
        return(wells_308[quadrant4])

# Loading the in the echo source plate map in grid format and
# reshape it into a more useful dataframe

# Todo write a load plate map function
    # will include the number of flies for that ligand in the plate map

##############################################################################
##############################################################################
######################## Writing functions ###################################
##############################################################################
##############################################################################


def write_pool_key_and_dict(plate_number, source_plate_map, well_layout,
    compression_schemes,drugs_in_order,positive_controls,num_positive_control_replicates,num_vehicle_controls, vehicle_name,
    num_media_controls, pools_path,key_path,pool_dict_path):

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
        num_vehicle_flies = num_drugs_per_pool * 2
        num_pos_control_vehicle_flies = num_vehicle_flies - 2



        if scheme_type == "random":
            pools = np.load(pools_path+scheme_type+"/"+str(compression_factor)+"x_"+str(num_replicates)+"r_"+scheme_type+".npy")
        elif scheme_type in ["randomPDAC1","randomPDAC2"]:
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
        for i in range(num_positive_control_replicates): # adding in five replicates of each positive control
            pool_source.append(positive_control)
            pool_index.append(np.nan)

    # Add in the vehicle control wells
    for i in range(num_vehicle_controls):
        pool_source.append(vehicle_name)
        pool_index.append(np.nan)


    # Add in the media only wells
    for i in range(num_media_controls):
        pool_source.append("media")
        pool_index.append(np.nan)


    # adding in the destination plate and well to the pool key
    destination_plate = [plate_number]*len(pool_source)
    # print("here")
    # print(len(pool_source))
    # print(len(pool_source[pool_source!="media"]))
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

    return pool_key, pool_dict,num_vehicle_flies,num_pos_control_vehicle_flies


def pool_key_dict_to_picklist(pool_key,pool_dict,num_vehicle_flies,num_pos_control_vehicle_flies,plate_number,source_plate_map, well_layout,drugs_in_order,positive_controls,
    num_vehicle_controls,vehicle_name,num_media_controls,picklist_path):

    # define lists that will become the columns of the picklist
    picklist_source_plate = []
    picklist_source_well = []
    picklist_destination_plate = []
    picklist_destination_well = []
    picklist_transfer_volume = []
    picklist_compound_name = []

    # Loop through the rows of the pool key
    for pool_key_index in range(pool_key.shape[0]):

        # What type of pool is this? Scheme or positive control or vehicle
        pool_source = pool_key.iloc[pool_key_index]['pool_source']
    
        # If that row is a positive control
        if pool_source in positive_controls:

            # find the wells with the least used ligand, randomly use one of those
            temp_plate_map = source_plate_map.loc[source_plate_map['Name']==pool_source,:]
            min_amount_used_indices = np.argwhere(temp_plate_map['amount_used'].values == np.amin(temp_plate_map['amount_used'].values)).flatten()
            use_well_index = np.random.choice(min_amount_used_indices)
            use_well = temp_plate_map['well'].iloc[use_well_index]
            # will need to add in the amount used


            # Add 5 uL of positive control
            picklist_source_plate.append(source_plate_map['plate_id'][source_plate_map['well']==use_well].values[0])
            picklist_source_well.append(use_well)
            picklist_destination_plate.append(pool_key.iloc[pool_key_index]['destination_plate'])
            picklist_destination_well.append(pool_key.iloc[pool_key_index]['destination_well'])
            
            # Figure out how much volume, add to picklist and track in source plate map
            amount_used = 2.5*source_plate_map['num_flies'][source_plate_map['well']==use_well].values[0]
            picklist_transfer_volume.append(amount_used)
            source_plate_map['amount_used'].loc[source_plate_map['well']==use_well] += amount_used
            picklist_compound_name.append(pool_source)

        # Removing the steps adding extra vehicle control for the ligand screen

            # # Choose a random source plate to add vehicle from
            # random_plate_id = np.random.choice(np.unique(source_plate_map['plate_id']))
            # # Choose random vehicle well in source plate
            # temp_plate_map = source_plate_map.loc[source_plate_map['plate_id']==random_plate_id,:]
            # vehicle_source_well = np.random.choice(temp_plate_map['well'][temp_plate_map['Name']=="vehicle"])

            # # add 2.5 uL * num_pos_control_vehicle_flies flies vehicle
            # picklist_source_plate.append(random_plate_id)
            # picklist_source_well.append(vehicle_source_well)
            # picklist_destination_plate.append(pool_key.iloc[pool_key_index]['destination_plate'])
            # picklist_destination_well.append(pool_key.iloc[pool_key_index]['destination_well'])
            # picklist_transfer_volume.append(2.5*num_pos_control_vehicle_flies)
            # picklist_compound_name.append("vehicle")


        # elif pool_source=="vehicle":
        #   # Choose a random source plate to add vehicle from
        #   random_plate_id = np.random.choice(np.unique(source_plate_map['plate_id']))
        #   # Choose random vehicle well in source plate
        #   temp_plate_map = source_plate_map.loc[source_plate_map['plate_id']==random_plate_id,:]
        #   vehicle_source_well = np.random.choice(temp_plate_map['well'][temp_plate_map['Name']=="vehicle"])
        #   # add 2.5 nL * num_vehicle_flies 
        #   picklist_source_plate.append(random_plate_id )
        #   picklist_source_well.append(vehicle_source_well)
        #   picklist_destination_plate.append(pool_key.iloc[pool_key_index]['destination_plate'])
        #   picklist_destination_well.append(pool_key.iloc[pool_key_index]['destination_well'])
        #   picklist_transfer_volume.append(2.5*num_vehicle_flies)
        #   picklist_compound_name.append("vehicle")
        # if we got a pool on our hands
        elif pool_source!='media':

            # Todo: set up the tracking of how much of the drugs we are using

            # get the pool matrix for the scheme at that row
            pool_matrix = pool_dict[pool_source]

            # Find the drugs in that pool
            matrix_row = pool_matrix[int(pool_key.iloc[pool_key_index]['pool_index']),:]
            non_zero_drug_indices = np.nonzero(matrix_row)[0]
            non_zero_drugs = drugs_in_order[non_zero_drug_indices]
            # Add those drugs to the picklistt
            for non_zero_drug in non_zero_drugs:

                # find the wells with the least used ligand, randomly use one of those
                temp_plate_map = source_plate_map.loc[source_plate_map['Name']==non_zero_drug,:]
                min_amount_used_indices = np.argwhere(temp_plate_map['amount_used'].values == np.amin(temp_plate_map['amount_used'].values)).flatten()
                use_well_index = np.random.choice(min_amount_used_indices)
                use_well = temp_plate_map['well'].iloc[use_well_index]

                picklist_source_plate.append(source_plate_map['plate_id'][source_plate_map['well']==use_well].values[0])
                picklist_source_well.append(use_well)
                picklist_destination_plate.append(pool_key.iloc[pool_key_index]['destination_plate'])
                picklist_destination_well.append(pool_key.iloc[pool_key_index]['destination_well'])
                
                amount_used = 2.5*source_plate_map['num_flies'][source_plate_map['well']==use_well].values[0]
                picklist_transfer_volume.append(amount_used)
                source_plate_map['amount_used'].loc[source_plate_map['well']==use_well] += amount_used
                picklist_compound_name.append(non_zero_drug)

    picklist = pd.DataFrame({"source_plate":picklist_source_plate,
                            "source_well":picklist_source_well,
                             "destination_plate":picklist_destination_plate,
                            "destination_well":picklist_destination_well,
                            "transfer_volume":picklist_transfer_volume,
                            "compound_name":picklist_compound_name})
    picklist = picklist.sort_values(by=['source_plate'])
    picklist.to_csv(picklist_path+plate_number+".csv")
    return(picklist, source_plate_map)
    
def write_well_drug_dict(pool_key,picklist,plate_number,well_drug_path):

    # Go through the picklist and add the drugs to the well drug dictionary
    final_dict = {}
    for i in range(picklist.shape[0]):
        dest_well = picklist['destination_well'].iloc[i]
        if dest_well not in final_dict:
            final_dict[dest_well] = [picklist['compound_name'].iloc[i]]
        else:
            final_dict[dest_well].append(picklist['compound_name'].iloc[i])
    # add in the media only wells from the pool_key
    media_only_key = pool_key.loc[pool_key['pool_source']=="media",:]
    for i in range(media_only_key.shape[0]):
        final_dict[media_only_key['destination_well'].iloc[i]] = []

    with open(well_drug_path+plate_number+'.pickle', 'wb') as handle:
        pickle.dump(final_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)
    return(final_dict)

def write_metadata(pool_key,plate_number,positive_controls,metadata_path,vehicle_name):

    # Metadata_Plate = []
    Metadata_Well = []
    Metadata_perturbation = []
    Metadata_compression = []
    Metadata_replicates = []

    for i in range(pool_key.shape[0]):
        # Metadata_Plate.append(pool_key['destination_plate'].iloc[i])

        well = pool_key['destination_well'].iloc[i]

        # add leading zero single digit numbers
        if len(well)==2:
            well = well[0] +"0"+well[1]

        Metadata_Well.append(well)

        pool_source = pool_key['pool_source'].iloc[i]
        if not pool_source in positive_controls and not pool_source == "media" and not pool_source==vehicle_name:
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
            # If all of the wells are media then just manually set the compression
            if pool_key.loc[pool_key['pool_source']=="media"].shape[0] == pool_key.shape[0]:
                compression = 0
                replicates = 0

            Metadata_perturbation.append(pool_source)
            Metadata_compression.append(compression)
            Metadata_replicates.append(replicates)


    metadata = pd.DataFrame({"Metadata_Well":Metadata_Well,
                            "Metadata_perturbation":Metadata_perturbation,
                            "Metadata_compression":Metadata_compression,
                            "Metadata_replicates":Metadata_replicates})
    metadata.to_csv(metadata_path+plate_number+".csv")
    return(metadata)


# Tester function that makes sure that everything is in order
def check_plate(plate_number, pool_key,pool_dict,picklist,well_drug_dict,metadata,
    compression_schemes,num_vehicle_controls,vehicle_name,num_media_controls,positive_controls,num_positive_control_replicates,drugs_in_order):

    print("Testing "+plate_number)

    # Calculate the total number of needed pools
    num_pools = 0
    num_schemes = 0
    for compression_scheme in compression_schemes:
        num_schemes +=1
        scheme_type = compression_schemes[compression_scheme]['scheme_type']
        compression_factor = compression_schemes[compression_scheme]['compression_factor']
        num_replicates = compression_schemes[compression_scheme]['num_replicates']
        if compression_factor !=0:
            num_pools += len(drugs_in_order)/(compression_factor/2)

            # also calculate the useful number of drugs per pool
            num_drugs_per_pool = [(float(compression_factor)/2.0)*float(num_replicates)]
            num_per_pool_if_constant_float = len(drugs_in_order) / (compression_factor/2.0)
            if num_per_pool_if_constant_float != math.ceil(num_per_pool_if_constant_float):
                num_drugs_per_pool.append(num_drugs_per_pool[0] - 1)

    if compression_factor != 0:
        num_pools += num_vehicle_controls
        num_pools += num_positive_control_replicates*len(positive_controls)
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

        if plate_number == "randomHEPG2":
            if int(num_replicates) != (num_reps_in_picklist - num_positive_control_replicates)  / num_schemes:
                print(num_reps_in_picklist)
                print(num_schemes)
                print("replicate numbers messed up")
                return(False)
        elif int(num_replicates) != (num_reps_in_picklist)  / num_schemes:
            print(num_reps_in_picklist)
            print(num_schemes)
            print("replicate numbers messed up")
            return(False)

        # calculate the needed volumes and flies for things
        # if len(num_drugs_per_pool)==1:
        #   total_volume = [num_drugs_per_pool[0] * 5.0]
        #   print("Expected total volume: "+str(total_volume[0]))
        # else:
        #   total_volume = [num_drugs_per_pool[0] * 5.0,num_drugs_per_pool[1] * 5.0]

        # volumes = picklist.groupby(['destination_well']).sum()['transfer_volume']
        # print(volumes.value_counts())
        # # if not np.all(volumes == total_volume):
        # if not np.sum(np.isin(volumes,total_volume)) == len(volumes):
        #   print("total volumes are wrong")
        #   return(False)


        # total_flies = num_schemes*(len(drugs_in_order)*num_replicates) + 2*(len(positive_controls)*num_positive_control_replicates) + num_vehicle_controls
        # print("Expected total number of echo flies: "+str(total_flies))
        # if total_flies != picklist.shape[0]:
        #   print("total fly number wrong")
        #   return(False)

        # checking to make sure no repeated wells in pool key
        if len(np.unique(pool_key['destination_well']))!=len(pool_key['destination_well']):
            print("repeated wells in pool key")
            print(len(np.unique(pool_key['destination_well'])))
            print(len(pool_key['destination_well']))
            return(False)
        
        # Checking that all landmarks exist
        # each landmark should have num_positive_control_replicates + num_replicates instances
        for positive_control in positive_controls:
            num_for_positive = num_positive_control_replicates+ num_replicates*num_schemes
            actual_num = picklist.loc[picklist['compound_name']==positive_control].shape[0]
            if actual_num != num_for_positive:
                print("Wrong number of positive controls")
                return(False)

        # Checking that vehicle exists
        if num_vehicle_controls == 0:
            num_for_vehicle = 0
        else:
            num_for_vehicle = num_vehicle_controls + num_positive_control_replicates*len(positive_controls)
        actual_vehicle = picklist.loc[picklist['compound_name']==vehicle_name].shape[0]
        if num_for_vehicle != actual_vehicle:
            print("wrong number of vehicle")
            return(False)

        # checking that the media only wells exist in the pool key
        media_only_key = pool_key.loc[pool_key['pool_source']=="media"]
        if media_only_key.shape[0] != num_media_controls:
            print(media_only_key.shape[0])
            print(num_media_controls)
            print("wrong number of media only wells")
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
    # Number of vehicle controls (integer, number of vehicle only control wells)

def write_plate(plate_number, source_plate_map, well_layout,
    compression_schemes,drugs_in_order,positive_controls,num_positive_control_replicates,num_vehicle_controls,vehicle_name,num_media_controls,
    pools_path,key_path,pool_dict_path,picklist_path,well_drug_path,metadata_path,source_plate_map_save_path):

    pool_key, pool_dict,num_vehicle_flies,num_pos_control_vehicle_flies = write_pool_key_and_dict(plate_number=plate_number, source_plate_map=source_plate_map, 
        well_layout=well_layout,compression_schemes=compression_schemes,
        drugs_in_order=drugs_in_order,positive_controls=positive_controls,num_positive_control_replicates=num_positive_control_replicates,
        num_vehicle_controls=num_vehicle_controls,vehicle_name=vehicle_name,num_media_controls=num_media_controls,
        pools_path=pools_path,key_path=key_path,pool_dict_path=pool_dict_path)

    picklist, source_plate_map_with_amount_used = pool_key_dict_to_picklist(pool_key=pool_key,pool_dict=pool_dict,plate_number=plate_number,
        num_vehicle_flies=num_vehicle_flies,num_pos_control_vehicle_flies=num_pos_control_vehicle_flies,
        source_plate_map=source_plate_map, well_layout=well_layout,drugs_in_order=drugs_in_order,
        positive_controls=positive_controls,num_vehicle_controls=num_vehicle_controls,
        vehicle_name=vehicle_name,num_media_controls=num_media_controls,
        picklist_path=picklist_path)

    source_plate_map_with_amount_used.to_csv(source_plate_map_save_path)

    well_drug_dict = write_well_drug_dict(pool_key=pool_key,picklist=picklist,plate_number=plate_number,
        well_drug_path=well_drug_path)

    #check that the wells in the well drug dictionary line up with the pool keys

    metadata = write_metadata(pool_key=pool_key,plate_number=plate_number,positive_controls=positive_controls,
        metadata_path=metadata_path,vehicle_name=vehicle_name)

    plate_is_good = check_plate(plate_number=plate_number,pool_key=pool_key,pool_dict=pool_dict,
        picklist=picklist,well_drug_dict=well_drug_dict,metadata=metadata,
        compression_schemes=compression_schemes,num_vehicle_controls=num_vehicle_controls,
        vehicle_name=vehicle_name,num_media_controls=num_media_controls,
        positive_controls=positive_controls,num_positive_control_replicates=num_positive_control_replicates,
        drugs_in_order=drugs_in_order)

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

    # Specifying the name for the ligand source plate
    final_source_plate_id = "ligand_october2021"

    # Loading in the number of flies that is needed for each screen
    fly_design_for_screens = pd.read_csv("fly_design_october2021.csv")
    print(fly_design_for_screens)
    # ligands = load the plate map (will be different bc num flies different for mac/hep vs pdac/glio)
    # Columns that wil be needed
        # num_files - the number of 2.5nL flies
        # amount_used - 0 for all 

    ##########################################################################
    ###############################     PDAC V2    ###########################
    ##########################################################################

    source_plate_layout = pd.read_csv("source_plate_ligands_oct2021.csv",index_col=0)
    source_plate_wells = get_source_plate_wells_for_drugs(source_plate_layout)

    make_source_plate_map(fly_design=fly_design_for_screens, screen="PDAC_V2",source_plate_id = final_source_plate_id,
        source_plate_wells=source_plate_wells)
    
    # quadrant1
    wells = generate_wells_384_quadrant1_picklist()
    positive_controls = ["TGFB1","IFNG","TNF"]
    reps_for_each_positive_control = 4

    ligands_plate_map= pd.read_csv("source_plate_maps/source_plate_map_PDAC_V2.csv")
    ligand_map_save_path = "source_plate_maps/source_plate_map_PDAC_V2_used.csv"

    drugs_in_order = np.unique(ligands_plate_map['Name'])
    drugs_in_order = drugs_in_order[drugs_in_order!="EMPTY"]
    pools_read_path = "pools/"
    media_control_replicate_number = 12

    plate_schemes = {'randomPDAC1':{'random1_4x5r':{"scheme_type":'randomPDAC1',
                                'compression_factor':4,
                                'num_replicates':5}}}

    for plate in plate_schemes:
        write_plate(plate_number=plate,source_plate_map=ligands_plate_map,
            well_layout=wells, compression_schemes=plate_schemes[plate],
            drugs_in_order=drugs_in_order,
            positive_controls=positive_controls,num_positive_control_replicates = reps_for_each_positive_control,
            num_vehicle_controls=0,vehicle_name='NA',num_media_controls=media_control_replicate_number,
            pools_path=pools_read_path,key_path="pool_keys/",pool_dict_path='pool_dicts/',
            picklist_path="picklists/",well_drug_path="well_drug_dicts/",metadata_path="metadata/",
            source_plate_map_save_path = ligand_map_save_path)


    # quadrant 2
    wells = generate_wells_384_quadrant2_picklist()
    positive_controls = ["TGFB1","IFNG","TNF"]
    reps_for_each_positive_control = 4

    ligands_plate_map = pd.read_csv("source_plate_maps/source_plate_map_PDAC_V2_used.csv")
    ligand_map_save_path = "source_plate_maps/source_plate_map_PDAC_used2.csv"

    drugs_in_order = np.unique(ligands_plate_map['Name'])
    drugs_in_order = drugs_in_order[drugs_in_order!="EMPTY"]
    pools_read_path = "pools/"
    media_control_replicate_number = 12
    plate_schemes = {'randomPDAC2':{'random2_4x5r':{"scheme_type":'randomPDAC2',
                                'compression_factor':4,
                                'num_replicates':5}}}
    for plate in plate_schemes:
        write_plate(plate_number=plate,source_plate_map=ligands_plate_map,
            well_layout=wells, compression_schemes=plate_schemes[plate],
            drugs_in_order=drugs_in_order,
            positive_controls=positive_controls,num_positive_control_replicates = reps_for_each_positive_control,
            num_vehicle_controls=0,vehicle_name='NA',num_media_controls=media_control_replicate_number,
            pools_path=pools_read_path,key_path="pool_keys/",pool_dict_path='pool_dicts/',
            picklist_path="picklists/",well_drug_path="well_drug_dicts/",metadata_path="metadata/",
            source_plate_map_save_path = ligand_map_save_path)

    # concatenating the two parts together
    final_name = "randomPDAC"

    # picklist
    picklist1 = pd.read_csv("picklists/"+final_name+"1.csv",index_col=0)
    picklist2 = pd.read_csv("picklists/"+final_name+"2.csv",index_col=0)
    picklist = pd.concat([picklist1,picklist2])
    picklist['destination_plate'] = final_name
    picklist.to_csv('picklists/'+final_name+'_final.csv')

    # pool_key
    poolkey1 = pd.read_csv("pool_keys/"+final_name+"1.csv",index_col=0)
    poolkey2 = pd.read_csv("pool_keys/"+final_name+"2.csv",index_col=0)
    poolkey = pd.concat([poolkey1,poolkey2])
    poolkey['final_destination_plate'] = final_name
    poolkey.to_csv('pool_keys/'+final_name+'_final.csv')

    # metadata
    metadata1 = pd.read_csv("metadata/"+final_name+"1.csv",index_col=0)
    metadata1['quadrant'] = 1
    metadata2 = pd.read_csv("metadata/"+final_name+"2.csv",index_col=0)
    metadata2['quadrant'] = 2
    metadata = pd.concat([metadata1,metadata2])
    metadata.to_csv("metadata/"+final_name+"_final.csv")

    # well_drug_dictionary
    with open('well_drug_dicts/'+final_name+'1.pickle', 'rb') as handle:
        wd1 = pickle.load(handle)
    with open('well_drug_dicts/'+final_name+'2.pickle', 'rb') as handle:
        wd2 = pickle.load(handle)
    final_well_drug_dict = {}
    final_well_drug_dict['quadrant1'] = wd1
    final_well_drug_dict['quadrant2'] = wd2
    with open('well_drug_dicts/'+final_name+'_final.pickle', 'wb') as handle:
        pickle.dump(final_well_drug_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)











if __name__ == '__main__':
    main()