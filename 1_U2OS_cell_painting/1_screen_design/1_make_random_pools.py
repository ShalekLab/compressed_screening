#!/usr/bin/env python3

# Comp Screening - Example code for making random drug pools
## For designing compressed screens
#### CK Nov. 15, 2022_

### What this does: Example function for making random pools
## - Defines initialize_drugs function for making random pools
### - Example use case

import math
import numpy as np
import pandas as pd

# Helper function to make random pools of drugs
def initialize_drugs(n_per_pool, n_drugs, n_replicates):
    
    n_pools = int(math.ceil(n_drugs * n_replicates / n_per_pool))
    
    # Creating a random order to load the drugs into the pools
    drug_sequence = np.repeat(np.arange(n_drugs),n_replicates)
    np.random.shuffle(drug_sequence)

    # Loading the drugs in
    pools = np.zeros((n_pools,n_drugs)) #initial emptry matrix of zeros
    for j in range(len(drug_sequence)):
        no_drug_in_pool = np.where(pools[:,drug_sequence[j]]==0)[0] # finds indices of pools without drug drug_sequence[j] in them
        pool_sums = pools[no_drug_in_pool,:].sum(axis=1)
        minimum_rowsum_indices = np.where(pool_sums==pool_sums.min())[0] #of pools with no drug, find all of the most empty pools
        possible_pool_indices = no_drug_in_pool[minimum_rowsum_indices] # assign those indices as possible pool locations
        pools[np.random.choice(possible_pool_indices),drug_sequence[j]] = 1 # randomly assign drug to one of possible choices
    return pools


# Example: 
total_num_drugs = 316
compression = 6 # relative to 6 replicates GT
n_reps = 3

random1_2x3r = initialize_drugs((compression/6)*n_reps,total_num_drugs,n_reps)
random2_2x3r = initialize_drugs((compression/6)*n_reps,total_num_drugs,n_reps)

np.save("../Optimization/pools/random1/2x_3r_random1.npy",random1_2x3r)
np.save("../Optimization/pools/random2/2x_3r_random2.npy",random2_2x3r)