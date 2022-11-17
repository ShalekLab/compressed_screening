This directory contains example code for designing a compressed screen. It is currently bespoke to the screening setup at the Koch Institute of MIT, but the general structure of our code, as well as many of the functions we created may be helpful to a user looking to conduct a compressed screen.

**1_making_CS_pool_design_matrices.ipynb**: Script for making the design matrices for the CS pools
    - Initialize drug function makes the pools
        - Input the compressed scheme design
        - Outputs .npy file containing numpy array of the pool design

**2_picklist_design.py**: Designing the PDAC compressed screens, makes filesfor running a compressed screens with an Echo acoustic liquid handler.With an Echo, need to tell it to flies drops from source plates to destination plates
- Inputs
    - Drug library plate layout csv
    - Npy files of the random pool layour from make_random_pools.py
- Outputs for each desintation plate:
    - picklist: csv input to Echo, tells it what drops come from  which source plates and go to which wells in destination plate
    - Metadata: csv that records for each well in the desintation plate what the dose, time, comrpession, replicates, etc are
    - Well_drug_dict: dictionary describing wells in destation plate, keys are wells, values are lists of the perturbations in that well
    - Pool_key: intermediate file used for building above, not used in downstream analysis
