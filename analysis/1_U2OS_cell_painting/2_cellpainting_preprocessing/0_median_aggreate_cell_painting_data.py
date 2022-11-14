### Script for median aggregating cell painting features across all cells in a given well

import os
import sqlite3
import pandas as pd
import numpy as np
import multiprocessing as mp
from itertools import repeat

def Median_by_Well(well, table, string, connect_path):
    
    well_query = ' and '+table+'_Metadata_Well = "'+well+'"'
    query = string+well_query
    
    feature_table = pd.read_sql(query, con=sqlite3.connect(connect_path))
    
    n = len(feature_table)
    well_data = feature_table.median()

    well_data = pd.DataFrame([well_data], columns= well_data.index.tolist())
    well_data.insert(0, "Number_of_Cells", n, True)
    well_data.insert(0, "Metadata_Well", well, True)  

    return well_data

def Build_QCfiltered_median_feature_table(table, QC_fail, plate, path_to_data):

    # Step 1: Load features for compartment table
    string = "pragma table_info("+table+")"
    
    connect_path = path_to_data+plate+'.sqlite'
    
    all_features = list(pd.read_sql(string, con=sqlite3.connect(connect_path))['name'])

    features = []
    metadata = []
    for f in all_features:
        if  table in f and \
            ('_Number_' not in f) and \
            ('_Parent_' not in f) and \
            ('_Metadata_' not in f) and \
            ('_Center_' not in f) and \
            ('_Location_' not in f):
            features.append(f)
        if 'Metadata' in f and 'QC' not in f:
            metadata.append(f)

    features.sort()
    metadata.sort()

    # Step 2: Get QC failed images
    fail_site = list(QC_fail.loc[QC_fail.Metadata_Plate == plate, 'ImageNumber_Original'].values)
    
    # Step 3: Build data table of features per object in non-failed wells
    string = 'select ' + ','.join(metadata) + ',' + ','.join(features) + \
             ' from ' + table + \
             ' where ImageNumber not in (' + ','.join([str(i) for i in fail_site]) + ')'
    
    feature_table = pd.read_sql(string, con=sqlite3.connect(connect_path))

    # Step 4: Median aggregate over cells and all features per well, count total cells / well
    wells = list(set(feature_table[table + '_Metadata_Well']))
    median_feature_table = pd.DataFrame()
    
    # run in parallel over wells    
    if len(wells) > (mp.cpu_count()-1):
        p_size = (mp.cpu_count()-1)
    else:
        p_size = len(wells)
    
    pool = mp.Pool(p_size)
    median_feature_table = pd.concat(pool.starmap(Median_by_Well, zip(wells, 
                                                                      repeat(table), 
                                                                      repeat(string),
                                                                      repeat(connect_path))), ignore_index=True)
    pool.close()
    pool.join()

    del median_feature_table[table + '_Metadata_Site']    
    median_feature_table.insert(0, "Metadata_Plate", plate, True)
    
    return median_feature_table
        
def Aggregate(path_to_data, plate, QC_fail, platemaps, platemap_dir, batch):
    
    # Step 1: Load QC, platemaps    
    platemap = platemaps['Plate_Map_Name'][platemaps.Assay_Plate_Barcode == plate].iloc[0]
    metadata_table = pd.read_csv(path_to_data + platemap_dir + platemap + '.csv')
    
    # Step 2: Aggregate by compartment, merge with metadata
    Cells_table = Build_QCfiltered_median_feature_table('Cells', QC_fail, plate, path_to_data)
    Cyto_table = Build_QCfiltered_median_feature_table('Cytoplasm', QC_fail, plate, path_to_data)
    
    merge_feature_table = Cells_table.merge(Cyto_table, on=['Metadata_Plate','Metadata_Well','Number_of_Cells'])
    
    Nuc_table = Build_QCfiltered_median_feature_table('Nuclei', QC_fail, plate, path_to_data)
    
    merge_feature_table = merge_feature_table.merge(Nuc_table, on=['Metadata_Plate','Metadata_Well','Number_of_Cells'])
    
    merge_feature_table = merge_feature_table.merge(metadata_table, on=['Metadata_Well'])
    
    # Step 3: Clean up order of features in final DF
    metadata_cols = [col for col in merge_feature_table.columns if 'Metadata' in col]
    
    feature_cols = [col for col in \
                    merge_feature_table.columns if 'Metadata' not in col and \
                    'Number_of_Cells' not in col]
    
    final_col_order = metadata_cols+['Number_of_Cells']+feature_cols
    
    merge_feature_table = merge_feature_table[final_col_order]
    
    # Save median'd values
    print('Saving %s plate to Output/median_plates' %plate)
    save = 'Output/median_plates/QC_median_aggragated_'+batch+'_'+plate+'_feature_table.gz'
    merge_feature_table.to_csv(save, index=False, compression='gzip')

def Batch_Aggregate(path_to_data, QC_fail, platemaps, platemap_dir, batch, date):
    
    ## expects all inputs as strings (file names w/o extension)
    
    ## expects following dir structure:
    # -- Working_dir
    #   -- Aggregate.py script (this file)
    #   -- Output (where final tables saved)
    #       -- median_plates
    # -- 0_RAW_data
    #   -- QC_fail.csv
    #   -- platemaps.csv
    #   -- n x plate.sqlite
    #   -- platemap
    #     -- n x platemap.csv 
    
    # Step 1: Load QC, platemaps to aggregate by batch
    QC_fail = pd.read_csv(path_to_data+QC_fail+'.csv')
    platemaps = pd.read_csv(path_to_data+platemaps+'.csv')
    
    # Step 2: Aggregate each plate in batch and append to final data table
    plates = list(set(platemaps['Assay_Plate_Barcode']))

    print('Aggregating %d plates ' %len(plates) + 'for %s' %batch)
    
    for p in plates:
        Aggregate(path_to_data, p, QC_fail, platemaps, platemap_dir, batch)
        
    print('Merging %d plates ' %len(plates) + 'for %s' %batch)
    batch_feature_table = pd.DataFrame()     
    
    for p in plates:
        file = 'Output/median_plates/QC_median_aggragated_'+batch+'_'+p+'_feature_table.gz'
        plate = pd.read_csv(file,low_memory=False)
        batch_feature_table = batch_feature_table.append(plate, ignore_index = True)
    
    batch_feature_table.to_csv('Output/'+date+'_QCfiltered_median_aggragated_'+batch+'_feature_table.gz', \
                               index=False, compression='gzip')
    
    print('DONE! - tables saved in "Output"')

#### Main Script starts here
if __name__ == '__main__':
    
    path_to_data = '../0_RAW_data/'
    QC_fail = '06152021_Image+WellQC_fail_compressed_plates_run3_merge'
    platemaps = 'compressed_screen_run3'
    platemap_dir = 'platemap_run3/'

    # UPDATE for a given run
    run = 'compressed_screen_run3'

    # UPDATE to the date you run this
    date = '06232021'

    Batch_Aggregate(path_to_data, QC_fail, platemaps, platemap_dir, run, date)
