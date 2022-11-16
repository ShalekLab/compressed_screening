To preprocess the scRNA-seq data, we ran the following scripts on each scRNA-seq dataset that we gathered

**1_qc_filtering.R**:
Runs basic qc filtering on the scrna-seq data, removing cells with two few UMIs, to high of a percent of reads that are 
mitochondrial etc

**2_cellhashing_demultiplexing.R**:
Demultiplexes samples based on antibody hashing calls **Preprocess_QC_BEM_v1.5_HASHING.R**

**3_quality_metric_visualization.R**:
Visualizes scRNA-seq quality metrics for samples after demultiplexing

