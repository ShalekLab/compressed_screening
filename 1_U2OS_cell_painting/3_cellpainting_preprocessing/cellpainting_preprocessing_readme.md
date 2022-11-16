Preprocessing scripts for cell painting features output by distributed cell profiler analysis

**0_median_aggreate_cell_painting_data.py**: By well median aggregation
- Median aggregates cell painting data across cells in a given well

**1_scale_features.ipynb**: (Python)
Robust scaling of cell painting features
- On a per batch basis
    - remove EMPTY wells and DMSO, TREAT wells < 50 cells
    - drop uninformative features (same value, 0s)
    - integrate all batches into single consistent full dataset
- On a per plate basis
    - scale each feature with sklearn RobustScaler
    - OR- Perform transformation of features with sklearn PowerTransformer (Yeo-Johnson)

**2_CV_Effect_filtering.ipynb**: (Python)
Together with next script, removes redundant/correlated cell painting features by selecting for informative features 
- Aim is to select features with minimal noise (control variance) and maximal information (treat variance) regardless of imaging batch.
- Divide the dataset into the 2 analyses:
    - Dose-time GT data of 316 compounds (DT)
    - 1uM & 24hr GT & CS of 316 compounds (OG316)
- On a per-dataset basis:
    - Calculate per plate ratio of non-DMSO to DMSO well dispersion (MAD)
    - Record features with dispersion ratio > 66th pctl
    - Retain features which pass above condition in > 50% of plates
    

**3_mRMR_batch.py**: (Python)
Apply maximum relevancy - minimum redundancy feature selection to cell painting features
- Run mRMR feature selection (retain top 150)
- https://github.com/fbrundu/pymrmr

**4_PCA_harmony.ipynb**: (Python)
Performs PCA reductions (keeping PCs accounting for 95% variance exp.), and also saves a harmony batch-corrected PCA space
- Harmony: https://github.com/immunogenomics/harmony

