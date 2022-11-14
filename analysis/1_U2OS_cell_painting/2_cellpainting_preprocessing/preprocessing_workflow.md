**0_median_aggreate_cell_painting_data.py**: By well median aggregation
- Median aggregates cell painting data across cells in a given well

**1_scale_features.ipynb**: Robust scaling of cell painting features
- On a per batch basis
    - remove EMPTY wells and DMSO, TREAT wells < 50 cells
    - drop uninformative features (same value, 0s)
    - integrate all batches into single consistent full dataset
- On a per plate basis
    - scale each feature with sklearn RobustScaler
    - OR- Perform transformation of features with sklearn PowerTransformer (Yeo-Johnson)

**2_CV_Effect_filtering.ipynb**: Removes redundant/correlated cell painting features by selecting for informative features 
- Aim is to select features with minimal noise (control variance) and maximal information (treat variance) regardless of imaging batch.
- Divide the dataset into the 2 analyses:
    - Dose-time GT data of 316 compounds (DT)
    - 1uM & 24hr GT & CS of 316 compounds (OG316)
- On a per-dataset basis:
    - Calculate per plate ratio of non-DMSO to DMSO well dispersion (MAD)
    - Record features with dispersion ratio > 66th pctl
    - Retain features which pass above condition in > 50% of plates
    - Run mRMR feature selection (retain top 150)
