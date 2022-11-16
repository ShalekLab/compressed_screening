Workflow for analyzing the ground truth cell painting data

**1_dose_time_coef_of_variatoin_calculations.ipynb**: (Python)
For choosing dose and time point for cell painting experiments
- Analysis of the 2 replicates 3 doses (0.1,1, 10 uM) X 3 time points (6 hours, 24 hours, 48 hours) dataset
- For each time point
    - calculates the effect of compounds relative to DMSO (measured by Mahalanobis distance
    - calcualtes the coeficient of variation of the effects across compounds
    - Largest coef of variation at 24 hour timepoint
- For each dose in the 24 hour timepoint
    - Calcualte coef of variation of mahalanobis distances to access best dose to screen at

**2_GT_screen_analaysis.ipynb**: (Python)
Analysis of GT screen, 316 perturbations, 6 replicates each; DMSO negative controls; 8 positive control single compounds
- Loads in anndata file of harmonized GT data, adds in appropriate metadata
- Calculates UMAP 
- Calculates Mahalanobis distance of each perturbation replicate from DMSO
    - Calls significant perturbations as those that have Mahalan. distance > mean(Mahala. distances) + std(Mahala. distances)
- Identifying clusters of perturbations
    - Runs louivan clustering over the samples with resolution scaninng to choose best clustering resolution
        - This defines GT cellular phenotypes
        - Runs enrichment tests to identify types of cell painting features enriched in each phenotype
    - For each drug,
        - Runs enrichment tests to see if the replicates from the drug are enriched in any given phenotype
    - Clusters drugs over the enrichment scores identified in the above step
        - This defines the GT drug clusters