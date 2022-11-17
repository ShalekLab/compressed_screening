The following notebooks contain our PDAC validation experiment analaysis

**1_preprocessing_anndata_with_cNMF_results.ipynb**: (Python)
- Makes a single anndata file containing the scRNA-seq data and the outputs of running cNMF on the data
- Visualizes the validation cNMF modules and identifies those that are highly correlated with the highly variable CS cNMF modules

**2_ligand_effects_with_batch_regressed_out.ipynb**: (R)
- Defines function, run_module_lm_by_ligand_and_day: function to identify ligand effects on cNMF modules with batch effects regressed out
- Runs on all single-ligands for the validation cNMF modules and saves results in a csv file

**3_PDAC_validation_experiment_analysis.ipynb**: (Python)
- Loads in anndata from 1_preprocessing_anndata_with_cNMF_results.ipynb
- Runs PCA, harmony, UMAP for data visualization
- Runs Decoupler to score cells on PROGENy and MsigDB genesets
    - Correlates signature scores with cNMF module usage
- Loads in ligand effects from 2_ligand_effects_with_batch_regressed_out.ipynb and visualizes them