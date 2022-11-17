The followingg notebooks contain our analysis of publically available PDAC tumor datasets

**1_preprocessing_TCGA_data.ipynb**: (Python)
- Loads in TCGA PDAC dataset
- Adds module scores to the data for
    - Prognostic PDAC signatures (e.g. Moffit et al Classical/Basal scores)
    - cNMF modules from validation experiments
    - MsigDB genesets related to the ligands in validation experiments

**2_preprocessing_Raghavan_scRNAseq.ipynb**: (Python)
- Loads in Raghavan et al PDAC tumor scRNA-seq dataset and does standard scRNAseq preprocessing
- Adds module scores to the data for
    - Prognostic PDAC signatures (e.g. Moffit et al Classical/Basal scores)
    - cNMF modules from validation experiments
    - MsigDB genesets related to the ligands in validation experiments

**3_prognostic_signature_vs_cNMF_module_correlation_TCGA.ipynb** (Python)
- Examining the correlation between expression of Basal/Classical axis and the expression of day 7 validation cNMF modules in TGCA PAAD samples 
 - Loads in TCGA metadata from 1_preprocessing_TCGA_data.ipynb
 - Plots the correlations between cNMF modules and Basal/classical scores

**4_prognostic_signature_vs_cNMF_module_correlation_raghavan_single_cell.ipynb**: (Python)
- Examining the correlation between expression of Basal/Classical axis and the expression of day 7 validation cNMF modules in Raghavan et al single cells
- Loads in Raghvan et al anndata from 2_preprocessing_Raghvan_scRNAseq.ipynb
- Plots the correlations between cNMF modules and Basal/classical scores

**5_scatterplots_of_t2I_ligand_genesets_vs_classical_score.ipynb** (R)
- Plotting type 2 immunity (t2I) gene set expression vs Moffit et al classical score in TCGA and Raghavan datasets
- Loads in metadata csvs of TCGA and Raghavan datasets with module scores
- Makes scatterplots

**6_macrophage_expression_analysis.ipynb**: (Python)
- Loads in non malignant cells from Raghavan et al and subsets down to macrophages
- Runs differential expression analysis between macrophage subtypes
- Plots IL4I1 expression across macrophage subtypes