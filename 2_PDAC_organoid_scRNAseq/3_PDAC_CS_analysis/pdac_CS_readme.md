The following notebooks contain our analysis for the PDAC organoid compressed screens

**1_adding_cNMF_output_and_design_matrices_to_anndata.ipynb**: (Python)
- Input data:
    - Adata files with Q1Q2 (plate 1+2) post QC & hash demultiplex scRNA-seq dataset
    - cNMF gene spectra scores, module usage, and highly variable genes from running cNMF wdl file
    - Design matrices indicating which drugs went into which wells in the compressed screens
- Loads in all inputs into anndata file
- Identifies cNMF modules that are highly variable across cells

**2_visualization_geneset_scoring_and_hit_deconvolution.ipynb**: (Python)
- Loads in anndata from 1_adding_cNMF_output_and_design_matrices_to_anndata.ipynb
- Runs PCA, harmony, UMAP for data visualization
- Runs Decoupler to score cells on PROGENy and MsigDB genesets
    - Correlates signature scores with cNMF module usage
- Runs linear module deconvolution
    - Using functions from **CompressedScreen.py** and **linear_deconvolution.py**

**3_key_geneset_intersections.ipynb**: (Python)
- For Type 2 cytokine signaling and IFNgamma signaling genesets
- Makes a venn diagram between published genesets and the top genes in the corresponding cNMF module in our data
