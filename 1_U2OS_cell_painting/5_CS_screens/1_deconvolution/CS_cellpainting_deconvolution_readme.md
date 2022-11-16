The following scripts contain our workflow for formatting cell painting data for our linear regression deconvolution, running the deconvolution, and running a simpler alternative approach to deconvolution as a benchmark for our regression method.

**1_format_cell_painting_data_for_deconvolution.py**: Script for ingesting cell painting. Writes 3 csvs
- Features: the cell painting features
- Metadata: Screen, well, perturbation, etc 
- Design matrix: Binary matrix (0/1) indicating which perturbations went into which wells

**2_run_deconvolution_save_pvals.py**: Loops through all CS cell painting screens and runs deconvolution code. Calls
- **CompressedScreen.py**: script containing python CompressedScreen class that is built by importing the 3 csvs above. Structures data nicely for running deconvolution
- **linear_deconvolution.py**: Script containing the functions for running regression and permutatoin testing

**3_saving_regression_results_at_multiple_permutation_thresholds.ipynb**: Loads in the results from the deconvolution and permutation testing and saves a coeficient matrix at multiple permutation test p value thresholds for latter examination of the impact of the permutation test stringency on the deconvolution results.

**4_CS_deconvolution_benchmark_calculating_Mahalanobis_distances_and_bootstrap_distrubtions.ipynb**:
Conducting two simpler approaches to running deconvolution on the compressed screen data. So that we can benchmark our regression method against these (spolier: our regression method did a lot better than these approaches)
1) The simplest approach we try is that for each perturbation in the compressed dataset, we calculate the average of the Mahalanobis distances of the wells that the perturbation occurred in. 
2) We then try a more sophisticated approach where rather than simply finding the average of the Mahalanobis distances of the wells a perturbation, we do a bootstrap resampling of these that we can later use to try to find a bootstrapped Mahalanobis distance that is less influenced by outlier wells for a given perturbation. (i.e. say perturbation A has little effect, and occurs in 5 wells, in one of which perturbation B (which has a strong effect) is present as well, causing that one well to have a large Mahalanobis distance. In this case, the boostrap resampling may help us to exlcude that effect from perturbation B).
