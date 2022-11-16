These scripts describe our analysis workflow and generate the figures presented in the manuscript

**1_summarizing_deconvolution_effects.ipynb**: (Python)
- Loads in the GT effects data
    - Including which compounds were "significant" GT hits (Mahalanobis distance > mean(Mahala. dist) + std(Mahala. dist))
- For all CS screens
    - Reads in the results of decon at many permutation testing stringencies
    - Saves the scaled L1 norm for each perturbation as a measure of CS effect
- For each CS screen at each stringency level of permutation testing
    - Calulcate Pearson & Spearman correlation between GT effects and CS inferred effects
    - Calculate TPR and FPR for finding significant hits from GT
        - True positive = Hit in CS is a significant hit in GT
        - False positive = Hit in CS is not aa significant hit in GT
        - True negative = Not a hit in CS and is not a significant hit in GT
        - False negative = Not a hit in CS and is a significant hit in GT

**2_CS_effect_vs_GT_effect_scatterplots.ipynb**: (R)
Reads in data from 1_summarizing_deconvolution_effects.ipynb and for each screen and three permutation test pvalues, makes a scatterplot of CS effect (scaled L1 norm) vs GT effect (Mahalanobis distane) with drugs that are a hit in both colored by the GT drug cluster

**3_CS_hits_dotplot.ipynb**: (R)
Makes the dot plot summarizing which drugs were consistent hits across CS screens

**4_CS_vs_GT_pearson_by_permutation_test_pvalue_and_TPR_vs_FPR_plots.ipynb**: (R)
- Loads in results from applying different permutation testing pvalues (0.01 to 1 by steps of 0.01) to deconvolution of each screen
- Makes Pearson and Spearman correlation plots for pvalue thresholds at 0.001, 0.01, and 0.05
- Makes TPR vs FPR plots

**5_collating_CS_decon_results_and_calculating_correlation_with_GT_effects.ipynb**: (Python)
- Loads in data across all CS screens on the pearson correlation between GT & CS for
    - Linear model deconvolution
    - Mahalanobis distance deconvolution
    - Boostrapped mahalanobis distance deconvolution
- Formats data for plotting
- Generates a plot comparing the methods

**6_plotting_deconvolution_approach_correlations_with_GT.ipynb**: (R)
- Makes a plot showing the correlation between GT & CS for the three different deconvolution scripts