# ROC Analysis

This directory stores two model fits of ROC-SEMPPR as implemented in the AnaCoDa R package relevant to this paper and the rebuttal.

- 00_Ecoli_refseq_ena_K_1:
    - The primary model fit used for the analyses presented in the main manuscript currently associated with this repo.
- 01_Ecoli_refseq_no_yp_ena_outliers_different_CSP
    - Performed allowing for variation in the codon-specific parameter (CSP) mutation bias and selection for outlier genes as identified in the file `../Ecoli/clara_clustering_for_outliers.tsv`
    - The parameters from this model fit are used for the simulations to allow variation in codon usage across protein secondary structures using selection parameters taken from a previous publication (see https://doi.org/10.1186/s12864-022-08635-0) which excluded these supposed outlier genes.
