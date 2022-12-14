# Codon_usage_prot_structure_angles
A recent article in Nature Communications by [Rosenberg *et al.*](https://doi.org/10.1038/s41467-022-30390-9) found correlations between synonymous codon usage and protein dihedral bond angles in *E. coli*. Although they corrected for potential differences in selection on codon usage across protein secondary structures, they did not account for differences in synonymous codon usage across genes. Notably, their test statistic for determining if the dihedral bond angles significantly differed between synonymous codons was correlated with differences in codon-specific elongation rates. Although this could be indicative of a mechanism, it could also indicate they are actually detecting signals related to selection for translation efficiency, which might occur if not properly controlling for gene expression in their analysis.

A reanalysis of the results from Rosenberg et al. Nat. Comm. 2022 using simulated data. Figures in Matters Arising can be re-created using the `analysis.Rmd` file. `analysis.pdf` is a compiled version briefly describing our methods and results.

A preprint of the article currently submitted to Nat. Comm. as a Matters Arising can be found [here](https://doi.org/10.1101/2022.10.26.513858).

# Installing AnaCoDa

We recommend installing `AnaCoDa` from https://github.com/acope3/RibModelFramework/tree/develop. The develop branch is currently more up-to-data. `AnaCoDa` can be installed locally using the commands:

```
R CMD build RibModelFramework .
R CMD INSTALL <tar.gz file produced from build>
```

