# Codon_usage_prot_structure_angles
A recent article in Nature Communications by [Rosenberg *et al.*](https://doi.org/10.1038/s41467-022-30390-9) found correlations between synonymous codon usage and protein dihedral bond angles in *E. coli*. Although they corrected for potential differences in selection on codon usage across protein secondary structures, they did not account for differences in synonymous codon usage across genes. Notably, their test statistic for determining if the dihedral bond angles significantly differed between synonymous codons was correlated with differences in codon-specific elongation rates. Although this could be indicative of a mechanism, it could also indicate they are actually detecting signals related to selection for translation efficiency, which might occur if not properly controlling for gene expression in their analysis.

`Akeju_and_Cope_Response_to_Rosenberg_etal_Nat_Comm_2022.pdf` describes our methods and findings. 

A reanalysis of the results from Rosenberg et al. Nat. Comm. 2022 using simulated data. Figures in Matters Arising can be re-created using the `analysis.Rmd` file. 

# Installing AnaCoDa

We recommend installing `AnaCoDa` from https://github.com/acope3/RibModelFramework/tree/develop. The develop branch is currently more up-to-data. `AnaCoDa` can be installed locally using the commands:

```
R CMD build RibModelFramework .
R CMD INSTALL <tar.gz file produced from build>
```

