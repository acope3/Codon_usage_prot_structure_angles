Model was run with bash command
```
Rscript --vanilla Scripts/rocAnalysis.R -i /data2/cope/Codon_usage_prot_structure_angles/Ecoli/Ecoli_refseq_no_yp_pseudo_and_ena.fasta -o /data2/cope/Codon_usage_prot_structure_angles/ROC_analysis/2024-01-04_Ecoli_refseq_no_yp_ena_K3_1_2_vs_3/ -d 50 -s 20000 -a 20 -t 5 -n 8 --est_csp --est_phi --est_hyp --mixture_assignment Ecoli/clustering/K_2/Ecoli_refseq_no_yp_pseudo_and_ena_k3_1_2_vs_3.tsv --mixture_definition allUnique --max_num_runs 2 --codon_table 1 --development /home/cope/AnaCoDa_installs/develop/
```
