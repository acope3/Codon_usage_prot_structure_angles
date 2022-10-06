#!/bin/bash
INPUT="/data2/cope/Codon_usage_prot_structure_angles/Ecoli/Ecoli_refseq_and_ena.fasta"
OUTPUT="/data2/cope/Codon_usage_prot_structure_angles/Results/Ecoli_refseq_ena_K_1/"
nohup Rscript --vanilla Scripts/rocAnalysis.R -i "$INPUT" -o "$OUTPUT" -d 50 -s 20000 -a 20 -t 5 -n 16 --est_csp --est_phi --est_hyp --max_num_runs 2 --codon_table 1 --development /home/cope/AnaCoDa_installs/develop/ &> /data2/cope/Codon_usage_prot_structure_angles/Results/Ecoli_refseq_ena_K1.Rout &
