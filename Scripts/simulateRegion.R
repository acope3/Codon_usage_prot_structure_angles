## Author: Alex Cope
## Code for simulating genomes.

library(AnaCoDa,lib.loc = "~/AnaCoDa_installs/ignore_bad_codons_simulation/")

simulate <- function(genome.file,sel.file,mut.file,phi.file,output.dir,output.fasta)
{
  dir.create(output.dir) ## Creates directory 
  genome <- initializeGenomeObject(genome.file) ## Load genome from genome.file
  phi <- read.table(phi.file,sep=",",header=T) ## read in gene expression data, should be comma separated
  phi <- phi$Mean ## get posterior mean estimates 
  
  ## Set up parameter object
  parameter <- initializeParameterObject(genome,
                                         sphi=c(0.01),
                                         num.mixtures = 1,
                                         gene.assignment = rep(1,length(genome)),
                                         mixture.definition = "allUnique",
                                         initial.expression.values = phi, ## specify phi (gene expression) values from previous ROC run
                                         split.serine = TRUE)
  parameter$initSelectionCategories(sel.file,1,F) ## Initialize selection values from previous runs, from output of AnaCoDa getCSPEstimates, does not expect reference codon
  parameter$initMutationCategories(mut.file,1,F) ## Initialize mutation values from previous runs, from output of AnaCoDa getCSPEstimates, does not expect reference codon

  
  model <- initializeModelObject(parameter,model="ROC",with.phi=FALSE)
  
  ## Parameters are initialized at values from previous ROC fit, then simulate genome and output to file
  model$simulateGenome(genome)
  genome$writeFasta(paste(output.dir,output.fasta,sep="/"),simulated=T)
  
}
sel.file <- "/data2/cope/Codon_usage_prot_structure_angles/Results/Ecoli_refseq_ena_K_1/run_1/Parameter_est/Ecoli_refseq_and_ena.fasta_no_ref_Selection.csv"
mut.file <- "/data2/cope/Codon_usage_prot_structure_angles/Results/Ecoli_refseq_ena_K_1/run_1/Parameter_est/Ecoli_refseq_and_ena.fasta_no_ref_Mutation.csv"
genome.file <- "/data2/cope/Codon_usage_prot_structure_angles/Ecoli/Ecoli_refseq_and_ena.fasta"
phi.file <- "/data2/cope/Codon_usage_prot_structure_angles/Results/Ecoli_refseq_ena_K_1/run_1/Parameter_est/gene_expression.txt"
output.dir <- "/data2/cope/Codon_usage_prot_structure_angles/Ecoli/"
output.fasta <- "simulatedEcoli_w_ena.fasta"
simulate(genome.file,sel.file,mut.file,phi.file,output.dir,output.fasta)

