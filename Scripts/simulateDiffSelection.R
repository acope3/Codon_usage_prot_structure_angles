library(AnaCoDa)
library(tidyverse)
library(parallel)

## This is based on groupings used to train PsiPred, which was the tool used to predict 
## protein secondary structures for our previous work looking at the relationship between codon
## usage and protein secondary structures.
real.df <- read_csv("../Rosenberg_analysis/dataverse_files/data-precs.csv") %>%
  mutate(Structure = case_when(
    secondary == "H" | secondary == "G" | secondary == "I" ~ "H",
    secondary == "E" | secondary == "B" ~ "E",
    secondary == "T" | secondary == "S" | secondary == "-" ~ "C"
  ))
ena.unp.map <- read_csv("../Rosenberg_analysis/dataverse_files/meta-structs_filtered.csv") %>%
  dplyr::select(unp_id,ena_id) %>%
  distinct()
id.map <- read_tsv("../ecoli_id_map.tsv") %>%
  dplyr::select(Locus,uniprot_id) 
id.map <- id.map %>% 
  left_join(ena.unp.map,by=c("uniprot_id"="unp_id")) %>%
  mutate(Locus = ifelse(!is.na(ena_id),ena_id,Locus))
cluster <- read_tsv("../Ecoli/clara_clustering_for_outliers.tsv") %>%
  left_join(id.map,by=c("Gene"="Locus"))


unique.seq <- real.df %>%
  distinct(unp_id,unp_idx,.keep_all = T) %>% 
  mutate(name = ifelse(codon %in% c("AGC","AGT"),"Z",name)) %>%
  left_join(cluster,c("unp_id"="uniprot_id")) %>%
  filter(Cluster != 3)

real.filt.df <- real.df %>%
  left_join(cluster,c("unp_id"="uniprot_id")) %>%
  filter(Cluster != 3) %>%
  dplyr::select(-Cluster,-Structure,-Gene,-ena_id)

write_csv(real.filt.df,"../Rosenberg_analysis/dataverse_files_remove_possible_exogenous/data-precs.csv")


sel.helix <- read_csv("../Ecoli/CUB_across_secondary_structures/Helix_Selection.csv") %>%
  mutate(Structure = "H") %>%
  dplyr::rename(Selection = Mean) %>%
  dplyr::select(AA,Codon,Selection,Structure)
sel.sheet <- read_csv("../Ecoli/CUB_across_secondary_structures/Sheet_Selection.csv")  %>%
  mutate(Structure = "E")  %>%
  dplyr::rename(Selection = Mean) %>%
  dplyr::select(AA,Codon,Selection,Structure)
sel.coil <- read_csv("../Ecoli/CUB_across_secondary_structures/Coil_Selection.csv")  %>%
  mutate(Structure = "C")  %>%
  dplyr::rename(Selection = Mean) %>%
  dplyr::select(AA,Codon,Selection,Structure)

sel.coef <- list(sel.helix,
                 sel.sheet,
                 sel.coil) %>%
  bind_rows()

## Assume mutation bias does not vary across protein secondary structures
mut.bias <- read_csv("../ROC_analysis/01_Ecoli_refseq_no_yp_ena_outliers_different_CSP/run_2/Parameter_est/1_Mutation.csv") %>%
  dplyr::select(AA,Codon,Mean) %>%
  dplyr::rename(Mutation = Mean)
phi <- read_csv("../ROC_analysis/01_Ecoli_refseq_no_yp_ena_outliers_different_CSP/run_2/Parameter_est/gene_expression.txt")  %>%
  dplyr::select(GeneID,Mean) %>%
  dplyr::rename(Phi = Mean)

codon.prob.by.gene <- sel.coef %>% 
  left_join(mut.bias) %>%
  left_join(phi,by=character()) %>% 
  mutate(Numerator = exp(-Mutation - Selection * Phi))

codon.total.prob.by.gene <- sel.coef %>% 
  left_join(mut.bias) %>%
  left_join(phi,by=character()) %>% 
  group_by(AA,Structure,GeneID) %>%
  summarize(Denominator = sum(exp(-Mutation - Selection * Phi)))

codon.prob.by.gene <- codon.prob.by.gene %>%
  left_join(codon.total.prob.by.gene,by=c("AA","Structure","GeneID")) %>%
  mutate(Codon.Probability = Numerator/Denominator)

unique.seq.split <- unique.seq %>%
  filter(codon != "---") %>% 
  group_by(unp_id,Locus,name,Structure) %>%
  group_split()

simCodons <- function(n,codon.prob.df)
{
  codon.order <- codon.prob.df$Codon 
  prob.vec <- codon.prob.df$Codon.Probability
  names(prob.vec) <- codon.order
  codon.counts <- rmultinom(1,size=n,prob = prob.vec)
  simulated.codons <- sample(rep(rownames(codon.counts),codon.counts))
  return(simulated.codons)
}


simulated.codons <- mclapply(unique.seq.split,function(df) {
                 if (!df$name[1] %in% c("M","W","X"))
                 {
                   codon.prob.vec <- codon.prob.by.gene %>%
                     filter(AA %in% df$name & Structure %in% df$Structure & GeneID %in% df$Locus)
                   df$Sim.Codons <- simCodons(nrow(df),codon.prob.vec)
                 } else{
                   df$Sim.Codons <- df$codon
                 }
                 df
            },
            mc.cores = 12
) %>% bind_rows()

real.count <- simulated.codons %>% 
  group_by(unp_id,codon) %>%
  summarize(Real=n()) 
sim.count <- simulated.codons %>%
  group_by(unp_id,Sim.Codons) %>%
  summarize(Sim=n())
both.counts <- real.count %>%
  left_join(sim.count,by=c("unp_id","codon"="Sim.Codons"))
cor.test(both.counts$Real,both.counts$Sim)

to.replace <- simulated.codons %>%
  dplyr::select(unp_id,unp_idx,Sim.Codons)

new.df <- real.df %>% 
  left_join(to.replace,by=c("unp_id","unp_idx")) %>%
  filter(!is.na(Sim.Codons)) %>%
  mutate(codon = Sim.Codons,codon_opts=Sim.Codons) %>%
  dplyr::select(-Sim.Codons,-Structure)

write_csv(new.df,"../Rosenberg_analysis/dataverse_files_sim_w_ena_remove_possible_exogenous_different_cub_across_structures/data-precs.csv")
