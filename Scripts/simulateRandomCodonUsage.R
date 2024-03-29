library(AnaCoDa,lib.loc="~/AnaCoDa_installs/ignore_bad_codons_simulation/")
library(tidyverse)

real.data <- read_csv("Rosenberg_analysis/dataverse_files/data-precs.csv",col_types = cols()) %>%
  filter(codon_score == 1 & codon != "---" & secondary != "-") 

aa.used <- real.data %>%
  filter(str_detect(codon,"[AGCT]{3}")) %>%
  distinct(unp_id,unp_idx,.keep_all = T) %>%
  dplyr::select(unp_id,unp_idx,name) %>%
  rowwise() %>%
  mutate(new.codon = case_when(
    name == "S" ~ dqsample::sample(c(AAToCodon(name),AAToCodon("Z")),1),
    T ~ dqsample::sample(AAToCodon(name),1)
  )
  )

sim.data <- real.data %>% 
  left_join(aa.used,by=c("unp_id","unp_idx","name")) %>%
  mutate(codon = ifelse(str_detect(codon,"[AGCT]{3}"),new.codon,codon),
         codon_opts = ifelse(str_detect(codon,"[AGCT]{3}"),new.codon,codon)) %>%
  dplyr::select(-new.codon)


## Uncomment line below to save data to file
write_csv(sim.data,"Rosenberg_analysis/dataverse_files_sim_w_ena_random/data-precs.csv")


```