library(tidyverse)
library(here)
library(speedyseq)

bl.phy <- readRDS(here('data/cleaned/blphy_filt_7sams.rds') )
hamm <- readRDS(here('data/analysis/hammdist.rds'))
abprev <- readRDS(here('data/analysis/abund_prev.rds'))
tax <- as(tax_table(bl.phy), 'matrix') %>% as_tibble(rownames = 'asv')

asvs.sel <- abprev %>%  filter(prevalence >= 0.3) %>% pull(asv)

asv.sel.hamm5 <- hamm %>% 
  filter(hammdist <= 5) %>% 
  filter(asv %in% asvs.sel, asv2 %in% asvs.sel)

# selecting based in being in the most typical genus 
genus_more4_nichepart <- tax %>% 
  filter(asv %in% unique(c(asv.sel.hamm5$asv, asv.sel.hamm5$asv2))) %>% 
  group_by(genus) %>% 
  filter(n() >= 4) %>% 
  filter(!is.na(genus)) 

env <- sample_data(bl.phy) %>% as_tibble() %>%
  select(Sample_name, Temperature, Chla_total,
         Secchi,
         PO4, NH4, NO2, NO3, Si,
         PNF_Micro, HNF_Micro)


bl.phy.gen <- subset_taxa(bl.phy,
                          taxa_names(bl.phy) %in% genus_more4_nichepart$asv )

samples.sel <- env$Sample_name
bl.phy.gen.filt <- subset_samples(bl.phy.gen,
                                  sample_names(bl.phy.gen) %in% samples.sel)


### HERE I wanted to check how the coefficients change with the scaling 
samdat.scaled <- env %>% 
  mutate_at(~(scale(.) %>% as.vector), .vars = vars(-Sample_name))  %>% 
  column_to_rownames(var = 'Sample_name') %>% 
  sample_data()
 
sample_data(bl.phy.gen.filt) <- samdat.scaled


### HERE we wanted to see how the outliers affect the distribution 
# samdata_corrected <- sample_data(bl.phy.gen.filt) %>%
#   as_tibble(rownames = 'sams')
# 
# samdata_corrected[26,'PO4'] <- NaN
# samdata_corrected[2,'NO3'] <- NaN
# samdata_corrected[26,'NO3'] <- NaN
# samdata_corrected[10,'NO3'] <- NaN
# samdata_corrected[12,'PNF_Micro'] <- NaN
# 
# new <- column_to_rownames(samdata_corrected, var = 'sams') %>% 
#   sample_data()
#                      
# sample_data(bl.phy.gen.filt) <- new
# 
