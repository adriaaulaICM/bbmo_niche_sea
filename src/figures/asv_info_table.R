library(tidyverse)
library(speedyseq)

# PHyloseqs
bl.phy <- readRDS("data/cleaned/blphyloseq.rds") 
tax <- as(tax_table(bl.phy), 'matrix') %>% as_tibble(rownames = 'asv') %>% 
  select(asv, domain, phylum, class, order, family, genus)

abunprev <- readRDS('data/analysis/abund_prev.rds') %>% 
  select(asv, presence, behavior, crt) 
  

lomb02 <- readRDS('data/analysis/lombsea.rds')

table <- tax %>% 
  left_join(abunprev, by = 'asv') %>% 
  mutate( seasonal = ifelse(asv  %in% lomb02$asv, TRUE, FALSE)) %>% 
  rename('distribution' = behavior)



table %>% 
  write_csv('results/tables/table2_asv_info.csv') 
