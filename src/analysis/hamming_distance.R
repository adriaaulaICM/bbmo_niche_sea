library(tidyverse)
library(dada2)
library(phyloseq)

bl.phy <- readRDS('data/cleaned/blphyloseq.rds') 

tax <- as(tax_table(bl.phy), 'matrix') %>% 
  data.frame() %>% 
  rownames_to_column(var = 'asv') %>% 
  mutate_if(.predicate = is.factor,.funs = as.character)

comparison_hamm <- function(chr){
  chr %>%
    cross2(., ., .filter = `==`) %>% 
    map(setNames, c("seq1", "seq2")) %>% 
    bind_rows() %>% 
    mutate( hammdist = nwhamming(seq1,seq2))
}

genus.hamm <- tax %>% group_by(genus) %>% 
  summarise(n = n()) %>% 
  filter( n > 1) %>% 
  pull(genus)

hammdist.all <- tax %>%
  filter(genus  %in% genus.hamm) %>%
  split(.$genus) %>%
  map('seq') %>%
  map(~comparison_hamm(.)) 


# HORRENDOUS WAY OF doing things
doubled <- bind_rows(hammdist.all, .id = 'genus')  %>% 
  left_join(tax, by = c('seq1' = 'seq')) %>% 
  left_join(tax, by = c('seq2' = 'seq')) %>% 
  select(genus,asv.x,asv.y,hammdist) %>% 
  gather( 'compar', 'asv', asv.x, asv.y) %>% 
  filter(compar == 'asv.y') %>% 
  distinct(genus,asv,hammdist) %>% 
  pull(asv)


#I know, I know..
hammrds  <- bind_rows(hammdist.all, .id = 'genus')  %>% 
  left_join(tax, by = c('seq1' = 'seq')) %>% 
  left_join(tax, by = c('seq2' = 'seq')) %>% 
  select(genus,asv.x,asv.y,hammdist) %>% 
  rename( asv = 'asv.x', asv2 = 'asv.y') 

  # Simmetrical comparisons, we need only one of them
  # gather( 'compar', 'asv', asv.x, asv.y) %>% 
  # distinct(otu,asv,hammdist) %>% 
  # mutate(asv2 = doubled)


write_rds(hammrds, 'data/analysis/hammdist.rds')
