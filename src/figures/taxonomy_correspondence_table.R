library(tidyverse)
library(speedyseq)

bl.phy <- readRDS('data/cleaned/blphy_filt_7sams.rds')
sea <- readRDS('data/analysis/lombsea.rds')
abunprev <- readRDS('data/analysis/abund_prev.rds')

tax <- as_tibble(as(tax_table(bl.phy), 'matrix'),
                            rownames = 'asv')

tax.silva <- readRDS('data/cleaned/dada2/03_taxonomy/bbmo_silva/bbmo_silva_tax_assignation.rds') %>% 
  as_tibble(rownames = 'seq') %>% 
  mutate(asv = str_c('asv', 1:nrow(.)),
         `order + family (SILVA)` = str_c(order,family, sep = ', ')) %>% 
  rename( `genus (SILVA)` = 'genus') %>% 
  select(asv, `genus (SILVA)`, `order + family (SILVA)`)

totals <- data.frame( asv = taxa_names(bl.phy),
                      total.read = taxa_sums(bl.phy),
                      genus = tax_table(bl.phy)[,'genus']) %>% 
  as_tibble() 

totals.all <- taxa_sums(readRDS('data/cleaned/blphyloseq.rds')) %>% sum()

relab.genus.gtdb <- totals %>% 
  group_by(genus) %>% 
  summarize( relab = ((sum(total.read) * 100) / totals.all) %>%
               round(., 2) %>%
               str_c(. , '%'))
  
genus.interest <- tax %>% 
  mutate(genus = as.character(genus)) %>% 
  group_by(genus) %>% 
  summarize(nasvs = n()) %>% 
  filter(nasvs >= 5, !is.na(genus)) 

n.seasonal.asvs <- sea  %>% 
  left_join(tax, by='asv') %>% 
  filter(genus %in% genus.interest$genus) %>% 
  group_by(genus) %>% 
  summarize(nsea = n())

gtdb.silva.correspondence <- tax %>% 
  mutate(`order + family (GTDB)` = str_c(order,family, sep = ', ')) %>% 
  left_join(tax.silva, by = 'asv') %>% 
  select(genus,
         `order + family (GTDB)`,
         `genus (SILVA)`,
         `order + family (SILVA)`) %>% 
  distinct()

tibble( genus = genus.interest$genus) %>% 
  left_join(gtdb.silva.correspondence, by = 'genus') %>% 
  left_join(n.seasonal.asvs, by = 'genus') %>% 
  left_join(genus.interest, by = 'genus') %>% 
  left_join(relab.genus.gtdb, by = 'genus') %>% 
  rename(`genus (GTDB)` = genus,
         `n. seasonal` = nsea,
          `n. tested ASVs` = nasvs,
          `total relative abundance` = relab) %>% 
  mutate(`General information genus` = NA,
         `References` = NA) %>%   
  write_tsv('results/supp_files/description_main_groups.tsv')
  

