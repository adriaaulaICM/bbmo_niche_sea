library(tidyverse)
library(speedyseq)
library(patchwork)
library(lubridate)
library(lomb)

source("src/util/params-graphs.R") 
source("src/util/main_functions.R")

bl.phy <- readRDS('data/cleaned/blphy_filt_7sams.rds') 
bl.phy <- transform_sample_counts(bl.phy, function(x) x / sum(x))
lombsea <- read_rds('data/analysis/lombsea.rds')
iter <- 300 

# Functions ---------------------------------------------------------------

prepare_dataset <- function(physeq){

 tsdf.0.2 <- physeq %>%
   psmelt_dplyr() %>%
   mutate(decimaldat = decimal_date(Date))

 df.ts <- tsdf.0.2 %>%
   arrange(Date) %>%
   select(Date, decimaldat, Abundance, OTU) %>%
   split(.$OTU, drop = T)

 return(df.ts)
}

check_lomb <- function(dfts){

 lomb.02  <-  dfts %>%
   map(~randlsp( x =.x$Abundance,
                 times = .x$decimaldat,
                 type = 'period',
                 plot = F,trace = FALSE))

 return(lomb.02)

}

extract_lomb_res <- function(lomb){

 lomb.df <- tibble( asv = names(lomb),
                    pval = map_dbl(lomb, ~.x$p.value),
                    peak = map_dbl(lomb, ~.x$peak),
                    interval  = map(lomb, ~.x$peak.at), 
                    int.min = map_dbl(interval, ~.[[2]]),
                    int.max = map_dbl(interval, ~.[[1]]))

 return(lomb.df)
}

calculate_histogram_distro <- function(physeq, keepasvs){
  bl.phy.rem <- prune_taxa(taxa = keepasvs, physeq)
  bl.agg <- tax_glom(bl.phy.rem, taxrank = 'class')
  # bl.agg.asinh <- transform_sample_counts(bl.agg, function(x) asinh(x))
  
  prepare_dataset(bl.agg) %>% 
    check_lomb() %>%
    extract_lomb_res() %>% 
    return()
}

check_at_rank_level <- function(rank.sel, rank.level){
  
  print(rank.sel)
  taxa <- as(tax_table(bl.phy), 'matrix') %>% as_tibble(rownames = 'asv') %>% 
    pivot_longer(names_to = 'rank', values_to = 'value',  cols = -asv)
  asv.sel <- taxa %>%
    filter( value == rank.sel, rank == rank.level) %>%
    pull(asv) %>% unique()
  bl.phy.rank <- prune_taxa(taxa = asv.sel, x = bl.phy)
  
  iterations <- 1:iter %>% 
    map(~calculate_histogram_distro(bl.phy.rank,
                                    sample(taxa_names(bl.phy.rank),
                                           size = ntaxa(bl.phy.rank) * 0.8 ) 
    ))
  
  return(bind_rows(iterations))
}


check_normal <- function(rank.sel, rank.level){
  
  taxa <- as(tax_table(bl.phy), 'matrix') %>% as_tibble(rownames = 'asv') %>% 
    pivot_longer(names_to = 'rank', values_to = 'value',  cols = -asv)
  asv.sel <- taxa %>%
    filter( value == rank.sel, rank == rank.level) %>%
    pull(asv) %>% unique()
  bl.phy.rank <- prune_taxa(taxa = asv.sel, x = bl.phy)
  
  normal <- calculate_histogram_distro(bl.phy.rank,
                                       taxa_names(bl.phy.rank)) %>% 
    mutate( nasvs = length(asv.sel))
  
  return(normal)
}

check_situation_max <- function(rank.sel, rank.level){
  
  taxa <- as(tax_table(bl.phy), 'matrix') %>% as_tibble(rownames = 'asv') %>% 
    pivot_longer(names_to = 'rank', values_to = 'value',  cols = -asv)
  asv.sel <- taxa %>%
    filter( value == rank.sel, rank == rank.level) %>%
    pull(asv) %>% unique()
  bl.phy.rank <- prune_taxa(taxa = asv.sel, x = bl.phy)
  
  theones <- taxa_sums(bl.phy.rank) %>%
    sort(decreasing = TRUE) %>%
    .[1:ntaxa(bl.phy.rank) * 0.2] %>% names()
  
  thealt <- calculate_histogram_distro(bl.phy.rank,
                                       setdiff(taxa_names(bl.phy.rank),
                                               theones))
  
  
  
  return(thealt)
}


rawdat <- as(tax_table(bl.phy), 'matrix') %>%  
  as_tibble(rownames = 'asv') %>% 
  filter(!is.na(class)) %>% select(-seq) %>%
  filter(class %in% c('Alphaproteobacteria',
                      'Gammaproteobacteria',
                      'Bacteroidia')) %>% 
  pivot_longer(names_to = 'rank', values_to = 'value', cols = -asv) %>%
  filter(rank != 'domain', rank != 'phylum') %>%
  filter(!is.na(value))

seasonallevels <- rawdat %>% 
  filter(asv %in% lombsea$asv)

ranks_more10 <-  rawdat  %>%
  group_by(rank,value) %>%
  filter(n() >= 10) %>%
  filter(value %in% seasonallevels$value)  %>% 
  distinct(rank,value)

distributions <- ranks_more10 %>% 
  mutate( values = map2(.x = value, .y = rank,
                        ~check_at_rank_level(.x, .y)))

without20max <- ranks_more10 %>% 
  mutate( values = map2(.x = value, .y = rank,
                        ~check_situation_max(.x, .y)))

normal <- ranks_more10 %>% 
  mutate( values = map2(.x = value, .y = rank,
                        ~check_normal(.x, .y)))

saveRDS(distributions,'data/analysis/cohesivity_dist.rds')
saveRDS(without20max, 'data/analysis/cohesivity_max.rds')
saveRDS(normal, 'data/analysis/cohesivity_normal.rds')
  