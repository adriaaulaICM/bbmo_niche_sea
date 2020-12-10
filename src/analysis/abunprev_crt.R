library(phyloseq) 
library(tidyverse) 
library(ggthemes)
library(janitor)

bl.phy <- readRDS('data/cleaned/blphyloseq.rds')
bl.phy.relab <- transform_sample_counts(bl.phy, function(x) x / sum(x))

calculate_abund_prevalen.df <- function(physeq){
  
  abundance.prevalence.df <-   as(otu_table(physeq), 'matrix') %>%
    data.frame() %>%
    gather(key = 'asv', value = 'Abundance') %>%
    group_by(asv) %>%
    summarize( ocurrence = sum(Abundance > 0),
               ocurrence_more1 = sum(Abundance >= 0.01),
               mean.reads = mean(Abundance)) %>%
    mutate(prevalence = ocurrence / nsamples(physeq),
           presence = ifelse( ocurrence_more1 > 0, 'abundant',
                              'rare'),
           behavior = case_when(prevalence >= 0.75 ~ 'Broad',
                                prevalence <= 0.10 ~ 'Narrow',
                                TRUE ~ 'Intermediate'))
  
  return(abundance.prevalence.df)
}

abunprev <- calculate_abund_prevalen.df(bl.phy.relab)




# Incorporate CRT  --------------------------------------------------------

library(vegan)
library(TSA)
source("src/util/crt_functions.R")


incorporate_crt_results <- function(prev.df, physeq, string){
  
  taxa.hyp.crt <- prev.df %>% 
    filter( ocurrence_more1 > 0, behavior %in% c('Narrow', 'Intermediate')) %>% 
    pull(asv)
  
  table.crt <- otu_table(physeq)[,taxa.hyp.crt] %>%  as(. , 'matrix') %>% 
    t() %>% 
    data.frame()
  
  
  crt.df <- SimpleRareToPrev.f(
    otu_fp = table.crt,
    outitle = string,
    abund_thresh = 0.005,
    abund_thresh_ALL = FALSE,
    b_thresh = 0.90,
    rdp_lastcol = FALSE
  )
  
  prev.df <- prev.df %>% 
    mutate( crt = ifelse(asv %in% crt.df$OTUID, TRUE, FALSE)) 
  
  print(prev.df %>% 
          select(presence,behavior,crt) %>% 
          table())
  
  return(prev.df)
}  

abunprev.crt <- incorporate_crt_results(abunprev, bl.phy, string = 'fl') %>% 
  mutate( behavior = ifelse(crt, 'CRT', behavior), 
          presence = factor(presence),
          behavior = factor(behavior,
                            levels = c('CRT', 'Narrow', 'Intermediate', 'Broad')))

saveRDS(abunprev.crt, 'data/analysis/abund_prev.rds')
