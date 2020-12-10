library(tidyverse)
library(speedyseq)
library(propr)

args <- commandArgs(trailingOnly = TRUE)

phy <- args[1]
name <- args[2]
  
bl.phy <- readRDS(phy) 

otu <- as(otu_table(bl.phy), 'matrix')

pr <- propr(otu, # rows as samples, like it should be
            metric = "rho", # or "phi", "phs", "cor", "vlr"
            ivar = "clr", # or can use "iqlr" instead
            alpha = NA, # use to handle zeros
            p = 999) # used by updateCutoffs

updat.pr <- updateCutoffs(pr,
              cutoff = seq(0, 1, .1)) # cutoffs at which to estimate FDR




# With 3% we get a FDR of 2% which is good enough? We will also filter by 
#the rho value 
res.propr <- getResults(updat.pr, cutoff = 0.3) %>% 
  as_tibble()

res.propr.07 <- res.propr %>% 
  filter(abs(propr) >= 0.7)

saveRDS(res.propr, str_c('data/analysis/propr_res',name,'.rds'))
saveRDS(res.propr.07, str_c('data/analysis/propr_res07',name,'.rds'))