library(tidyverse)
library(speedyseq)
library(DivNet)


bl.phy <- readRDS('data/cleaned/blphyloseq.rds')
bl.phy.genus <- tax_glom(bl.phy, taxrank = 'genus')

divnet_genus.sea <- bl.phy.genus %>% divnet(X = 'season', ncores = 4)
saveRDS(divnet_genus.sea,'data/analysis/divnet_genussea.rds')


divnet_genus.mon <- bl.phy.genus %>% divnet(X = 'month', ncores = 4)
saveRDS(divnet_genus.mon,'data/analysis/divnet_genusmon.rds')

# divnet_asv <- bl.phy %>% divnet(ncores = 3, base = 'asv1')
# saveRDS(divnet_asv,'data/analysis/divnet_asv.rds')
