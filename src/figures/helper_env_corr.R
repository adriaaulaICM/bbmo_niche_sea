library(tidyverse)
library(corrr)
library(speedyseq)

bl.phy <- readRDS('data/cleaned/blphy_filt_7sams.rds')

scale2 <- function(x, na.rm = TRUE) (x - mean(x, na.rm = na.rm)) / sd(x, na.rm)
env <- sample_data(bl.phy) %>% as_tibble() %>%
  select(Sample_name, Temperature, Chla_total, PO4, 
         NO2, NO3, PNF_Micro, HNF_Micro) %>% 
  mutate_if(is.numeric, scale2) 

env %>% 
  naniar::vis_miss()
# we have removed salinity due to a lot of NAs, we can always process her
# separately 
env.filt <- env  %>% 
  na.omit()


corrdf <- env.filt %>% 
  select(-Sample_name) %>% 
  correlate()  %>% 
  rearrange()  %>% 
  shave() 

cor.pvals <- outer(env.filt[,-1], env.filt[,-1], function(X, Y){
  mapply(function(...) cor.test(..., na.action = "na.exclude")$p.value,
         X, Y)
})

cor.pvals %>% 
  data.frame() %>% 
  rownames_to_column(var = 'parameter') %>% 
  pivot_longer(names_to = 'par', values_to = 'pval', cols = -parameter) %>% 
  filter(pval <= 0.01)

env.filt
fashion(corrdf)
corrplot <- rplot(corrdf,print_cor = T)
ggsave(corrplot,filename = 'results/figures/env_corr.pdf')

network_plot(env.filt %>% 
  select(-Sample_name) %>% 
  correlate(), min_cor = .5)

# smaller version 

ggplot(env.filt, aes(Temperature, NO2)) + 
         geom_point() + 
  geom_smooth() 

       