library(tidyverse)
library(speedyseq)
library(gt)
library(gtsummary)

# PHyloseqs
bl.phy <- readRDS("data/cleaned/blphyloseq.rds") 

abunprev <- readRDS('data/analysis/abund_prev.rds')

lomb02 <- readRDS('data/analysis/lombsea.rds')

table_preval <- function(df, string,samnum){
  
  df %>% 
    select(OTU, Abundance, seasonal, prevalence, presence,crt, behavior) %>% 
    group_by(presence, behavior) %>% 
    summarize( `Count ASVs` = length(unique(OTU)),
               `Count CRT` = sum(crt) / samnum, 
               `Seasonal ASVs` = sum(seasonal) / samnum,
               total = sum(Abundance),
               `Median ocurrence` = median(prevalence)) %>% 
    ungroup() %>% 
    mutate( `Relative abundance` = total / sum(total),
            presence = str_to_title(presence)) %>% 
    rename( 'Distribution' = behavior) %>% 
    select(-total) %>% 
    group_by(presence) %>% 
    gt() %>% 
    fmt_percent( columns = vars(`Relative abundance`,
                                `Median ocurrence`), decimals = 1) %>%   
    tab_footnote(
      footnote = "Broad = in >75% of samples, Narrow = in <10% samples, Intermediate = middle",
      locations = cells_column_labels(columns = vars(Distribution))
    )  %>% 
    tab_footnote(
      footnote = "Seasonality based in lombscargle test. qval<0.01, PN>=10",
      locations = cells_column_labels(columns = vars(`Seasonal ASVs`))
    )  %>% 
    tab_header( title = md(paste0("**Ocurrence and relative abundance in ",
                                  string, "**")) )
  
}


dat.bl <- psmelt(bl.phy) %>% 
  left_join(abunprev, 
            by = c('OTU' = 'asv')) %>% 
  mutate( seasonal = ifelse(OTU  %in% lomb02$asv, TRUE, FALSE),
          behavior = factor(behavior,
                            levels = c('Broad', 'Intermediate',
                                       'Narrow', 'CRT')))

gtable.obj <- dat.bl %>% 
  table_preval(string = ' the BBMO dataset', samnum = 131)


gtsave(data = gtable.obj, filename = 'results/tables/table1.pdf')
