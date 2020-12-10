library(tidyverse)
library(speedyseq)
library(breakaway)

bl.phy <- readRDS('data/cleaned/blphyloseq.rds') 


# General -----------------------------------------------------------------
  
rich.res <- breakaway(bl.phy)
rich.df <- data.frame(summary(rich.res),
                      sample_data(bl.phy)) %>% 
  as_tibble()

bt <- betta(summary(rich.res)$estimate,
            summary(rich.res)$error,
            make_design_matrix(bl.phy, "month"))

bt.sea <- betta(summary(rich.res)$estimate,
            summary(rich.res)$error,
            make_design_matrix(bl.phy, "season"))

write_rds(rich.df, 'data/analysis/richness_break_all.rds')
write_rds(bt, 'data/analysis/betta_rich_month.rds')
write_tsv(as_tibble(bt$table, rownames = 'predictor'),
          'results/tests/alphadiv_betta_months.tsv')

write_rds(bt.sea, 'data/analysis/betta_rich_season.rds')
write_tsv(as_tibble(bt.sea$table, rownames = 'predictor'),
          'results/tests/alphadiv_betta_sea.tsv')
  


# Testing with diferent baselines  ----------------------------------------

iterate_predictors <- function(predictors, pred){
  predictors <- predictors %>% as.factor %>% relevel(ref = pred)
  mX <- model.matrix(~predictors) 
  betta(summary(rich.res)$estimate,
        summary(rich.res)$error,
        mX)$table %>% 
    as_tibble(rownames = 'predictor')
}

betta_wposthoc <- function(physeq, var){
  
  predictors <- bl.phy %>% sample_data %>% get_variable(var)
  vec.levels <- unique(predictors)
  
  iterations.list <- map(vec.levels, ~iterate_predictors(predictors, .x)) %>% 
    set_names(nm = vec.levels) 
  
  iterations.list %>% 
    bind_rows( .id = 'level') %>% 
    return() 
  
}

betta.month  <- betta_wposthoc(bl.phy, 'month')
betta.season <- betta_wposthoc(bl.phy, 'season')



