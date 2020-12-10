library(tidyverse)
library(speedyseq)
library(corncob)

source('src/analysis/helper_prepare_nichepart_data.R')


corncob_test <- function(par = 'Temperature'){

  selsamples <- data.frame(sam = sample_names(bl.phy.gen.filt),
                           var = get_variable(bl.phy.gen.filt, par)) %>% 
    filter(!is.na(var)) %>% 
    pull(sam) %>% 
    as.character()
  
  thephyseq <- prune_samples(samples = selsamples, bl.phy.gen.filt)
  
  testtemp <- differentialTest(formula = as.formula(str_c('~', par)),
                               formula_null = ~1,
                               phi.formula = ~1, 
                               phi.formula_null =  ~1,
                               data = thephyseq,
                               test = 'Wald',
                               boot = FALSE,
                               fdr_cutoff = 0.05,
                               na.rm = TRUE)

  return(testtemp)

}

params <- colnames(env)[-1]
estimates.par <-  params %>%
  set_names(x =  map(., ~corncob_test(.x)) ,
            nm = .)


pars.results <- estimates.par %>%
  map_dbl(~.x$significant_taxa %>%
            length()) %>% .[ . > 0] %>% names()


obtain_coeffs <- function(obj){
  data.frame(plot(obj, level = 'genus', data_only = T) ,
             asv = obj$significant_taxa)
}


coeffs.corncob <- pars.results %>%
  set_names(x = map(., ~obtain_coeffs(estimates.par[[.x]])),
            nm = .) %>%
  bind_rows(.id = 'par') %>%
  as_tibble()  %>%
  left_join(tax, by = 'asv')  %>%
  rename( 'estimate' = x,
          'upper' = xmax,
          'lower' = xmin)
  # Only strong coefficients
  # filter(abs(estimate) > 0.01)  %>% 
  # Only estimates with the confidence interval without taking the zero value
  # filter(sign(upper) == sign(lower))


write_rds(coeffs.corncob, 'data/analysis/difftest_par.rds')

difft <- differentialTest(formula = ~ Temperature + Chla_total + 
                            PO4+ NH4+ NO2+ NO3+ Si+ PNF_Micro+ HNF_Micro,
                          formula_null = ~1,
                          phi.formula = ~1,
                          phi.formula_null = ~1,
                          data = bl.phy.gen.filt,
                          test = 'Wald', 
                          boot = FALSE,
                          fdr_cutoff = 0.05, na.rm = T)

coefs.all <- obtain_coeffs(difft) %>% 
  left_join(tax, by = 'asv')  %>% 
  rename( 'estimate' = x, 
          'upper' = xmax,
          'lower' = xmin,
          'par' = variable) %>% 
  mutate(par = str_replace(par,
                           '\nDifferential Abundance', ''))  %>% 
  # Only strong coefficients
  filter(abs(estimate) > 0.01)  %>% 
  # Only estimates with the confidence interval without taking the zero value
  filter(sign(upper) == sign(lower)) %>% 
  as_tibble()

write_rds(coefs.all, 'data/analysis/niche_par_corncob.rds')

compadrecompare <- list(  sep = coeffs.corncob %>%
                            mutate( whole = str_c(asv,par)) %>% 
                            pull(whole),
                          together = coefs.all %>%
                            mutate( whole = str_c(asv,par)) %>% 
                            pull(whole))

library(UpSetR)
  
upset(fromList(compadrecompare), order.by = "freq") 

coeffs.corncob <- readRDS('data/analysis/difftest_par.rds')
coeffs.month <- readRDS('data/analysis/difftest_par_wmonth.rds')

compadrecompare <- list(  nomonth = coeffs.corncob %>%
                            mutate( whole = str_c(asv,par)) %>% 
                            pull(whole),
                          month = coeffs.month %>%
                            mutate( whole = str_c(asv,par)) %>% 
                            pull(whole))

library(UpSetR)
  
upset(fromList(compadrecompare), order.by = "freq") 
