library(tidyverse)
library(glue)
library(gt)
library(mgcv)
library(broom)
library(ggbeeswarm)
library(speedyseq)
library(ggtext)
library(patchwork)

source("src/util/backbone_functions.R") 
source("src/util/backbone_params-graphs.R") 
source("src/util/params-graphs.R") 

bl.phy <- readRDS('data/cleaned/blphy_filt_7sams.rds')
rho <- readRDS('data/analysis/propr_resblphy.rds') 
hamm <- readRDS('data/analysis/hammdist.rds')
sea <- readRDS('data/analysis/lombsea.rds')

dist <- readRDS('data/analysis/bbmodech_dist_min95.rds') 
  # regrettably the names are switched 
  

tax <- as(tax_table(bl.phy), 'matrix') %>% as_tibble(rownames = 'asv')


rho_hamm_genus3asvs <- rho %>% 
  rename('asv' = Partner, 'asv2' = Pair) %>% 
  left_join( hamm, by = c('asv', 'asv2'))  %>% 
  left_join( dist, by = c('asv', 'asv2'))  %>% 
  filter(hammdist <= 5) %>% 
  group_by(genus) %>% 
  filter(n() > 10)

rho_hamm_genus3asvs %>% 
  ggplot(aes(propr, genus)) +
  geom_quasirandom(alpha = 0.8, groupOnX = FALSE)  

# distance version of the plot
ggplot(rho_hamm_genus3asvs, aes(as.character(round(distance,3)), propr)) + 
  geom_jitter(width = 0.2, alpha = 0.8)  + 
  geom_point(stat = "summary", color = 'red' ) + 
  geom_smooth(method = 'lm', se = T) + 
  facet_wrap(~genus, scales = 'free_x') + 
  ylab('Rho proportionality') + 
  xlab('Nucleotide divergence') + 
  theme_classic(base_size = 16)
  

# Considering all the taxa with only a hammdist >= 0
rho %>% 
  rename('asv' = Partner, 'asv2' = Pair) %>% 
  left_join( hamm, by = c('asv', 'asv2'))  %>% 
  filter(hammdist <= 10) %>% 
  ggplot(aes(hammdist, propr)) + 
  geom_jitter(width = 0.2, alpha = 0.8)  + 
  geom_point(stat = "summary", fun.y = "mean", color = 'red' ) + 
  geom_smooth(method = 'lm', se = T) + 
  scale_x_continuous(breaks = c(1,2,3,4,5,6,7,8,9,10)) + 
  ylab('Rho proportionality') + 
  xlab('Nucleotide divergence') + 
  ggtitle('Correlation tendency along nucleotide differences') + 
  theme_classic(base_size = 16)

# taking out pelagibacter 
rho %>% 
  rename('asv' = Partner, 'asv2' = Pair) %>% 
  left_join( hamm, by = c('asv', 'asv2'))  %>% 
  filter(genus != 'Pelagibacter') %>%
  filter(hammdist <= 10) %>% 
  ggplot(aes(hammdist, propr)) + 
  geom_quasirandom(alpha = 0.8)  + 
  geom_point(stat = "summary", fun.y = "mean",
             color = 'red', shape = 21, size = 3) + 
  geom_smooth(method = 'lm', se = T) + 
  scale_x_continuous(breaks = c(1,2,3,4,5,6,7,8,9,10)) +
  ylab('Rho proportionality') + 
  xlab('Nucleotide divergence') + 
  ggtitle('Correlation tendency along nucleotide differences',
          subtitle = 'Pelagibacter genus is extracted from this comparison') + 
  theme_classic(base_size = 16)

ggsave('results/figures/rho_vs_hammdist_wo_pelag.pdf',
       width = 8, height = 6)

# Applying GAMs to check distribution -------------------------------------


rho_hamm <- rho %>% 
  rename('asv' = Partner, 'asv2' = Pair) %>% 
  left_join( hamm, by = c('asv', 'asv2'))  %>% 
  left_join( dist, by = c('asv', 'asv2'))  %>% 
  filter(hammdist <= 5) %>% 
  mutate(hammdist = as.character(hammdist))


gammodels <- rho_hamm %>% 
  filter(genus %in% c('Pelagibacter', 'Pelagibacter_A', 'Synechococcus_C',
                      'AG-337-I02', 'SAR86A', 'Luminiphilus', 'Litoricola')) %>%
  group_by(genus) %>% 
  nest() %>% 
  mutate(model = map(data, ~gam(propr ~distance,
                                data = .x)),
         results = map(model, glance),
         pval.term = map_dbl(model, ~ summary(.)$p.pv['distance']) ,
         R.square = map_dbl(model, ~ summary(.)$r.sq)) 


gammos <- gammodels %>% 
  select(genus, results, pval.term, R.square) %>% 
  mutate( pval.term = ifelse(pval.term < 0.0001,
                             '<0.0001',
                             as.character(round(pval.term, 3)))) %>% 
  mutate(genus = ifelse(genus %in% genus.italy,
                             str_c('*', genus, '*'),
                             as.character(genus))) %>% 
  unnest(c(results)) %>% 
  ungroup() %>% 
  gt() %>% 
  fmt_markdown(columns = vars(genus)) %>% 
  fmt_number(columns = 3:6, decimals = 1) %>% 
  fmt_number(columns = 3:6, decimals = 1) %>% 
  fmt_number(columns = vars(R.square), decimals = 3)
  
gtsave(gammos, "gamm_rho_hamm.pdf")

# I have to move the table manually since this shitty program cannot move files

file.remove("results/tables/gam_rho_hamm.pdf")
file.copy("gamm_rho_hamm.pdf", "results/tables/gam_rho_hamm.pdf")
file.remove("gamm_rho_hamm.pdf")

# final plot

significantmodels <- gammodels %>% 
  filter(pval.term <= 0.055)

rho_hamm_genus3asvs <- rho_hamm_genus3asvs %>% 
  mutate( significant = ifelse(genus %in% significantmodels$genus,
                               TRUE, FALSE),
          genus.ital = ifelse(genus %in% genus.italy,
                              str_c('*', genus, '*'),
                              as.character(genus))) %>% 
  mutate( genus = factor(genus,
                         levels = unique(genus),
                         labels = unique(genus.ital)))


textos <- significantmodels %>% 
  mutate(pval.term = round(pval.term, digits = 2),
         R.square = round(R.square, digits = 3),
         R.square = str_sub(R.square, 2)) %>% 
  ungroup() %>% 
  mutate( pval.term = ifelse(pval.term < 0.01, 0.01, pval.term),
          text = glue("*p*<{pval.term}<br>R^(2): {R.square}"),
          hammdist = 5, 
          propr = 0.9) %>% 
  mutate(text = ifelse(genus == 'SAR86A', str_replace(text, '<', '='),
                       text),
         genus.ital = ifelse(genus %in% genus.italy,
                             str_c('*', genus, '*'),
                             as.character(genus))) %>% 
  mutate( genus = factor(genus,
                         levels = unique(genus),
                         labels = unique(genus.ital)))

syne.dummy <- data.frame( genus = 'Synechococcus_C',
                          x = 5, 
                          y = 0.847)

plot <- ggplot(rho_hamm_genus3asvs, aes(hammdist, propr)) +
  # geom_hline(data = syne.dummy, aes(yintercept = y), color = 'red', linetype = 2) + 
  # geom_vline(data = syne.dummy, aes(xintercept = x), color = 'red', linetype = 2) + 
  geom_quasirandom(alpha = 0.8, groupOnX = TRUE)  +
  geom_richtext(data = textos, aes(label = text), fill = NA, label.color = NA) +
  # geom_point(stat = "summary", fun.y = "mean_se", color = 'red' ) +
  stat_smooth(aes(group = genus,
                  color = significant),
              method = "gam", formula = y ~ x, show.legend = FALSE) + 
  facet_wrap(~genus) +
  scale_x_continuous(breaks = c(1,2,3,4,5)) +
  scale_color_manual(values = c('grey', 'dodgerblue')) +
  ylab('Rho proportionality') +
  xlab('Nucleotide divergence')  + 
  lil.strip + 
  theme(strip.text.x = element_markdown()) 


dummy.gg <- data.frame( pattern = c('Habitat filtering',
                        'Competitive exclusion',
                        'Random pattern'),
            x.ini = c( 0, 0, 0),
            x.end = c(10, 10, 10),
            y.ini = c( 10, 1, 5.5),
            y.end = c(1, 10, 5.5), 
            x.label = c(3,4,8.8),
            y.label = c(10, 0.5, 5)) %>% 
  ggplot(aes(x = x.ini, y = y.ini, color = pattern)) +
  geom_segment(aes(xend = x.end,
                   yend = y.end),
               show.legend = F, 
               size = 1.3) + 
  geom_text(aes(label = pattern, x = x.label, y = y.label), 
            show.legend = F) + 
  xlab('Nucleotide divergence') +
  ylab('Niche similarity') +
  scale_color_manual( values = c('lightcoral', 'dodgerblue','grey')) + 
  scale_x_continuous(limits = c(-1,11)) + 
  scale_y_continuous(limits = c(-1,11)) + 
  theme_void()


plot.f <-  plot + inset_element(p = dummy.gg,
                                left =  0.7,
                                bottom = 0.07 ,
                                right =  1,
                                top = 0.3,
                                align_to = 'full')

ggsave(plot = plot.f, 'results/figures/rho_vs_hammdist.pdf',
       width = 10, height = 9)

# If we only keep the ones non-significant, is it significant? 
# NOPE
rho %>% 
  rename('asv' = Partner, 'asv2' = Pair) %>% 
  left_join( hamm, by = c('asv', 'asv2'))  %>% 
  filter(!genus %in% c('Pelagibacter', 'Pelagibacter_A', 'SAR86A')) %>%
  filter(hammdist <= 5) %>% 
  ggplot(aes(hammdist, propr)) + 
  geom_quasirandom(alpha = 0.8)  + 
  geom_point(stat = "summary", fun.y = "mean",
             color = 'red', shape = 21, size = 3) + 
  geom_smooth(method = 'lm', se = T) + 
  scale_x_continuous(breaks = c(1,2,3,4,5)) +
  ylab('Rho proportionality') + 
  xlab('Nucleotide divergence') + 
  ggtitle('No tendency in case of erasing the signigicant ones!') + 
  theme_classic(base_size = 16)

ggsave('results/figures/rho_vs_hammdist_wo_significantones.pdf',
       width = 8, height = 6)



# Some nambers!  ----------------------------------------------------------

# span of numbers 
rho_hamm_genus3asvs %>% 
  group_by(genus) %>% 
  summarise(min = min(propr), 
            mean = mean(propr),
            max = max(propr))

# number of ASVs per groupie 
rho_hamm_genus3asvs %>% 
  select(asv,asv2,genus) %>% 
  pivot_longer(names_to = 'pos', values_to = 'asv', cols = -genus) %>% 
  distinct(genus,asv) %>% 
  group_by(genus) %>% 
  summarise( counts = n()) %>% 
  arrange(-counts)

# Poster version

theme_set(theme_bw(base_size = 25))

plot.p <- ggplot(rho_hamm_genus3asvs, aes(hammdist, propr)) +
  geom_hline(data = syne.dummy, aes(yintercept = y), color = 'red', linetype = 2) + 
  geom_vline(data = syne.dummy, aes(xintercept = x), color = 'red', linetype = 2) + 
  geom_quasirandom(alpha = 0.8, groupOnX = TRUE)  +
  geom_richtext(data = textos, aes(label = text), fill = NA, label.color = NA) +
  # geom_point(stat = "summary", fun.y = "mean_se", color = 'red' ) +
  stat_smooth(aes(group = genus,
                  color = significant),
              method = "gam", formula = y ~ x, show.legend = FALSE) + 
  facet_wrap(~genus, ncol = 4) +
  scale_x_continuous(breaks = c(1,2,3,4,5)) +
  scale_color_manual(values = c('grey', 'dodgerblue')) +
  ylab('Rho proportionality') +
  xlab('Nucleotide divergence')  + 
  lil.strip + 
  theme(strip.text.x = element_markdown()) 

plot.p

ggsave(plot = plot.p, 'results/figures/rho_vs_hammdist_poster.pdf',
       width = 15, height = 8)
