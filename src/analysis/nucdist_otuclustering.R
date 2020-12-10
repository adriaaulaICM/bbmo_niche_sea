library(tidyverse)
# Packages that are required but not loaded:
library(speedyseq)
library(DECIPHER)
library(Biostrings)

nproc <- 4 # set to number of cpus/processors to use for the clustering


args <- commandArgs(trailingOnly = TRUE)

seq.loc <- args[1]
name <- args[2]

seqtab <- readRDS(seq.loc)
asv.names <- tibble(asv = str_c('asv', 1:ncol(seqtab)),
                    seq = colnames(seqtab))

asv_sequences <- colnames(seqtab)
sample_names <- rownames(seqtab)
dna <- Biostrings::DNAStringSet(asv_sequences)

## Find clusters of ASVs to form the new OTUs
# alignment and distance 
aln <- DECIPHER::AlignSeqs(dna, processors = nproc)
d <- DECIPHER::DistanceMatrix(aln, processors = nproc)

# The clusters 
clusters97 <- DECIPHER::IdClusters(
  d, 
  method = "complete",
  cutoff = 0.03, # use `cutoff = 0.03` for a 97% OTU 
  processors = nproc)

clusters99 <- DECIPHER::IdClusters(
  d, 
  method = "complete",
  cutoff = 0.01, # use `cutoff = 0.03` for a 97% OTU 
  processors = nproc)

clusters99.5 <- DECIPHER::IdClusters(
  d, 
  method = "complete",
  cutoff = 0.005, # use `cutoff = 0.03` for a 97% OTU 
  processors = nproc)

## Use dplyr to merge the columns of the seqtab matrix for ASVs in the same OTU
# prep by adding sequences to the `clusters` data frame
# I have changed the structure this can have consequences in some scripts 
# bewareeee
clust.tibb <- tibble( cluster97 = clusters97,
        cluster99 = clusters99,
        cluster995 = clusters99.5,
        seq = asv_sequences,
        asv = asv.names$asv )

write_rds(clust.tibb,
          str_c('data/analysis/', name, '_otus_decipher.rds'))

d.tibb <- d

d.tibb[lower.tri(d.tibb, diag = TRUE)]  <- NA

colnames(d.tibb) <- asv.names$asv
dist.df <- as_tibble(d.tibb) %>% 
  mutate(asv2 = colnames(d.tibb)) %>% 
  pivot_longer(names_to = 'asv', values_to = 'distance', cols = -asv2) %>% 
  filter(!is.na(distance)) %>% 
  filter(distance <= 0.05)

write_rds(dist.df, str_c('data/analysis/',name,'dech_dist_min95.rds'))

