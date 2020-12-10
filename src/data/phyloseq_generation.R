library(tidyverse)
library(phyloseq)

fold <- 'data/cleaned/dada2/02_nochimera_mergeruns/bbmo_fl_20032013/'
name <- 'bbmo_fl_20032013_seqtab_final.rds'
otu <- readRDS(str_c(fold,name))

# One sample has a problem in the amplification
# removing it 
otu <- otu[rowSums(otu) > 0,]
  

conversion <- tibble::tribble(
  ~oldname, ~realname,
  "BL1",	"BL030128",
  "BL2",  "BL030325",
  "BL3", "BL030422",
  "BL4", "BL030513",
  "BL5", "BL030714",
  "BL6", "BL030804",
  "BL7", "BL030916",
  "BL8", "BL031216"
)

rownames(otu) <- str_replace_all(rownames(otu),
                                 c( "BL060704" = "BL060705",
                                    "BL061010" = "BL061009",
                                    "BL080311" = "BL080312",
                                    "BL090521" = "BL090512",
                                    "BL100413" = "BL100414",
                                    "BL110704" = "BL110705",
                                    "BL120518" = "BL120511",
                                    "BL131204" = "BL131215"))

rownames(otu) <- c(rownames(otu)[1:123], conversion$realname)
otu <- otu[sort(rownames(otu)),]

tax <- readRDS('data/cleaned/dada2/03_taxonomy/bbmo_gtdb/bbmo_gtdb_tax_assignation.rds') %>% 
  as.data.frame()

envdata <- read_tsv('data/raw/metadata/metadata_Blanes_compact_may2017.tsv') %>% 
  filter(Sample_name %in% rownames(otu)) %>% 
  distinct(Sample_name, .keep_all = T) 
rownames(envdata) <- envdata$Sample_name



# Taxonomy changes  -------------------------------------------------------

tax$seq <- rownames(tax)

colnames(otu) <- str_c('asv', 1:ncol(otu))
rownames(tax) <- str_c('asv', 1:nrow(tax))

# Erasing all Chloroplasts from the tax and otu datasets
tax.final <- subset(tax, genus != 'Chloroplast'| is.na(genus)) 
otu.final <- otu[,rownames(tax.final)]


# Creating phyloseq object ------------------------------------------------

OTU <- otu_table(otu.final, taxa_are_rows = F)
TAX <- tax_table(as.matrix(tax.final))
DAT <- sample_data(envdata)

bl.phy <- phyloseq(OTU, TAX, DAT)
write_rds(bl.phy,'data/cleaned/blphyloseq.rds')


# Creating secondary dataset filtering tax for at least >= 7 samples ----------
# In this dataset we will also apply zCompositions to try to predict when a 
# 0 is not a real 0
library(zCompositions)

bl.phy.filttax <- filter_taxa(bl.phy, function(x) sum(x > 0) >= 7, TRUE)
otu.filt <- as(otu_table(bl.phy.filttax), 'matrix')
otu.zcomp <- cmultRepl(otu.filt, method = "CZM", output = "p-counts") 
# 7% of the values have changed from 0 to 1

bl.phy.zcomp <-  merge_phyloseq(bl.phy.filttax,
                                otu_table(otu.zcomp, taxa_are_rows = TRUE))


write_rds(bl.phy.filttax, 'data/cleaned/blphy_filt_7sams.rds')
write_rds(bl.phy.zcomp, 'data/cleaned/blphy_filt_zcomp.rds')
