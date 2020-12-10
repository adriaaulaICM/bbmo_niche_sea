library(dada2)
library(tidyverse)

cat(paste0('\n',"You are using DADA2 version ", packageVersion('dada2'),'\n'))

cat('################################\n\n')

args <- commandArgs(trailingOnly = TRUE)

seqtab.nochim <- args[1]
output <- args[2]
name <- args[3]
tax_db <- args[4]
tax_db_sp <- args[5]


dir.create(file.path(output, "03_taxonomy"), showWarnings = FALSE)
dir.create(file.path(output, "03_taxonomy", name), showWarnings = FALSE)

output <- paste0(output,"/03_taxonomy/",name,"/")

# Assign taxonomy (general)
set.seed(42) #random  generator necessary for reproducibility

seqtab <- readRDS(seqtab.nochim)

tax <- assignTaxonomy(seqtab,
                      tax_db,
                      multithread=TRUE,minBoot=80)
# minboot: N of idntical bootstraps to gen. identification

print("Taxonomy assigned, to Genus level")

head(unname(tax))

# Assign species
tax.sp <- addSpecies(tax, tax_db_sp, verbose=TRUE, allowMultiple=3)


print("Taxonomy assigned, to Species level!")

# Write to disk
saveRDS(tax.sp, paste0(output, name, "_tax_assignation.rds"))


