library(tidyverse)

args = commandArgs(trailingOnly=TRUE)
species = args[1]

ncbi <- read_delim(paste0("./resources/",species,"/ncbi.txt"), "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
names(ncbi) <- c("id","entrez","symbol")
ncbi <- ncbi %>% select(symbol, entrez)
ncbi <- ncbi %>% filter(!(ncbi %>% duplicated()))

gencode <- read_delim(paste0("./resources/",species,"/gencode.txt") , "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
names(gencode) <- c("ensembl","symbol","havana")

ensembl <- read_delim(paste0("./resources/",species,"/ensembl.txt"), "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
names(ensembl) <- c("ensembl","symbol","havana")

join1 <- full_join(gencode,ncbi, by = c("symbol"))

join3 <- full_join(join1, ensembl, by = c("ensembl"))

kk <- join3  %>% select(-symbol.y, -havana.y)
names(kk) <- c("ensembl","symbol","havana","entrez")

saveRDS(kk, paste0("./resources/",species,"/geneAnnotation.RDS") )