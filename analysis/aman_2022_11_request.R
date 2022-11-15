library(tidyverse)
library(SummarizedExperiment)

load('data/rse_gene_03.Rdata')

base_dir <- '~/data/metaRPE_RSE/'
gff <- rtracklayer::import.gff(paste0(base_dir, '/human/annotations/gene_sums/human.gene_sums.G029.gtf.gz')) %>% as_tibble()

samples <- c('KB3','15_RS','1_RS')


matrix <- assay(rse_gene)[, samples]
colnames(matrix) <- samples


matrix %>% as_tibble(rownames = 'gene_id') %>% left_join(gff %>% dplyr::select(type, gene_id, gene_type, gene_name)) %>% relocate(type, gene_id, gene_type, gene_name) %>% 
  write_csv(file = 'data/aman_request_counts_2022_11_14.csv')
