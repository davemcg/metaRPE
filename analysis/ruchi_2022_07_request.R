library(tidyverse)
library(SummarizedExperiment)

load('data/rse_gene_03.Rdata')

base_dir <- '~/data/metaRPE_RSE/'
gff <- rtracklayer::import.gff(paste0(base_dir, '/human/annotations/gene_sums/human.gene_sums.G029.gtf.gz')) %>% as_tibble()



colData(rse_gene)$column_field <- paste(colData(rse_gene)$Owner, colData(rse_gene)$Sample, colData(rse_gene)$CellType, colData(rse_gene)$Class, colData(rse_gene)$cell_line, sep = '__')


samples <- row.names(colData(rse_gene) %>% data.frame() %>% dplyr::select(Sample:DrugSimple) %>% filter(DrugSimple != 'Drug' | !grepl('^HS',serum), Class == 'Control' ))
samples <- c(samples, 'AG4','AG5', 'AG2', 'AG1') # ruchi request these specifically
sample_names <- colData(rse_gene)[samples,]  %>% data.frame() %>% pull(column_field)


matrix <- assay(rse_gene)[, samples]
colnames(matrix) <- sample_names


matrix %>% as_tibble(rownames = 'gene_id') %>% left_join(gff %>% dplyr::select(type, gene_id, gene_type, gene_name)) %>% relocate(type, gene_id, gene_type, gene_name) %>% 
  write_csv(file = 'data/ruchi_control_counts_04.csv')
