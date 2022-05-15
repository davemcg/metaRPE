# build trackDb

library(tidyverse)
library(glue)
library(DESeq2)
base_dir <- '~/data/metaRPE_RSE/'
load('../data/rse_gene_select.Rdata')

# set XY ratio
colData(rse_gene_select)$XYratio <- (colData(rse_gene_select)$"recount_qc.aligned_reads%.chrx" + 0.01) / (colData(rse_gene_select)$"recount_qc.aligned_reads%.chry" + 0.01)

# make merged cell line / owner field
colData(rse_gene_select)$CL_Owner <- paste(colData(rse_gene_select)$cell_line, colData(rse_gene_select)$Owner, sep = '_')

# bigwig files
## in biowulf2
## attached as SMB....so this crawl is slow
bw_files <- list.files('/Volumes/data/projects/nei/bharti/metaRPE/pump_output/', full.names = TRUE, recursive = TRUE)
bw_files <- grep('local.all.bw', bw_files, value = TRUE)

bw_meta <- bw_files %>% enframe(value = 'bw_path') %>% 
  mutate(base_path = basename(bw_path),
         external_id = basename(bw_path) %>% 
           gsub('\\!.*','',.)) %>% 
  left_join(colData(rse_gene_select) %>% as_tibble(), by = 'external_id') %>% 
  filter(!is.na(rail_id))


# build track DB
bw_meta %>% mutate(CL_Owner = gsub(' ', '_', CL_Owner),
                   trackDB = case_when(Class == 'Control' ~ 
                                         glue("track {CL_Owner}_{lane_sample}
parent Control
bigDataUrl {base_path}
shortLabel {CL_Owner}_{lane_sample}
longLabel {CL_Owner}_{lane_sample}
color 100,50,50
type bigWig\n\n"),
                   TRUE ~ glue("track {CL_Owner}_{lane_sample}
parent Degeneration
bigDataUrl {base_path}
shortLabel {CL_Owner}_{lane_sample}
longLabel {CL_Owner}_{lane_sample}
color 50,100,100
type bigWig\n\n"))) %>% 
  pull(trackDB) %>% 
  writeLines()

