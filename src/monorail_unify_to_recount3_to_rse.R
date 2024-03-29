# import in monorail unify output for OGVFB Bharti RNA-seq samples
library(tidyverse)
library(recount3)
library(DESeq2)
options(recount3_url = '~/data/metaRPE_RSE')

hp <- available_projects()
hp

# Create RSE
## Also get the counts returned because the `raw_counts` matrix is  the sum of the coverage across a feature
## load this custom create_rse script which fixes a little bug I found involving a mistmatched GTF and count matrix
source('src/create_rse_manual.R')
rse_list <- list()
for (i in hp$project){
  print(i)
  rse_list[[i]] <- create_rse2(hp[hp$project == i, ],  type = 'gene', annotation = "gencode_v29")
  assay(rse_list[[i]], "counts") <- transform_counts(rse_list[[i]])
  assay(rse_list[[i]], "tpm") <- recount::getTPM(rse_list[[i]], 
                                                 length_var = 'score',
                                                 mapped_var = 'recount_qc.star.all_mapped_reads') %>% 
    data.frame()
  assay(rse_list[[i]], "counts") <- compute_read_counts(rse_list[[i]])  # fill in UNSCALED read counts (`transform_counts` applies a 40e6 total scaling)
}

# create one rse for all samples
rse_gene <- do.call(cbind, rse_list)



# Import in my metadata and add to colData

lane_info <- read_tsv('data/ogvfb_working_meta_seq_lane.tsv')
meta_info <- read_tsv('data/ogvfb_full_meta.tsv') %>% filter(!Owner %in% c('kelcy','Mitra'))
kelcy <- read_tsv('data/ogvfb_full_meta_02.tsv') %>% filter(Owner == 'kelcy')
mitra <-  read_tsv('data/ogvfb_full_meta_02.tsv') %>% filter(Owner == 'Mitra')
info <- lane_info %>% 
  mutate(Class = case_when(is.na(Class) ~ 'Development', 
                           TRUE ~ Class)) %>% 
  left_join(meta_info %>% 
              select(-cell_line, -CellType, -Treatment, -Owner, -Class, -Repo), 
            by = 'Sample') %>% 
  bind_rows(kelcy %>% dplyr::rename(lane_sample = fastq_file)) %>%  # add kelcy 2022 05
  bind_rows(mitra %>% dplyr::mutate(lane_sample = gsub('\\/quant.sf','',file_path))) %>%  # add mitra 2022 07
  janitor::remove_empty()  # janitor removes columns where all NA


info <- data.frame(info)
row.names(info) <- info$lane_sample
info <- info[row.names(colData(rse_gene)),]
if (!identical(row.names(info), row.names(colData(rse_gene)))){
  print('Sample <-> Fastq Name Issue!!!')
  stop()
}

colData(rse_gene) <- cbind(colData(rse_gene), info)



# Collapse technical (lane) replicates
rse_gene <- collapseReplicates(rse_gene, groupby = colData(rse_gene)$Sample)
save(rse_gene, file = 'data/rse_gene_04.Rdata') # updated to "02" for the kelcy data add on 2022-05-19
                                                # updated to "03" for the mitra add on 2022-07-08
                                                # updated to "04" to add in TPM data on 2022-11-15

########### NOT UPDATED WITH THE KELCY or MITRA DATA UPDATE ##############

# Select Control and Degeneration Samples
## Only retain samples that Kapil labelled as controls or degeneration. I make one "audible" and include Karla's sham shRNA treatment as a control as I'd rather have more studies where was have both controls and degeneration samples. 

keep_samples <- colData(rse_gene) %>% as_tibble() %>% filter(Class %in% c('Degeneration','Control')) %>% pull(Sample)
# add in Karla mock shRNA as control
karla_mock_shRNA <-  colData(rse_gene) %>% as_tibble() %>% filter(Owner == 'karla', Condition == 'Ctrl shRNA Monolayer') %>% pull(Sample)
keep_samples <- c(keep_samples,
                  karla_mock_shRNA)
rse_gene_select <- rse_gene[, keep_samples]
# edit colData Karla mock shRNA to control
colData(rse_gene_select)[karla_mock_shRNA,'Class'] <- c('Control','Control')
# merge cellLine/owner
colData(rse_gene_select)$CL_owner <- paste(colData(rse_gene_select)$cell_line, colData(rse_gene_select)$Owner, sep = '_' )
save(rse_gene_select, file = 'data/rse_gene_select.Rdata')
#################################################################
