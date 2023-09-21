# 2023 09 20
# Add three studies to metaRPE
# zander time series iRPE
# mitra metformin
# dominik ciliopathy  
library(tidyverse)
load('data/rse_gene_04.Rdata')
# rse_gene_05 and rse_gene_06 are the same as 04, just with some assay tweaking to accomodate geyser

# ugh, so I'm swapping the quant method from rail-rna/recount3/Snakerail to EiaD 2023's salmon based
# approach, which is way simpler
# but that means re-quantifying everything
# also this approach needs the sample sheet to specify the sample run type (e.g. 3C_D0    D3C_D0_1__HJLT7DSX3_19270591_S73_L002   paired )
# and the metadata in the colData slot doesn't consistently have that....so
# I had to sort of hand rebuild that bit here:
# biowulf2:/data/mcgaugheyd/projects/nei/bharti/metaRPE_eiad/sample_table_2023_09_20.tsv
# most used the existing metadata and searched for their individual (git repo usually) lane info tables (usually from NISC)
# also had to convert sample names to have a leading "X" if they started with a digit (because R)
sample_table <- read_tsv("data/sample_table_2023_09_20.tsv")

existing_meta <- colData(rse_gene) %>% as_tibble() %>% 
  mutate(Sample2=case_when(grepl("^\\d",Sample) ~ paste0("X",Sample), 
                           TRUE ~ Sample),
         Sample2 = gsub("-","_",Sample2))

# yay, all are still there
table(existing_meta$Sample2 %in% sample_table$sample_accession)

# now the new three (?)
## dominik ciliopathy cohort
dominik_meta <- readxl::read_xlsx('data/dominik_ciliopathy_2023/Experiment_layout_metadata.xlsx') %>% 
  dplyr::rename(Sample2 = `Sample name`) %>% 
  mutate( Sample2 = case_when(str_extract(Sample2, '\\d+') %>% as.integer() <= 24 ~ paste0(Sample2, '_2'),
                              TRUE ~ Sample2) ,
          study = 'dominik_ciliopathy',
          Owner = 'dominik',
          CellType = 'iRPE',
          Class = case_when(Disease == 'Healthy Control' ~ 'Control',
                            TRUE ~ "Ciliopathy"),
          Repo = 'zander_bharti_iRPE_differentiation',
          Note = case_when(Sample2 %in% c('DR19','DR20','DR21', 'DR25','DR26','DR27') ~ 
                             'Hand remove, strange phenotype in iRPE'))

## mitra metformin/abca4 KO (the iRPE samples, not the mouse scRNA/bulk)
mitra_meta <- readxl::read_xlsx('data/mitra_metformin_bulkRNA/Bharti Metformin Sample Submission Form (Human RPE cells)-2020.xlsx', 
                                sheet = 'Submission form', skip = 21) %>% 
  dplyr::rename(sample = `Sample ID***`) %>% 
  mutate(sample = gsub("-","_",sample),
         treatment = gsub("^.* ","", `...10`, ) %>% gsub('-','',.),
         experiment = case_when(grepl("LORD", `...10`) ~ "LORD",
                                TRUE ~ "Stargardt_ABCA4"),
         genotype = case_when(grepl("ontrol", `...10`) ~ 'Control',
                              grepl("patient", `...10`) ~ "Patient",
                              grepl("KO", `...10`) ~ 'KO'),
         RIN = gsub("RIN: |RIN:", "", RIN) %>% as.numeric())

mitra_meta_n <- 
  tibble::tribble(
    ~`Sample.ID***`, ~Species,             ~Library.Type, ~`Enrichment.Method.for.RNA-Seq`, ~`Average.Fragment.Size (for.ChIP-Seq.or.other.pre-sheared.samples)`, ~`Sample.Vol.(ul)`, ~`Sample.Conc.(ng/ul)`, ~`Comments/Special.Instructions`,      ~RIN,  ~`Sample.description `,       ~Tissue,
    "MF- RNA 1", "Human ", "Strand-specific RNA-seq",                         "Ploy-A",                                                                   NA,                10L,                   148L,                               NA, "RIN:9.5",  "LORD control POS+MET",   "RPE Cells",
    "MF- RNA 2", "Human ", "Strand-specific RNA-seq",                         "Ploy-A",                                                                   NA,                10L,                   188L,                               NA, "RIN:9.3",   "ABCA4 KO C2 POS+MET",   "RPE Cells",
    "MF- RNA 3", "Human ", "Strand-specific RNA-seq",                         "Ploy-A",                                                                   NA,                10L,                   219L,                               NA, "RIN:9.3",   "ABCA4 KO C2 POS+MET",   "RPE Cells",
    "MF- RNA 4", "Human ", "Strand-specific RNA-seq",                         "Ploy-A",                                                                   NA,                10L,                   150L,                               NA, "RIN:9.3", "Stargardt Control POS",   "RPE Cells"
  ) %>% 
  dplyr::rename(sample = `Sample.ID***`) %>% 
  mutate(sample = gsub("-","_",sample) %>% gsub(" ", "", .),
         treatment = gsub("^.* ","", `Sample.description `, ) %>% gsub('-','',.),
         experiment = case_when(grepl("LORD", `Sample.description `) ~ "LORD",
                                TRUE ~ "Stargardt_ABCA4"),
         genotype = case_when(grepl("ontrol", `Sample.description `) ~ 'Control',
                              grepl("patient", `Sample.description `) ~ "Patient",
                              grepl("KO", `Sample.description `) ~ 'KO'),
         RIN = gsub("RIN: |RIN:", "", RIN) %>% as.numeric()) 

mitra_meta_bind <- bind_rows(mitra_meta %>% mutate(Batch = 'One'),
                             mitra_meta_n %>% mutate(Batch = 'Two')) %>% 
  mutate(Sample2 = case_when(grepl("RNA", sample) ~ paste0(sample, 'metforminRedo'),
                             TRUE ~ paste0(sample, 'metformin')),
         CellType = 'iRPE',
         study = 'mitra_metformin',
         Owner = 'mitra',
         Class = case_when(genotype == 'Control' ~ 'Control',
                           TRUE ~ 'Perturbed'),
         Repo = 'mitra_metformin_bulkRNA')

## zander iRPE differentiation series
zander_samples <- (readxl::read_xlsx('data/zander_bharti_iRPE_differentiation/NISC Submission Form.xlsx',skip=21, sheet = 'Submission Form') %>% 
                     pull(1))[1:36]

zander_meta <- zander_samples %>% enframe(value = 'Sample2') %>% 
  select(-name) %>% 
  separate(Sample2, into = c('Cell.line','time','biological_replicate'), sep = '_', remove = FALSE) %>% 
  mutate(age_days = str_extract(time, '\\d+') %>% as.integer(),
         CellType = 'iRPE',
         study = 'zander',
         Owner = 'zander',
         Class = 'Development',
         Repo = 'zander_bharti_iRPE_differentiation')


new_three <- bind_rows(dominik_meta, mitra_meta_bind, zander_meta) %>% 
  dplyr::rename(Sample = Sample2)

all_together <- bind_rows(existing_meta %>% dplyr::mutate(Sample = Sample2) %>% select(-Sample2),
                          new_three) %>% 
  select(-contains("recount"), -contains("BigW"), -rail_id, -external_id)

new_full_meta <- left_join(sample_table, 
                           all_together,
                           by = c('sample_accession' = 'Sample')) %>% 
  mutate(cell_line2 = case_when(is.na(cell_line) ~ Cell.line,
                               TRUE ~ cell_line)) %>% 
  mutate(cell_line = cell_line2) %>% 
  select(-cell_line2) %>% 
  dplyr::rename(Sample = sample_accession)

# OK!
# time to load in the new counts
# cd ~/data/metaRPE_2023_09_20/
# rsync -Prav h2:/data/mcgaugheyd/projects/nei/bharti/metaRPE_eiad/counts/* .
gene_annot <- data.table::fread('~/data/metaRPE_2023_09_20/gene_anno.csv.gz')
counts <- data.table::fread('~/data/metaRPE_2023_09_20/gene_counts.csv.gz')
tpm <- data.table::fread('~/data/metaRPE_2023_09_20/gene_tpm.csv.gz')

# add human readable info
gene_data <- counts %>% select(Gene) %>% 
  left_join(gene_annot %>% 
              select(gene_id, gene_name, type) %>% 
              group_by(gene_id, gene_name) %>% 
              summarise(type = paste0(unique(type), collapse = ' ,')), 
            by = c("Gene" = "gene_id")) %>% 
  mutate(ngene = paste0(gene_name, ' (', Gene, ')')) 
counts <- counts[,-1] %>% as.matrix()
row.names(counts) <- gene_data$ngene

tpm <- tpm[,-1] %>% as.matrix()
row.names(tpm) <- gene_data$ngene

# ensure metadata is aligned with col names
rse_meta <- colnames(counts) %>% enframe(value = 'Sample') %>% select(-name) %>% left_join(new_full_meta %>% select(-run_accession) %>% unique())

table(rse_meta$Sample == colnames(counts))

# set up gene info
gene_data <- gene_data %>% data.frame()
row.names(gene_data) <- gene_data$ngene

# build it
rse_big <- SummarizedExperiment(assays = list(counts = counts, tpm = tpm),
                            metadata = rse_meta,
                            rowData = gene_data)
rse <- SummarizedExperiment(assays = list(counts = counts), 
                            metadata = rse_meta,
                            rowData = gene_data)
save(rse, file = 'data/rse_gene_07.2023_09_20.Rdata')
save(rse_big, file = 'data/rse_gene_07.2023_09_20.big.Rdata')
write_tsv(new_full_meta, file = 'data/rse_gene_07.full_meta.tsv.gz')