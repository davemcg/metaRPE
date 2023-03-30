library(tidyverse)
library(SummarizedExperiment)

load('data/rse_gene_03.Rdata')

base_dir <- '~/data/metaRPE_RSE/'
gff <- rtracklayer::import.gff(paste0(base_dir, '/human/annotations/gene_sums/human.gene_sums.G029.gtf.gz')) %>% as_tibble()



colData(rse_gene)$column_field <- paste(colData(rse_gene)$Owner, colData(rse_gene)$Sample, colData(rse_gene)$Group, colData(rse_gene)$Class, colData(rse_gene)$cell_line, sep = '__')


samples <- row.names(colData(rse_gene) %>% data.frame() %>% filter(study == 'mitra'))
sample_names <- colData(rse_gene) %>% data.frame()%>% filter(study == 'mitra') %>% pull(column_field)


matrix <- assay(rse_gene)[, samples]
colnames(matrix) <- sample_names

Group 1
Abca4 ko including  (C35, C26, C46)- unfed VS control (best4c)- unfed  
Group 2
Abca4 ko including  (C35, C26, C46)- pos fed VS control (best4c)- fed  
Group 3
Stargardt patient (step)- unfed vs uneffaced control (F24FC and D3C)- unfed
Group 4
Stargardt patient (stap)- posfed vs uneffaced control (F24FC and D3C)- posfed
Best

matrix %>% 
  as_tibble(rownames = 'gene_id') %>% 
  left_join(gff %>% dplyr::select(type, gene_id, gene_type, gene_name)) %>% 
  relocate(type, gene_id, gene_type, gene_name) %>% 
  pivot_longer(cols = contains("MF"), names_to = 'Sample', values_to = 'counts') %>% 
  filter(gene_name == 'PIKFYVE') %>% separate(Sample, c("scientist","sample",'group','class','line'), '__') %>% 
  separate(group, c('fed','genotype'), '_') %>% 
  mutate(fed = case_when(fed == 'Control' ~ 'unfed', TRUE ~ 'fed')) %>% 
  mutate(Genotype = case_when(grepl("^C\\d+|BEST", line) ~ 'ABCA4KO',
                                 TRUE ~ 'Stargardt')) %>% 
  mutate(condition = case_when(line == 'BEST4C' ~ 'Control',
                             grepl("^C\\d+", line) ~ 'Experiment',
                             line == 'STAP' ~ 'Patient',
                             TRUE ~ 'Control')) %>% 
  ggplot(aes(x=condition, 
             y = log2(counts), 
             group = condition, label = sample, color = line)) + 
  geom_boxplot() + 
  geom_point() +
  ggrepel::geom_label_repel() + 
  scale_color_manual(values = pals::glasbey()) + cowplot::theme_cowplot() +
  facet_wrap(~interaction(fed, Genotype), scales = 'free_x')
