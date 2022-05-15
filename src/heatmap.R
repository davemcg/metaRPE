#BiocManager::install('ComplexHeatmap', 'org.Hs.eg.db')
#install.packages('pals')
library(ComplexHeatmap)
library(tidyverse)
library(org.Hs.eg.db)

gene_swapper <- function(gene_id){
  # dds uses the ENSGENE gene id
  ## this returns a tibble
  ## with the gene name
  gene_id %>% enframe(value = 'gene_id') %>% left_join(gff %>% as_tibble())
}

meaner <- function(matrix){
  # normalizes values by the row
  # takes the mean
  # and subtracts the mean from all values
  matrixRM <- rowMeans(matrix)
  matrix - matrixRM
}

make_hm <- function(dds, genes,
                    aggregate = TRUE,
                    mean_correction = 'all',
                    title = '',
                    ...){
  extra_args <- list(...)
  # colData(dds)
  # sort columns (samples) to put the same class (then owner) next to each other
  samples <- colData(dds) %>% as_tibble(rownames = 'sample') %>% arrange(Class, Owner) %>% pull(Sample)
  # reorder dds
  dds <- dds[,samples]



  ha_column = HeatmapAnnotation(df = data.frame(Class = colData(dds)$Class,
                                                Owner = colData(dds)$Owner
                                                ),
                                col = list(Owner = c("karla" = pals::alphabet2(20)[3] %>% unname(),
                                                     "davide" = pals::alphabet2(20)[6] %>% unname(), 
                                                     "qin" = pals::alphabet2(20)[9] %>% unname(),
                                                     "ruchi" = pals::alphabet2(20)[12] %>% unname(),
                                                     "dominik" = pals::alphabet2(20)[16] %>% unname()),
                                           Class = c("Control" = viridis::magma(10)[3], 
                                                     "Degeneration" = viridis::magma(10)[7])))


  heat_mat <- assay(dds, slot = 'counts') %>% data.frame()

  
  # if aggregate is TRUE
  # then merge multiple transcripts into one gene by summing
  if (aggregate){
    heat_mat$Gene <- gene_swapper(row.names(heat_mat))$gene_name
    heat_mat <- heat_mat %>% filter(Gene %in% genes) 
    
    heat_mat <- heat_mat %>% group_by(Gene) %>%
      summarise(across(where(is.numeric), sum)) %>% 
      data.frame()
    
    row.names(heat_mat) <- heat_mat$Gene
  } else {
    heat_mat$Gene <- gene_swapper(row.names(heat_mat))$gene_name
    heat_mat <- heat_mat %>% filter(Gene %in% genes) %>% data.frame()
    row.names(heat_mat) <- heat_mat$Gene
  }

  # remove gene column and turn into matrix (ComplexHeatmap requires this)
  heat_mat <- subset(heat_mat, select = -c(Gene)) %>% as.matrix()

  # two ways to plot
  ## either take the actual expression values and log scale
  ## or take the mean of the values, then subtract all from this
  if (mean_correction == 'all'){
    heat_mat <- meaner(log2(heat_mat+1))
    name = 'Relative Expression\n(log2(Counts))'
    color = pals::coolwarm(20)
  } else if (mean_correction == 'none') {
    heat_mat <- log2(heat_mat+1)
    name = 'log2(Counts)'
    color = viridisLite::viridis(20)
  } else {
    print(paste("Invalid `mean_correction` input", mean_correction ))
    stop()
  }

  # remove "outliers" with using given max val
  if ('max_val' %in% names(extra_args)){
    max_val <- extra_args$max_val
    heat_mat[heat_mat > max_val] <- max_val
  }
  # remove "outliers" with using given min val
  if ('min_val' %in% names(extra_args)){
    min_val <- extra_args$min_val
    heat_mat[heat_mat < min_val] <- min_val
  }
  # finally plot the damn heatmap
  Heatmap(heat_mat,
          column_title = title,
          col=color,
          top_annotation = ha_column,
          name = name,
          show_column_names = FALSE,
          cluster_columns = FALSE)
}
