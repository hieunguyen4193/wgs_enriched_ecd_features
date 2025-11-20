gc()
rm(list = ls())
srcdir <- "/home/hieunguyen/src/ecd_wgs_enriched_and_kmer_features"
path.to.main.src <- "/home/hieunguyen/src/ecd_wgs_enriched_and_kmer_features/TSS_Feature_pipeline"
source(file.path(path.to.main.src, "helper_functions.R"))
library(comprehenr)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

new.bioc.pkgs <- c("liftOver", "biomaRt")
for (pkg in new.bioc.pkgs){
  if (pkg %in% installed.packages() == FALSE){
    install.packages("dplyr")
    BiocManager::install(pkg, update = FALSE)  
  }  
}


library(comprehenr)
library(biomaRt)
library(tidyverse)
library(dplyr)

# input.list.version <- "v0.2"
for (i in seq(1, 10)){
  input.list.version <- sprintf("RANDOM_gene_lists_v%s", i)
  if (input.list.version == "v0.1"){
    # get the manual gene list 
    source(file.path("/home/hieunguyen/src/ecd_wgs_enriched_and_kmer_features/TSS_Feature_pipeline/gene_lists_v0.1.R"))  
  } else if (input.list.version == "v0.2"){
    # generate gene.list from TCGA bulk data output. 
    gene.list <- readRDS("/home/hieunguyen/src/ecd_wgs_enriched_and_kmer_features/TSS_Feature_pipeline/gene_lists_v0.2.rds")
  } else if (grepl("RANDOM", input.list.version) == TRUE){
    gene.list <- readRDS(sprintf("/home/hieunguyen/src/ecd_wgs_enriched_and_kmer_features/TSS_Feature_pipeline/RANDOM_gene_lists_v%s.rds", i))
  }
  
  path.to.save.biomart <- "/media/HNSD01/storage/biomart"
  utrdf.hg19 <- read.table(file.path(path.to.save.biomart, "5prime_UTR.hg19.bed"))[, c("V1", "V2", "V3", "V6")]
  colnames(utrdf.hg19) <- c("chrom", "start", "end", "gene")
  
  for (cancer.type in names(gene.list)){
    for (list.version in names(gene.list[[cancer.type]])){
      full.genes <- c(gene.list[[cancer.type]][[list.version]]$up, gene.list[[cancer.type]][[list.version]]$down)
      up.genes <- gene.list[[cancer.type]][[list.version]]$up
      down.genes <- gene.list[[cancer.type]][[list.version]]$down
      
      utrdf.up <- subset(utrdf.hg19, utrdf.hg19$gene %in% up.genes)
      utrdf.down <- subset(utrdf.hg19, utrdf.hg19$gene %in% down.genes)
      utrdf.full <- subset(utrdf.hg19, utrdf.hg19$gene %in% full.genes)
      
      path.to.save.output <- file.path(srcdir, 
                                       "assets", 
                                       "UTR5p_Feature_pipeline", 
                                       "up_down_regulated_genes", 
                                       cancer.type, 
                                       list.version, 
                                       "biomart")
      dir.create(path.to.save.output, showWarnings = FALSE, recursive = TRUE)
      dir.create(file.path(path.to.save.output, "up_genes"), showWarnings = FALSE, recursive = TRUE)
      dir.create(file.path(path.to.save.output, "down_genes"), showWarnings = FALSE, recursive = TRUE)
      dir.create(file.path(path.to.save.output, "all_genes"), showWarnings = FALSE, recursive = TRUE)
      
      write.table(utrdf.up, file.path(path.to.save.output, "up_genes", "5prime_UTR_up_genes.hg19.bed"), quote = FALSE,
                  sep = "\t", row.names = FALSE, col.names = FALSE)
      write.table(utrdf.down, file.path(path.to.save.output, "down_genes", "5prime_UTR_down_genes.hg19.bed"), quote = FALSE,
                  sep = "\t", row.names = FALSE, col.names = FALSE)
      write.table(utrdf.full, file.path(path.to.save.output, "all_genes", "5prime_UTR_all_genes.hg19.bed"), quote = FALSE,
                  sep = "\t", row.names = FALSE, col.names = FALSE)
    }
  }
}