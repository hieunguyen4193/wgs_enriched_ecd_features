gc()
rm(list = ls())

library(tidyverse)
library(dplyr)
library(ggplot2)
library(stringr)

path.to.main.src <- "/home/hieunguyen/src/wgs_enriched_ecd_features"

path.to.main.output <- file.path(path.to.main.src, "assets", "TSS_UTR5p_regions")
path.to.01.output <- file.path(path.to.main.output, "01_output")
path.to.02.output <- file.path(path.to.main.output, "02_output")
path.to.03.output <- file.path(path.to.main.output, "03_output")

for (input.list.version in c("Tumor", "Normal", "TUMOR_minus_NORMAL")){
  path.to.04.output <- file.path(path.to.main.output, "04_output", input.list.version)
  
  source(file.path(path.to.main.src, "TSS_UTR5p_Feature_pipeline" ,"TSS_UTR5p_helper_functions.R"))
  source(file.path(path.to.main.src, "TSS_UTR5p_Feature_pipeline", "helper_functions.R"))
  
  dir.create(file.path(path.to.04.output, "TSS_beds"), showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(path.to.04.output, "UTR5_beds"), showWarnings = FALSE, recursive = TRUE)
  
  # ***** manual curated gene lists *****
  
  if (input.list.version == "v0.1"){
    source(file.path(path.to.main.src, "TSS_UTR5p_Feature_pipeline", "gene_lists_v0.1.R"))  
  } else if (input.list.version == "v0.2"){
    gene.list <- readRDS(file.path(path.to.main.src, "TSS_UTR5p_Feature_pipeline", "gene_lists_v0.2.rds"))
  } else if (grepl("RANDOM", input.list.version) == TRUE){
    gene.list <- readRDS(file.path(path.to.main.src, "TSS_UTR5p_Feature_pipeline", sprintf("RANDOM_gene_lists_v%s.rds", i)))
  } else if (grepl("Tumor", input.list.version) == TRUE){
    path.to.TOO.TCGA.output <- "/home/hieunguyen/src/wgs_enriched_ecd_features/assets/TCGA_bulkRNAseq/03_output/all_TOO_gene_lists.rds"
    gene.list <- readRDS(path.to.TOO.TCGA.output)[["Tumor"]]
  } else if (grepl("Normal", input.list.version) == TRUE){
    path.to.TOO.TCGA.output <- "/home/hieunguyen/src/wgs_enriched_ecd_features/assets/TCGA_bulkRNAseq/03_output/all_TOO_gene_lists.rds"
    gene.list <- readRDS(path.to.TOO.TCGA.output)[["Normal"]]
  }
  
  # reading in saved TSS dataframe. 
  tssdf <- read.csv(file.path(path.to.01.output, "tssdf.hg19.csv")) %>%
    subset(select = -c(X))
  utr5pdf <- read.csv(file.path(path.to.03.output, "5prime_UTR.hg19.bed"), sep = "\t", header = FALSE)[, c("V1", "V2", "V3", "V6")]
  colnames(utr5pdf) <- c("chrom", "start", "end", "gene")
  
  for (cancer.type in names(gene.list)){
    for (list.version in names(gene.list[[cancer.type]])){
      full.genes <- c(gene.list[[cancer.type]][[list.version]]$up, gene.list[[cancer.type]][[list.version]]$down)
      up.genes <- gene.list[[cancer.type]][[list.version]]$up
      down.genes <- gene.list[[cancer.type]][[list.version]]$down
      
      # generate TSS regions of selected genes. 
      generate_tss_regions_for_input_genes(tssdf = tssdf, 
                                           input.genes = full.genes, 
                                           outputdir = file.path(path.to.04.output, 
                                                                 "TSS_beds", 
                                                                 cancer.type, 
                                                                 list.version, 
                                                                 "full_genes"), 
                                           up.flank = 1000, 
                                           down.flank = 1000)
      generate_tss_regions_for_input_genes(tssdf = tssdf, 
                                           input.genes = up.genes, 
                                           outputdir = file.path(path.to.04.output, 
                                                                 "TSS_beds", 
                                                                 cancer.type, 
                                                                 list.version, 
                                                                 "up_genes"), 
                                           up.flank = 1000, 
                                           down.flank = 1000)
      generate_tss_regions_for_input_genes(tssdf = tssdf, 
                                           input.genes = down.genes, 
                                           outputdir = file.path(path.to.04.output, 
                                                                 "TSS_beds", 
                                                                 cancer.type, 
                                                                 list.version, 
                                                                 "down_genes"), 
                                           up.flank = 1000, 
                                           down.flank = 1000)
      
      # generate UTR5p regions of selected genes
      utrdf.up <- subset(utr5pdf, utr5pdf$gene %in% up.genes)
      utrdf.down <- subset(utr5pdf, utr5pdf$gene %in% down.genes)
      utrdf.full <- subset(utr5pdf, utr5pdf$gene %in% full.genes)
      
      dir.create(file.path(path.to.04.output, "UTR5_beds", cancer.type, list.version, "all_genes"), showWarnings = FALSE, recursive = TRUE)
      dir.create(file.path(path.to.04.output, "UTR5_beds", cancer.type, list.version, "up_genes"), showWarnings = FALSE, recursive = TRUE)
      dir.create(file.path(path.to.04.output, "UTR5_beds", cancer.type, list.version, "down_genes"), showWarnings = FALSE, recursive = TRUE)
      
      write.table(utrdf.up, file.path(path.to.04.output, "UTR5_beds", cancer.type, list.version, "up_genes", "5prime_UTR_up_genes.hg19.bed"), quote = FALSE,
                  sep = "\t", row.names = FALSE, col.names = FALSE)
      write.table(utrdf.down, file.path(path.to.04.output, "UTR5_beds", cancer.type, list.version, "down_genes", "5prime_UTR_down_genes.hg19.bed"), quote = FALSE,
                  sep = "\t", row.names = FALSE, col.names = FALSE)
      write.table(utrdf.full, file.path(path.to.04.output, "UTR5_beds", cancer.type, list.version, "all_genes", "5prime_UTR_all_genes.hg19.bed"), quote = FALSE,
                  sep = "\t", row.names = FALSE, col.names = FALSE)
    }
  }
}



