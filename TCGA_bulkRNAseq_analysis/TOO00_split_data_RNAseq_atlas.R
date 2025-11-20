gc()
rm(list = ls())

rerun <- TRUE
if (rerun == TRUE){
  ##### prepare full information
  library(comprehenr)
  library(DESeq2)
  library(rjson)
  
  
  new.pkgs <- c("tidyverse")
  if ("tidyverse" %in% installed.packages() == FALSE){
    install.packages("tidyverse")
  }
  
  new.bioc.pkgs <- c("biomaRt")
  for (pkg in new.bioc.pkgs){
    if (pkg %in% installed.packages() == FALSE){
      BiocManager::install(pkg, update = FALSE)    
    }
  }
  library(tidyverse)
  library(dplyr)
  library(biomaRt)
  library(ggrepel)
  path.to.main.src <- "/home/hieunguyen/src/ecd_tcga_bulk_rnaseq_data_analysis"
  source(file.path(path.to.main.src, "helper_functions.R"))
  
  manifestdir <- file.path(path.to.main.src, "rna_manifests")
  
  all.manifests <- Sys.glob(file.path(manifestdir, "*.txt"))
  names(all.manifests) <- to_vec(
    for (item in basename(all.manifests)) str_replace(item, ".RNAseq.txt", ""))
  
  mdf <- data.frame()
  for (i in seq(1, length(all.manifests))){
    manifest.name <- names(all.manifests)[[i]]
    tmp.mdf <- read.csv(file.path(all.manifests[[i]]), sep = "\t")
    tmp.mdf$cohort <- manifest.name  
    mdf <- rbind(mdf, tmp.mdf)
  }
  
  main.data.dir <- "/media/HNSD03/raw_data/ecd_tcga_bulk_rnaseq_data_analysis/raw_from_tcga"
  glob.files <- Sys.glob(file.path(main.data.dir, "/*/*/*.tsv"))
  
  dataset.metadata <- data.frame(
    path = glob.files,
    filename = basename(glob.files)
  )
  
  full.dataset.metadata <- merge(mdf, dataset.metadata, by.x = "filename", by.y = "filename")  %>%
    rowwise() %>%
    mutate(sampleid = str_replace(filename, ".rna_seq.augmented_star_gene_counts.tsv", ""))
  
  sample.metadata <- read.csv(file.path(path.to.main.src, "gdc_sample_sheet.2025-10-13.tsv"), sep = "\t")
  
  full.dataset.metadata <- merge(full.dataset.metadata, sample.metadata, 
                                 by.x = "filename", by.y = "File.Name")
  all.cohorts <- unique(full.dataset.metadata$cohort)
  
  cutoff.adjp <- 0.05
  # ***** MAIN RUN*****
  outdir <- "/media/HNSD01/outdir"
  path.to.main.output <- file.path(outdir, "ecd_tcga_bulk_rnaseq_data_analysis_TOO")
  
  
  # chck number of normal samples
  # table(subset(full.dataset.metadata, full.dataset.metadata$Tissue.Type == "Normal")$cohort)
  
  comparedf <- data.frame(combn(all.cohorts, 2) %>% t())
  colnames(comparedf) <- c("condition1", "condition2")
  
  path.to.00.output <- file.path(path.to.main.output, "00_output")
  dir.create(path.to.00.output, showWarnings = FALSE, recursive = TRUE)
  
  countdf <- table(full.dataset.metadata$Tissue.Type, full.dataset.metadata$cohort) %>%
    t() %>% data.frame() %>%
    pivot_wider(names_from = "Var1", values_from = "Freq")
  writexl::write_xlsx(countdf, file.path(path.to.00.output, "count_sample.xlsx"))
  
  train.dataset <- list()
  val.dataset <- list()
  
  rate <- 0.2
  for (input.tissue.type in unique(full.dataset.metadata$Tissue.Type)){
    train.dataset[[input.tissue.type]] <- list()
    val.dataset[[input.tissue.type]] <- list()
    for (input.cohort in all.cohorts){
      all.samples <- subset(full.dataset.metadata, 
                            full.dataset.metadata$cohort == input.cohort &
                              full.dataset.metadata$Tissue.Type == input.tissue.type)$Sample.ID
      if (input.tissue.type == "Normal"){
        val.samples <- sample(all.samples, floor(0.2 * length(all.samples)))
        train.samples <- setdiff(all.samples, val.samples)      
      } else {
        train.samples <- sample(all.samples, 300)
        val.samples <- setdiff(all.samples, train.samples)
      }
      train.dataset[[input.tissue.type]][[input.cohort]] <- train.samples
      val.dataset[[input.tissue.type]][[input.cohort]] <- val.samples
      
      print(sprintf("Cohort %s, tissue type %s, Num Train = %s",
                    input.cohort, input.tissue.type, length(train.samples)))
      print(sprintf("Cohort %s, tissue type %s, Num Val = %s",
                    input.cohort, input.tissue.type, length(val.samples)))
    }
  }
  saveRDS(train.dataset, file.path(path.to.00.output, "train_dataset.rds"))
  saveRDS(val.dataset, file.path(path.to.00.output, "val_dataset.rds"))
}
