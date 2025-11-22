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
path.to.04.output <- file.path(path.to.main.output, "04_output")
path.to.05.output <- file.path(path.to.main.output, "05_output")
dir.create(path.to.05.output, showWarnings = FALSE, recursive = TRUE)

path.to.TOO.TCGA.output <- "/home/hieunguyen/src/wgs_enriched_ecd_features/assets/TCGA_bulkRNAseq/03_output/all_TOO_gene_lists.rds"
gene.list <- readRDS(path.to.TOO.TCGA.output)

all.bed.files <- Sys.glob(file.path(path.to.04.output, "*", "*", "*", "*", "*", "*.bed"))

path.to.hpc.assets <- "/mnt/nvme/DATA_HIEUNGUYEN/resource/assets"
path.to.current.assets <- "/home/hieunguyen/src/wgs_enriched_ecd_features/assets"

maindf <- data.frame(inputbed = all.bed.files) %>%
  rowwise() %>%
  mutate(up_or_down = basename(dirname(inputbed))) %>%
  mutate(version.id = str_split(inputbed, "/")[[1]][[12]]) %>%
  mutate(cancer.type = str_split(inputbed, "/")[[1]][[11]]) %>%
  mutate(tss_or_utr5p = str_split(inputbed, "/")[[1]][[10]]) %>%
  mutate(tumor_or_normal = str_split(inputbed, "/")[[1]][[9]]) %>%
  mutate(run_name = sprintf("%s_%s_%s_%s_%s",
                            up_or_down,
                            version.id,
                            cancer.type,
                            tss_or_utr5p,
                            tumor_or_normal)) %>%
  mutate(inputbed = str_replace(inputbed, path.to.current.assets, path.to.hpc.assets))

write.table(maindf[, c("inputbed", "run_name")], sep = "\t", row.names = FALSE, col.names = TRUE, file.path(path.to.main.src, "bed_list.tsv"), quote = FALSE)
