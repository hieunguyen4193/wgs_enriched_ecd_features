gc()
rm(list = ls())

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
names(all.manifests) <- to_vec(for (item in basename(all.manifests)) str_replace(item, ".RNAseq.txt", ""))

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

# ***** MAIN RUN*****
# outdir <- "/media/HNSD01/outdir"
# path.to.main.output <- file.path(outdir, "ecd_tcga_bulk_rnaseq_data_analysis_TOO")

outdir <- "/home/hieunguyen/src/wgs_enriched_ecd_features/assets"
path.to.main.output <- file.path(outdir, "TCGA_bulkRNAseq")
path.to.00.output <- file.path(path.to.main.output, "00_output")
path.to.01.output <- file.path(path.to.main.output, "01_output")
path.to.02.output <- file.path(path.to.main.output, "02_output")
path.to.03.output <- file.path(path.to.main.output, "03_output")
dir.create(path.to.03.output, showWarnings = FALSE, recursive = TRUE)

# input.tissue.type <- "Normal"
input.tissue.type <- "Tumor"
cutoff.adjp <- 0.05
topN <- 500

