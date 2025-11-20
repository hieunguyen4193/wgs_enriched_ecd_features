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

# input.cohort <- "TCGA-BRCA"
# input.cohort <- "TCGA-LIHC"
# input.cohort <- "TCGA-STAD"
# input.cohort <- "TCGA-LUSC"
# input.cohort <- "TCGA-COAD"
# input.cohort <- "TCGA-LUAD"

convert.cohort.names <- list(
  "TCGA-BRCA" = "Breast",
  "TCGA-STAD" = "Gastric",
  "TCGA-LUSC" = "Lung",
  "TCGA-LUAD" = "Lung",
  "TCGA-COAD" = "Colorectal",
  "TCGA-LIHC" = "Liver"
)
all.cohorts <- c("TCGA-BRCA", "TCGA-LIHC", "TCGA-STAD", "TCGA-LUSC", "TCGA-LUAD", "TCGA-COAD")

outdir <- "/media/HNSD01/outdir/ecd_tcga_bulk_rnaseq_data_analysis"

list.version <- "v0.2"
path.to.save.output <- "/home/hieunguyen/src/ecd_wgs_enriched_and_kmer_features"

if (file.exists(file.path(path.to.save.output, "TSS_Feature_pipeline", "gene_lists_v0.2.rds")) == FALSE){
  gene.list <- list()
  
  for (input.cohort in all.cohorts){
    path.to.main.output <- file.path(outdir, input.cohort)
    path.to.01.output <- file.path(path.to.main.output, "01_output", sprintf("run_%s", "*"))
    path.to.02.output <- file.path(path.to.main.output, "02_output")
    
    up.genedf <- readxl::read_excel(file.path(path.to.02.output, "up_genedf.xlsx"))
    down.genedf <- readxl::read_excel(file.path(path.to.02.output, "down_genedf.xlsx"))
    cancer.type <- convert.cohort.names[[input.cohort]]
    if (cancer.type %in% names(gene.list) == FALSE){
      gene.list[[cancer.type]] <- list()    
    }
    if (input.cohort %in% c("TCGA-LUSC", "TCGA-LUAD")){
      gene.list[[cancer.type]][[sprintf("%s_%s", list.version, str_replace(input.cohort, "TCGA-", ""))]] <- list(
        up = up.genedf$gene,
        down = down.genedf$gene
      )
    } else {
      gene.list[[cancer.type]][[list.version]] <- list(
        up = up.genedf$gene,
        down = down.genedf$gene
      )    
    }
  }
  
  # save output to the directory where we generate enriched features
  saveRDS(gene.list, file.path(path.to.save.output, "TSS_Feature_pipeline", "gene_lists_v0.2.rds"))
} else {
  print("gene list v0.2 exists, not generate a new one")
}

#####----------------------------------------------------------------------#####
# ***** generate list of random genes *****
#####----------------------------------------------------------------------#####

for (i in seq(1, 10)){
  list.version <- sprintf("random%s", i)
  random.genelist <- list()
  
  for (input.cohort in all.cohorts){
    path.to.main.output <- file.path(outdir, input.cohort)
    path.to.01.output <- file.path(path.to.main.output, "01_output", sprintf("run_%s", "*"))
    path.to.02.output <- file.path(path.to.main.output, "02_output")
    
    up.genedf <- readxl::read_excel(file.path(outdir, "full_genelist.xlsx"))
    down.genedf <- readxl::read_excel(file.path(outdir, "full_genelist.xlsx"))
    
    random.up.genes <- sample(up.genedf$gene, 100)
    random.down.genes <- sample(down.genedf$gene, 100)
    
    up.genedf <- subset(up.genedf, up.genedf$gene %in% random.up.genes)
    down.genedf <- subset(down.genedf, down.genedf$gene %in% random.down.genes)
    
    up.genedf <- up.genedf[, c("gene", "logFC", "padj")]
    down.genedf <- down.genedf[, c("gene", "logFC", "padj")]
    
    cancer.type <- convert.cohort.names[[input.cohort]]
    if (cancer.type %in% names(random.genelist) == FALSE){
      random.genelist[[cancer.type]] <- list()    
    }
    if (input.cohort %in% c("TCGA-LUSC", "TCGA-LUAD")){
      random.genelist[[cancer.type]][[sprintf("%s_%s", list.version, str_replace(input.cohort, "TCGA-", ""))]] <- list(
        up = up.genedf$gene,
        down = down.genedf$gene
      )
    } else {
      random.genelist[[cancer.type]][[list.version]] <- list(
        up = up.genedf$gene,
        down = down.genedf$gene
      )    
    }
  }
  print(random.genelist$Breast %>% head())
  saveRDS(random.genelist, file.path(path.to.save.output, "TSS_Feature_pipeline", sprintf("RANDOM_gene_lists_v%s.rds", i)))
}

