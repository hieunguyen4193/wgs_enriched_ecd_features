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
dir.create(path.to.02.output, showWarnings = FALSE, recursive = TRUE)

# input.tissue.type <- "Normal"
input.tissue.type <- "Tumor"
cutoff.adjp <- 0.05
topN <- 500

if (file.exists(file.path(path.to.02.output, sprintf("%s_topN_%s_up_genes.rds", input.tissue.type, topN))) == FALSE){
  train.dataset <- readRDS(file.path(path.to.00.output, "train_dataset.rds"))
  val.dataset <- readRDS(file.path(path.to.00.output, "val_dataset.rds"))
  
  to.run.metadata <- subset(full.dataset.metadata, full.dataset.metadata$Tissue.Type == input.tissue.type &
                              full.dataset.metadata$Sample.ID %in% unlist(train.dataset[[input.tissue.type]]))
  to.run.metadata <- to.run.metadata[!duplicated(to.run.metadata$Sample.ID), ]
  to.run.metadata <- to.run.metadata %>% arrange(desc(cohort))
  
  #------------------------------------------------------------------------------#
  # ***** read in the data of gene expression of all training samples *****
  #------------------------------------------------------------------------------#
  i <- 1
  input.path <- to.run.metadata$path[[i]]
  mat <- read.table(input.path, 
                    sep = "\t", 
                    header = TRUE)
  tx2gene <- mat[, c("gene_name")]
  sample.id <- to.run.metadata$Sample.ID[[i]]
  mat <- tail(mat, nrow(mat) - 4)
  mat <- subset(mat, select = c(gene_name, unstranded))
  mat <- mat %>% group_by(gene_name) %>% summarise(count = sum(unstranded))
  colnames(mat) <- c("gene_name", sample.id)  
  mat <- subset(mat, mat[[sample.id]] >= 100)
  
  for (i in seq(2, nrow(to.run.metadata))){
    print(i)
    input.path <- to.run.metadata$path[[i]]
    tmp.mat <- read.table(input.path, 
                          sep = "\t", 
                          header = TRUE)
    tx2gene <- tmp.mat[, c("gene_name")]
    sample.id <- to.run.metadata$Sample.ID[[i]]
    tmp.mat <- tail(tmp.mat, nrow(tmp.mat) - 4)
    tmp.mat <- subset(tmp.mat, select = c(gene_name, unstranded))
    tmp.mat <- tmp.mat %>% group_by(gene_name) %>% summarise(count = sum(unstranded))
    colnames(tmp.mat) <- c("gene_name", sample.id)  
    tmp.mat <- subset(tmp.mat, tmp.mat[[sample.id]] >= 100)
    mat <- merge(mat, tmp.mat, by.x = "gene_name", by.y = "gene_name")
  }
  
  row.names(mat) <- NULL
  mat <- mat %>% column_to_rownames("gene_name")
  
  #------------------------------------------------------------------------------#
  # ***** get list of DE genes of each pairwise comparison *****
  #------------------------------------------------------------------------------#
  all.res.files <- Sys.glob(file.path(path.to.01.output, input.tissue.type, "*", "*sig*"))
  cohort.up.genes <- list()
  cohort.down.genes <- list()
  
  for (input.cohort in all.cohorts){
    if (input.cohort == "TCGA-LUSC"){
      tmp.all.res.files <- all.res.files[grepl("LUAD", all.res.files) == FALSE]
    } else if (input.cohort == "TCGA-LUAD"){
      tmp.all.res.files <- all.res.files[grepl("LUSC", all.res.files) == FALSE]
    } else {
      tmp.all.res.files <- all.res.files
    }
    selected.res <- tmp.all.res.files[grepl(input.cohort, tmp.all.res.files)]
    resdf <- list()
    up.genes <- list()
    down.genes <- list()
    
    for (file in selected.res){
      filenames <- str_split(basename(dirname(file)),"_vs_")[[1]]
      tmpdf <- readxl::read_excel(file)
      tmpdf$abs_log2FoldChange <- abs(tmpdf$log2FoldChange)
      if (input.cohort == filenames[[2]]){
        tmpdf$log2FoldChange <- -1 * tmpdf$log2FoldChange
      }
      tmpdf <- tmpdf %>% arrange(desc(abs_log2FoldChange))
      resdf[[basename(dirname(file))]] <- tmpdf
      up.genes[[basename(dirname(file))]] <- subset(tmpdf, tmpdf$log2FoldChange < 0 ) %>% head(topN) %>% pull(gene_name)
      down.genes[[basename(dirname(file))]] <- subset(tmpdf, tmpdf$log2FoldChange > 0 ) %>% head(topN) %>% pull(gene_name)
    }
    
    cohort.up.genes[[input.cohort]] <- Reduce(intersect, lapply(up.genes, na.omit))
    cohort.down.genes[[input.cohort]] <- Reduce(intersect, lapply(down.genes, na.omit))
  }
  
  plot.genes <- c()
  for (input.cohort in unique(to.run.metadata$cohort)){
    plot.genes <- c(plot.genes, cohort.up.genes[[input.cohort]])
    # plot.genes <- c(plot.genes, cohort.down.genes[[input.cohort]])
  }
  
  # generate dds object for input matrix. 
  dds.metadata <- data.frame(sample = colnames(mat))
  dds.metadata <- merge(dds.metadata, to.run.metadata[, c("Sample.ID", "cohort")], by.x = "sample", by.y = "Sample.ID")
  dds <- DESeqDataSetFromMatrix(countData = mat,
                                colData = dds.metadata,
                                design = ~ cohort)
  dds <- dds[rowSums(counts(dds)) > 100, ]
  vst.dds <- vst(dds)
  
  de.mat <- assay(vst.dds)
  de.mat <- de.mat[intersect(plot.genes, row.names(mat)), to.run.metadata$Sample.ID]
  
  library(viridis)
  library(RColorBrewer)
  pal <- brewer.pal(9, "YlOrRd") 
  color_gradient <- colorRampPalette(pal)(100)
  pheatmap::pheatmap(de.mat, 
                     cluster_rows = FALSE, 
                     cluster_cols = FALSE, 
                     show_rownames = FALSE, 
                     show_colnames = FALSE,
                     color = rev(color_gradient), 
                     scale = "none")
  
  saveRDS(cohort.down.genes, file.path(path.to.02.output, sprintf("%s_topN_%s_down_genes.rds", input.tissue.type, topN)))
  saveRDS(cohort.up.genes, file.path(path.to.02.output, sprintf("%s_topN_%s_up_genes.rds", input.tissue.type, topN)))
} else {
  print("File exists")
}
