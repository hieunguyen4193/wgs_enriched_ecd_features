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

cutoff.adjp <- 0.05
# ***** MAIN RUN*****
outdir <- "/media/HNSD01/outdir"
path.to.main.output <- file.path(outdir, "ecd_tcga_bulk_rnaseq_data_analysis_TOO")
path.to.00.output <- file.path(path.to.main.output, "00_output")
dir.create(path.to.00.output, showWarnings = FALSE, recursive = TRUE)
train.dataset <- readRDS(file.path(path.to.00.output, "train_dataset.rds"))
val.dataset <- readRDS(file.path(path.to.00.output, "val_dataset.rds"))

# chck number of normal samples
# table(subset(full.dataset.metadata, full.dataset.metadata$Tissue.Type == "Normal")$cohort)
input.tissue.type <- "Tumor"

comparedf <- data.frame(combn(all.cohorts, 2) %>% t())
colnames(comparedf) <- c("condition1", "condition2")

for (row_i in seq(1, nrow(comparedf))){
  condition1 <- comparedf[row_i, ]$condition1
  condition2 <- comparedf[row_i, ]$condition2
  print(sprintf("--------------------------------------------------"))
  print(sprintf("Working on condition %s vs condition %s", condition1, condition2))
  print(sprintf("--------------------------------------------------"))
  
  path.to.01.output <- file.path(path.to.main.output, 
                                 "01_output", 
                                 input.tissue.type, 
                                 sprintf("%s_vs_%s", condition1, condition2))
  dir.create(path.to.01.output, showWarnings = FALSE, recursive = TRUE)
  
  meta.data <- subset(full.dataset.metadata, 
                      full.dataset.metadata$Tissue.Type == input.tissue.type &
                      full.dataset.metadata$Sample.ID %in% c(train.dataset[[input.tissue.type]][[condition1]],
                                                             train.dataset[[input.tissue.type]][[condition2]]))
  
  to.run.metadata <- subset(meta.data, meta.data$cohort %in% c(condition1, condition2))
  
  print(table(to.run.metadata$cohort))
  
  to.run.metadata <- to.run.metadata[!duplicated(to.run.metadata$Sample.ID), ]
  
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
  
  deseq.metadata <- to.run.metadata[, c("Sample.ID", "cohort")]
  colnames(deseq.metadata) <- c("sample", "condition")
  deseq.metadata$condition <- factor(deseq.metadata$condition, levels = c(condition1, condition2))
  
  dds <- DESeqDataSetFromMatrix(countData = mat,
                                colData = deseq.metadata,
                                design = ~ condition)
  deseq.output <- run_DESeq2_and_preprocess(deseq.dataset = dds, 
                                            thresh.pval = 0.05)
  input.df <- deseq.output$all.resultdf
  input.df <- input.df %>% mutate(abs.log2FoldChange = abs(log2FoldChange))
  input.df <- input.df %>% rowwise() %>%
    mutate(show.gene.name = ifelse(padj < cutoff.adjp, gene_name, NA))
  
  volcano.plot <- ggplot(data=input.df, 
                         aes(x=log2FoldChange, y=-log10(padj), col=sig, label=show.gene.name)) + 
    geom_point() + 
    scale_color_manual(values=c("#c0d2f0", "#f28095")) +
    theme_minimal() +
    geom_vline(xintercept=c(-1, 1), col="#9a9fa6", linetype='dotted') +
    geom_hline(yintercept=-log10(cutoff.adjp), col="#9a9fa6", linetype='dotted') +
    geom_text_repel() +
    ggtitle(sprintf("Volcano plot, fold-change = %s / %s (right/left)", condition2, condition1)) +
    theme_bw() + 
    theme(plot.title = element_text(hjust=0.5, face="bold", size = 12)) +
    xlim(c(-max(input.df$abs.log2FoldChange), max(input.df$abs.log2FoldChange)))
  ggsave(plot = volcano.plot, 
         filename = "volcano_plot.svg", 
         path = path.to.01.output, 
         dpi = 300, width = 14, height = 10)
  
  saveRDS(deseq.output, file.path(path.to.01.output, "deseq_output.rds"))
  writexl::write_xlsx(deseq.metadata, file.path(path.to.01.output, "metadata.xlsx"))
  writexl::write_xlsx(deseq.output$resultdf.sig, file.path(path.to.01.output, "sig.xlsx"))
  
  # prepare and save output for further downstream analysis, i.e. pathway analysis
  input.df <- input.df[, c("gene_name", "pvalue", "log2FoldChange", "padj", "abs.log2FoldChange")]
  colnames(input.df) <- c("gene", "pval", "logFC", "padj", "absLogFC")
  writexl::write_xlsx(input.df, file.path(path.to.01.output, "full_genelist.xlsx"))
}

