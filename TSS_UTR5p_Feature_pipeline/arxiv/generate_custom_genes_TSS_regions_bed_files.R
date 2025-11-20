gc()
rm(list = ls())

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

##### note
# downgrade dbplyr if you get the error
# Error in `collect()`: ! Failed to collect lazy table. Caused by error in `db_collect()`: 
# ! Arguments in `...` must be used. 
# Problematic argument: • ..1 = Inf 
# Did you misspell an argument name? Run `rlang::last_trace()` to see where the error occurred.
# install.packages('https://cran.r-project.org/src/contrib/Archive/dbplyr/dbplyr_2.3.4.tar.gz', repos = NULL)

ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl", mirror = "www")
tss_data <- getBM(
  attributes = c("chromosome_name", "transcription_start_site", 
                 "strand", "external_gene_name", "ensembl_gene_id"),
  mart = ensembl
)

bed <- tss_data
bed$chromosome_name <- paste0("chr", bed$chromosome_name)
bed$start <- bed$transcription_start_site - 1  # BED is 0-based
bed$end <- bed$transcription_start_site

bed_formatted <- bed[, c("chromosome_name", "start", "end", "external_gene_name", "ensembl_gene_id", "strand")]
colnames(bed_formatted) <- c("chr", "start", "end", "name", "gene_id", "strand")

# tssdf.hg38 <- read.csv(path.to.TSS.file, sep = "\t", header = FALSE)[c("V2", "V3", "V4", "V7", "V12")]
tssdf.hg38 <- bed_formatted[, c("chr", "start", "end", "strand", "name")]
colnames(tssdf.hg38) <- c("chrom", "start", "end", "strand", "gene")
tssdf.hg38$strand <- to_vec(
  for (item in tssdf.hg38$strand){
    if (item == 1){
      return("+")
    } else {
      return("-")
    }
  }
)
tssdf.hg38 <- tssdf.hg38 %>% rowwise() %>%
  mutate(tssName = sprintf("%s_%s_%s", chrom, start, end))
all.chroms <- to_vec(for (i in seq(1,22)) sprintf("chr%s", i))
tssdf.hg38 <- subset(tssdf.hg38, tssdf.hg38$chrom %in% all.chroms) %>% arrange(chrom)

# finished generating tssdf.hg38

# liftOver tssdf.hg38 from hg38 to hg19
library(liftOver)
library(GenomicRanges)
path.to.chain.file <-  "/home/hieunguyen/src/ecd_wgs_enriched_and_kmer_features/resources/hg38ToHg19.over.chain"

tssdf.hg38.grange <- makeGRangesFromDataFrame(df = tssdf.hg38, seqnames.field = "chrom", start.field = "start", end.field = "end", keep.extra.columns = TRUE)
chain <- import.chain(path.to.chain.file)
tssdf.hg19 <- liftOver(tssdf.hg38.grange, chain) %>% as.data.frame() %>%
  subset(select = c(seqnames, start, end, strand, gene, tssName))
tssdf.hg19 <- tssdf.hg19 %>% rowwise() %>% 
  mutate(tssName = sprintf("%s_%s_%s", seqnames, start, end))
colnames(tssdf.hg19) <- c("chrom",
                          "start",
                          "end",
                          "strand",
                          "gene",
                          "tssName")

# cancer.type <- "Lung"
# list.version <- "v0.1"
# list.version <- "v0.2"
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
  
  for (cancer.type in names(gene.list)){
    for (list.version in names(gene.list[[cancer.type]])){
      path.to.assets <- "/home/hieunguyen/src/ecd_wgs_enriched_and_kmer_features/assets/TSS_Feature_pipeline"
      path.to.save.output <- file.path(path.to.assets, "up_down_regulated_genes", cancer.type, list.version)
      dir.create(path.to.save.output, showWarnings = FALSE, recursive = TRUE)
      
      full.genes <- c(gene.list[[cancer.type]][[list.version]]$up, gene.list[[cancer.type]][[list.version]]$down)
      up.genes <- gene.list[[cancer.type]][[list.version]]$up
      down.genes <- gene.list[[cancer.type]][[list.version]]$down
      
      # choose tssdf.hg19 for generating tssdf at custom selected genes. 
      tssdf.up <- subset(tssdf.hg19, tssdf.hg19$gene %in% up.genes)
      tssdf.down <- subset(tssdf.hg19, tssdf.hg19$gene %in% down.genes)
      tssdf.full <- subset(tssdf.hg19, tssdf.hg19$gene %in% full.genes)
      
      promoterdf <- define_promoter_regions(tssdf = tssdf.up, 
                                            up.flank = 1000, 
                                            down.flank = 1000, 
                                            outputdir = file.path(path.to.save.output, "biomart", "up_genes"))
      promoterdf <- merge(promoterdf, tssdf.up, by.x = "promoter_tssName", by.y = "tssName")
      write.table(promoterdf, file.path(path.to.save.output, "biomart", "up_genes", "full_info.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
      
      promoterdf <- define_promoter_regions(tssdf = tssdf.down, 
                                            up.flank = 1000, 
                                            down.flank = 1000, 
                                            outputdir = file.path(path.to.save.output, "biomart", "down_genes"))
      promoterdf <- merge(promoterdf, tssdf.down, by.x = "promoter_tssName", by.y = "tssName")
      write.table(promoterdf, file.path(path.to.save.output, "biomart", "down_genes", "full_info.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
      
      promoterdf <- define_promoter_regions(tssdf = tssdf.full, 
                                            up.flank = 1000, 
                                            down.flank = 1000, 
                                            outputdir = file.path(path.to.save.output, "biomart", "all_genes"))
      promoterdf <- merge(promoterdf, tssdf.full, by.x = "promoter_tssName", by.y = "tssName")
      write.table(promoterdf, file.path(path.to.save.output, "biomart", "all_genes", "full_info.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
      print(nrow(promoterdf))
    }
  }
  # EOF
}
