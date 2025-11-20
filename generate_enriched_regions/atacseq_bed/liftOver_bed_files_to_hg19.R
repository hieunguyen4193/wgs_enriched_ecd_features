gc()
rm(list = ls())
if ("liftOver" %in% installed.packages() == FALSE){
  BiocManager::install("liftOver", update = FALSE)  
}
library(liftOver)
library(tidyverse)
library(dplyr)

assets.dir <- "/home/hieunguyen/src/ecd_wgs_enriched_and_kmer_features/assets/enrichedFeature_pipeline/atacseq_bed"
path.to.chain.file <- "/home/hieunguyen/src/ecd_wgs_enriched_and_kmer_features/resources/hg38ToHg19.over.chain"
dir.create(file.path(assets.dir, "hg19"), showWarnings = FALSE, recursive = TRUE)

all.bed.files <- Sys.glob(file.path(assets.dir, "hg38", "*.bed"))

for (file in all.bed.files){
  bedfile <- read.csv(file, sep = "\t", header = FALSE)
  colnames(bedfile) <- c("chr", "pos", "end")
  bedfile <- makeGRangesFromDataFrame(df = bedfile, seqnames.field = "chr", start.field = "pos", end.field = "end", keep.extra.columns = TRUE)
  chain <- import.chain(path.to.chain.file)
  bedfile.hg19 <- liftOver(bedfile, chain) %>% as.data.frame()  %>% subset(select = c(seqnames, start, end))
  
  write.table(bedfile.hg19, sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE, 
              file = file.path(assets.dir, "hg19", str_replace(basename(file), "[.]bed", ".hg19.bed")))
}

