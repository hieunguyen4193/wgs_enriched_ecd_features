gc()
rm(list = ls())

path.to.main.src <- "/home/hieunguyen/src/wgs_enriched_ecd_features"

path.to.main.output <- file.path(path.to.main.src, "assets", "TSS_UTR5p_regions")
path.to.02.output <- file.path(path.to.main.output, "02_output")
dir.create(path.to.02.output, showWarnings = FALSE, recursive = TRUE)

source(file.path(path.to.main.src, "TSS_UTR5p_Feature_pipeline", "helper_functions.R"))

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

library(liftOver)
library(GenomicRanges)
path.to.chain.file <-  file.path(path.to.main.src, "resources/hg38ToHg19.over.chain")

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

# check all datasets in ensembl.
# all.ensembl.dataset <- listDatasets(ensembl) %>% data.frame()
# subset(all.ensembl.dataset, all.ensembl.dataset$dataset == "hsapiens_gene_ensembl") # only hg38 exists, need liftover

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
write.csv(tssdf.hg19, file.path(path.to.02.output, "tssdf.hg19.csv"))

customRegiondf <- define_customRegion_regions(tssdf = tssdf.hg19, 
                                              from.flank = 250, 
                                              to.flank = 500, 
                                              direction = "upstream", 
                                              outputdir = file.path(path.to.02.output, "biomart"),
                                              rerun = TRUE)
customRegiondf <- define_customRegion_regions(tssdf = tssdf.hg19, 
                                              from.flank = 500, 
                                              to.flank = 1500, 
                                              direction = "upstream", 
                                              outputdir = file.path(path.to.02.output, "biomart"),
                                              rerun = TRUE)
customRegiondf <- define_customRegion_regions(tssdf = tssdf.hg19, 
                                              from.flank = 500, 
                                              to.flank = 1000, 
                                              direction = "upstream", 
                                              outputdir = file.path(path.to.02.output, "biomart"),
                                              rerun = TRUE)

customRegiondf <- define_customRegion_regions(tssdf = tssdf.hg19, 
                                              from.flank = 250, 
                                              to.flank = 500, 
                                              direction = "downstream", 
                                              outputdir = file.path(path.to.02.output, "biomart"),
                                              rerun = TRUE)
customRegiondf <- define_customRegion_regions(tssdf = tssdf.hg19, 
                                              from.flank = 500, 
                                              to.flank = 1500, 
                                              direction = "downstream", 
                                              outputdir = file.path(path.to.02.output, "biomart"),
                                              rerun = TRUE)
customRegiondf <- define_customRegion_regions(tssdf = tssdf.hg19, 
                                              from.flank = 500, 
                                              to.flank = 1000, 
                                              direction = "downstream", 
                                              outputdir = file.path(path.to.02.output, "biomart"),
                                              rerun = TRUE)