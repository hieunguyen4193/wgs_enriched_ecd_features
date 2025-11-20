gc()
rm(list = ls())

path.to.main.src <- "/home/hieunguyen/src/wgs_enriched_ecd_features"
path.to.main.output <- file.path(path.to.main.src, "assets", "TSS_UTR5p_regions")
path.to.03.output <- file.path(path.to.main.output, "03_output")
dir.create(path.to.03.output, showWarnings = FALSE, recursive = TRUE)

source(file.path(path.to.main.src, "TSS_UTR5p_Feature_pipeline", "helper_functions.R"))

library(comprehenr)
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

# get ensembl mart data
ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")

if (file.exists(file.path(path.to.03.output, "TSS.hg38.rds")) == FALSE){
  # get tss data. 
  tss_data <- getBM(
    attributes = c("chromosome_name", "transcription_start_site", 
                   "strand", "external_gene_name", "ensembl_gene_id"),
    mart = ensembl
  ) %>% data.frame()
  saveRDS(tss_data, file.path(path.to.03.output, "TSS.hg38.rds"))
} else {
  print("reading in tss data for hg38 ...")
  tss_data <- readRDS(file.path(path.to.03.output, "TSS.hg38.rds"))
}

all.chroms <- to_vec(for (i in seq(1,22)) sprintf("%s", i))

tss_data <- subset(tss_data, tss_data$external_gene_name != "")
tss_data <- subset(tss_data, grepl("MT-", tss_data$external_gene_name) == FALSE)
tss_data <- subset(tss_data, tss_data$chromosome_name %in% all.chroms)
# input.gene <- head(tss_data$external_gene_name)
# input.gene <- sample(tss_data$external_gene_name, 10)

if (file.exists(file.path(path.to.03.output, "5prime_UTR.hg38.rds")) == FALSE){
  N <- 100
  split.batch.genes <- split(tss_data$external_gene_name, ceiling(seq_along(tss_data$external_gene_name) / N))
  utrdf <- data.frame()
  for (i in names(split.batch.genes)){
    print(sprintf("working on step %s", i))
    batch.gene <- split.batch.genes[[i]]
    query.bm <- getBM(attributes=c("chromosome_name", 
                                   "external_gene_name", 
                                   "5_utr_start",
                                   "5_utr_end"),
                      filters = 'external_gene_name', values = batch.gene, mart = ensembl, 
                      verbose = TRUE) %>%
      data.frame() %>%
      subset(is.na(`X5_utr_start`) == FALSE) %>%
      subset(is.na(`X5_utr_end`) == FALSE) 
    colnames(query.bm) <- c("chrom", "gene", "5_prime_utr_start", "5_prime_utr_end")
    query.bm$width <- query.bm$`5_prime_utr_end` - query.bm$`5_prime_utr_start`
    query.bm$mid.point <- 0.5*(query.bm$`5_prime_utr_end` - query.bm$`5_prime_utr_start`) + query.bm$`5_prime_utr_start`
    query.bm$start <- query.bm$mid.point - 1000
    query.bm$end <- query.bm$mid.point + 1000
    tmp.utrdf <- query.bm[, c("chrom", "start", "end", "gene", "5_prime_utr_start", "5_prime_utr_end", "width", "mid.point")]  
    utrdf <- rbind(utrdf, tmp.utrdf)
  }
  saveRDS(utrdf, file.path(path.to.03.output, "5prime_UTR.hg38.rds"))
} else {
  print("File 5prime_UTR.hg38.rds exsits.")
  print("reading in saved 5' UTR for hg38 ...")
  utrdf <- readRDS(file.path(path.to.03.output, "5prime_UTR.hg38.rds"))
  utrdf <- subset(utrdf, utrdf$chrom %in% all.chroms)
  utrdf <- utrdf %>% rowwise() %>%
    mutate(chrom = sprintf("chr%s", chrom))
}

re.save.bed.file <- FALSE
if (re.save.bed.file == TRUE){
  write.table(utrdf, file.path(path.to.03.output, "5prime_UTR.hg38.bed"), quote = FALSE,
              sep = "\t", row.names = FALSE, col.names = FALSE)
}

###### lift over to hg19
library(liftOver)
library(GenomicRanges)
path.to.chain.file <-  "/home/hieunguyen/src/ecd_wgs_enriched_and_kmer_features/resources/hg38ToHg19.over.chain"

utrdf.grange <- makeGRangesFromDataFrame(df = utrdf, 
                                         seqnames.field = "chrom",
                                         start.field = "start", 
                                         end.field = "end", 
                                         keep.extra.columns = TRUE)
chain <- import.chain(path.to.chain.file)
utrdf.hg19 <- liftOver(utrdf.grange, chain) %>% as.data.frame() 

utrdf.hg19 <- utrdf.hg19[, c("seqnames", 
                             "start",
                             "end", 
                             "width",
                             "strand",
                             "gene",
                             "X5_prime_utr_start",
                             "X5_prime_utr_end",
                             "mid.point")]
colnames(utrdf.hg19) <- c("chrom", 
                          "start",
                          "end", 
                          "width",
                          "strand",
                          "gene",
                          "5_prime_utr_start",
                          "5_prime_utr_end",
                          "mid.point")
write.table(utrdf.hg19, file.path(path.to.03.output, "5prime_UTR.hg19.bed"), quote = FALSE,
            sep = "\t", row.names = FALSE, col.names = FALSE)

##### overlap utrdf.hg19 with nucleosome map, 

if (file.exists(file.path(path.to.03.output, "5prime_UTR_overlap_nucleosomeMap.hg19.bed")) == FALSE){
  # read in nucleosome map
  path.to.nucleosome.map <- "/media/HNSD01/storage/resources/rpr_map_EXP0779.sorted.bed"
  nucleosome.map <- read.csv(path.to.nucleosome.map, header = FALSE, sep = "\t")
  colnames(nucleosome.map) <- c("chrom", "start", "end", "region", "V5", "strand", "mid.point.start", "mid.point.end")
  nucleosome.map.grange <- makeGRangesFromDataFrame(df = nucleosome.map, 
                                                    seqnames.field = "chrom", 
                                                    start.field = "start",
                                                    end.field = "end", 
                                                    keep.extra.columns = TRUE)
  
  # main overlap scripts
  utrdf.hg19.grange <- makeGRangesFromDataFrame(df = utrdf.hg19, 
                                                seqnames.field = "chrom", 
                                                start.field = "start",
                                                end.field = "end", 
                                                keep.extra.columns = TRUE)
  
  overlap.res <- findOverlaps(query = utrdf.hg19.grange, subject = nucleosome.map.grange)
  
  firstdf <- utrdf.hg19.grange[queryHits(overlap.res)] %>% data.frame()
  colnames(firstdf) <- to_vec(
    for (item in colnames(firstdf)){
      sprintf("utrdf.hg19_%s", item)
    }
  )
  
  seconddf <- nucleosome.map.grange[subjectHits(overlap.res)] %>% data.frame()
  colnames(seconddf) <- to_vec(
    for (item in colnames(seconddf)){
      sprintf("NuMAP_%s", item)
    }
  )
  
  overlapdf <- cbind(firstdf, seconddf) %>%
    rowwise() %>%
    mutate(chrom = utrdf.hg19_seqnames) %>%
    mutate(new.start = min(utrdf.hg19_start, NuMAP_start)) %>%
    mutate(new.end = max(utrdf.hg19_end, NuMAP_end)) %>%
    subset(select = c(chrom, new.start, new.end)) %>%
    mutate(region = sprintf("%s:%s-%s", chrom, new.start, new.end)) %>%
    mutate(region.len = new.end - new.start)
  
  overlapdf <- overlapdf[!duplicated(overlapdf$region),]
  colnames(overlapdf) <- c("chrom", "start", "end", "region", "region.len")
  write.table(overlapdf, file.path(path.to.03.output, "5prime_UTR_overlap_nucleosomeMap.hg19.bed"), quote = FALSE,
              sep = "\t", row.names = FALSE, col.names = FALSE)
  
} else {
  print(sprintf("File exists at %s", 
                file.path(path.to.03.output, "5prime_UTR_overlap_nucleosomeMap.hg19.bed")))
}

# EOF
