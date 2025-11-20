gc()
rm(list = ls())

library(dplyr)
library(tidyverse)
library(ggplot2)
library(vroom)
library(limma)
library(GenomicRanges)
library(comprehenr)


path.to.assets <- "/home/hieunguyen/src/ecd_wgs_enriched_and_kmer_features/assets"
path.to.save.output <- file.path(path.to.assets, "enrichedFeature_pipeline", "panel_7_8_vs_TSS")
dir.create(path.to.save.output, showWarnings = FALSE, recursive = TRUE)

all.beds <- Sys.glob(file.path(path.to.assets, "enrichedFeature_pipeline", "atacseq_bed", "hg19","*.hg19.bed"))
names(all.beds) <- to_vec(
  for (item in all.beds){
    str_replace(basename(item), ".bed", "")
  }
)

# read in TSS map from healthy samples
# we use the tss files generated from the k-mer-feature repo. 
tss.dir <- "/home/hieunguyen/src/ecd_wgs_enriched_and_kmer_features/assets/TSS_Feature_pipeline/tss_beds"
tss.source <- "biomart"

for (region.len in c(250, 500, 1000, 1500)){
  # up.flanking.size <- 1000
  # down.flanking.size <- 1000
  
  up.flanking.size <- region.len
  down.flanking.size <- region.len
  
  tss.input.file <- sprintf(file.path(tss.dir, 
                                      sprintf("%s/promoter_regions_up_%s_down_%s.bed", 
                                              tss.source, 
                                              up.flanking.size, 
                                              down.flanking.size)))
  tssdf <- read.csv(tss.input.file, sep = "\t", header =  FALSE)
  colnames(tssdf) <- c("chrom", "start", "end", "region")
  
  tss.grange <- makeGRangesFromDataFrame( df = tssdf, 
                                          seqnames.field = "chrom", 
                                          start.field = "start",
                                          end.field = "end", 
                                          keep.extra.columns = TRUE)
  
  for (bed.name in names(all.beds)){
    print(sprintf("Working on bed file %s", bed.name))
    atacseq.beddf <- read.csv(all.beds[[bed.name]], sep = "\t", header = TRUE)
    colnames(atacseq.beddf) <- c("chrom", "start", "end")
    atacseq.beddf$bed.name <- bed.name
    
    atacseq.grange <- makeGRangesFromDataFrame(df = atacseq.beddf, 
                                               seqnames.field = "chrom", 
                                               start.field = "start",
                                               end.field = "end", 
                                               keep.extra.columns = TRUE)
    
    overlap.res <- findOverlaps(atacseq.grange, tss.grange)
    
    firstdf <- atacseq.grange[queryHits(overlap.res)] %>% data.frame()
    colnames(firstdf) <- to_vec(
      for (item in colnames(firstdf)){
        sprintf("ATACseq_%s", item)
      }
    )
    
    seconddf <- tss.grange[subjectHits(overlap.res)] %>% data.frame()
    colnames(seconddf) <- to_vec(
      for (item in colnames(seconddf)){
        sprintf("NuMAP_%s", item)
      }
    )
    
    overlapdf <- cbind(firstdf, seconddf) %>%
      rowwise() %>%
      mutate(chrom = ATACseq_seqnames) %>%
      mutate(new.start = min(ATACseq_start, NuMAP_start)) %>%
      mutate(new.end = max(ATACseq_end, NuMAP_end)) %>%
      subset(select = c(chrom, new.start, new.end)) %>%
      mutate(region = sprintf("%s:%s-%s", chrom, new.start, new.end)) %>%
      mutate(region.len = new.end - new.start)
    
    overlapdf <- overlapdf[!duplicated(overlapdf$region),]
    colnames(overlapdf) <- c("chrom", "start", "end", "region", "region.len")
    
    write.table(overlapdf, file.path(path.to.save.output, sprintf("%s_overlap_TSS_up_%s_down_%s_%s.bed", bed.name, up.flanking.size, down.flanking.size, tss.source)), 
                row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)
  }
}

