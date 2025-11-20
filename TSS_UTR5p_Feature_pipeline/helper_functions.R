#####
##### DEFINE NORMLA PROMOTER REGIONS
#####
define_promoter_regions <- function(tssdf, up.flank, down.flank, outputdir, rerun = TRUE){
  dir.create(outputdir, showWarnings = FALSE, recursive = TRUE)
  if (file.exists(file.path(outputdir, sprintf("promoter_regions_up_%s_down_%s.csv", up.flank, down.flank))) == FALSE | rerun == TRUE){
    library(GenomicRanges)
    library(dplyr)
    library(tidyverse)
    library(comprehenr)
    
    path.to.main.src.tmp <- "/media/HNSD01/storage/resources/tss_beds"
    path.to.TSS.file <- file.path(path.to.main.src.tmp, "UCSC_TSS.bed")
    path.to.refseq <- file.path(path.to.main.src.tmp, "hg19_refseq_all.bed")
    
    path.to.TSS.annotated.file <- file.path(path.to.main.src.tmp, "UCSC_TSS.annotated.bed")
    
    hg19.refseq <- read.csv(file.path(path.to.refseq), sep = "\t")[,c("chrom", "txStart", "txEnd", "name", "name2", "strand")]
    all.valid.chroms <- to_vec( for(item in seq(1,22)) sprintf("chr%s", item))
    hg19.refseq <- subset(hg19.refseq, hg19.refseq$chrom %in% all.valid.chroms)
    colnames(hg19.refseq) <- c("chrom", "start", "end", "tx_name", "gene", "strand")
    hg19.refseq.grange <- makeGRangesFromDataFrame(df = hg19.refseq, start.field = "start", end.field = "end", seqnames.field = "chrom", keep.extra.columns = TRUE, strand.field = "strand")
    
    tss.grange <- makeGRangesFromDataFrame(df = tssdf, 
                                           start.field = "start", 
                                           end.field = 'end', 
                                           seqnames.field = "chrom", 
                                           keep.extra.columns = TRUE)
    
    
    if (file.exists(path.to.TSS.annotated.file) == FALSE){
      add_suffix1 <- "TSS"
      add_suffix2 <- "hg19"
      
      overlap.idxs <-  findOverlaps(tss.grange, hg19.refseq.grange)
      tmpdf1 <- tss.grange[queryHits(overlap.idxs), ] %>% data.frame() %>%
        rowwise() %>%
        mutate(name = sprintf("%s.%s.%s", seqnames, start, end))
      colnames(tmpdf1) <- to_vec ( for (item in colnames(tmpdf1)) sprintf("%s_%s", add_suffix1, item))
      tmpdf2 <- hg19.refseq.grange[subjectHits(overlap.idxs), ] %>% data.frame() %>%
        rowwise() %>%
        mutate(promoter.name = sprintf("%s.%s.%s", seqnames, start, end))
      colnames(tmpdf2) <- to_vec ( for (item in colnames(tmpdf2)) sprintf("%s_%s", add_suffix2, item))
      
      tss.annotated.df <- cbind(tmpdf1, tmpdf2)
      write.csv(tss.annotated.df, path.to.TSS.annotated.file)
    }
    
    promoterdf <- tssdf %>% rowwise() %>% 
      mutate(promoter.start = ifelse(strand == "+", start - up.flank, end - down.flank)) %>%
      mutate(promoter.end = ifelse(strand == "+", start + down.flank, end + up.flank))
    
    promoter.grange <- makeGRangesFromDataFrame(df = subset(promoterdf, select = c(chrom, promoter.start, promoter.end, strand, tssName)), 
                                                seqnames.field = "chrom",
                                                start.field = "promoter.start",
                                                end.field = "promoter.end",
                                                strand.field = "strand",
                                                keep.extra.columns = TRUE)
    
    add_suffix1 <- "promoter"
    add_suffix2 <- "hg19"
    
    overlap.idxs <-  findOverlaps(promoter.grange, hg19.refseq.grange)
    
    tmpdf1 <- promoter.grange[queryHits(overlap.idxs), ] %>% data.frame() %>%
      rowwise() %>%
      mutate(name = sprintf("%s.%s.%s", seqnames, start, end))
    colnames(tmpdf1) <- to_vec ( for (item in colnames(tmpdf1)) sprintf("%s_%s", add_suffix1, item))
    
    tmpdf2 <- hg19.refseq.grange[subjectHits(overlap.idxs), ] %>% data.frame() 
    colnames(tmpdf2) <- to_vec ( for (item in colnames(tmpdf2)) sprintf("%s_%s", add_suffix2, item))
    
    promoterdf <- cbind(tmpdf1, tmpdf2)
    promoterdf <- promoterdf[!duplicated(promoterdf$promoter_name), ]
    
    promoterdf <- merge(promoterdf, tssdf, by.x = "promoter_tssName", by.y = "tssName")
    promoterdf <- promoterdf %>%
      rowwise() %>%
      mutate(promoter_fullname = sprintf("%s_%s", promoter_tssName, gene))
    
    write.csv(promoterdf, file.path(outputdir, sprintf("promoter_regions_up_%s_down_%s.csv", up.flank, down.flank)))
  } else {
    promoterdf <- read.csv(file.path(outputdir, sprintf("promoter_regions_up_%s_down_%s.csv", up.flank, down.flank))) %>%
      subset(select = -c(X))
  }
  write.table(promoterdf[, c("promoter_seqnames", "promoter_start", "promoter_end", "promoter_fullname")],
              file.path(outputdir, sprintf("promoter_regions_up_%s_down_%s.bed", up.flank, down.flank)),
              row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)
  return(promoterdf)
}

#####
##### DEFINE CUSTOM PROMOTER REGIONS
#####
define_customRegion_regions <- function(tssdf, from.flank, to.flank, direction, outputdir, rerun = TRUE){
  dir.create(outputdir, showWarnings = FALSE, recursive = TRUE)
  if (file.exists(file.path(outputdir, sprintf("customRegion_regions_from_%s_to_%s_%s.csv", from.flank, to.flank, direction))) == FALSE | rerun == TRUE){
    library(GenomicRanges)
    library(dplyr)
    library(tidyverse)
    library(comprehenr)
    
    path.to.main.src.tmp <- "/media/HNSD01/storage/resources/tss_beds"
    path.to.TSS.file <- file.path(path.to.main.src.tmp, "UCSC_TSS.bed")
    path.to.refseq <- file.path(path.to.main.src.tmp, "hg19_refseq_all.bed")
    
    path.to.TSS.annotated.file <- file.path(path.to.main.src.tmp, "UCSC_TSS.annotated.bed")
    
    hg19.refseq <- read.csv(file.path(path.to.refseq), sep = "\t")[,c("chrom", "txStart", "txEnd", "name", "name2", "strand")]
    all.valid.chroms <- to_vec( for(item in seq(1,22)) sprintf("chr%s", item))
    hg19.refseq <- subset(hg19.refseq, hg19.refseq$chrom %in% all.valid.chroms)
    colnames(hg19.refseq) <- c("chrom", "start", "end", "tx_name", "gene", "strand")
    hg19.refseq.grange <- makeGRangesFromDataFrame(df = hg19.refseq, start.field = "start", end.field = "end", seqnames.field = "chrom", keep.extra.columns = TRUE, strand.field = "strand")
    
    tss.grange <- makeGRangesFromDataFrame(df = tssdf, 
                                           start.field = "start", 
                                           end.field = 'end', 
                                           seqnames.field = "chrom", 
                                           keep.extra.columns = TRUE)
    if (file.exists(path.to.TSS.annotated.file) == FALSE){
      add_suffix1 <- "TSS"
      add_suffix2 <- "hg19"
      
      overlap.idxs <-  findOverlaps(tss.grange, hg19.refseq.grange)
      tmpdf1 <- tss.grange[queryHits(overlap.idxs), ] %>% data.frame() %>%
        rowwise() %>%
        mutate(name = sprintf("%s.%s.%s", seqnames, start, end))
      colnames(tmpdf1) <- to_vec ( for (item in colnames(tmpdf1)) sprintf("%s_%s", add_suffix1, item))
      tmpdf2 <- hg19.refseq.grange[subjectHits(overlap.idxs), ] %>% data.frame() %>%
        rowwise() %>%
        mutate(customRegion.name = sprintf("%s.%s.%s", seqnames, start, end))
      colnames(tmpdf2) <- to_vec ( for (item in colnames(tmpdf2)) sprintf("%s_%s", add_suffix2, item))
      
      tss.annotated.df <- cbind(tmpdf1, tmpdf2)
      write.csv(tss.annotated.df, path.to.TSS.annotated.file)
    }    
    
    if (direction == "upstream"){
      customRegiondf <- tssdf %>% rowwise() %>% 
        mutate(customRegion.start = ifelse(strand == "+", start - to.flank, end - to.flank)) %>%
        mutate(customRegion.end = ifelse(strand == "+", start - from.flank, end - from.flank)) %>%
        subset(customRegion.start > 0)      
    } else if (direction == "downstream"){
      customRegiondf <- tssdf %>% rowwise() %>% 
        mutate(customRegion.start = ifelse(strand == "+", start + from.flank, start + from.flank)) %>%
        mutate(customRegion.end = ifelse(strand == "+", start + to.flank, start + to.flank)) %>%
        subset(customRegion.start > 0)
    }
    
    customRegion.grange <- makeGRangesFromDataFrame(df = subset(customRegiondf, select = c(chrom, customRegion.start, customRegion.end, strand, tssName)), 
                                                    seqnames.field = "chrom",
                                                    start.field = "customRegion.start",
                                                    end.field = "customRegion.end",
                                                    strand.field = "strand",
                                                    keep.extra.columns = TRUE)
    
    add_suffix1 <- "customRegion"
    add_suffix2 <- "hg19"
    
    overlap.idxs <-  findOverlaps(customRegion.grange, hg19.refseq.grange)
    
    tmpdf1 <- customRegion.grange[queryHits(overlap.idxs), ] %>% data.frame() %>%
      rowwise() %>%
      mutate(name = sprintf("%s.%s.%s", seqnames, start, end))
    colnames(tmpdf1) <- to_vec ( for (item in colnames(tmpdf1)) sprintf("%s_%s", add_suffix1, item))
    
    tmpdf2 <- hg19.refseq.grange[subjectHits(overlap.idxs), ] %>% data.frame() 
    colnames(tmpdf2) <- to_vec ( for (item in colnames(tmpdf2)) sprintf("%s_%s", add_suffix2, item))
    
    customRegiondf <- cbind(tmpdf1, tmpdf2)
    write.csv(customRegiondf, file.path(outputdir, sprintf("customRegion_regions_from_%s_to_%s_%s.csv", from.flank, to.flank, direction)))
  } else {
    customRegiondf <- read.csv(file.path(outputdir, sprintf("customRegion_regions_from_%s_to_%s_%s.csv", from.flank, to.flank, direction))) 
    if ("X" %in% colnames(customRegiondf)){
      customRegiondf <- customRegiondf %>% subset(select = -c(X))
    }
  }
  write.table(customRegiondf[, c("customRegion_seqnames", "customRegion_start", "customRegion_end", "customRegion_tssName")],
              file.path(outputdir, sprintf("customRegion_regions_from_%s_to_%s_%s.bed", from.flank, to.flank, direction)),
              row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)
  return(customRegiondf)
}