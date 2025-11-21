generate_tss_regions_for_input_genes <- function(tssdf, input.genes, outputdir, up.flank = 1000, down.flank = 1000){
  dir.create(outputdir, showWarnings = FALSE, recursive = TRUE)
  tssdf.filter <- subset(tssdf, tssdf$gene %in% input.genes)
  promoterdf <- define_promoter_regions(tssdf = tssdf.filter, 
                                        up.flank = up.flank, 
                                        down.flank = down.flank, 
                                        outputdir = outputdir)
  promoterdf <- merge(promoterdf, tssdf.filter, by.x = "promoter_tssName", by.y = "tssName")
  write.table(promoterdf, file.path(outputdir, "full_info.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
}

generate_utr5p_regions_for_input_genes <- function(utr5pdf, input.genes, outputdir){
  utr5pdf.filter <- subset(utr5pdf, utr5pdf$gene %in% input.genes)
  write.table(utr5pdf.filter, file.path(outputdir, "5prime_UTR_up_genes.hg19.bed"), quote = FALSE,
              sep = "\t", row.names = FALSE, col.names = FALSE)
}