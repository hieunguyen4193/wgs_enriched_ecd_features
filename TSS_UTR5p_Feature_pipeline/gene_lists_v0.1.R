gene.list <- list(
  Lung = list(
    v0.1 = list(
      up = c("TOP2A", "TPX2", "ASPM", "CDK2", "NCAPG", "DLGAP5", "SLC2A1", "EGFR", "KRAS", "HDAC2", "VEGFA"),
      down = c("BTNL9", "CPED1", "FABP4", "AGER", "MGP", "PECAM1")
    )
    # ,
    # v0.2 = list(
    #   up = c("TOP2A", "CCNB1", "CCNA2", "CDK1", "UBE2C", "AURKA", "KIF20A", "ASPM", "CDC20", "MAD2L1", "EGFR", "KRAS", "MYC", "VEGFA", "TPX2", "DLGAP5", "SLC2A1", "HDAC2", "NCAPG", "MET", "CCND1", "MKI67", "MMP9"),
    #   down = c("TP53", "CDKN2A", "RB1", "PTEN", "FHIT", "RASSF1A", "SEMA3B", "BAP1", "BTNL9", "CPED1", "FABP4", "AGER", "MGP", "PECAM1", "STK11", "NKX2-1", "CDH1")
    # )
  ),
  Gastric = list(
    v0.1 = list(
      up = c("MYC", "CCND1", "PIK3CA", "VEGFA", "MMP2", "MMP9", "S100A4", "KRT7", "PI3", "E2F1", "TIMP1", "S100A9", "SERPINE1", "CCL3", "CCL20", "CXCL8", "IL1B", "INHBA", "COL1A1", "COL1A2", "CDH17", "SPP1", "REG4"),
      down = c("CDH1", "TP53", "ARID1A", "CLDN18", "ATP4A", "ATP4B", "GKN1", "GKN2", "LIPF", "GATA6", "SMAD4", "APC", "CDKN1A", "CDKN2A", "TFF1")
    )
  ),
  Colorectal = list(
    v0.1 = list(
      up = c("MYC", "CCND1", "KRAS", "PIK3CA", "TOP2A", "CCNB1", "CDC20", "AURKA", "UBE2C", "VEGFA", "MMP7", "MMP2", "CDK4", "EGFR", "LGR5", "TCF7L2", "AXIN2", "ASCL2", "EPHB2"),
      down = c("APC", "TP53", "SMAD4", "FBXW7", "PIK3R1", "PTEN", "DCC", "MLH1", "MSH2", "CDKN2A", "CDKN1A", "CDH1", "TGFBR2", "BAX")
    )
  )
)

##### run once after changing the gene list above only. 

# library(tidyverse)
# library(dplyr)
# 
# path.to.main.src <- "/home/hieunguyen/src/k-mer-feature"
# path.to.save.output <- file.path(path.to.main.src, "genelist")
# dir.create(path.to.save.output, showWarnings = FALSE, recursive = TRUE)
# 
# for (cancer.type in names(gene.list)){
#   for (list.version in names(gene.list[[cancer.type]])){
#     for (up.or.down in c("up", "down")){
#       list.of.genes <- gene.list[[cancer.type]][[list.version]][[up.or.down]]
#       tmpdf <- data.frame(Gene = list.of.genes)
#       write.table(tmpdf, file.path(path.to.save.output, sprintf("%s_%s_%s.csv", 
#                                                               cancer.type, list.version, up.or.down)),
#                   sep = ",", quote = FALSE, row.names = FALSE)
#     }
#   }
# }
