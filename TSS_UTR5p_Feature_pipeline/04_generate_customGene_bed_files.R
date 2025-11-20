gc()
rm(list = ls())

path.to.main.src <- "/home/hieunguyen/src/wgs_enriched_ecd_features"

path.to.main.output <- file.path(path.to.main.src, "assets", "TSS_UTR5p_regions")
path.to.01.output <- file.path(path.to.main.output, "01_output")
path.to.02.output <- file.path(path.to.main.output, "02_output")
path.to.03.output <- file.path(path.to.main.output, "03_output")
path.to.04.output <- file.path(path.to.main.output, "04_output")

source(file.path(path.to.main.src, "TSS_UTR5p_Feature_pipeline" ,"TSS_UTR5p_helper_functions.R"))
source(file.path(path.to.main.src, "TSS_UTR5p_Feature_pipeline", "helper_functions.R"))

dir.create(file.path(path.to.04.output, "TSS_beds"))
dir.create(file.path(path.to.04.output, "UTR5_beds"))

# ***** manual curated gene lists *****
input.list.version <- "Tumor"
if (input.list.version == "v0.1"){
  source(file.path(path.to.main.src, "TSS_UTR5p_Feature_pipeline", "gene_lists_v0.1.R"))  
} else if (input.list.version == "v0.2"){
  gene.list <- readRDS(file.path(path.to.main.src, "TSS_UTR5p_Feature_pipeline", "gene_lists_v0.2.rds"))
} else if (grepl("RANDOM", input.list.version) == TRUE){
  gene.list <- readRDS(file.path(path.to.main.src, "TSS_UTR5p_Feature_pipeline", sprintf("RANDOM_gene_lists_v%s.rds", i)))
} else if (grepl("Tumor", input.list.version) == TRUE){
  path.to.TOO.TCGA.output <- "/home/hieunguyen/src/wgs_enriched_ecd_features/assets/TCGA_bulkRNAseq/03_output/all_TOO_gene_lists.rds"
  gene.list <- readRDS(path.to.TOO.TCGA.output)[["Tumor"]]
} else if (grepl("Normal", input.list.version) == TRUE){
  path.to.TOO.TCGA.output <- "/home/hieunguyen/src/wgs_enriched_ecd_features/assets/TCGA_bulkRNAseq/03_output/all_TOO_gene_lists.rds"
  gene.list <- readRDS(path.to.TOO.TCGA.output)[["Normal"]]
}