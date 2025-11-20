gc()
rm(list = ls())

outdir <- "/home/hieunguyen/src/wgs_enriched_ecd_features/assets"
path.to.main.output <- file.path(outdir, "TCGA_bulkRNAseq")
path.to.00.output <- file.path(path.to.main.output, "00_output")
path.to.01.output <- file.path(path.to.main.output, "01_output")
path.to.02.output <- file.path(path.to.main.output, "02_output")
path.to.03.output <- file.path(path.to.main.output, "03_output")
dir.create(path.to.03.output, showWarnings = FALSE, recursive = TRUE)

topN <- 500
list.version <- "v0.1"

all.organs <- c("Liver",
                "Colorectal",
                "Lung",
                "Gastric",
                "Breast")

all.gene.lists <- list()

for (input.tissue.type in c("Tumor", "Normal")){
  down.genes <- list()
  up.genes <- list()
  
  cohort.down.genes <- readRDS(file.path(path.to.02.output, sprintf("%s_topN_%s_down_genes.rds", input.tissue.type, topN)))
  cohort.up.genes <- readRDS(file.path(path.to.02.output, sprintf("%s_topN_%s_up_genes.rds", input.tissue.type, topN)))
  
  down.genes$Liver <- cohort.down.genes$`TCGA-LIHC`
  down.genes$Colorectal <- cohort.down.genes$`TCGA-COAD`
  down.genes$Lung<- unique(c(cohort.down.genes$`TCGA-LUAD`, cohort.down.genes$`TCGA-LUSC`))
  down.genes$Gastric <- cohort.down.genes$`TCGA-STAD`
  down.genes$Breast <- cohort.down.genes$`TCGA-BRCA`
  
  up.genes$Liver <- cohort.up.genes$`TCGA-LIHC`
  up.genes$Colorectal <- cohort.up.genes$`TCGA-COAD`
  up.genes$Lung <- unique(c(cohort.up.genes$`TCGA-LUAD`, cohort.up.genes$`TCGA-LUSC`))
  up.genes$Gastric<- cohort.up.genes$`TCGA-STAD`
  up.genes$Breast <- cohort.up.genes$`TCGA-BRCA`
  
  all.gene.lists[[input.tissue.type]] <- list()
  
  for (input.cancer in names(down.genes)){
    all.gene.lists[[input.tissue.type]][[input.cancer]][[list.version]][["up"]] <- up.genes[[input.cancer]]
    all.gene.lists[[input.tissue.type]][[input.cancer]][[list.version]][["down"]] <- down.genes[[input.cancer]]
  }
}

up.genes <- list()
down.genes <- list()
for (input.organ in all.organs){
   up.genes[[input.organ]] <- unlist(setdiff(all.gene.lists$Tumor[[input.organ]][[list.version]]$up,  
                                                             all.gene.lists$Normal[[input.organ]][[list.version]]$up))
   down.genes[[input.organ]] <- unlist(setdiff(all.gene.lists$Tumor[[input.organ]][[list.version]]$down,  
                                                               all.gene.lists$Normal[[input.organ]][[list.version]]$down))
}

for (input.cancer in names(down.genes)){
  all.gene.lists[["TUMOR_minus_NORMAL"]][[input.cancer]][[list.version]][["up"]] <- up.genes[[input.cancer]]
  all.gene.lists[["TUMOR_minus_NORMAL"]][[input.cancer]][[list.version]][["down"]] <- down.genes[[input.cancer]]
}

saveRDS(all.gene.lists, file.path(path.to.03.output, "all_TOO_gene_lists.rds"))
