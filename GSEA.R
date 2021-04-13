library(pbapply)
library(fgsea)
.PAR <- par(no.readonly=T)
pboptions(type="timer")

# Create pathways lists ----
# http://www.baderlab.org/GeneSets
# ^ ALL ----
temp_file <- "~/Dropbox/GDB_archive/Human_GOBP_AllPathways_no_GO_iea_March_01_2021_entrezgene.gmt" 
temp_lines <- strsplit(readLines(temp_file), "\t")
pathways_ALL <- lapply(temp_lines,tail,-2)
pathwayNames_ALL <- sapply(temp_lines,function(X) X[2])
names(pathwayNames_ALL) <- names(pathways_ALL) <- sapply(temp_lines,head,1)

# ^ GO biological process ----
temp_file <- "~/Dropbox/GDB_archive/Human_GO_bp_no_GO_iea_entrezgene.gmt" # http://www.baderlab.org/GeneSets
temp_lines <- strsplit(readLines(temp_file), "\t")
pathways_GOBP <- lapply(temp_lines,tail,-2)
pathwayNames_GOBP <- sapply(temp_lines,function(X) X[2])
names(pathwayNames_GOBP) <- names(pathways_GOBP) <- sapply(temp_lines,head,1)

# ^ Reactome ----
temp_file <- "~/Dropbox/GDB_archive/Human_Reactome_March_01_2021_Entrezgene.gmt" # http://www.baderlab.org/GeneSets
temp_lines <- strsplit(readLines(temp_file), "\t")
pathways_REACT <- lapply(temp_lines,tail,-2)
pathwayNames_REACT <- sapply(temp_lines,function(X) X[2])
names(pathwayNames_REACT) <- names(pathways_REACT) <- sapply(temp_lines,head,1)

rm(list=grep("^temp",ls(),value=T))


# Lvl5 ----
if (!file.exists("~/Dropbox/GDB_archive/CMapCorr_files/GSEA_lvl5.RData")) {
  load("~/Dropbox/GDB_archive/CMapCorr_files/lvl5_inputs.RData")
  
  gsea_nes_GOBP <- pbapply(lvl5_data@mat,2,function(X) 
    fgseaSimple(pathways=pathways_GOBP,
                stats=X,
                minSize=10,
                maxSize=200,
                nperm=1e3)$NES)
  rownames(gsea_nes_GOBP) <- fgseaSimple(pathways=pathways_GOBP,
                                         stats=lvl5_data@mat[,1],
                                         minSize=10,
                                         maxSize=200,
                                         nperm=1e3)$pathway
  
  gsea_nes_REACT <- pbapply(lvl5_data@mat,2,function(X) 
    fgseaSimple(pathways=pathways_REACT,
                stats=X,
                minSize=10,
                maxSize=200,
                nperm=1e3)$NES)
  rownames(gsea_nes_REACT) <- fgseaSimple(pathways=pathways_REACT,
                                          stats=lvl5_data@mat[,1],
                                          minSize=10,
                                          maxSize=200,
                                          nperm=1e3)$pathway
  
  gsea_nes_ALL <- pbapply(lvl5_data@mat,2,function(X) 
    fgseaSimple(pathways=pathways_ALL,
                stats=X,
                minSize=10,
                maxSize=200,
                nperm=1e3)$NES)
  rownames(gsea_nes_ALL) <- fgseaSimple(pathways=pathways_ALL,
                                        stats=lvl5_data@mat[,1],
                                        minSize=10,
                                        maxSize=200,
                                        nperm=1e3)$pathway
  
  save(list=grep("^gsea",ls(),value=T),file="~/Dropbox/GDB_archive/CMapCorr_files/GSEA_lvl5.RData")
}



# Lvl4 ----
if (!file.exists("~/Dropbox/GDB_archive/CMapCorr_files/GSEA_lvl4.RData")) {
  load("~/Dropbox/GDB_archive/CMapCorr_files/lvl4_inputs.RData")
  temp_ligcountsperct <- sapply(unique(lvl4_data@cdesc$cell_id),function(CT)
    length(unique(lvl4_data@cdesc[lvl4_data@cdesc$cell_id == CT,"pert_iname"])))
  ct9 <- names(temp_ligcountsperct)[temp_ligcountsperct > 100]
  names(ct9) <- sapply(ct9,function(X) names(ct14)[ct14 == X])
  
  temp_ctcountsperlig <- sapply(unique(lvl4_data@cdesc$pert_iname),function(LIG)
    length(unique(lvl4_data@cdesc[lvl4_data@cdesc$pert_iname == LIG &
                                    lvl4_data@cdesc$cell_id %in% ct9,"cell_id"])))
  lig295 <- names(temp_ctcountsperlig)[temp_ctcountsperlig == 9]
  lig295 <- sort(lig295)
  
  IDlig295 <- rownames(lvl4_data@cdesc)[lvl4_data@cdesc$pert_iname %in% lig295 & 
                                          lvl4_data@cdesc$cell_id %in% ct9]
  lvl4_data@mat <- lvl4_data@mat[,IDlig295]
  lvl4_data@cdesc <- lvl4_data@cdesc[IDlig295,]
  
  gsea_nes_GOBP <- pbapply(lvl4_data@mat,2,function(X) 
    fgseaSimple(pathways=pathways_GOBP,
                stats=X,
                minSize=10,
                maxSize=200,
                nperm=1e3)$NES)
  rownames(gsea_nes_GOBP) <- fgseaSimple(pathways=pathways_GOBP,
                                         stats=lvl4_data@mat[,1],
                                         minSize=10,
                                         maxSize=200,
                                         nperm=1e3)$pathway
  
  gsea_nes_REACT <- pbapply(lvl4_data@mat,2,function(X) 
    fgseaSimple(pathways=pathways_REACT,
                stats=X,
                minSize=10,
                maxSize=200,
                nperm=1e3)$NES)
  rownames(gsea_nes_REACT) <- fgseaSimple(pathways=pathways_REACT,
                                          stats=lvl4_data@mat[,1],
                                          minSize=10,
                                          maxSize=200,
                                          nperm=1e3)$pathway
  
  gsea_nes_ALL <- pbapply(lvl4_data@mat,2,function(X) 
    fgseaSimple(pathways=pathways_ALL,
                stats=X,
                minSize=10,
                maxSize=200,
                nperm=1e3)$NES)
  rownames(gsea_nes_ALL) <- fgseaSimple(pathways=pathways_ALL,
                                        stats=lvl4_data@mat[,1],
                                        minSize=10,
                                        maxSize=200,
                                        nperm=1e3)$pathway
  
  save(list=grep("^gsea",ls(),value=T),file="~/Dropbox/GDB_archive/CMapCorr_files/GSEA_lvl4.RData")
}