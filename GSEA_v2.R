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


# Lvl4 ----
# ^ load data ----
load("~/Dropbox/GDB_archive/CMapCorr_files/lvl4_inputs_allgenes.RData")
lvl4_data <- lvl4_data_all
rm(lvl4_data_all)
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
rm(list=c("lig16","ct14",grep("^temp",ls(),value=T)))


# ^ filter pathways ----
pathways_REACT <- sapply(pathways_REACT,function(X) 
  X[X %in% lvl4_data@rid],simplify=F)
pathwayNames_REACT <- pathwayNames_REACT[sapply(pathways_REACT,length) >= 10 & 
                                           sapply(pathways_REACT,length) <= 200]
pathways_REACT <- pathways_REACT[sapply(pathways_REACT,length) >= 10 & 
                                   sapply(pathways_REACT,length) <= 200]

pathways_GOBP <- sapply(pathways_GOBP,function(X) 
  X[X %in% lvl4_data@rid],simplify=F)
pathwayNames_GOBP <- pathwayNames_GOBP[sapply(pathways_GOBP,length) >= 10 & 
                                         sapply(pathways_GOBP,length) <= 200]
pathways_GOBP <- pathways_GOBP[sapply(pathways_GOBP,length) >= 10 & 
                                 sapply(pathways_GOBP,length) <= 200]

pathways_ALL <- sapply(pathways_ALL,function(X) 
  X[X %in% lvl4_data@rid],simplify=F)
pathwayNames_ALL <- pathwayNames_ALL[sapply(pathways_ALL,length) >= 10 & 
                                       sapply(pathways_ALL,length) <= 200]
pathways_ALL <- pathways_ALL[sapply(pathways_ALL,length) >= 10 & 
                               sapply(pathways_ALL,length) <= 200]


# ^ do GSEA ----
gsea_nes_REACT <- pbapply(lvl4_data@mat,2,function(X) 
  fgseaSimple(pathways=pathways_REACT,stats=X,nperm=1e4)[,c("NES","padj")],
  cl=4)
gsea_nes_REACT <- append(gsea_nes_REACT,list(PATHWAYS=pathwayNames_REACT),0)
save(gsea_nes_REACT,file="~/Dropbox/GDB_archive/CMapCorr_files/GSEA_lvl4_REACT.RData")

gsea_nes_GOBP <- pbapply(lvl4_data@mat,2,function(X) 
  fgseaSimple(pathways=pathways_GOBP,stats=X,nperm=1e4)[,c("NES","padj")],
  cl=4)
gsea_nes_GOBP <- append(gsea_nes_GOBP,list(PATHWAYS=pathwayNames_GOBP),0)
save(gsea_nes_GOBP,file="~/Dropbox/GDB_archive/CMapCorr_files/GSEA_lvl4_GOBP.RData")

gsea_nes_ALL <- pbapply(lvl4_data@mat,2,function(X) 
  fgseaSimple(pathways=pathways_ALL,stats=X,nperm=1e4)[,c("NES","padj")],
  cl=4)
gsea_nes_ALL <- append(gsea_nes_ALL,list(PATHWAYS=pathwayNames_ALL),0)
save(gsea_nes_ALL,file="~/Dropbox/GDB_archive/CMapCorr_files/GSEA_lvl4_ALL.RData")
