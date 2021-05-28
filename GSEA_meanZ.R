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


# load meanZ ----
for (FN in list.files("~/Dropbox/GDB_archive/CMapCorr_files","lig295_DE_[a-z]+_FDR\\.RData",
                      full.names=T)) {
  load(FN)
}


# filter pathways ----
pathways_REACT <- sapply(pathways_REACT,function(X) 
  X[X %in% rownames(meanZ_ct)],simplify=F)
pathwayNames_REACT <- pathwayNames_REACT[sapply(pathways_REACT,length) >= 10 & 
                                           sapply(pathways_REACT,length) <= 200]
pathways_REACT <- pathways_REACT[sapply(pathways_REACT,length) >= 10 & 
                                   sapply(pathways_REACT,length) <= 200]

pathways_GOBP <- sapply(pathways_GOBP,function(X) 
  X[X %in% rownames(meanZ_ct)],simplify=F)
pathwayNames_GOBP <- pathwayNames_GOBP[sapply(pathways_GOBP,length) >= 10 & 
                                         sapply(pathways_GOBP,length) <= 200]
pathways_GOBP <- pathways_GOBP[sapply(pathways_GOBP,length) >= 10 & 
                                 sapply(pathways_GOBP,length) <= 200]

pathways_ALL <- sapply(pathways_ALL,function(X) 
  X[X %in% rownames(meanZ_ct)],simplify=F)
pathwayNames_ALL <- pathwayNames_ALL[sapply(pathways_ALL,length) >= 10 & 
                                       sapply(pathways_ALL,length) <= 200]
pathways_ALL <- pathways_ALL[sapply(pathways_ALL,length) >= 10 & 
                               sapply(pathways_ALL,length) <= 200]


# ^ do GSEA ----
for (PATH in grep("^pathways",ls(),value=T)) {
  GSEA <- strsplit(PATH,"_")[[1]][2]
  print(paste(GSEA,"----------"))
  GSEA_MZ <- list()
  for (MZ in grep("^meanZ",ls(),value=T)) {
    L <- strsplit(MZ,"_")[[1]][2]
    print(paste("^",L,"-----"))
    GSEA_MZ[[L]] <- pbapply(get(MZ),2,function(X) 
      as.data.frame(fgseaSimple(pathways=get(PATH),
                                stats=X,
                                nperm=1e4)[,c("NES","padj")]),
      cl=4)
  }
  save(GSEA_MZ,file=paste0("~/Dropbox/GDB_archive/CMapCorr_files/GSEA_meanZ_",GSEA,".RData"))
}
