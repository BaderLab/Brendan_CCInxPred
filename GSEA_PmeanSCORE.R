library(cmapR)
library(pbapply)
pboptions(type="timer")

# Load lvl4_data ----
load("~/Dropbox/GDB_archive/CMapCorr_files/lvl4_inputs.RData")
rm(lvl4_data_ctl)

temp_ligcountsperct <- sapply(unique(lvl4_data@cdesc$cell_id),function(CT)
  length(unique(lvl4_data@cdesc[lvl4_data@cdesc$cell_id == CT,"pert_iname"])))
ct9 <- names(temp_ligcountsperct)[temp_ligcountsperct > 100]
names(ct9) <- sapply(ct9,function(X) names(ct14)[ct14 == X])

temp_ctcountsperlig <- sapply(unique(lvl4_data@cdesc$pert_iname),function(LIG)
  length(unique(lvl4_data@cdesc[lvl4_data@cdesc$pert_iname == LIG &
                                  lvl4_data@cdesc$cell_id %in% ct9,"cell_id"])))
lig295 <- names(temp_ctcountsperlig)[temp_ctcountsperlig == 9]
lig295 <- sort(lig295)
temp_id <- lvl4_data@cid[lvl4_data@cdesc$cell_id %in% ct9 &
                           lvl4_data@cdesc$pert_iname %in% lig295]
lvl4_data@mat <- lvl4_data@mat[,temp_id]
CDESC <- lvl4_data@cdesc[temp_id,]
CID <- temp_id

rm(list=c("lvl4_data","ct14","lig16",grep("^temp",ls(),value=T)))


# Load GSEA results ----

minP <- 1e4

for (FN in list.files("~/Dropbox/GDB_archive/CMapCorr_files/","^GSEA_lvl4_[A-Z]")) {
  temp <- load(paste0("~/Dropbox/GDB_archive/CMapCorr_files/",FN))
  GSEA <- strsplit(temp,"_")[[1]][3]
  if (file.exists(paste0("~/Dropbox/GDB_archive/CMapCorr_files/pfdrGSEA_",GSEA,".RData"))) { 
    rm(list=temp)
    next 
  }

  print(paste(GSEA,"----"))
  
  temp_NES <- sapply(get(temp)[-1],"[[","NES")
  rownames(temp_NES) <- as.vector(get(temp)$PATHWAYS)
  
  temp_padj <- sapply(get(temp)[-1],"[[","padj")
  rownames(temp_padj) <- as.vector(get(temp)$PATHWAYS)
  
  GSEA_score <- -log10(temp_padj) * sign(temp_NES)
  rownames(GSEA_score) <- as.vector(get(temp)$PATHWAYS)
  
  rm(list=c(temp,grep("^temp",ls(),value=T)))


  # Mean signed -log10(FDR) ----
  # ^ CT ----
  print("^ CT ----")
  if (file.exists(paste0("~/Dropbox/GDB_archive/CMapCorr_files/pfdrGSEA_",GSEA,"_BKGD_CT.RData"))) {
    load(paste0("~/Dropbox/GDB_archive/CMapCorr_files/pfdrGSEA_",GSEA,"_BKGD_CT.RData"))
  } else {
    print("^^ Generating null distribution ----")
    temp_counts_uniq <- sort(unique(
      sapply(ct9,function(CT) sum(CDESC$cell_id == CT)),
    ),decreasing=T)
    BKGD_CT <- pbsapply(temp_counts_uniq,function(N) {
      replicate(minP,{
        rowMeans(GSEA_score[,sample(colnames(GSEA_score),N)])
      })
    },simplify=F)
    names(BKGD_CT) <- temp_counts_uniq
    save(BKGD_CT,
         file=paste0("~/Dropbox/GDB_archive/CMapCorr_files/pfdrGSEA_",GSEA,"_BKGD_CT.RData"))
  }
  
  print("^^ Calculating ----")
  temp_counts <- sapply(ct9,function(CT) sum(CDESC$cell_id == CT))
  temp_mFDR <- sapply(ct9,function(CT)
    rowMeans(GSEA_score[,CID[CDESC$cell_id == CT]]))
  Pscore_CT <- pbsapply(names(temp_counts),function(I) {
    temp <- rowSums(apply(BKGD_CT[[as.character(temp_counts[I])]],2,
                          function(X) abs(X) >= abs(temp_mFDR[,I])))
    temp <- p.adjust(temp / ncol(BKGD_CT[[as.character(temp_counts[I])]]),method="fdr")
    temp[temp <= 0] <- 1 / ncol(BKGD_CT[[as.character(temp_counts[I])]])
    return(-log10(temp) * sign(temp_mFDR[,I]))
  })
  rm(list=c("BKGD_CT",grep("^temp",ls(),value=T)))  
  
  # ^ Lig ----
  print("^ Lig ----")
  if (file.exists(paste0("~/Dropbox/GDB_archive/CMapCorr_files/pfdrGSEA_",GSEA,"_BKGD_LIG.RData"))) {
    load(paste0("~/Dropbox/GDB_archive/CMapCorr_files/pfdrGSEA_",GSEA,"_BKGD_LIG.RData"))
  } else {
    print("^^ Generating null distribution ----")
    temp_counts_uniq <- sort(unique(
      sapply(lig295,function(LIG) sum(CDESC$pert_iname == LIG)),
    ),decreasing=T)
    BKGD_LIG <- pbsapply(temp_counts_uniq,function(N) {
      replicate(minP,{
        rowMeans(GSEA_score[,sample(colnames(GSEA_score),N)])
      })
    },simplify=F)
    names(BKGD_LIG) <- temp_counts_uniq
    save(BKGD_LIG,
         file=paste0("~/Dropbox/GDB_archive/CMapCorr_files/pfdrGSEA_",GSEA,"_BKGD_LIG.RData"))
  }
  
  print("^^ Calculating ----")
  temp_counts <- sapply(lig295,function(LIG) sum(CDESC$pert_iname == LIG))
  temp_mFDR <- sapply(lig295,function(LIG)
    rowMeans(GSEA_score[,CID[CDESC$pert_iname == LIG]]))
  Pscore_Lig <- pbsapply(names(temp_counts),function(I) {
    temp <- rowSums(apply(BKGD_LIG[[as.character(temp_counts[I])]],2,
                          function(X) abs(X) >= abs(temp_mFDR[,I])))
    temp <- p.adjust(temp / ncol(BKGD_LIG[[as.character(temp_counts[I])]]),method="fdr")
    temp[temp <= 0] <- 1 / ncol(BKGD_LIG[[as.character(temp_counts[I])]])
    return(-log10(temp) * sign(temp_mFDR[,I]))
  })
  rm(list=c("BKGD_LIG",grep("^temp",ls(),value=T)))  

  # ^ LigCT ----
  print("^ LigCT ----")
  if (file.exists(paste0("~/Dropbox/GDB_archive/CMapCorr_files/pfdrGSEA_",GSEA,"_BKGD_LIGCT.RData"))) {
    load(paste0("~/Dropbox/GDB_archive/CMapCorr_files/pfdrGSEA_",GSEA,"_BKGD_LIGCT.RData"))
  } else {
    print("^^ Generating null distribution ----")
    temp_counts_uniq <- sort(unique(as.vector(
      sapply(ct9,function(CT)
        sapply(lig295,function(LIG) 
          sum(CDESC$pert_iname == LIG & CDESC$cell_id == CT))),
    )),decreasing=T)
    BKGD_LIGCT <- pbsapply(temp_counts_uniq,function(N) {
      replicate(minP,{
        rowMeans(GSEA_score[,sample(colnames(GSEA_score),N)])
      })
    },simplify=F)
    names(BKGD_LIGCT) <- temp_counts_uniq
    save(BKGD_LIGCT,
         file=paste0("~/Dropbox/GDB_archive/CMapCorr_files/pfdrGSEA_",GSEA,"_BKGD_LIGCT.RData"))
  }
  
  print("^^ Calculating ----")
  temp_counts <- as.vector(sapply(ct9,function(CT)
    sapply(lig295,function(LIG) 
      sum(CDESC$pert_iname == LIG & CDESC$cell_id == CT))))
  names(temp_counts) <- as.vector(sapply(ct9,function(X) paste(X,lig295,sep="_")))
  
  temp_mFDR <- sapply(ct9,function(CT)
    sapply(lig295,function(LIG)
      rowMeans(GSEA_score[,CID[CDESC$pert_iname == LIG &
                                                 CDESC$cell_id == CT]])),
    simplify=F)
  temp_mFDR <- do.call(cbind,temp_mFDR)
  colnames(temp_mFDR) <- names(temp_counts)
  
  Pscore_LigCT <- pbsapply(names(temp_counts),function(I) {
    temp <- rowSums(apply(BKGD_LIGCT[[as.character(temp_counts[I])]],2,
                          function(X) abs(X) >= abs(temp_mFDR[,I])))
    temp <- p.adjust(temp / ncol(BKGD_LIGCT[[as.character(temp_counts[I])]]),method="fdr")
    temp[temp <= 0] <- 1 / ncol(BKGD_LIGCT[[as.character(temp_counts[I])]])
    return(-log10(temp) * sign(temp_mFDR[,I]))
  })
  rm(list=c("BKGD_LIGCT",grep("^temp",ls(),value=T)))
  save(list=grep("^Pscore",ls(),value=T),
       file=paste0("~/Dropbox/GDB_archive/CMapCorr_files/pfdrGSEA_",GSEA,".RData"))
  rm(list=grep("^Pscore",ls(),value=T))
}
