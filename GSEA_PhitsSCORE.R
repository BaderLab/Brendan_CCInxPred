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


minP <- 1e4
allCuts <- seq(0.1,2,by=0.1)
names(allCuts) <- allCuts

for (FN in list.files("~/Dropbox/GDB_archive/CMapCorr_files/","^pfdrGSEA_[A-Z]+\\.RData")) {
  temp <- load(paste0("~/Dropbox/GDB_archive/CMapCorr_files/",FN))
  GSEA <- gsub("(pfdrGSEA_)|(\\.RData)","",FN)
  if (file.exists(paste0("~/Dropbox/GDB_archive/CMapCorr_files/phitsGSEA_",GSEA,".RData"))) { 
    rm(list=temp)
    next 
  }

  print(paste(GSEA,"----"))
  
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
  Mscore_CT <- sapply(ct9,function(CT)
    rowMeans(GSEA_score[,CID[CDESC$cell_id == CT]]))
  
  Nhits_CT <- sapply(allCuts,function(CUT) colSums(abs(Mscore_CT) >= CUT))
  
  Phits_CT <- pbsapply(allCuts,function(CUT) {
    temp_hit_BKGD <- sapply(BKGD_CT,function(X) colSums(abs(X) >= CUT))
    OUT <- sapply(names(temp_counts),function(L)
      mean(temp_hit_BKGD[,as.character(temp_counts[L])] >= Nhits_CT[L,as.character(CUT)]))
    OUT[OUT <= 0] <- 1 / ncol(BKGD_CT[[1]])
    return(OUT)
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
  Mscore_LIG <- sapply(lig295,function(LIG)
    rowMeans(GSEA_score[,CID[CDESC$pert_iname == LIG]]))
  
  Nhits_LIG <- sapply(allCuts,function(CUT) colSums(abs(Mscore_LIG) >= CUT))
  
  Phits_LIG <- pbsapply(allCuts,function(CUT) {
    temp_hit_BKGD <- sapply(BKGD_LIG,function(X) colSums(abs(X) >= CUT))
    OUT <- sapply(names(temp_counts),function(L)
      mean(temp_hit_BKGD[,as.character(temp_counts[L])] >= Nhits_LIG[L,as.character(CUT)]))
    OUT[OUT <= 0] <- 1 / ncol(BKGD_LIG[[1]])
    return(OUT)
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
  
  Mscore_LIGCT <- sapply(ct9,function(CT)
    sapply(lig295,function(LIG)
      rowMeans(GSEA_score[,CID[CDESC$pert_iname == LIG &
                                                 CDESC$cell_id == CT]])),
    simplify=F)
  Mscore_LIGCT <- do.call(cbind,Mscore_LIGCT)
  colnames(Mscore_LIGCT) <- names(temp_counts)
  
  Nhits_LIGCT <- sapply(allCuts,function(CUT) colSums(abs(Mscore_LIGCT) >= CUT))
  
  Phits_LIGCT <- pbsapply(allCuts,function(CUT) {
    temp_hit_BKGD <- sapply(BKGD_LIGCT,function(X) colSums(abs(X) >= CUT))
    OUT <- sapply(names(temp_counts),function(L)
      mean(temp_hit_BKGD[,as.character(temp_counts[L])] >= Nhits_LIGCT[L,as.character(CUT)]))
    OUT[OUT <= 0] <- 1 / ncol(BKGD_LIGCT[[1]])
    return(OUT)
  })
  rm(list=c("BKGD_LIGCT",grep("^temp",ls(),value=T)))
  
  save(list=c(grep("^Phits",ls(),value=T),
              grep("^Nhits",ls(),value=T),
              grep("^Mscore",ls(),value=T)),
       file=paste0("~/Dropbox/GDB_archive/CMapCorr_files/phitsGSEA_",GSEA,".RData"))
  rm(list=c(grep("^Phits",ls(),value=T),
            grep("^Nhits",ls(),value=T),
            grep("^Mscore",ls(),value=T)))
}
