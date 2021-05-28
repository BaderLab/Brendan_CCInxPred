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

for (FN in list.files("~/Dropbox/GDB_archive/CMapCorr_files/","^GSEA_lvl4_[A-Z]")) {
  temp <- load(paste0("~/Dropbox/GDB_archive/CMapCorr_files/",FN))
  GSEA <- strsplit(temp,"_")[[1]][3]
  if (file.exists(paste0("~/Dropbox/GDB_archive/CMapCorr_files/pnesGSEA_",GSEA,".RData"))) { 
    rm(list=temp)
    next 
  }

  print(paste(GSEA,"----"))
  
  GSEA_NES <- sapply(get(temp)[-1],"[[","NES")
  rownames(GSEA_NES) <- as.vector(get(temp)$PATHWAYS)

  # make NES into a normal(ish) distribution:
  # GSEA_NES[abs(GSEA_NES) < 1] <- 0
  # GSEA_NES[GSEA_NES < 0] <- GSEA_NES[GSEA_NES < 0] + 1
  # GSEA_NES[GSEA_NES > 0] <- GSEA_NES[GSEA_NES > 0] - 1  
  # |NES values| <= 1 have at best an FDR around 0.5, so fuck them.

  # Make NES => |NES| and do it one-sided?
  GSEA_NES <- abs(GSEA_NES)
    
  rm(list=c(temp,grep("^temp",ls(),value=T)))


  # Mean signed -log10(FDR) ----
  # ^ CT ----
  print("^ CT ----")
  if (file.exists(paste0("~/Dropbox/GDB_archive/CMapCorr_files/pnesGSEA_",GSEA,"_BKGD_CT.RData"))) {
    load(paste0("~/Dropbox/GDB_archive/CMapCorr_files/pnesGSEA_",GSEA,"_BKGD_CT.RData"))
  } else {
    print("^^ Generating null distribution ----")
    temp_counts_uniq <- sort(unique(
      sapply(ct9,function(CT) sum(CDESC$cell_id == CT)),
    ),decreasing=T)
    BKGD_CT <- pbsapply(temp_counts_uniq,function(N) {
      replicate(minP,{
        rowMeans(GSEA_NES[,sample(colnames(GSEA_NES),N)])
      })
    },simplify=F)
    names(BKGD_CT) <- temp_counts_uniq
    save(BKGD_CT,
         file=paste0("~/Dropbox/GDB_archive/CMapCorr_files/pnesGSEA_",GSEA,"_BKGD_CT.RData"))
  }
  
  print("^^ Calculating ----")
  temp_counts <- sapply(ct9,function(CT) sum(CDESC$cell_id == CT))
  Mnes_CT <- sapply(ct9,function(CT)
    rowMeans(GSEA_NES[,CID[CDESC$cell_id == CT]]))
  
  Pnes_CT <- pbsapply(colnames(Mnes_CT),function(X) 
    rowMeans(abs(Mnes_CT[,X]) < abs(BKGD_CT[[as.character(temp_counts[X])]])))
  Pnes_CT[Pnes_CT <= 0] <- 1 / ncol(BKGD_CT[[1]])
  Qnes_CT <- apply(Pnes_CT,2,p.adjust,method="fdr")

  rm(list=c("BKGD_CT",grep("^temp",ls(),value=T)))  
  
  # ^ Lig ----
  print("^ Lig ----")
  if (file.exists(paste0("~/Dropbox/GDB_archive/CMapCorr_files/pnesGSEA_",GSEA,"_BKGD_LIG.RData"))) {
    load(paste0("~/Dropbox/GDB_archive/CMapCorr_files/pnesGSEA_",GSEA,"_BKGD_LIG.RData"))
  } else {
    print("^^ Generating null distribution ----")
    temp_counts_uniq <- sort(unique(
      sapply(lig295,function(LIG) sum(CDESC$pert_iname == LIG)),
    ),decreasing=T)
    BKGD_LIG <- pbsapply(temp_counts_uniq,function(N) {
      replicate(minP,{
        rowMeans(GSEA_NES[,sample(colnames(GSEA_NES),N)])
      })
    },simplify=F)
    names(BKGD_LIG) <- temp_counts_uniq
    save(BKGD_LIG,
         file=paste0("~/Dropbox/GDB_archive/CMapCorr_files/pnesGSEA_",GSEA,"_BKGD_LIG.RData"))
  }
  
  print("^^ Calculating ----")
  temp_counts <- sapply(lig295,function(LIG) sum(CDESC$pert_iname == LIG))
  Mnes_LIG <- sapply(lig295,function(LIG)
    rowMeans(GSEA_NES[,CID[CDESC$pert_iname == LIG]]))
  
  Pnes_LIG <- pbsapply(colnames(Mnes_LIG),function(X) 
    rowMeans(abs(Mnes_LIG[,X]) < abs(BKGD_LIG[[as.character(temp_counts[X])]])))
  Pnes_LIG[Pnes_LIG <= 0] <- 1 / ncol(BKGD_LIG[[1]])
  Qnes_LIG <- apply(Pnes_LIG,2,p.adjust,method="fdr")

  rm(list=c("BKGD_LIG",grep("^temp",ls(),value=T)))  

  # ^ LigCT ----
  print("^ LigCT ----")
  if (file.exists(paste0("~/Dropbox/GDB_archive/CMapCorr_files/pnesGSEA_",GSEA,"_BKGD_LIGCT.RData"))) {
    load(paste0("~/Dropbox/GDB_archive/CMapCorr_files/pnesGSEA_",GSEA,"_BKGD_LIGCT.RData"))
  } else {
    print("^^ Generating null distribution ----")
    temp_counts_uniq <- sort(unique(as.vector(
      sapply(ct9,function(CT)
        sapply(lig295,function(LIG) 
          sum(CDESC$pert_iname == LIG & CDESC$cell_id == CT))),
    )),decreasing=T)
    BKGD_LIGCT <- pbsapply(temp_counts_uniq,function(N) {
      replicate(minP,{
        rowMeans(GSEA_NES[,sample(colnames(GSEA_NES),N)])
      })
    },simplify=F)
    names(BKGD_LIGCT) <- temp_counts_uniq
    save(BKGD_LIGCT,
         file=paste0("~/Dropbox/GDB_archive/CMapCorr_files/pnesGSEA_",GSEA,"_BKGD_LIGCT.RData"))
  }
  
  print("^^ Calculating ----")
  temp_counts <- as.vector(sapply(ct9,function(CT)
    sapply(lig295,function(LIG) 
      sum(CDESC$pert_iname == LIG & CDESC$cell_id == CT))))
  names(temp_counts) <- as.vector(sapply(ct9,function(X) paste(X,lig295,sep="_")))
  
  Mnes_LIGCT <- sapply(ct9,function(CT)
    sapply(lig295,function(LIG)
      rowMeans(GSEA_NES[,CID[CDESC$pert_iname == LIG &
                                                 CDESC$cell_id == CT]])),
    simplify=F)
  Mnes_LIGCT <- do.call(cbind,Mnes_LIGCT)
  colnames(Mnes_LIGCT) <- names(temp_counts)
  
  Pnes_LIGCT <- pbsapply(colnames(Mnes_LIGCT),function(X) 
    rowMeans(abs(Mnes_LIGCT[,X]) < abs(BKGD_LIGCT[[as.character(temp_counts[X])]])))
  Pnes_LIGCT[Pnes_LIGCT <= 0] <- 1 / ncol(BKGD_LIGCT[[1]])
  Qnes_LIGCT <- apply(Pnes_LIGCT,2,p.adjust,method="fdr")
  
  rm(list=c("BKGD_LIGCT",grep("^temp",ls(),value=T)))
  
  save(list=grep("nes",ls(),value=T),
       file=paste0("~/Dropbox/GDB_archive/CMapCorr_files/pnesGSEA_",GSEA,".RData"))
  rm(list=grep("nes",ls(),value=T))
}
