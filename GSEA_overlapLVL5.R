library(pbapply)
pboptions(type="timer")

# load lvl5 ----
load("~/Dropbox/GDB_archive/CMapCorr_files/lvl5_inputs_allgenes.RData")
temp_ligcountsperct <- sapply(unique(lvl5_data@cdesc$cell_id),function(CT)
  length(unique(lvl5_data@cdesc[lvl5_data@cdesc$cell_id == CT,"pert_iname"])))
ct9 <- names(temp_ligcountsperct)[temp_ligcountsperct > 100]
names(ct9) <- sapply(ct9,function(X) names(ct14)[ct14 == X])

temp_ctcountsperlig <- sapply(unique(lvl5_data@cdesc$pert_iname),function(LIG)
  length(unique(lvl5_data@cdesc[lvl5_data@cdesc$pert_iname == LIG &
                                  lvl5_data@cdesc$cell_id %in% ct9,"cell_id"])))
lig295 <- names(temp_ctcountsperlig)[temp_ctcountsperlig == 9]
lig295 <- sort(lig295)

temp_id <- rownames(lvl5_data@cdesc)[lvl5_data@cdesc$pert_iname %in% lig295 & 
                                       lvl5_data@cdesc$cell_id %in% ct9]
lvl5_data@mat <- lvl5_data@mat[,temp_id]
lvl5_data@cdesc <- lvl5_data@cdesc[temp_id,]
lvl5_data@cid <- temp_id
rm(list=c("lig16","ct14",grep("^temp",ls(),value=T)))

cutFDRlist <- c(0.1,0.05,0.01)

# load GSEA ----
GSEA_FILENAMES <- rev(list.files("~/Dropbox/GDB_archive/CMapCorr_files/","^GSEA_lvl5",full.names=T))
for (GSEA_FN in GSEA_FILENAMES) {
  GSEA <- load(GSEA_FN)
  if (file.exists(
    paste0("~/Dropbox/GDB_archive/CMapCorr_files/povlpGSEA_",
           sub("gsea_nes_","",GSEA),".RData")
  )) { rm(list=GSEA); next }
  print(paste(GSEA,"----"))
  
  # calc background ----
  print("^ Background ----")
  
  # ^ CT ----
  print("^^ CT ----")
  CT_counts <- sort(unique(sapply(ct9,function(CT) sum(lvl5_data@cdesc$cell_id == CT))))
  CT_bkgd <- pbsapply(CT_counts,function(N) {
    temp <- sapply(1:1e4,function(X) {
      id <- sample(names(get(GSEA)[-1]),N)
      apply(sapply(get(GSEA)[id],"[[","padj"),1,function(X) 
        sapply(cutFDRlist,function(Y) sum(X <= Y)))
    })
    return(sapply(c("1"=1,"05"=2,"01"=3),function(X) 
      temp[seq(X,nrow(temp),by=3),],simplify=F))
  },cl=3)


  # ^ Lig ----
  print("^^ Lig ----")
  LIG_counts <- sort(unique(sapply(lig295,function(LIG) sum(lvl5_data@cdesc$pert_iname == LIG))))
  LIG_bkgd <- pbsapply(LIG_counts,function(N) {
    temp <- sapply(1:1e4,function(X) {
      id <- sample(names(get(GSEA)[-1]),N)
      apply(sapply(get(GSEA)[id],"[[","padj"),1,function(X) 
        sapply(cutFDRlist,function(Y) sum(X <= Y)))
    })
    return(sapply(c("1"=1,"05"=2,"01"=3),function(X) 
      temp[seq(X,nrow(temp),by=3),],simplify=F))
  },cl=4)
  

  # ^ LigCT ----
  print("^^ LigCT ----")
  LIGCT_counts <- sort(unique(as.vector(
    sapply(lig295,function(LIG) 
      sapply(ct9,function(CT)
        sum(lvl5_data@cdesc$pert_iname == LIG &
              lvl5_data@cdesc$cell_id == CT)))
  )))
  
  LIGCT_bkgd <- pbsapply(LIGCT_counts,function(N) {
    temp <- sapply(1:1e4,function(X) {
      id <- sample(names(get(GSEA)[-1]),N)
      apply(sapply(get(GSEA)[id],"[[","padj"),1,function(X) 
        sapply(cutFDRlist,function(Y) sum(X <= Y)))
    })
    return(sapply(c("1"=1,"05"=2,"01"=3),function(X) 
      temp[seq(X,nrow(temp),by=3),],simplify=F))
  },cl=5)
  colnames(LIGCT_bkgd) <- LIGCT_counts
  
  LIGCT_names <- as.vector(sapply(lig295,function(X) paste(X,ct9,sep="_")))
  LIGCT_counts <- sapply(lig295,function(LIG) 
    sapply(ct9,function(CT) 
      sum(lvl5_data@cdesc$pert_iname == LIG &
            lvl5_data@cdesc$cell_id == CT)
    ),simplify=F)
  LIGCT_counts <- do.call(c,LIGCT_counts)
  names(LIGCT_counts) <- LIGCT_names
  

  # calc p-val ----  
  # probability of this number of samples sharing each 
  # enriched (at `cutFDR`) gene set happening by chance
  
  print("^ Calc p-val ----")  
  Poverlap_ct <- Poverlap_lig <- Poverlap_ligct <- list()
  Noverlap_ct <- Noverlap_lig <- Noverlap_ligct <- list()
  for (cutFDR in cutFDRlist) {
    print(paste("^^ FDR",cutFDR,"----"))
    L <- sub("0.","",cutFDR,fixed=T)

    # ^ CT ----
    print("^^^ CT ----")
    Noverlap_ct[[L]] <- sapply(ct9,function(CT) {
      id <- lvl5_data@cid[lvl5_data@cdesc$cell_id == CT]
      apply(sapply(get(GSEA)[id],"[[","padj"),1,function(X) sum(X <= cutFDR))
    })
    Poverlap_ct[[L]] <- pbsapply(colnames(Noverlap_ct[[L]]),function(CT) {
      rowMeans(apply(
        CT_bkgd[[L,which(CT_counts == sum(lvl5_data@cdesc$cell_id == ct9[CT]))]],
        MARGIN=2,function(X) X >= Noverlap_ct[[L]][,CT]))
    })
    Poverlap_ct[[L]][Poverlap_ct[[L]] %in% 0] <- 1 / ncol(CT_bkgd[[1]])
    
    # ^ Lig ----
    print("^^^ Lig ----")
    Noverlap_lig[[L]] <- sapply(lig295,function(LIG) {
      id <- lvl5_data@cid[lvl5_data@cdesc$pert_iname == LIG]
      apply(sapply(get(GSEA)[id],"[[","padj"),1,function(X) sum(X <= cutFDR))
    })
    Poverlap_lig[[L]] <- pbsapply(colnames(Noverlap_lig[[L]]),function(LIG) {
      rowMeans(apply(
        LIG_bkgd[[L,which(LIG_counts == sum(lvl5_data@cdesc$pert_iname == LIG))]],
        MARGIN=2,function(X) X >= Noverlap_lig[[L]][,LIG]))
    })
    Poverlap_lig[[L]][Poverlap_lig[[L]] %in% 0] <- 1 / ncol(LIG_bkgd[[1]])
    
    # ^ LigCT ----
    print("^^^ LigCT ----")
    Noverlap_ligct[[L]] <- sapply(lig295,function(LIG) 
      sapply(ct9,function(CT) {
        id <- lvl5_data@cid[lvl5_data@cdesc$pert_iname == LIG &
                              lvl5_data@cdesc$cell_id == CT]
        apply(sapply(get(GSEA)[id],"[[","padj"),1,function(X) sum(X <= cutFDR))
      }),simplify=F)
    Noverlap_ligct[[L]] <- do.call(cbind,Noverlap_ligct[[L]])
    colnames(Noverlap_ligct[[L]]) <- LIGCT_names
    
    Poverlap_ligct[[L]] <- pbsapply(colnames(Noverlap_ligct[[L]]),function(LIGCT) {
      rowMeans(apply(
        LIGCT_bkgd[[L,as.character(LIGCT_counts[LIGCT])]],
        MARGIN=2,function(X) X >= Noverlap_ligct[[L]][,LIGCT]))
    })
    Poverlap_ligct[[L]][Poverlap_ligct[[L]] %in% 0] <- 1 / ncol(LIGCT_bkgd[[1]])
  }
  
  save(list=grep("overlap",ls(),value=T),
       file=paste0("~/Dropbox/GDB_archive/CMapCorr_files/povlpGSEA_",
                   sub("gsea_nes_","",GSEA),".RData"))
  rm(list=c("L",GSEA))
}



