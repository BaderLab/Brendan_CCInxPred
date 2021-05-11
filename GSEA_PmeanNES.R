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
lvl4_data@cdesc <- lvl4_data@cdesc[temp_id,]
lvl4_data@cid <- temp_id

rm(list=c("ct14","lig16",grep("^temp",ls(),value=T)))


# Load GSEA results ----
GSEA_nes <- list() 
# GSEA_padj <- list()
# GSEA_score <- list()
for (FN in list.files("~/Dropbox/GDB_archive/CMapCorr_files/","^GSEA_lvl4_[A-Z]")) {
  temp <- load(paste0("~/Dropbox/GDB_archive/CMapCorr_files/",FN))
  temp_name <- strsplit(temp,"_")[[1]][3]
  
  GSEA_nes[[temp_name]] <- sapply(get(temp)[-1],"[[","NES")
  rownames(GSEA_nes[[temp_name]]) <- as.vector(get(temp)$PATHWAYS)
  
  # GSEA_padj[[temp_name]] <- sapply(get(temp)[-1],"[[","padj")
  # rownames(GSEA_padj[[temp_name]]) <- as.vector(get(temp)$PATHWAYS)
  
  # GSEA_score[[temp_name]] <- -log10(GSEA_padj[[temp_name]]) * sign(GSEA_nes[[temp_name]])
  # rownames(GSEA_score[[temp_name]]) <- as.vector(get(temp)$PATHWAYS)
  
  rm(list=temp)
}
rm(list=c("FN",grep("^temp",ls(),value=T)))

minP <- 1e4

for (GSEA in rev(names(GSEA_nes))) {
  if (file.exists(paste0("~/Dropbox/GDB_archive/CMapCorr_files/pnesGSEA_",GSEA,".RData"))) { next }
  print(paste(GSEA,"----"))
  # Calc BKGD ----
  print("^ BKGD ----")
  COUNTS <- sort(unique(c(
    sapply(ct9,function(CT) sum(lvl4_data@cdesc$cell_id == CT)),
    sapply(lig295,function(LIG) sum(lvl4_data@cdesc$pert_iname == LIG)),
    sapply(ct9,function(CT) 
      sapply(lig295,function(LIG) 
        sum(lvl4_data@cdesc$cell_id == CT & lvl4_data@cdesc$pert_iname == LIG)))
  )),decreasing=T)
  BKGD <- pbsapply(COUNTS,function(N) {
    replicate(minP,{
      rowMeans(GSEA_nes[[GSEA]][,sample(colnames(GSEA_nes[[GSEA]]),N)])
    })
  },simplify=F)
  
  # Overlap ----
  # ^ CT ----
  print("^ CT ----")
  temp_counts <- sapply(ct9,function(CT) sum(lvl4_data@cdesc$cell_id == CT))
  temp_mNES <- sapply(ct9,function(CT)
    rowMeans(GSEA_nes[[GSEA]][,lvl4_data@cid[lvl4_data@cdesc$cell_id == CT]]))
  Pscore_CT <- pbsapply(names(temp_counts),function(I) {
    temp <- rowSums(apply(BKGD[[which(COUNTS == temp_counts[I])]],2,
                          function(X) abs(X) >= abs(temp_mNES[,I])))
    temp <- p.adjust(temp / minP,method="fdr")
    temp[temp <= 0] <- 1 / minP
    return(-log10(temp) * sign(temp_mNES[,I]))
  })
  rm(list=grep("^temp",ls(),value=T))
  
  # ^ Lig ----
  print("^ Lig ----")
  temp_counts <- sapply(lig295,function(LIG) sum(lvl4_data@cdesc$pert_iname == LIG))
  temp_mNES <- sapply(lig295,function(LIG)
    rowMeans(GSEA_nes[[GSEA]][,lvl4_data@cid[lvl4_data@cdesc$pert_iname == LIG]]))
  Pscore_Lig <- pbsapply(names(temp_counts),function(I) {
    temp <- rowSums(apply(BKGD[[which(COUNTS == temp_counts[I])]],2,
                          function(X) abs(X) >= abs(temp_mNES[,I])))
    temp <- p.adjust(temp / minP,method="fdr")
    temp[temp <= 0] <- 1 / minP
    return(-log10(temp) * sign(temp_mNES[,I]))
  })
  rm(list=grep("^temp",ls(),value=T))
  
  # ^ LigCT ----
  print("^ LigCT ----")
  temp_counts <- as.vector(sapply(ct9,function(CT)
    sapply(lig295,function(LIG) 
      sum(lvl4_data@cdesc$pert_iname == LIG & lvl4_data@cdesc$cell_id == CT))))
  names(temp_counts) <- as.vector(sapply(ct9,function(X) paste(X,lig295,sep="_")))
  
  temp_mNES <- sapply(ct9,function(CT)
    sapply(lig295,function(LIG)
      rowMeans(GSEA_nes[[GSEA]][,lvl4_data@cid[lvl4_data@cdesc$pert_iname == LIG &
                                                 lvl4_data@cdesc$cell_id == CT]])),
    simplify=F)
  temp_mNES <- do.call(cbind,temp_mNES)
  colnames(temp_mNES) <- names(temp_counts)
  
  Pscore_LigCT <- pbsapply(names(temp_counts),function(I) {
    temp <- rowSums(apply(BKGD[[which(COUNTS == temp_counts[I])]],2,
                          function(X) abs(X) >= abs(temp_mNES[,I])))
    temp <- p.adjust(temp / minP,method="fdr")
    temp[temp <= 0] <- 1 / minP
    return(-log10(temp) * sign(temp_mNES[,I]))
  })
  rm(list=grep("^temp",ls(),value=T))
  
  
  save(list=grep("^Pscore",ls(),value=T),
       file=paste0("~/Dropbox/GDB_archive/CMapCorr_files/pnesGSEA_",GSEA,".RData"))
  rm(list=grep("^Pscore",ls(),value=T))
}
