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

# Calc BKGD ----


# Overlap ----

# ^ ct NOT EDITTED! ----
for (GSEA in names(GSEA_nes)) {
  pSCORE_ct <- pbsapply(ct9,function(CT) {
    id <- lvl4_data@cid[lvl4_data@cdesc$cell_id == CT]
    meanSCORE <- rowMeans(GSEA_score[[GSEA]][,id,drop=F])
    bkgd <- sapply(1:minP,function(L) 
      rowMeans(GSEA_score[[GSEA]][,sample(ncol(GSEA_score[[GSEA]]),length(id))]))
    pval <- rowMeans(apply(bkgd,2,function(X) abs(X) >= abs(meanSCORE)))
    pval[pval == 0] <- 1/ncol(bkgd)
    return(pval)
  },cl=5)
  save(pSCORE_ct,
       file=paste0("~/Dropbox/GDB_archive/CMapCorr_files/pfdrGSEA_",
                   GSEA,"_ct.RData"))
  rm(pSCORE_ct)
}


# ^ lig ----
for (GSEA in names(GSEA_nes)) {
  pSCORE_lig <- pbsapply(lig295,function(LIG) {
    id <- lvl4_data@cid[lvl4_data@cdesc$pert_iname == LIG]
    meanSCORE <- rowMeans(GSEA_nes[[GSEA]][,id,drop=F])
    bkgd <- sapply(1:minP,function(L) 
      rowMeans(GSEA_score[[GSEA]][,sample(ncol(GSEA_score[[GSEA]]),length(id))]))
    pval <- sapply(seq_along(meanSCORE),function(X) 
      sum(abs(meanSCORE[X]) <= abs(bkgd[X,])) / minP)
    pval[pval <= 0] <- 1 / minP

    return(pval)
  },cl=5)
  save(pSCORE_lig,
       file=paste0("~/Dropbox/GDB_archive/CMapCorr_files/pfdrGSEA_",
                   GSEA,"_lig.RData"))
  rm(pSCORE_lig)
}


# permute gene sets?





# ^ ligct ----
for (GSEA in names(GSEA_nes)) {
  pSCORE_ligct <- pbsapply(lig295,function(LIG) {
    sapply(ct9,function(CT) {
      id <- lvl4_data@cid[lvl4_data@cdesc$pert_iname == LIG &
                            lvl4_data@cdesc$cell_id == CT]
      meanSCORE <- rowMeans(GSEA_score[[GSEA]][,id,drop=F])
      bkgd <- sapply(1:minP,function(L) 
        rowMeans(GSEA_score[[GSEA]][,sample(ncol(GSEA_score[[GSEA]]),length(id))]))
      pval <- sapply(seq_along(meanSCORE),function(X) 
        sum(abs(meanSCORE[X]) <= abs(bkgd[X,])) / minP)
      return(pval)
    })
  },simplify=F,cl=5)
  save(pSCORE_ligct,
       file=paste0("~/Dropbox/GDB_archive/CMapCorr_files/pfdrGSEA_",
                   GSEA,"_ligct.RData"))
  rm(pSCORE_ligct)
}


# ^ rep ----
for (GSEA in names(GSEA_nes)) {
  temp_tx <- unique(lvl4_data@cdesc[,c("pert_iname","cell_id","pert_dose","pert_time")])
  rownames(temp_tx) <- apply(temp_tx,1,paste,collapse="_")
  
  pSCORE_rep <- pbapply(temp_tx,1,function(X) {
    id <- lvl4_data@cid[lvl4_data@cdesc$pert_iname == X[1] &
                          lvl4_data@cdesc$cell_id == X[2] &
                          lvl4_data@cdesc$pert_dose == X[3] &
                          lvl4_data@cdesc$pert_time == X[4]]
    meanSCORE <- rowMeans(GSEA_score[[GSEA]][,id,drop=F])
    bkgd <- sapply(1:minP,function(L) 
      rowMeans(GSEA_score[[GSEA]][,sample(ncol(GSEA_score[[GSEA]]),length(id))]))
    pval <- sapply(seq_along(meanSCORE),function(X) 
      sum(abs(meanSCORE[X]) <= abs(bkgd[X,])) / minP)
    return(pval)
  },cl=5)
  save(pSCORE_rep,
       file=paste0("~/Dropbox/GDB_archive/CMapCorr_files/pfdrGSEA_",
                   GSEA,"_rep.RData"))
  rm(pSCORE_rep)
}
