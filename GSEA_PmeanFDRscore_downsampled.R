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
DS <- 5
allCuts <- seq(log10(minP),0,by=-0.01)

# for (FN in rev(list.files("~/Dropbox/GDB_archive/CMapCorr_files/","^GSEA_lvl4_[A-Z]"))) {
for (FN in rev(list.files("~/Dropbox/GDB_archive/CMapCorr_files/","^GSEA_lvl4_[A-Z]"))[2]) {
  temp <- load(paste0("~/Dropbox/GDB_archive/CMapCorr_files/",FN))
  GSEA <- strsplit(temp,"_")[[1]][3]
  if (file.exists(paste0("~/Dropbox/GDB_archive/CMapCorr_files/pfdrGSEA_",GSEA,"_DS.RData"))) { 
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

  
  # ^ BKGD ----
  # if (file.exists(paste0("~/Dropbox/GDB_archive/CMapCorr_files/pfdrGSEA_",GSEA,"_DS_BKGD.RData"))) {
  #   load(paste0("~/Dropbox/GDB_archive/CMapCorr_files/pfdrGSEA_",GSEA,"_DS_BKGD.RData"))
  # } else {
    print("^ Generating null distribution ----")
    BKGD <- pbreplicate(minP,{
      rowMeans(GSEA_score[,sample(colnames(GSEA_score),DS)])
    })
  #   save(BKGD,
  #        file=paste0("~/Dropbox/GDB_archive/CMapCorr_files/pfdrGSEA_",GSEA,"_DS_BKGD.RData"))
  # }
  

  # ^ CT ----
  print("^ CT ----")
  
  temp_counts <- sapply(ct9,function(CT) sum(CDESC$cell_id == CT))
  
  Mscore_CT <- pbsapply(ct9,function(CT) 
    replicate(ceiling(max(temp_counts / DS)),{
      rowMeans(GSEA_score[,sample(CID[CDESC$cell_id == CT],DS)])
    }),simplify=F)
  names(Mscore_CT) <- ct9
  
  Pscore_CT <- pbsapply(Mscore_CT,function(X) {
    temp <- apply(X,2,function(Y) rowMeans(abs(BKGD) >= abs(Y)))
    temp[temp <= 0] <- 1 / ncol(BKGD)
    return(temp)
  },simplify=F)
  
  NPhits_CT <- pbsapply(allCuts,function(CUT) 
    sapply(Pscore_CT,function(X) mean(-log10(X) >= CUT)))
  
  Qscore_CT <- sapply(Pscore_CT,p.adjust,method="fdr",simplify=F)
  Qscore_CT <- sapply(Qscore_CT,matrix,
                      nrow=nrow(Pscore_CT[[1]]),
                      dimnames=list(rownames(Pscore_CT[[1]])),
                      simplify=F)
  
  NQhits_CT <- pbsapply(allCuts,function(CUT) 
    sapply(Qscore_CT,function(X) mean(-log10(X) >= CUT)))
  
  rm(list=grep("^temp",ls(),value=T))
  

  # ^ Lig ----
  print("^ Lig ----")
  
  temp_counts <- sapply(lig295,function(LIG) sum(CDESC$pert_iname == LIG))
  
  Mscore_LIG <- pbsapply(lig295,function(LIG) 
    replicate(ceiling(max(temp_counts / DS)),{
      rowMeans(GSEA_score[,sample(CID[CDESC$pert_iname == LIG],DS)])
    }),simplify=F)
  
  Pscore_LIG <- pbsapply(Mscore_LIG,function(X) {
    temp <- apply(X,2,function(Y) rowMeans(abs(BKGD) >= abs(Y)))
    temp[temp <= 0] <- 1 / ncol(BKGD)
    return(temp)
  },simplify=F)
  
  NPhits_LIG <- pbsapply(allCuts,function(CUT) 
    sapply(Pscore_LIG,function(X) mean(-log10(X) >= CUT)))

  Qscore_LIG <- sapply(Pscore_LIG,p.adjust,method="fdr",simplify=F)
  Qscore_LIG <- sapply(Qscore_LIG,matrix,
                       nrow=nrow(Pscore_LIG[[1]]),
                       dimnames=list(rownames(Pscore_LIG[[1]])),
                       simplify=F)
  
  NQhits_LIG <- pbsapply(allCuts,function(CUT) 
    sapply(Qscore_LIG,function(X) mean(-log10(X) >= CUT)))
  
  rm(list=grep("^temp",ls(),value=T))
  
  
  # ^ LigCT ----
  print("^ LigCT ----")
  temp_counts <- as.vector(sapply(ct9,function(CT)
    sapply(lig295,function(LIG) 
      sum(CDESC$pert_iname == LIG & CDESC$cell_id == CT))))
  names(temp_counts) <- as.vector(sapply(ct9,function(X) paste(X,lig295,sep="_")))
  
  temp_ids <- sapply(ct9,function(CT)
    sapply(lig295,function(LIG)
      CID[CDESC$pert_iname == LIG & CDESC$cell_id == CT],
      simplify=F),
    simplify=F)
  temp_ids <- unlist(temp_ids,recursive=F,use.names=F)
  names(temp_ids) <- names(temp_counts)
  temp_ids <- temp_ids[temp_counts >= 5]
  
  Mscore_LIGCT <- sapply(temp_ids,function(ID)
    replicate(ceiling(max(temp_counts / DS)),{
      rowMeans(GSEA_score[,sample(ID,DS)])
    }),simplify=F)
  
  Pscore_LIGCT <- pbsapply(Mscore_LIGCT,function(X) {
    temp <- apply(X,2,function(Y) rowMeans(abs(BKGD) >= abs(Y)))
    temp[temp <= 0] <- 1 / ncol(BKGD)
    return(temp)
  },simplify=F)

  NPhits_LIGCT <- pbsapply(allCuts,function(CUT) 
    sapply(Pscore_LIGCT,function(X) mean(-log10(X) >= CUT)))
  
  Qscore_LIGCT <- sapply(Pscore_LIGCT,p.adjust,method="fdr",simplify=F)
  Qscore_LIGCT <- sapply(Qscore_LIGCT,matrix,
                         nrow=nrow(Pscore_LIGCT[[1]]),
                         dimnames=list(rownames(Pscore_LIGCT[[1]])),
                         simplify=F)
  
  NQhits_LIGCT <- pbsapply(allCuts,function(CUT) 
    sapply(Qscore_LIGCT,function(X) mean(-log10(X) >= CUT)))
  
  rm(list=grep("^temp",ls(),value=T))
  
    
  save(list=c(grep("score",ls(),value=T),grep("N[PQ]hits",ls(),value=T)),
       file=paste0("~/Dropbox/GDB_archive/CMapCorr_files/pfdrGSEA_",GSEA,"_DS.RData"))
  rm(list=c("BKGD",grep("score",ls(),value=T),grep("N[PQ]hits",ls(),value=T)))
}
