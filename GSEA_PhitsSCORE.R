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


allCuts <- seq(4,0,by=-0.01)

# par(mfrow=c(1,3),mar=c(3,3,2,1),mgp=2:0)
# require(colorspace)
# 
# for (FN in list.files("~/Dropbox/GDB_archive/CMapCorr_files/","^pfdrGSEA_[A-Z]+\\.RData")) {
#   temp <- load(paste0("~/Dropbox/GDB_archive/CMapCorr_files/",FN))
#   GSEA <- gsub("(pfdrGSEA_)|(\\.RData)","",FN)
#   
#   Nhits_CT <- pbsapply(allCuts,function(CUT) mean(-log10(Qscore_CT) >= CUT))
#   Nhits_Lig <- pbsapply(allCuts,function(CUT) mean(-log10(Qscore_LIG) >= CUT))
#   Nhits_LigCT <- pbsapply(allCuts,function(CUT) mean(-log10(Qscore_LIGCT) >= CUT))
#   
#   plot(NA,NA,xlim=c(4,0),ylim=c(0,1),xaxt="n",main=GSEA,
#        ylab="Proportion of all |mean signed FDR| \u2264 x",xlab=NA)
#   axis(side=1,labels=c(.0001,.001,.01,.05,.1,.5,1),
#        at=-log10(c(.0001,.001,.01,.05,.1,.5,1)))
#   mtext("Gene set |mean signed FDR| threshold",side=1,line=2)
#   lines(allCuts,Nhits_CT,lty=1,lwd=2,
#         col=qualitative_hcl(3,palette="dark3")[1])
#   lines(allCuts,Nhits_Lig,lty=1,lwd=2,
#         col=qualitative_hcl(3,palette="dark3")[2])
#   lines(allCuts,Nhits_LigCT,lty=1,lwd=2,
#         col=qualitative_hcl(3,palette="dark3")[3])
#   legend("topleft",legend=c("Cell line","Ligand","Ligand per cell line"),
#          bty="n",col=qualitative_hcl(3,palette="dark3"),lty=1,lwd=2)
#   
#   save(list=grep("^Nhits",ls(),value=T),
#        file=paste0("~/Dropbox/GDB_archive/CMapCorr_files/nhitsGSEA_",GSEA,".RData"))
# }  




# Old and sketchy? ----
# The sketchy part is Nhits / Phits - we should be making ECDFs from the Q values (NQhits),
# not the averages scores (Nhits).  Taking probabilities of Nhits <= chance isn't terribly meaningful.

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
  load(paste0("~/Dropbox/GDB_archive/CMapCorr_files/pfdrGSEA_",GSEA,"_BKGD_CT.RData"))
  temp_counts <- sapply(ct9,function(CT) sum(CDESC$cell_id == CT))
  Nhits_CT <- pbsapply(allCuts,function(CUT) colMeans(abs(Mscore_CT) >= CUT))
  Phits_CT <- pbsapply(seq_along(allCuts),function(Z) {
    temp_hit_BKGD <- sapply(BKGD_CT,function(X) colMeans(abs(X) >= allCuts[Z]))
    OUT <- sapply(names(temp_counts),function(L)
      mean(temp_hit_BKGD[,as.character(temp_counts[L])] >= Nhits_CT[L,Z]))
    OUT[OUT <= 0] <- 1 / ncol(BKGD_CT[[1]])
    return(OUT)
  })
  NQhits_CT <- pbsapply(allCuts,function(CUT) colMeans(-log10(Qscore_CT) >= CUT))
  rm(list=c("BKGD_CT",grep("^temp",ls(),value=T)))  
  
  # ^ Lig ----
  print("^ Lig ----")
  load(paste0("~/Dropbox/GDB_archive/CMapCorr_files/pfdrGSEA_",GSEA,"_BKGD_LIG.RData"))
  temp_counts <- sapply(lig295,function(LIG) sum(CDESC$pert_iname == LIG))
  Nhits_LIG <- pbsapply(allCuts,function(CUT) colMeans(abs(Mscore_LIG) >= CUT))
  Phits_LIG <- pbsapply(seq_along(allCuts),function(Z) {
    temp_hit_BKGD <- sapply(BKGD_LIG,function(X) colMeans(abs(X) >= allCuts[Z]))
    OUT <- sapply(names(temp_counts),function(L)
      mean(temp_hit_BKGD[,as.character(temp_counts[L])] >= Nhits_LIG[L,Z]))
    OUT[OUT <= 0] <- 1 / ncol(BKGD_LIG[[1]])
    return(OUT)
  })
  NQhits_LIG <- pbsapply(allCuts,function(CUT) colMeans(-log10(Qscore_LIG) >= CUT))
  rm(list=c("BKGD_LIG",grep("^temp",ls(),value=T)))  
  
  # ^ LigCT ----
  print("^ LigCT ----")
  load(paste0("~/Dropbox/GDB_archive/CMapCorr_files/pfdrGSEA_",GSEA,"_BKGD_LIGCT.RData"))
  temp_counts <- as.vector(sapply(ct9,function(CT)
    sapply(lig295,function(LIG) 
      sum(CDESC$pert_iname == LIG & CDESC$cell_id == CT))))
  names(temp_counts) <- as.vector(sapply(ct9,function(X) paste(X,lig295,sep="_")))
  Nhits_LIGCT <- pbsapply(allCuts,function(CUT) colMeans(abs(Mscore_LIGCT) >= CUT))
  Phits_LIGCT <- pbsapply(seq_along(allCuts),function(Z) {
    temp_hit_BKGD <- sapply(BKGD_LIGCT,function(X) colMeans(abs(X) >= allCuts[Z]))
    OUT <- sapply(names(temp_counts),function(L)
      mean(temp_hit_BKGD[,as.character(temp_counts[L])] >= Nhits_LIGCT[L,Z]))
    OUT[OUT <= 0] <- 1 / ncol(BKGD_LIGCT[[1]])
    return(OUT)
  })
  NQhits_LIGCT <- pbsapply(allCuts,function(CUT) colMeans(-log10(Qscore_LIGCT) >= CUT))
  rm(list=c("BKGD_LIGCT",grep("^temp",ls(),value=T)))
  
  save(list=grep("hits_",ls(),value=T),
       file=paste0("~/Dropbox/GDB_archive/CMapCorr_files/phitsGSEA_",GSEA,".RData"))
  rm(list=grep("hits_",ls(),value=T))
}
