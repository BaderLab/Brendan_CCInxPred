library(pbapply)
load("~/Dropbox/GDB_archive/CMapCorr_files/NN_ALLvCT_dat.RData")

# temp_lfc <- unlist(sapply(nn_lig_rep,function(X) as.vector(X$lfc)),use.names=F)
# temp_fdr <- unlist(sapply(nn_lig_rep,function(X) as.vector(X$qval)),use.names=F)
# 
# par(mar=c(3,3,1,1),mgp=2:0)
# plot(NA,NA,xlim=range(temp_lfc),ylim=range(-log10(temp_fdr)),
#      xlab="LogFC",ylab="-log10(FDR)")
# points(temp_lfc,-log10(temp_fdr),pch=".",col=scales::alpha("black",0.1))
# abline(h=-log10(0.1),col="red")
# abline(v=c(-1,1),col="red")



# cut_lfc <- log2(c(1.5,2,3,4))
# cut_fdr <- c(0.1,0.05,0.01,0.001)
cut_lfc <- log2(c(2,4))
cut_fdr <- c(0.1,0.01)


CUTOFFS <- t(do.call(cbind,lapply(cut_lfc,function(LFC) sapply(cut_fdr,function(FDR) c(LFC=LFC,FDR=FDR)))))

DEovlp_Pmat <- list()
for (L in 1:nrow(CUTOFFS)) {
  print(paste(L,"/",nrow(CUTOFFS)))
  # LFC <- CUTOFFS[16,"LFC"]; FDR <- CUTOFFS[16,"FDR"]
  LFC <- CUTOFFS[L,"LFC"]; FDR <- CUTOFFS[L,"FDR"] 
  # DE overlap background ----
  nn_DE <- list()
  temp_commongenes <- Reduce(intersect,sapply(nn_lig_rep,function(X) rownames(X$lfc)))
  for (LIG in nn_ligands) {
    nn_DE[[LIG]] <- mapply(intersect,
                           # apply(nn_lig_rep[[LIG]]$lfc,2,function(X) names(which(X >= 1))),
                           apply(nn_lig_rep[[LIG]]$lfc,2,function(X) names(which(abs(X) >= LFC))),
                           apply(nn_lig_rep[[LIG]]$qval,2,function(X) names(which(X <= FDR))),
                           SIMPLIFY=F)
    # nn_DE[[LIG]] <- sapply(nn_DE[[LIG]],function(X) X[X %in% temp_commongenes])
  }
  
  # par(mar=c(3,3,1,1),mgp=2:0)
  # hist(sapply(unlist(nn_DE,recursive=F),length),
  #      xlab="# DE per dataset (|logFC| >= 1, FDR <= 0.1)",main=NA)
  
  temp_DE <- unlist(nn_DE,recursive=F)
  DEbkgd <- pbsapply(min(sapply(DSinfo,nrow)):max(sapply(DSinfo,nrow)),function(X)
    sapply(1:1e6,function(Y)
      length(Reduce(intersect,sample(temp_DE,X)))
    ),cl=4)
  colnames(DEbkgd) <- min(sapply(DSinfo,nrow)):max(sapply(DSinfo,nrow))
  rm(temp_DE)
  
  
  # DE overlap calc ----
  de_all <- de_ct <- list()
  P_all <- P_ct <- list()
  for (LIG in nn_ligands) {
    temp_ct <- sapply(unique(DSinfo[[LIG]]$cell_type),function(X) 
      rownames(DSinfo[[LIG]])[DSinfo[[LIG]]$cell_type == X],
      simplify=F)
    temp_ds <- c()
    for (X in names(temp_ct)[sapply(temp_ct,length) > 1]) {
      temp_ds <- append(temp_ds,paste(temp_ct[[X]],collapse="."))
      de_ct[[LIG]][[X]] <- Reduce(intersect,nn_DE[[LIG]][temp_ct[[X]]])
      P_ct[[LIG]][[X]] <- sum(DEbkgd[,as.character(length(temp_ct[[X]]))] >= 
                                length(de_ct[[LIG]][[X]])) / nrow(DEbkgd)
      P_ct[[LIG]][[X]][P_ct[[LIG]][[X]] == 0] <- 0.1 / nrow(DEbkgd)
    }
    if (!paste(rownames(DSinfo[[LIG]]),collapse=".") %in% temp_ds) {
      de_all[[LIG]] <- Reduce(intersect,nn_DE[[LIG]])
      P_all[[LIG]] <- sum(DEbkgd[,as.character(length(nn_DE[[LIG]]))] >= 
                                 length(de_all[[LIG]])) / nrow(DEbkgd)
      P_all[[LIG]][P_all[[LIG]] == 0] <- 0.1 / nrow(DEbkgd)
    }
  }
  
  rm(list=c("X","LIG",grep("^temp",ls(),value=T)))
  save(nn_DE,de_all,de_ct,P_all,P_ct,
       file=paste0(
         "~/Dropbox/GDB_archive/CMapCorr_files/NN_ALLvCT_DEovlp_",
         "LFC",sub(".","",2^LFC,fixed=T),
         "FDR",sub("0.","",FDR,fixed=T),
         ".RData"))
  
  DEovlp_Pmat[[L]] <- list(ALL=unlist(P_all),CT=unlist(P_ct))
  rm(nn_DE,de_all,de_ct,P_all,P_ct)
}


Pmat <- sapply(DEovlp_Pmat,function(X) wilcox.test(X$ALL,X$CT)$p.value)
Pmat <- matrix(Pmat,nrow=4)
colnames(Pmat) <- paste("FC",unique(2^CUTOFFS[,"LFC"]))
rownames(Pmat) <- paste0("FDR ",unique(CUTOFFS[,"FDR"] * 100),"%")

Pmat

save(CUTOFFS,DEovlp_Pmat,Pmat,
     file="~/Dropbox/GDB_archive/CMapCorr_files/NN_ALLvCT_DEovlp_Pmatrix.RData")


# LFC 1,2
# FDR 10%, 1%