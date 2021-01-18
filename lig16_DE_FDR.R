library(cmapR)
library(colorspace)
library(pbapply)
pboptions(type="timer")
.PAR <- par(no.readonly=T)


load("~/Dropbox/GDB_archive/CMapCorr_files/lvl4_lig16_inputs.RData") 

ALPHA <- c(0.1,0.05,0.01,0.001)
names(ALPHA) <- c("10","5","1","01")




# LIG mean Z-score (unweighted) ----
meanZ_lig <- sapply(lig16,function(LIG)
  rowMeans(lvl4_data@mat[,lvl4_data@cdesc$pert_iname == LIG]))

FDR_lig <- apply(meanZ_lig,2,function(X) p.adjust(pnorm(-abs(X))))

samples_lig <- sapply(lig16,function(X) sum(lvl4_data@cdesc$pert_iname == X))

BKG_lig <- pbsapply(lig16,function(LIG) {
  temp <- sapply(1:1e4,function(L) 
    p.adjust(pnorm(-abs(
      rowMeans(lvl4_data@mat[,sample(ncol(lvl4_data@mat),samples_lig[LIG])])
    )))
  )
  return(
    sapply(ALPHA,function(Z) 
      apply(temp,2,function(Y) sum(Y <= Z)))
  )
},simplify=F,cl=8)

pDE_lig <- sapply(ALPHA,function(Z) {
  temp_count <- apply(FDR_lig,2,function(X) sum(X <= Z))
  sapply(lig16,function(LIG) {
    sum(BKG_lig[[LIG]][,names(which(ALPHA == Z))] >= temp_count[LIG]) / nrow(BKG_lig[[LIG]])
  })
})

save(meanZ_lig,FDR_lig,pDE_lig,
     file=paste0("~/Dropbox/GDB_archive/CMapCorr_files/lig16_DE_lig_FDR.RData"))
rm(list=grep("lig$",ls(),value=T))




# CT mean Z-score (unweighted) ----
meanZ_ct <- sapply(ct14,function(CT)
  rowMeans(lvl4_data@mat[,lvl4_data@cdesc$cell_id == CT]))

FDR_ct <- apply(meanZ_ct,2,function(X) p.adjust(pnorm(-abs(X))))

samples_ct <- sapply(ct14,function(X) sum(lvl4_data@cdesc$cell_id == X))

BKG_ct <- pbsapply(names(ct14),function(CT) {
  temp <- sapply(1:1e4,function(L) 
    p.adjust(pnorm(-abs(
      rowMeans(lvl4_data@mat[,sample(ncol(lvl4_data@mat),samples_ct[CT])])
    )))
  )
  return(
    sapply(ALPHA,function(Z) 
      apply(temp,2,function(Y) sum(Y <= Z)))
  )
},simplify=F,cl=8)

pDE_ct <- sapply(ALPHA,function(Z) {
  temp_count <- apply(FDR_ct,2,function(X) sum(X <= Z))
  sapply(names(ct14),function(CT) {
    sum(BKG_ct[[CT]][,names(which(ALPHA == Z))] >= temp_count[CT]) / nrow(BKG_ct[[CT]])
  })
})

save(meanZ_ct,FDR_ct,pDE_ct,
     file=paste0("~/Dropbox/GDB_archive/CMapCorr_files/lig16_DE_ct_FDR.RData"))
rm(list=grep("ct$",ls(),value=T))




# LIG / CT mean Z-score (unweighted) ----
meanZ_ligct <- sapply(ct14,function(CT) {
  sapply(lig16,function(LIG) {
    rowMeans(lvl4_data@mat[,lvl4_data@cdesc$cell_id == CT &
                             lvl4_data@cdesc$pert_iname == LIG])
  })
},simplify=F)
temp_colnames <- as.vector(sapply(names(ct14),function(X) paste(X,lig16,sep="_")))
meanZ_ligct <- do.call(cbind,meanZ_ligct)
colnames(meanZ_ligct) <- temp_colnames

FDR_ligct <- apply(meanZ_ligct,2,function(X) p.adjust(pnorm(-abs(X))))

samples_ligct <- sapply(ct14,function(CT) {
  sapply(lig16,function(LIG) {
    sum(lvl4_data@cdesc$cell_id == CT &
          lvl4_data@cdesc$pert_iname == LIG)
  })
})
samples_ligct <- as.vector(samples_ligct)
names(samples_ligct) <- colnames(FDR_ligct)

BKG_ligct <- pbsapply(samples_ligct,function(X) {
  temp <- sapply(1:1e4,function(L) 
    p.adjust(pnorm(-abs(
      rowMeans(lvl4_data@mat[,sample(ncol(lvl4_data@mat),X)])
    )))
  )
  return(
    sapply(ALPHA,function(Z) 
      apply(temp,2,function(Y) sum(Y <= Z)))
  )
},simplify=F,cl=8)

pDE_ligct <- sapply(ALPHA,function(Z) {
  temp_count <- apply(FDR_ligct,2,function(X) sum(X <= Z))
  sapply(colnames(FDR_ligct),function(Y) {
    sum(BKG_ligct[[Y]][,names(which(ALPHA == Z))] >= temp_count[Y]) / nrow(BKG_ligct[[Y]])
  })
})

save(meanZ_ligct,FDR_ligct,pDE_ligct,
     file=paste0("~/Dropbox/GDB_archive/CMapCorr_files/lig16_DE_ligct_FDR.RData"))
rm(list=grep("ligct$",ls(),value=T))
rm(list=grep("^temp",ls(),value=T))




# REP mean Z-score (unweighted) ----
temp_tx <- unique(lvl4_data@cdesc[,c("pert_iname","cell_id","pert_dose","pert_time")])
rownames(temp_tx) <- apply(temp_tx,1,paste,collapse="_")

meanZ_rep <- apply(temp_tx,1,function(X) {
  temp <- lvl4_data@cdesc$pert_iname == X[1] &
    lvl4_data@cdesc$cell_id == X[2] &
    lvl4_data@cdesc$pert_dose == X[3] &
    lvl4_data@cdesc$pert_time == X[4]
  return(
    switch((sum(temp) > 1) + 1,
           lvl4_data@mat[,temp],
           rowMeans(lvl4_data@mat[,temp]))
  )
})

FDR_rep <- apply(meanZ_rep,2,function(X) p.adjust(pnorm(-abs(X))))

samples_rep <- apply(temp_tx,1,function(X) 
  sum(lvl4_data@cdesc$pert_iname == X[1] &
        lvl4_data@cdesc$cell_id == X[2] &
        lvl4_data@cdesc$pert_dose == X[3] &
        lvl4_data@cdesc$pert_time == X[4]))

BKG_rep <- pbsapply(samples_rep,function(X) {
  temp <- sapply(1:1e4,function(L) 
    p.adjust(pnorm(-abs(switch(
      (X > 1) + 1,
      lvl4_data@mat[,sample(ncol(lvl4_data@mat),X)],
      rowMeans(lvl4_data@mat[,sample(ncol(lvl4_data@mat),X)])
    ))))
  )
  return(
    sapply(ALPHA,function(Z) 
      apply(temp,2,function(Y) sum(Y <= Z)))
  )
},simplify=F,cl=8)

pDE_rep <- sapply(ALPHA,function(Z) {
  temp_count <- apply(FDR_rep,2,function(X) sum(X <= Z))
  sapply(colnames(FDR_rep),function(Y) {
    sum(BKG_rep[[Y]][,names(which(ALPHA == Z))] >= temp_count[Y]) / nrow(BKG_rep[[Y]])
  })
})

save(meanZ_rep,FDR_rep,pDE_rep,
     file=paste0("~/Dropbox/GDB_archive/CMapCorr_files/lig16_DE_rep_FDR.RData"))
rm(list=grep("ligct$",ls(),value=T))
rm(list=grep("^temp",ls(),value=T))



