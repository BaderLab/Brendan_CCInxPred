library(cmapR)
library(colorspace)
library(pbapply)
pboptions(type="timer")
.PAR <- par(no.readonly=T)


load("~/Dropbox/GDB_archive/CMapCorr_files/lvl4_inputs_allgenes.RData")
lvl4_data <- lvl4_data_all
rm(lvl4_data_all)
lvl4_data@cdesc <- lvl4_data@cdesc[lvl4_data@cdesc$pert_iname %in% lig16,]
lvl4_data@mat <- lvl4_data@mat[,rownames(lvl4_data@cdesc)]


ALPHA <- c(0.1,0.05,0.01,0.001)
names(ALPHA) <- c("10","5","1","01")




# LIG mean Z-score (unweighted) ----
FDR_lig <- sapply(lig16,function(LIG)
  p.adjust(pnorm(-abs(
    rowMeans(lvl4_data@mat[,lvl4_data@cdesc$pert_iname == LIG])
  )))
)

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

save(FDR_lig,pDE_lig,
     file=paste0("~/Dropbox/GDB_archive/CMapCorr_files/lig16_DE_allgenes_lig_FDR.RData"))
rm(list=grep("lig$",ls(),value=T))




# CT mean Z-score (unweighted) ----
FDR_ct <- sapply(ct14,function(CT)
  p.adjust(pnorm(-abs(
    rowMeans(lvl4_data@mat[,lvl4_data@cdesc$cell_id == CT])
  )))
)

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

save(FDR_ct,pDE_ct,
     file=paste0("~/Dropbox/GDB_archive/CMapCorr_files/lig16_DE_allgenes_ct_FDR.RData"))
rm(list=grep("ct$",ls(),value=T))




# LIG / CT mean Z-score (unweighted) ----
FDR_ligct <- sapply(ct14,function(CT) {
  sapply(lig16,function(LIG) {
    p.adjust(pnorm(-abs(
      rowMeans(lvl4_data@mat[,lvl4_data@cdesc$cell_id == CT &
                               lvl4_data@cdesc$pert_iname == LIG])
    )))
  })
},simplify=F)
temp_colnames <- as.vector(sapply(names(ct14),function(X) paste(X,lig16,sep="_")))
FDR_ligct <- do.call(cbind,FDR_ligct)
colnames(FDR_ligct) <- temp_colnames

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

save(FDR_ligct,pDE_ligct,
     file=paste0("~/Dropbox/GDB_archive/CMapCorr_files/lig16_DE_allgenes_ligct_FDR.RData"))
rm(list=grep("ligct$",ls(),value=T))
rm(list=grep("^temp",ls(),value=T))




# REP mean Z-score (unweighted) ----
temp_tx <- unique(lvl4_data@cdesc[,c("pert_iname","cell_id","pert_dose","pert_time")])
rownames(temp_tx) <- apply(temp_tx,1,paste,collapse="_")

FDR_rep <- apply(temp_tx,1,function(X) {
  temp <- lvl4_data@cdesc$pert_iname == X[1] &
    lvl4_data@cdesc$cell_id == X[2] &
    lvl4_data@cdesc$pert_dose == X[3] &
    lvl4_data@cdesc$pert_time == X[4]
  return(
    p.adjust(pnorm(-abs(
      switch((sum(temp) > 1) + 1,
             lvl4_data@mat[,temp],
             rowMeans(lvl4_data@mat[,temp]))
    )))
  )
})

samples_rep <- apply(temp_tx,1,function(X) 
  sum(lvl4_data@cdesc$pert_iname == X[1] &
        lvl4_data@cdesc$cell_id == X[2] &
        lvl4_data@cdesc$pert_dose == X[3] &
        lvl4_data@cdesc$pert_time == X[4]))

BKG_rep <- pbsapply(samples_rep,function(X) {
  temp <- sapply(1:1e1,function(L) 
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

save(FDR_rep,pDE_rep,
     file=paste0("~/Dropbox/GDB_archive/CMapCorr_files/lig16_DE_allgenes_rep_FDR.RData"))
rm(list=grep("ligct$",ls(),value=T))
rm(list=grep("^temp",ls(),value=T))



