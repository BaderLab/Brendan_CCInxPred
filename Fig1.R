library(cmapR)
library(colorspace)
library(pbapply)
.PAR <- par(no.readonly=T)


if (exists("lvl4_data")) {
} else if (file.exists("~/Dropbox/GDB_archive/CMapCorr_files/lvl4_lig16_inputs.RData")) {
  load("~/Dropbox/GDB_archive/CMapCorr_files/lvl4_lig16_inputs.RData") 
} else {
  source("lvl4_lig16_inputs.R")
}


# LIG mean Z-score (unweighted) ----
meanZlig_all <- sapply(lig16,function(LIG)
  rowMeans(lvl4_data@mat[,lvl4_data@cdesc$pert_iname == LIG]))

countDElig <- apply(meanZlig_all,2,function(X) sum(X > 1.645))
samples_lig <- sapply(lig16,function(X) sum(lvl4_data@cdesc$pert_iname == X))

pDElig <- pbsapply(lig16,function(LIG) {
  temp <- sapply(1:1e4,function(Z) 
    sum(rowMeans(lvl4_data@mat[,sample(ncol(lvl4_data@mat),samples_lig[LIG])]) > 1.645))
  return(sum(temp >= countDElig[LIG]) / length(temp))
},cl=8)

save(meanZlig_all,countDElig,pDElig,
     file="~/Dropbox/GDB/CMapCorr/Fig1_lig.RData")


# CT mean Z-score (unweighted) ----

meanZct_all <- sapply(ct14,function(CT)
  rowMeans(lvl4_data@mat[,lvl4_data@cdesc$cell_id == CT]))

countDEct <- apply(meanZct_all,2,function(X) sum(X > 1.645))
samples_ct <- sapply(ct14,function(X) sum(lvl4_data@cdesc$cell_id == X))

pDEct <- pbsapply(names(ct14),function(CT) {
  temp <- sapply(1:1e4,function(Z) 
    sum(rowMeans(lvl4_data@mat[,sample(ncol(lvl4_data@mat),samples_ct[CT])]) > 1.645))
  return(sum(temp >= countDEct[CT]) / length(temp))
},cl=8)

save(meanZct_all,countDEct,pDEct,
     file="~/Dropbox/GDB/CMapCorr/Fig1_ct.RData")


# LIG / CT mean Z-score (unweighted) ----

meanZligct <- sapply(ct14,function(CT) {
  sapply(lig16,function(LIG) {
    rowMeans(lvl4_data@mat[,lvl4_data@cdesc$cell_id == CT &
                             lvl4_data@cdesc$pert_iname == LIG])
  })
},simplify=F)
temp_colnames <- as.vector(sapply(names(ct14),function(X) paste(X,lig16,sep="_")))
meanZligct_all <- do.call(cbind,meanZligct)
colnames(meanZligct_all) <- temp_colnames

countDEligct <- sapply(meanZligct,function(X)
  apply(X,2,function(Y) sum(Y > 1.645)))

samples_ligct <- sapply(ct14,function(CT) {
  sapply(lig16,function(LIG) {
    sum(lvl4_data@cdesc$cell_id == CT &
          lvl4_data@cdesc$pert_iname == LIG)
  })
})

pDEligct <- pbsapply(names(ct14),function(CT) {
  sapply(lig16,function(LIG) {
    temp <- sapply(1:1e4,function(Z) 
      sum(rowMeans(lvl4_data@mat[,sample(ncol(lvl4_data@mat),samples_ligct[LIG,CT])]) > 1.645))
    return(sum(temp >= countDEligct[LIG,CT]) / length(temp))
  })
},cl=8)

scoreDEligct <- pDEligct + 1e-4
scoreDEligct[scoreDEligct > 1] <- 1
scoreDEligct <- -log10(scoreDEligct)

save(meanZligct,countDEligct,pDEligct,scoreDEligct,
     file="~/Dropbox/GDB/CMapCorr/Fig1_ligct.RData")


# REP mean Z-score (unweighted) ----
rep_tx <- unique(lvl4_data@cdesc[,c("pert_iname","cell_id","pert_dose","pert_time")])
rownames(rep_tx) <- apply(rep_tx,1,paste,collapse="_")

meanZrep <- apply(rep_tx,1,function(X) {
  temp <- lvl4_data@cdesc$pert_iname == X[1] &
    lvl4_data@cdesc$cell_id == X[2] &
    lvl4_data@cdesc$pert_dose == X[3] &
    lvl4_data@cdesc$pert_time == X[4]
  if (sum(temp) > 1) {
    return(rowMeans(lvl4_data@mat[,temp]))
  } else {
    return(lvl4_data@mat[,temp])
  }
})

countDErep <- apply(meanZrep,2,function(X) sum(X > 1.645))
samples_rep <- apply(rep_tx,1,function(X) 
  sum(lvl4_data@cdesc$pert_iname == X[1] &
        lvl4_data@cdesc$cell_id == X[2] &
        lvl4_data@cdesc$pert_dose == X[3] &
        lvl4_data@cdesc$pert_time == X[4]))

pDErep <- pbsapply(names(samples_rep),function(X) {
  temp <- sapply(1:1e4,function(Z) {
    if (samples_rep[X] > 1) {
      return(sum(rowMeans(lvl4_data@mat[,sample(ncol(lvl4_data@mat),samples_rep[X])]) > 1.645))
    } else {
      return(sum(lvl4_data@mat[,sample(ncol(lvl4_data@mat),samples_rep[X])] > 1.645))
    }
  })
  return(sum(temp >= countDErep[X]) / length(temp))
},cl=8)

scoreDErep <- pDErep + 1e-4
scoreDErep[scoreDErep > 1] <- 1
scoreDErep <- -log10(scoreDErep)

save(meanZrep,countDErep,pDErep,scoreDErep,
     file="~/Dropbox/GDB/CMapCorr/Fig1_rep.RData")
