library(cmapR)
library(colorspace)
library(pbapply)
pboptions(type="timer")
.PAR <- par(no.readonly=T)


temp <- load("~/Dropbox/GDB_archive/CMapCorr_files/lvl4_inputs.RData") 

temp_ligcountsperct <- sapply(unique(lvl4_data@cdesc$cell_id),function(CT)
  length(unique(lvl4_data@cdesc[lvl4_data@cdesc$cell_id == CT,"pert_iname"])))
ct9 <- names(temp_ligcountsperct)[temp_ligcountsperct > 100]
names(ct9) <- sapply(ct9,function(X) names(ct14)[ct14 == X])

lvl4_data@cdesc <- lvl4_data@cdesc[lvl4_data@cdesc$cell_id %in% ct9,]
lvl4_data@mat <- lvl4_data@mat[,rownames(lvl4_data@cdesc)]

temp_ctcountsperlig <- sapply(unique(lvl4_data@cdesc$pert_iname),function(LIG)
  length(unique(lvl4_data@cdesc[lvl4_data@cdesc$pert_iname == LIG,"cell_id"])))

lig295 <- names(temp_ctcountsperlig)[temp_ctcountsperlig == 9]

lvl4_data@cdesc <- lvl4_data@cdesc[lvl4_data@cdesc$pert_iname %in% lig295,]
lvl4_data@mat <- lvl4_data@mat[,rownames(lvl4_data@cdesc)]

rm(list=temp[temp != "lvl4_data"])
rm(list=grep("^temp",ls(),value=T))


Zpercent <- c(1.282,1.645,2.326,3.09)
names(Zpercent) <- c("90","95","99","99.9")



# LIG mean Z-score (unweighted) ----
meanZlig_all <- sapply(lig295,function(LIG)
  rowMeans(lvl4_data@mat[,lvl4_data@cdesc$pert_iname == LIG]))

samples_lig <- sapply(lig295,function(X) sum(lvl4_data@cdesc$pert_iname == X))

for (Z in Zpercent) {
  countDElig <- apply(meanZlig_all,2,function(X) sum(X > Z))
  
  pDElig <- pbsapply(lig295,function(LIG) {
    temp <- sapply(1:1e4,function(L) 
      sum(rowMeans(lvl4_data@mat[,sample(ncol(lvl4_data@mat),samples_lig[LIG])]) > Z))
    return(sum(temp >= countDElig[LIG]) / length(temp))
  },cl=8)
  
  save(meanZlig_all,countDElig,pDElig,
       file=paste0("~/Dropbox/GDB_archive/CMapCorr_files/lig295_DE_lig_",
                   names(Zpercent)[Zpercent == Z],
                   ".RData"))
}



# CT mean Z-score (unweighted) ----

meanZct_all <- sapply(ct9,function(CT)
  rowMeans(lvl4_data@mat[,lvl4_data@cdesc$cell_id == CT]))

samples_ct <- sapply(ct9,function(X) sum(lvl4_data@cdesc$cell_id == X))

for (Z in Zpercent) {
  countDEct <- apply(meanZct_all,2,function(X) sum(X > Z))
  
  pDEct <- sapply(names(ct9),function(CT) {
    print(paste(which(CT == names(ct9)),"/",length(ct9)))
    temp <- pbsapply(1:1e4,function(L) 
      sum(rowMeans(lvl4_data@mat[,sample(ncol(lvl4_data@mat),samples_ct[CT])]) > Z),
      cl=8)
    return(sum(temp >= countDEct[CT]) / length(temp))
  })
  
  save(meanZct_all,countDEct,pDEct,
       file=paste0("~/Dropbox/GDB_archive/CMapCorr_files/lig295_DE_ct_",
                   names(Zpercent)[Zpercent == Z],
                   ".RData"))
}



# LIG / CT mean Z-score (unweighted) ----

meanZligct <- sapply(ct9,function(CT) {
  sapply(lig295,function(LIG) {
    rowMeans(lvl4_data@mat[,lvl4_data@cdesc$cell_id == CT &
                             lvl4_data@cdesc$pert_iname == LIG])
  })
},simplify=F)
temp_colnames <- as.vector(sapply(names(ct9),function(X) paste(X,lig295,sep="_")))
meanZligct_all <- do.call(cbind,meanZligct)
colnames(meanZligct_all) <- temp_colnames

samples_ligct <- sapply(ct9,function(CT) {
  sapply(lig295,function(LIG) {
    sum(lvl4_data@cdesc$cell_id == CT &
          lvl4_data@cdesc$pert_iname == LIG)
  })
})

for (Z in Zpercent) {
  countDEligct <- sapply(meanZligct,function(X)
    apply(X,2,function(Y) sum(Y > Z)))
  
  pDEligct <- pbsapply(names(ct9),function(CT) {
    sapply(lig295,function(LIG) {
      temp <- sapply(1:1e4,function(L) 
        sum(rowMeans(lvl4_data@mat[,sample(ncol(lvl4_data@mat),samples_ligct[LIG,CT])]) > Z))
      return(sum(temp >= countDEligct[LIG,CT]) / length(temp))
    })
  },cl=8)
  
  scoreDEligct <- pDEligct + 1e-4
  scoreDEligct[scoreDEligct > 1] <- 1
  scoreDEligct <- -log10(scoreDEligct)
  
  save(meanZligct,countDEligct,pDEligct,scoreDEligct,
       file=paste0("~/Dropbox/GDB_archive/CMapCorr_files/lig295_DE_ligct_",
                   names(Zpercent)[Zpercent == Z],
                   ".RData"))
}



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

samples_rep <- apply(rep_tx,1,function(X) 
  sum(lvl4_data@cdesc$pert_iname == X[1] &
        lvl4_data@cdesc$cell_id == X[2] &
        lvl4_data@cdesc$pert_dose == X[3] &
        lvl4_data@cdesc$pert_time == X[4]))

for (Z in Zpercent) {
  countDErep <- apply(meanZrep,2,function(X) sum(X > Z))
  
  pDErep <- pbsapply(names(samples_rep),function(X) {
    temp <- sapply(1:1e4,function(L) {
      if (samples_rep[X] > 1) {
        return(sum(rowMeans(lvl4_data@mat[,sample(ncol(lvl4_data@mat),samples_rep[X])]) > Z))
      } else {
        return(sum(lvl4_data@mat[,sample(ncol(lvl4_data@mat),samples_rep[X])] > Z))
      }
    })
    return(sum(temp >= countDErep[X]) / length(temp))
  },cl=8)
  
  scoreDErep <- pDErep + 1e-4
  scoreDErep[scoreDErep > 1] <- 1
  scoreDErep <- -log10(scoreDErep)
  
  save(meanZrep,samples_rep,countDErep,pDErep,scoreDErep,
       file=paste0("~/Dropbox/GDB_archive/CMapCorr_files/lig295_DE_rep_",
                   names(Zpercent)[Zpercent == Z],
                   ".RData"))
}

# pDE background ----
# pDEbackground <- pbsapply(2:max(table(lvl4_data@cdesc$cell_id)),function(X) {
#   sapply(1:1000,function(Y) sum(rowMeans(lvl4_data@mat[,sample(ncol(lvl4_data@mat),X)]) > 1.645))
# },cl=8)
# pDEbackground <- cbind(
#   sapply(1:1000,function(Y) sum(lvl4_data@mat[,sample(ncol(lvl4_data@mat),1)] > 1.645)),
#   pDEbackground
# )
# save(pDEbackground,file="~/Dropbox/GDB/CMapCorr/Fig3_pDEbackground.RData")
