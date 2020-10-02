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


same_lig <- sapply(lig16,function(X) rownames(lvl4_data@cdesc)[lvl4_data@cdesc$pert_iname == X])
same_lig_cor <- pbsapply(same_lig,function(X) {
  temp <- cor(lvl4_data@mat[,X])
  return(temp[upper.tri(temp)])
},simplify=F)

par(mfrow=c(2,1),mar=c(2,4,1,1),mgp=2:0)
boxplot(sapply(lig16,function(X) rowMeans(lvl4_data@mat[,lvl4_data@cdesc$pert_iname == X])),
        outline=F,xaxt="n",ylab=paste("Mean Z-score of","assayed genes",sep="\n"))
abline(h=0,col="red",lty=2)
par(mar=c(1,4,2,1))
boxplot(same_lig_cor,outline=F,xaxt="n",
        ylab=paste("Pearson correlation of","Z-scores across samples",sep="\n"))
abline(h=0,col="red",lty=2)
mtext(lig16,side=3,at=seq_along(lig16),las=2,xpd=NA,line=2,adj=0.5)

same_ct <- sapply(ct14,function(X) rownames(lvl4_data@cdesc)[lvl4_data@cdesc$cell_id == X])
same_ligct <- sapply(lig16,function(LIG)
  sapply(ct14,function(CT)
    rownames(lvl4_data@cdesc)[lvl4_data@cdesc$cell_id == CT & lvl4_data@cdesc$pert_iname == LIG]
  ),simplify=F)
same_ligct <- unlist(same_ligct,recursive=F)
same_reps <- apply(unique(lvl4_data@cdesc[,c("pert_iname","cell_id","pert_dose","pert_time")]),1,function(X)
  rownames(lvl4_data@cdesc)[lvl4_data@cdesc$pert_iname == X[1] &
                              lvl4_data@cdesc$cell_id == X[2] &
                              lvl4_data@cdesc$pert_dose == X[3] &
                              lvl4_data@cdesc$pert_time == X[4]])

temp1 <- sample(names(same_reps)[sapply(same_reps,length) > 1],1)
temp2 <- sample(same_reps[[temp1]],2)
plot(lvl4_data@mat[,temp2],pch=".")
