library(cmapR)
library(scales)
library(colorspace)
library(pbapply)
.PAR <- par(no.readonly=T)


if (exists("lvl4_data")) {
} else if (file.exists("~/Dropbox/GDB_archive/CMapCorr_files/lvl4_lig16_inputs.RData")) {
  load("~/Dropbox/GDB_archive/CMapCorr_files/lvl4_lig16_inputs.RData") 
} else {
  source("lvl4_lig16_inputs.R")
}
lig16 <- sort(lig16)

# LIG ----
lig_id <- sapply(lig16,function(X) rownames(lvl4_data@cdesc)[lvl4_data@cdesc$pert_iname == X])
lig_mean <- sapply(lig16,function(X) rowMeans(lvl4_data@mat[,lvl4_data@cdesc$pert_iname == X]))
lig_cor <- pbsapply(lig_id,function(X) {
  temp <- cor(lvl4_data@mat[,X])
  return(temp[upper.tri(temp)])
},simplify=F)
lig_col <- qualitative_hcl(16,palette="dark3")
names(lig_col) <- lig16

png("~/Dropbox/GDB/CMapCorr/Fig_corr_lig.png",
    width=5,height=5,units="in",res=300)
par(mfrow=c(2,1),mar=c(2,4,1,1),mgp=2:0)
boxplot(lig_mean,outline=T,pch=".",xaxt="n",
        border=lig_col,col=alpha(lig_col,0.5),
        ylab=paste("Mean Z-score of","assayed genes",sep="\n"))
abline(h=0,col="red",lty=2)
par(mar=c(1,4,2,1))
boxplot(lig_cor,outline=F,xaxt="n",
        border=lig_col,col=alpha(lig_col,0.5),
        ylab=paste("Pearson correlation of","Z-scores across samples",sep="\n"))
abline(h=0,col="red",lty=2)
mtext(lig16,side=3,at=seq_along(lig16),las=2,xpd=NA,line=2,adj=0.5,col=lig_col)
dev.off()
par(.PAR)



# CT ----
ct_id <- sapply(ct14,function(X) rownames(lvl4_data@cdesc)[lvl4_data@cdesc$cell_id == X])
ct_mean <- sapply(ct14,function(X) rowMeans(lvl4_data@mat[,lvl4_data@cdesc$cell_id == X]))
ct_cor <- pbsapply(ct_id,function(X) {
  temp <- cor(lvl4_data@mat[,X])
  return(temp[upper.tri(temp)])
},simplify=F)
ct_col <- qualitative_hcl(14,palette="dark3")

png("~/Dropbox/GDB/CMapCorr/Fig_corr_ct.png",
    width=5,height=5,units="in",res=300)
par(mfrow=c(2,1),mar=c(2,4,1,1),mgp=2:0)
boxplot(ct_mean,outline=T,pch=".",xaxt="n",
        border=ct_col,col=alpha(ct_col,0.5),
        ylab=paste("Mean Z-score of","assayed genes",sep="\n"))
abline(h=0,col="red",lty=2)
par(mar=c(1,4,2,1))
boxplot(ct_cor,outline=F,xaxt="n",
        border=ct_col,col=alpha(ct_col,0.5),
        ylab=paste("Pearson correlation of","Z-scores across samples",sep="\n"))
abline(h=0,col="red",lty=2)
mtext(paste(mapply(sub,pattern=paste0("_",ct14),x=names(ct14),replacement=""),ct14,sep="\n"),
      side=3,at=seq_along(ct14),las=2,xpd=NA,line=2,adj=0.5,cex=0.6,col=ct_col)
dev.off()
par(.PAR)



# LIG CT ----
ligct_id <- sapply(lig16,function(LIG)
  sapply(ct14,function(CT)
    rownames(lvl4_data@cdesc)[lvl4_data@cdesc$cell_id == CT & lvl4_data@cdesc$pert_iname == LIG]
  ),simplify=F)
ligct_id <- unlist(ligct_id,recursive=F)
ligct_mean <- sapply(ligct_id,function(X) rowMeans(lvl4_data@mat[,X]))
ligct_cor <- pbsapply(ligct_id,function(X) {
  temp <- cor(lvl4_data@mat[,X])
  return(temp[upper.tri(temp)])
},simplify=F)
ligct_col <- as.vector(sapply(lig_col,rep,times=14))

png("~/Dropbox/GDB/CMapCorr/Fig_corr_ligct.png",
    width=5,height=5,units="in",res=300)
par(mfrow=c(2,1),mar=c(2,4,1,1),mgp=2:0)
boxplot(ligct_mean,outline=T,pch=".",xaxt="n",
        border=ligct_col,col=alpha(ligct_col,0.5),
        ylab=paste("Mean Z-score of","assayed genes",sep="\n"))
abline(h=0,col="red",lty=2)
par(mar=c(1,4,2,1))
boxplot(ligct_cor,outline=F,xaxt="n",
        border=ligct_col,col=alpha(ligct_col,0.5),
        ylab=paste("Pearson correlation of","Z-scores across samples",sep="\n"))
abline(h=0,col="red",lty=2)
mtext(lig16,side=3,at=seq(7.5,by=14,length.out=16),
      las=2,xpd=NA,line=2,adj=0.5,col=lig_col)
dev.off()
par(.PAR)



# REP ----
rep_tx <- unique(lvl4_data@cdesc[,c("pert_iname","cell_id","pert_dose","pert_time")])
rownames(rep_tx) <- apply(rep_tx,1,paste,collapse="_")
rep_id <- apply(rep_tx,1,function(X)
  rownames(lvl4_data@cdesc)[lvl4_data@cdesc$pert_iname == X[1] &
                              lvl4_data@cdesc$cell_id == X[2] &
                              lvl4_data@cdesc$pert_dose == X[3] &
                              lvl4_data@cdesc$pert_time == X[4]])
rep_id <- rep_id[sort(names(rep_id))]
rep_id <- rep_id[sapply(rep_id,length) > 1]
rep_cor <- pbsapply(rep_id,function(X) {
  temp <- cor(lvl4_data@mat[,X])
  return(temp[upper.tri(temp)])
},simplify=F)
rep_col <- lig_col[sub("_.+$","",names(rep_id))]

png("~/Dropbox/GDB/CMapCorr/Fig_corr_rep.png",
    width=5,height=5,units="in",res=300)
par(mar=c(1,4,4,1),mgp=2:0)
boxplot(rep_cor,outline=F,xaxt="n",
        border=rep_col,col=alpha(rep_col,0.5),
        ylab=paste("Pearson correlation of","Z-scores across samples",sep="\n"))
abline(h=0,col="red",lty=2)
mtext(lig16,side=3,las=2,adj=-0.1,col=lig_col,
      at=sapply(lig16,function(X) mean(which(sub("_.+$","",names(rep_id)) == X))))
dev.off()
par(.PAR)


# Compare all ----
png("~/Dropbox/GDB/CMapCorr/Fig_corr_compareall.png",
    width=5,height=8,units="in",res=300)
par(mfrow=c(4,1),mar=c(4,3,1,1),mgp=2:0)
hist(unlist(ct_cor),xlim=c(-1,1),breaks=100,main=NA,
     xlab="PCC b/w Z-scores per cell line")
abline(v=median(unlist(ct_cor)),col="red")
legend("topleft",legend="Median PCC",lty=1,col="red",bty="n")

hist(unlist(lig_cor),xlim=c(-1,1),breaks=50,main=NA,
     xlab="PCC b/w Z-scores per ligand")
abline(v=median(unlist(lig_cor)),col="red")
legend("topleft",lty=c(1,NA),col=c("red",NA),bty="n",
       legend=c("Median PCC",paste("Difference =",
                                   signif(median(unlist(lig_cor)) - 
                                            median(unlist(ct_cor)),2))))

hist(unlist(ligct_cor),xlim=c(-1,1),breaks=50,main=NA,
     xlab="PCC b/w Z-scores per ligand / cell line")
abline(v=median(unlist(ligct_cor)),col="red")
legend("topleft",lty=c(1,NA),col=c("red",NA),bty="n",
       legend=c("Median PCC",paste("Difference =",
                                   signif(median(unlist(ligct_cor)) - 
                                            median(unlist(lig_cor)),2))))

hist(unlist(rep_cor),xlim=c(-1,1),breaks=50,main=NA,
     xlab="PCC b/w Z-scores per replicate")
abline(v=median(unlist(rep_cor)),col="red")
legend("topleft",lty=c(1,NA),col=c("red",NA),bty="n",
       legend=c("Median PCC",paste("Difference =",
                                   signif(median(unlist(rep_cor)) - 
                                            median(unlist(ligct_cor)),2))))
dev.off()
par(.PAR)



# Compare all per ligand----
AllCorrByLig <- pbsapply(lig16,function(X) {
  list(
    ByLig=lig_cor[[X]],
    ByLigCT=unlist(ligct_cor[sub("\\..+$","",names(ligct_cor)) == X]),
    ByRep=unlist(rep_cor[sub("_.+$","",names(rep_cor)) == X])
  )
},simplify=F)


png("~/Dropbox/GDB/CMapCorr/Fig_corr_compareallperligand.png",
    width=6,height=4,units="in",res=300)
par(mar=c(6,3,4,1),mgp=2:0)
plot(NA,NA,xlim=c(0,16*3+1),ylim=c(-1,1),xaxt="n",xaxs="i",yaxs="i",
     xlab=NA,ylab="Pearson correlation of Z-scores")
abline(h=0,lty=2,col="grey30")
boxplot(unlist(AllCorrByLig,recursive=F),pch=".",xaxt="n",add=T,
        # ylab="Pearson correlation of Z-scores",
        col=alpha(lig_col[as.vector(sapply(lig16,rep,3))],0.5),
        border=lig_col[as.vector(sapply(lig16,rep,3))])
mtext(lig16,las=2,side=3,at=seq(2,by=3,length.out=16),line=0.1,col=lig_col)
mtext(rep(c("By ligand","By ligand / cell line","By treatment condition")),
      las=2,side=1,at=seq(1,16*3),line=0.1,adj=1,cex=0.6)
dev.off()
par(.PAR)

