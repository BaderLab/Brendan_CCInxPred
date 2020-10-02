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

# load Crow et al. https://doi.org/10.1073/pnas.1802973116
DEprior <- read.table("~/Dropbox/GDB/DE_Prior.txt",header=T,
                      colClasses=c("integer","character","integer","numeric","character"))


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


# trim to genes with <10% probability of change by chance in at least one ligand
meanZlig <- meanZlig_all[apply(meanZlig_all,1,function(X) any(abs(X) > 1.645)),]

# order ligands
temp_hROW <- hclust(dist(t(meanZlig)))

# order genes
temp_hGENE <- hclust(dist(meanZlig))

# order genes by DE prior
# all(rownames(meanZlig) %in% DEprior$Gene_EntrezID)
# meanZlig <- meanZlig[DEprior$Gene_EntrezID[DEprior$Gene_EntrezID %in% rownames(meanZlig)],]

layout(cbind(2:1),heights=c(1,8))
par(mar=c(2.5,4,0,6),mgp=2:0,las=2)
image(z=meanZlig[temp_hGENE$order,temp_hROW$order] / max(abs(meanZlig)),
      x=1:nrow(meanZlig),y=1:ncol(meanZlig),
      col=diverging_hcl(99,palette="Tofino"),
      breaks=seq(-1,1,length.out=100),
      xaxt="n",yaxt="n",xlab=NA,ylab=NA)
mtext(colnames(meanZlig)[temp_hROW$order],side=2,
      at=1:ncol(meanZlig),adj=1.1)
temp_p <- paste0("(p=",pDElig,")")
temp_p[pDElig == 0] <- "(p < 1e-4)"
mtext(paste(countDElig,temp_p)[temp_hROW$order],
      side=4,at=1:ncol(meanZlig),adj=-0.1,cex=0.9)
mtext("# of genes with mean Z > 95th percentile",
      side=4,line=5,las=0,cex=0.9)
mtext(sapply(rownames(meanZlig),function(X) 
  DEprior$Gene_Name[DEprior$Gene_EntrezID == X])[temp_hGENE$order],
  side=1,at=1:nrow(meanZlig),adj=1.1,cex=0.5)

par(mar=c(0,4,0.5,6))
barplot(sapply(rownames(meanZlig),function(X) 
  DEprior$DE_Prior_Rank[DEprior$Gene_EntrezID == X])[temp_hGENE$order],
  ylim=0:1,xaxt="n",xaxs="i",yaxt="n",xlab=NA,ylab=NA,
  col="thistle4",border=NA)
box()
mtext("DE prior",side=4,at=0.5,adj=-0.1)

par(.PAR)
rm(list=grep("^temp",ls(),value=T))




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

# trim to genes with <10% probability of change by chance in at least one ligand
meanZct <- meanZct_all[apply(meanZct_all,1,function(X) any(abs(X) > 1.645)),]

# order ligands
temp_hROW <- hclust(dist(t(meanZct)))

# order genes
temp_hGENE <- hclust(dist(meanZct))

# order genes by DE prior
# all(rownames(meanZct) %in% DEprior$Gene_EntrezID)
# meanZct <- meanZct[DEprior$Gene_EntrezID[DEprior$Gene_EntrezID %in% rownames(meanZct)],]

layout(cbind(2:1),heights=c(1,8))
par(mar=c(2.5,7.5,0,6),mgp=2:0,las=2)
image(z=meanZct[temp_hGENE$order,temp_hROW$order] / max(abs(meanZct)),
      x=1:nrow(meanZct),y=1:ncol(meanZct),
      col=diverging_hcl(99,palette="Tofino"),
      breaks=seq(-1,1,length.out=100),
      xaxt="n",yaxt="n",xlab=NA,ylab=NA)
mtext(colnames(meanZct)[temp_hROW$order],side=2,
      at=1:ncol(meanZct),adj=1.01,cex=0.7)
temp_p <- paste0("(p=",pDEct,")")
temp_p[pDEct == 0] <- "(p < 1e-4)"
mtext(paste(countDEct,temp_p)[temp_hROW$order],
      side=4,at=1:ncol(meanZct),adj=-0.1,cex=0.9)
mtext("# of genes with mean Z > 95th percentile",
      side=4,line=5,las=0,cex=0.9)
mtext(sapply(rownames(meanZct),function(X) 
  DEprior$Gene_Name[DEprior$Gene_EntrezID == X])[temp_hGENE$order],
  side=1,at=1:nrow(meanZct),adj=1.1,cex=0.5)

par(mar=c(0,8,0.5,6))
temp <- sapply(rownames(meanZct),function(X)
  DEprior$DE_Prior_Rank[DEprior$Gene_EntrezID == X])
temp[sapply(temp,length) == 0] <- 0
barplot(unlist(temp)[temp_hGENE$order],
  ylim=0:1,xaxt="n",xaxs="i",yaxt="n",xlab=NA,ylab=NA,
  col="thistle4",border=NA)
box()
mtext("DE prior",side=4,at=0.5,adj=-0.1)


par(.PAR)
rm(list=grep("^temp",ls(),value=T))



# Correlation within ligand tx Z-scores ----
LIG <- "TNF"
zCor <- cor(lvl4_data@mat[,lvl4_data@cdesc$pert_iname == LIG],method="spearman")
hCOR <- hclust(dist(zCor))

layout(cbind(2:1),heights=c(1,8))
par(mar=c(1,1,0,1),mgp=2:0)
image(zCor[hCOR$order,hCOR$order],
      xaxt="n",yaxt="n",xlab=NA,ylab=NA)
par(mar=c(0,1,1,1))
barplot(rep(1,ncol(zCor)),space=0,border=NA,
        ylim=0:1,xaxt="n",xaxs="i",yaxt="n",xlab=NA,ylab=NA,
        col=rainbow(14)[as.factor(lvl4_data@cdesc[colnames(zCor),"cell_id"])[hCOR$order]])




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

# order ligands
temp_hLIG <- hclust(dist(scoreDEligct))

# order cell types
temp_hCT <- hclust(dist(t(scoreDEligct)))


par(mar=c(1,7.5,2.5,1),mgp=2:0,las=2)
image(scoreDEligct[temp_hLIG$order,temp_hCT$order],
      col=sequential_hcl(100,palette="inferno",rev=T),
      x=1:nrow(scoreDEligct),y=1:ncol(scoreDEligct),
      xaxt="n",yaxt="n",xlab=NA,ylab=NA)
mtext(colnames(scoreDEligct)[temp_hCT$order],side=2,
      at=1:ncol(scoreDEligct),adj=1.03,cex=0.7)
mtext(rownames(scoreDEligct)[temp_hLIG$order],side=3,
      at=1:nrow(scoreDEligct),adj=-0.1,cex=0.7)


mtext(colnames(meanZligct)[temp_hROW$order],side=2,
      at=1:ncol(meanZligct),adj=1.01,cex=0.7)
mtext(paste0(apply(meanZligct,2,function(X) sum(X > 1.645))," (",
             sapply(ct14,function(X) sum(lvl4_data@cdesc$cell_id == X)),
             ")")[temp_hROW$order],
      side=4,at=1:ncol(meanZligct),adj=-0.1,cex=0.9)
mtext("# of genes with mean Z scores in top 5% (# of samples)",
      side=4,line=3,las=0,cex=0.9)
mtext(sapply(rownames(meanZligct),function(X) 
  DEprior$Gene_Name[DEprior$Gene_EntrezID == X])[temp_hGENE$order],
  side=1,at=1:nrow(meanZligct),adj=1.1,cex=0.5)

par(mar=c(0,8,0.5,4))
temp <- sapply(rownames(meanZligct),function(X)
  DEprior$DE_Prior_Rank[DEprior$Gene_EntrezID == X])
temp[sapply(temp,length) == 0] <- 0
barplot(unlist(temp)[temp_hGENE$order],
        ylim=0:1,xaxt="n",xaxs="i",yaxt="n",xlab=NA,ylab=NA,
        col="thistle4",border=NA)
box()
mtext("DE prior",side=2,at=0.5,adj=1.1)


par(.PAR)
rm(list=grep("^temp",ls(),value=T))



# Compare Z-score distributions ----
temp_min <- min(c(meanZlig_all,meanZct_all,meanZligct_all))
temp_max <- max(c(meanZlig_all,meanZct_all,meanZligct_all))

par(mfrow=c(3,1),mar=c(3,3,1,1),mgp=2:0)
hist(as.vector(meanZlig_all),breaks=50,main=NA,
     xlab="Mean gene Z-score within ligand treatment",
     xlim=c(temp_min,temp_max),freq=F)
points(seq(temp_min,temp_max,length.out=1000),
       dnorm(seq(temp_min,temp_max,length.out=1000)),
       type="l")
hist(as.vector(meanZct_all),breaks=50,main=NA,
     xlab="Mean gene Z-score within cell type",
     xlim=c(temp_min,temp_max),freq=F)
points(seq(temp_min,temp_max,length.out=1000),
       dnorm(seq(temp_min,temp_max,length.out=1000)),
       type="l")
hist(as.vector(meanZligct_all),breaks=500,main=NA,
     xlab="Mean gene Z-score within ligand treatment per cell type",
     xlim=c(temp_min,temp_max),freq=F)
points(seq(temp_min,temp_max,length.out=1000),
       dnorm(seq(temp_min,temp_max,length.out=1000)),
       type="l")


