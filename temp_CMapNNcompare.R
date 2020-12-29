library(cmapR)
library(colorspace)
library(scales)
source("~/Dropbox/GDB/line2user.R")

LIG <- "TNF"


# CMap genes ----
load("~/Dropbox/GDB_archive/CMapCorr_files/lvl4_inputs_allgenes.RData")
CM <- lvl4_data_all
rm(lig16,lvl4_data_all)

CM_ogZS <- sapply(ct14,function(CT) 
  rowMeans(CM@mat[,CM@cdesc$pert_iname == LIG & CM@cdesc$cell_id == CT]))
CM_ogFDR <- apply(CM_ogZS,2,function(Z) p.adjust(pnorm(-abs(Z))))

CM_GENES <- unique(c(
  rownames(CM_ogFDR)[apply(CM_ogFDR,1,function(X) sum(X <= 0.1) >= 2)]
))


# NN genes ----
load("~/Dropbox/GDB_archive/CMapCorr_files/NN_all_dat.RData")
# NN_ogZS <- nn_lig_rep[[LIG]]$lfc
NN_ogZS <- nn_lig_rep[[LIG]]$lfc * -1 # because NN flipped the sign on their log-ratios.
NN_ogZS <- NN_ogZS[,!DSinfo[[LIG]][colnames(NN_ogZS),"time_series"]]
NN_ogFDR <- nn_lig_rep[[LIG]]$qval[,colnames(NN_ogZS)]
colnames(NN_ogFDR) <- colnames(NN_ogZS) <- paste(DSinfo[[LIG]][colnames(NN_ogZS),"CtAcc"],"*")

NN_GENES <- rownames(NN_ogZS)[apply(sapply(colnames(NN_ogZS),function(X)
  abs(NN_ogZS[,X]) >= 1 & NN_ogFDR[,X] <= 0.1),
  1,function(Y) sum(Y) >= 4)]
NN_GENES[!NN_GENES %in% CM@rdesc$pr_gene_symbol]
# "SELENOP" was once "SEPP1"


# Merge genes ----
temp_CMgenes <- GENES <- unique(c(NN_GENES,CM@rdesc[CM_GENES,"pr_gene_symbol"]))
print(temp_CMgenes[!temp_CMgenes %in% CM@rdesc$pr_gene_symbol])
temp_CMgenes[temp_CMgenes == "SELENOP"] <- "SEPP1"
temp_CMgenes <- sapply(temp_CMgenes,function(X) rownames(CM@rdesc)[CM@rdesc$pr_gene_symbol == X])
names(temp_CMgenes)[names(temp_CMgenes) == "SEPP1"] <- "SELENOP"

# Filter CMap ----
CM_ogZS <- CM_ogZS[temp_CMgenes,]
CM_ogFDR <- CM_ogFDR[temp_CMgenes,]
rownames(CM_ogZS) <- rownames(CM_ogFDR) <- GENES


# Filter NN ----
NN_ogZS <- NN_ogZS[GENES,]
NN_ogFDR <- NN_ogFDR[GENES,]


# Merge data ----
ogZS <- cbind(CM_ogZS,NN_ogZS)
ogFDR <- cbind(CM_ogFDR,NN_ogFDR)


# HClust ----
hS <- hclust(as.dist(1 - cor(ogZS,method="spearman")))
hG <- hclust(as.dist(1 - cor(t(ogZS),method="spearman")))
ogZS <- ogZS[hG$order,hS$order]
ogFDR <- ogFDR[hG$order,hS$order]



# Scale CMap ----
CM_cut <- 5
boxplot(CM_ogZS,ylim=c(-CM_cut,CM_cut))
CM_ZS <- CM_ogZS
CM_ZS[CM_ZS > 5] <- CM_cut; CM_ZS[CM_ZS < -5] <- -CM_cut
CM_ZS <- matrix(cut(c(-1,1,as.vector(CM_ogZS / max(abs(CM_ogZS)))),
                 breaks=1000,labels=F)[-(1:2)],nrow=nrow(CM_ogZS))
rownames(CM_ZS) <- rownames(CM_ogZS); colnames(CM_ZS) <- colnames(CM_ogZS)

CM_FDR <- -log10(CM_ogFDR)
CM_FDR[CM_FDR > 3] <- 3
CM_FDR <- CM_FDR / max(CM_FDR)
CM_FDR[CM_FDR <= 0] <- NA


# Scale NN ----
NN_cut <- 5
boxplot(NN_ogZS,ylim=c(-NN_cut,NN_cut))
NN_ZS <- NN_ogZS
NN_ZS[NN_ZS > 3] <- NN_cut; NN_ZS[NN_ZS < -3] <- -NN_cut
NN_ZS <- matrix(cut(c(-1,1,as.vector(NN_ogZS / max(abs(NN_ogZS)))),
                    breaks=1000,labels=F)[-(1:2)],nrow=nrow(NN_ogZS))
rownames(NN_ZS) <- rownames(NN_ogZS); colnames(NN_ZS) <- colnames(NN_ogZS)

NN_FDR <- -log10(NN_ogFDR)
NN_FDR[NN_FDR > 3] <- 3
NN_FDR <- NN_FDR / max(NN_FDR)
NN_FDR[NN_FDR <= 0] <- NA


# Merge plotting data ----
ZS <- t(cbind(CM_ZS,NN_ZS)[rownames(ogZS),colnames(ogZS)])
FDR <- t(cbind(CM_FDR,NN_FDR)[rownames(ogFDR),colnames(ogFDR)])


# PLOT ----
# png(paste0("docs/output_figs/CMapNN_",LIG,"_wrong.png"),
# png(paste0("docs/output_figs/CMapNN_",LIG,"_fixed.png"),
#     width=nrow(FDR) * 0.2 + 1.8,height=ncol(FDR) * 0.2 + 1.6,
#     units="in",res=120)

par(mar=c(7,8,1,1),mgp=2:0)
plot(x=NULL,y=NULL,xlim=c(0.5,nrow(FDR) + .5),ylim=c(0.5,ncol(FDR) + .5),
     xaxs="i",yaxs="i",xaxt="n",yaxt="n",xlab=NA,ylab=NA,bty="n",asp=1)
# abline(h=1:ncol(FDR),col="grey90")
abline(v=1:nrow(FDR),col="grey90")
symbols(x=rep(1:nrow(FDR),ncol(FDR)),
        y=as.vector(sapply(1:ncol(FDR),function(X) rep(X,nrow(FDR)))),
        circles=as.vector(FDR)/2,inches=F,add=T,xpd=NA,
        fg=diverging_hcl(1000,palette="Blue-Red3")[as.vector(ZS)],
        bg=diverging_hcl(1000,palette="Blue-Red3")[as.vector(ZS)])
mtext(ifelse(as.logical(as.integer(CM@rdesc[temp_CMgenes[colnames(FDR)],"pr_is_lm"])),
             yes=colnames(FDR),paste0(colnames(FDR),"^")),
      side=2,at=1:ncol(FDR),las=2,line=0.1,cex=0.8)
text(x=1:nrow(FDR),y=rep(line2user(0.3,1),nrow(FDR)),
     labels=rownames(FDR),xpd=NA,adj=1,srt=45,cex=0.8)

temp_yh <- par("usr")[4] - par("usr")[3]
segments(x0=rep(line2user(6,2),1000),x1=rep(line2user(5,2),1000),
         y0=seq(0.6 * temp_yh + par("usr")[3],
                0.9 * temp_yh + par("usr")[3],
                length.out=1000),
         y1=seq(0.6 * temp_yh + par("usr")[3],
                0.9 * temp_yh + par("usr")[3],
                length.out=1000),
         xpd=NA,col=diverging_hcl(1000,palette="Blue-Red3"))
mtext(c(paste0("\u2265",CM_cut),"Z-score",paste0("\u2264",-CM_cut)),
      side=2,line=6.1,adj=c(-0.1,0.5,1.1),
      at=c(0.9,0.75,0.6) * temp_yh + par("usr")[3])
mtext(c(paste0("\u2265",NN_cut),"*LogFC",paste0("\u2264",-NN_cut)),
      side=2,line=3.9,adj=c(-0.1,0.5,1.1),
      at=c(0.9,0.75,0.6) * temp_yh + par("usr")[3])


temp_q <- rev(c(0.2,0.1,0.05,0.01,0.001))
symbols(x=rep(line2user(5.5,2),length(temp_q)),
        y=seq(0.4 * temp_yh + par("usr")[3],
              by=-2,length.out=length(temp_q)),
        circles=-log10(temp_q)/4,inches=F,add=T,xpd=NA,
        fg="black")
mtext(c(paste0("\u2264",temp_q[1] * 100,"%"),
        paste0(temp_q[2:5] * 100,"%")),
      side=2,line=6.2,cex=0.8,
      at=seq(0.4 * temp_yh + par("usr")[3],
             by=-2,length.out=length(temp_q)))
mtext("FDR",side=2,line=7,
      at=0.4 * temp_yh + par("usr")[3] - 2 * 2)

mtext(LIG,side=3,line=-1,at=line2user(5.5,2),font=2,cex=1.2)
text(line2user(7,2),line2user(6.5,1),xpd=NA,adj=c(0,0.5),srt=90,cex=0.8,
     labels=paste("* microarray data","^ inferred in CMap",sep="\n"))
# text(line2user(1,4),line2user(6,1),xpd=NA,adj=c(1,0.5),srt=0,cex=0.8,
#      labels=paste("* microarray data","^ inferred in CMap",sep="\n"))

# dev.off()
