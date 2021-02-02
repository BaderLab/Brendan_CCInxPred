LIG <- "TNF"  

if (all(DSinfo[[LIG]]$time_series)) { next }
if (!LIG %in% lvl4_data@cdesc$pert_iname) { next }

# CMap genes ----
CM_ogZS <- sapply(unique(lvl4_data@cdesc$cell_id[lvl4_data@cdesc$pert_iname == LIG]),function(CT) 
  rowMeans(lvl4_data@mat[,lvl4_data@cdesc$pert_iname == LIG & lvl4_data@cdesc$cell_id == CT]))
CM_ogFDR <- apply(CM_ogZS,2,function(Z) p.adjust(pnorm(-abs(Z))))

CM_GENES <- unique(c(
  rownames(CM_ogFDR)[apply(CM_ogFDR,1,function(X) sum(X <= 0.1) >= 2)]
))


# NN genes ----
NN_ogZS <- nn_lig_rep[[LIG]]$lfc
# NN_ogZS <- NN_ogZS[,!DSinfo[[LIG]][colnames(NN_ogZS),"time_series"],drop=F]
# CAN'T DROP TIMESERIES
NN_ogFDR <- nn_lig_rep[[LIG]]$qval[,colnames(NN_ogZS),drop=F]
temp_ts <- DSinfo[[LIG]][colnames(NN_ogZS),"time_series"]
colnames(NN_ogFDR) <- colnames(NN_ogZS) <- paste(DSinfo[[LIG]][colnames(NN_ogZS),"CtAcc"],"*")
colnames(NN_ogFDR)[temp_ts] <- colnames(NN_ogZS)[temp_ts] <- paste0(colnames(NN_ogZS)[temp_ts],"*")
if (LIG == "VEGFA") {
  NN_GENES <- character()
} else {  
  temp_samples <- ncol(NN_ogFDR) + 1
  NN_GENES <- character()
  while (length(NN_GENES) < 1) {
    temp_samples <- temp_samples - 1
    if (temp_samples <= 1) { break }
    NN_GENES <- rownames(NN_ogZS)[apply(sapply(colnames(NN_ogZS),function(X)
      abs(NN_ogZS[,X]) >= 1 & NN_ogFDR[,X] <= 0.1),
      1,function(Y) sum(Y) >= temp_samples)]  
  }
}
if (any(any(!NN_GENES %in% lvl4_data@rdesc$pr_gene_symbol))) {
  # print(NN_GENES[!NN_GENES %in% lvl4_data@rdesc$pr_gene_symbol])
}


# Merge genes & fix gene names ----
temp_CMgenes <- GENES <- unique(c(NN_GENES,lvl4_data@rdesc[CM_GENES,"pr_gene_symbol"]))
# Fixing SELENOP (formerly SEPP1)
temp_CMgenes[temp_CMgenes == "SELENOP"] <- "SEPP1"
# Fixing missing genes
temp_CMgenes <- temp_CMgenes[temp_CMgenes %in% lvl4_data@rdesc$pr_gene_symbol]
GENES <- GENES[GENES %in% lvl4_data@rdesc$pr_gene_symbol]
temp_CMgenes <- temp_CMgenes[temp_CMgenes %in% rownames(NN_ogZS)]
GENES <- GENES[GENES %in% rownames(NN_ogZS)]

temp_CMgenes <- sapply(temp_CMgenes,function(X) 
  rownames(lvl4_data@rdesc)[lvl4_data@rdesc$pr_gene_symbol == X])
names(temp_CMgenes)[names(temp_CMgenes) == "SEPP1"] <- "SELENOP"


# Filter CMap ----
CM_ogZS <- CM_ogZS[temp_CMgenes,,drop=F]
CM_ogFDR <- CM_ogFDR[temp_CMgenes,,drop=F]
rownames(CM_ogZS) <- rownames(CM_ogFDR) <- GENES


# Filter NN ----
NN_ogZS <- NN_ogZS[GENES,,drop=F]
NN_ogFDR <- NN_ogFDR[GENES,,drop=F]


# Merge data ----
ogZS <- cbind(CM_ogZS,NN_ogZS)
ogFDR <- cbind(CM_ogFDR,NN_ogFDR)

# HClust ----
if (nrow(ogZS) > 1) {
  temp_mat <- -log10(ogFDR)
  temp_mat[temp_mat == 0] <- 1e-10
  temp_mat <- temp_mat * ogZS
  temp_corS <- cor(temp_mat,method="pearson")
  temp_corG <- cor(t(temp_mat),method="pearson")
  hS <- hclust(as.dist((1 - abs(temp_corS)) * sign(temp_corS)),method="ward.D2")
  hG <- hclust(as.dist((1 - abs(temp_corG)) * sign(temp_corG)),method="ward.D2")
  ogZS <- ogZS[hG$order,hS$order]
  ogFDR <- ogFDR[hG$order,hS$order]
}

# Scale CMap ----
CM_cut <- 5
# boxplot(CM_ogZS,ylim=c(-CM_cut,CM_cut))
CM_ZS <- CM_ogZS
CM_ZS[CM_ZS > CM_cut] <- CM_cut; CM_ZS[CM_ZS < -CM_cut] <- -CM_cut
CM_ZS <- matrix(cut(c(-1,1,as.vector(CM_ogZS / max(abs(CM_ogZS)))),
                    breaks=1000,labels=F)[-(1:2)],nrow=nrow(CM_ogZS))
rownames(CM_ZS) <- rownames(CM_ogZS); colnames(CM_ZS) <- colnames(CM_ogZS)

CM_FDR <- -log10(CM_ogFDR)
CM_FDR[CM_FDR > 3] <- 3
CM_FDR <- CM_FDR / max(CM_FDR)
CM_FDR[CM_FDR <= 0] <- NA


# Scale NN ----
NN_cut <- CM_cut # if you want these to be different, change the legend.
# boxplot(NN_ogZS,ylim=c(-NN_cut,NN_cut))
NN_ZS <- NN_ogZS
NN_ZS[NN_ZS > NN_cut] <- NN_cut; NN_ZS[NN_ZS < -NN_cut] <- -NN_cut
NN_ZS <- matrix(cut(c(-1,1,as.vector(NN_ogZS / max(abs(NN_ogZS)))),
                    breaks=1000,labels=F)[-(1:2)],nrow=nrow(NN_ogZS))
rownames(NN_ZS) <- rownames(NN_ogZS); colnames(NN_ZS) <- colnames(NN_ogZS)

NN_FDR <- -log10(NN_ogFDR)
NN_FDR[NN_FDR > 3] <- 3
NN_FDR <- NN_FDR / max(NN_FDR)
NN_FDR[NN_FDR <= 0] <- NA


# Merge plotting data ----
ZS <- cbind(CM_ZS,NN_ZS)[rownames(ogZS),colnames(ogZS)]
FDR <- cbind(CM_FDR,NN_FDR)[rownames(ogFDR),colnames(ogFDR)]


# PLOT ----
png(paste0("docs/output_figs/CMapNN_",LIG,".png"),
    width=nrow(FDR) * 0.2 + 2.45,height=ncol(FDR) * 0.2 + 1.2,
    units="in",res=120)

par(mar=c(4.5,11,1,1),mgp=2:0)
plot(x=NULL,y=NULL,xlim=c(0.5,nrow(FDR) + .5),ylim=c(0.5,ncol(FDR) + .5),
     xaxs="i",yaxs="i",xaxt="n",yaxt="n",xlab=NA,ylab=NA,bty="n",asp=1)
# abline(h=1:ncol(FDR),col="grey90")
abline(v=1:nrow(FDR),col="grey90")
abline(h=1:ncol(FDR),col="grey90")
symbols(x=rep(1:nrow(FDR),ncol(FDR)),
        y=as.vector(sapply(1:ncol(FDR),function(X) rep(X,nrow(FDR)))),
        circles=as.vector(FDR)/2,inches=F,add=T,xpd=NA,
        fg=diverging_hcl(1000,palette="Blue-Red3")[as.vector(ZS)],
        bg=diverging_hcl(1000,palette="Blue-Red3")[as.vector(ZS)])
text(1:nrow(FDR),rep(line2user(0.3,1),nrow(FDR)),
     ifelse(as.logical(as.integer(lvl4_data@rdesc[temp_CMgenes[rownames(FDR)],"pr_is_lm"])),
            yes=rownames(FDR),paste0(rownames(FDR),"^")),
     xpd=NA,adj=1,srt=45,cex=0.8)
text(x=rep(line2user(0.3,2),ncol(FDR)),y=1:ncol(FDR),
     labels=colnames(FDR),xpd=NA,adj=1,srt=45,cex=0.8)

mtext(LIG,side=3,line=0,at=line2user(8.4,2),font=2,cex=1.2)

temp_yh <- par("usr")[4] - par("usr")[3]
mtext("Z-score /",side=2,line=9.9,at=0.75 * temp_yh + par("usr")[3])
mtext("LogFC *",side=2,line=9.1,at=0.75 * temp_yh + par("usr")[3])
segments(x0=rep(line2user(9,2),1000),x1=rep(line2user(7.8,2),1000),
         y0=seq(0.5 * temp_yh + par("usr")[3],
                temp_yh + par("usr")[3],
                length.out=1000),
         y1=seq(0.5 * temp_yh + par("usr")[3],
                temp_yh + par("usr")[3],
                length.out=1000),
         xpd=NA,col=diverging_hcl(1000,palette="Blue-Red3"))
text(rep(line2user(8.4,2),2),c(1,0.5) * temp_yh + par("usr")[3],
     labels=c(paste0("\u2265",CM_cut),paste0("\u2264",-CM_cut)),
     col="grey80",pos=c(1,3),offset=0.1,xpd=NA)

mtext("FDR",side=2,line=9,las=2,
      at=1 + 0.35 * temp_yh + par("usr")[3])
temp_q <- rev(c(0.2,0.1,0.05,0.01,0.001))
symbols(x=rep(line2user(8.4,2),length(temp_q)),
        y=seq(0.35 * temp_yh + par("usr")[3],
              by=-1,length.out=length(temp_q)),
        circles=-log10(temp_q) / (max(-log10(temp_q)) *2),
        inches=F,add=T,xpd=NA,fg="black")
mtext(c(paste0("\u2264",temp_q[1] * 100,"%"),
        paste0(temp_q[2:5] * 100,"%")),
      side=2,line=9,cex=0.8,las=2,
      at=seq(0.35 * temp_yh + par("usr")[3],
             by=-1,length.out=length(temp_q)))

mtext(paste("* microarray data",
            "** timeseries data",
            "^ inferred in CMap",
            sep="; "),
      side=1,line=3.5,cex=0.8,at=line2user(0.5,4),adj=1)

dev.off()
