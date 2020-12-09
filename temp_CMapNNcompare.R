load("~/Dropbox/GDB_archive/CMapCorr_files/lvl5_inputs_allgenes.RData")
rm(lig16)

Y <- rownames(lvl5_data@cdesc)[lvl5_data@cdesc$pert_iname == "TNF"]
X <- rownames(lvl5_data@mat)[unique(apply(lvl5_data@mat[,Y],2,which.max))]




FDR <- ogFDR[,colnames(ogFDR) %in% names(mean_lfc)]
FDR <- FDR[apply(FDR,1,function(X) any(X <= 0.1)),,drop=F]
# temp_hROW <- hclust(dist(t(FDR)),method="ward.D2") 
temp_hROW <- list(order=order(colnames(FDR),decreasing=T))
temp_hGENE <- hclust(dist(FDR),method="ward.D2")
FDR <- FDR[temp_hGENE$order,temp_hROW$order]
FDR <- -log10(FDR)
FDR[FDR > 2] <- 2
FDR <- FDR / max(FDR)
FDR[FDR == 0] <- NA

ogZS <- sapply(colnames(FDR),function(LIG)
  rowMeans(lvl4_data@mat[rownames(FDR),lvl4_data@cdesc$pert_iname == LIG]))
ZS <- matrix(cut(c(-1,1,as.vector(ogZS / max(abs(ogZS)))),
                 breaks=1000,labels=F)[-(1:2)],nrow=nrow(ogZS))


par(mar=c(5,7,2,1),mgp=2:0)
plot(x=NULL,y=NULL,xlim=c(0.5,nrow(FDR) + .5),ylim=c(0.5,ncol(FDR) + .5),
     xaxs="i",yaxs="i",xaxt="n",yaxt="n",xlab=NA,ylab=NA,bty="n",asp=1)
abline(h=1:ncol(FDR),col="grey90")
symbols(x=rep(1:nrow(FDR),ncol(FDR)),
        y=as.vector(sapply(1:ncol(FDR),function(X) rep(X,nrow(FDR)))),
        circles=as.vector(FDR)/2,inches=F,add=T,xpd=NA,
        fg=diverging_hcl(1000,palette="Blue-Red")[as.vector(ZS)],
        bg=diverging_hcl(1000,palette="Blue-Red")[as.vector(ZS)])
mtext(colnames(FDR),side=2,at=seq_along(colnames(FDR)),las=2,line=0.1,cex=0.8)
mtext("Ligands",side=2,line=2.5)
mtext(lvl4_data@rdesc[rownames(FDR),"pr_gene_symbol"],
      side=1,las=2,at=1:nrow(FDR),line=0.1,cex=0.8)
mtext("Genes differentially expressed",side=1,line=3.5)

temp_yh <- par("usr")[4] - par("usr")[3]
segments(x0=rep(line2user(5,2),1000),x1=rep(line2user(4,2),1000),
         y0=seq(0.6 * temp_yh + par("usr")[3],
                0.9 * temp_yh + par("usr")[3],
                length.out=1000),
         y1=seq(0.6 * temp_yh + par("usr")[3],
                0.9 * temp_yh + par("usr")[3],
                length.out=1000),
         xpd=NA,col=diverging_hcl(1000,palette="Blue-Red"))
mtext(c(signif(max(abs(ogZS)),2),"Z-score",-1 * signif(max(abs(ogZS)),2)),
      side=2,line=c(4,5.1,4),adj=c(-0.1,0.5,1.1),
      at=c(0.9,0.75,0.6) * temp_yh + par("usr")[3])

temp_q <- rev(c(0.5,0.2,0.1,0.05,0.01))
symbols(x=rep(line2user(4.5,2),length(temp_q)),
        y=seq(0.4 * temp_yh + par("usr")[3],
              by=-1.5,length.out=length(temp_q)),
        circles=-log10(temp_q)/4,inches=F,add=T,xpd=NA,
        fg="black")
mtext(c(paste0("\u2264",temp_q[1] * 100,"%"),
        paste0(temp_q[2:5] * 100,"%")),
      side=2,line=5,cex=0.8,
      at=seq(0.4 * temp_yh + par("usr")[3],
             by=-1.5,length.out=length(temp_q)))
mtext("FDR",side=2,line=5.8,
      at=0.4 * temp_yh + par("usr")[3] - 3)

mtext(MAIN,side=3,line=0.5,font=2,cex=1.2)

