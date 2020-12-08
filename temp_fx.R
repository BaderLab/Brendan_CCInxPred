temp_DS <- rownames(DSinfo[[LIG]])
temp_DS <- temp_DS[!temp_DS %in% exclude]

temp_ord_ds <- hclust(dist(t(nn_lig_rep[[LIG]]$lfc[GENES,temp_DS])),"ward.D2")$order
temp_ord_gene <- hclust(dist(nn_lig_rep[[LIG]]$lfc[GENES,temp_DS]),"ward.D2")$order

CTcol <- qualitative_hcl(length(unique(DSinfo[[LIG]]$cell_type)),palette="dark3")
names(CTcol) <- unique(DSinfo[[LIG]]$cell_type)

ogLFC <- nn_lig_rep[[LIG]]$lfc[GENES,temp_DS]
ogLFC <- ogLFC[temp_ord_gene,temp_ord_ds]
temp_dsname <- DSinfo[[LIG]][colnames(ogLFC),"CtAcc"]
temp_dsname[DSinfo[[LIG]][colnames(ogLFC),"time_series"]] <-
  paste(temp_dsname[DSinfo[[LIG]][colnames(ogLFC),"time_series"]],"*")
temp_labelcol <- CTcol[DSinfo[[LIG]][colnames(ogLFC),"cell_type"]]
names(temp_labelcol) <- colnames(ogLFC) <- temp_dsname
if (!is.null(mean_lfc[[LIG]])) {
  LFC <- cbind(ogLFC,rep(0,nrow(ogLFC)))
  temp_labelcol <- append(temp_labelcol,"black")
  if (!is.null(mean_lfc[[LIG]]$nts_lfc)) {
    LFC <- cbind(LFC,"Mean (no *)"=mean_lfc[[LIG]][GENES[temp_ord_gene],"nts_lfc"])
    temp_labelcol <- append(temp_labelcol,"black")
  }
  if (!is.null(mean_lfc[[LIG]]$ts_lfc)) {
    LFC <- cbind(LFC,"Mean (only *)"=mean_lfc[[LIG]][GENES[temp_ord_gene],"ts_lfc"])
    temp_labelcol <- append(temp_labelcol,"black")
  }
  if (!is.null(mean_lfc_CT[[LIG]])){
    LFC <- cbind(LFC,sapply(mean_lfc_CT[[LIG]],function(X) X[GENES[temp_ord_gene],"lfc"]))
    temp_labelcol <- append(temp_labelcol,CTcol[names(mean_lfc_CT[[LIG]])])
  }
  if (!is.null(mean_lfc_DS[[LIG]])){
    LFC <- cbind(LFC,sapply(mean_lfc_DS[[LIG]],function(X) X[GENES[temp_ord_gene],"lfc"]))
    temp_labelcol <- append(temp_labelcol,temp_labelcol[names(mean_lfc_DS[[LIG]])])
  }
} else {
  LFC <- ogLFC
}
LFC <- t(matrix(
  cut(c(-1,1,as.vector(LFC / max(LFC))),
      breaks=1000,labels=F)[-1*(1:2)],
  nrow(LFC)))


ogFDR <- nn_lig_rep[[LIG]]$qval[GENES,temp_DS]
ogFDR <- ogFDR[temp_ord_gene,temp_ord_ds]
colnames(ogFDR) <- temp_dsname
if (!is.null(mean_lfc[[LIG]])) {
  FDR <- cbind(ogFDR,rep(1,nrow(ogFDR)))
  if (!is.null(mean_lfc[[LIG]]$nts_lfc)) {
    FDR <- cbind(FDR,"Mean (no *)"=mean_lfc[[LIG]][GENES[temp_ord_gene],"nts_fdr"])
  }
  if (!is.null(mean_lfc[[LIG]]$ts_lfc)) {
    FDR <- cbind(FDR,"Mean (only *)"=mean_lfc[[LIG]][GENES[temp_ord_gene],"ts_fdr"])
  }
  if (!is.null(mean_lfc_CT[[LIG]])){
    FDR <- cbind(FDR,sapply(mean_lfc_CT[[LIG]],function(X) X[GENES[temp_ord_gene],"fdr"]))
  }
  if (!is.null(mean_lfc_DS[[LIG]])){
    FDR <- cbind(FDR,sapply(mean_lfc_DS[[LIG]],function(X) X[GENES[temp_ord_gene],"fdr"]))
  }
} else {
  FDR <- ogFDR
}
FDR <- -log10(FDR)
FDR[FDR > 2] <- 2
FDR <- t(FDR / max(FDR))
FDR[FDR == 0] <- NA


DE <- t(ogLFC >= 1 & ogFDR <= 0.1)
if (!is.null(mean_lfc[[LIG]])) {
  DE <- rbind(DE,rep(F,ncol(DE)))
  if (!is.null(mean_lfc[[LIG]]$nts_fdr)) {
    DE <- rbind(DE,mean_lfc[[LIG]][GENES[temp_ord_gene],"nts_lfc"] >= 1 & 
                  mean_lfc[[LIG]][GENES[temp_ord_gene],"nts_fdr"] <= 0.1)
  }
  if (!is.null(mean_lfc[[LIG]]$ts_fdr)) {
    DE <- rbind(DE,mean_lfc[[LIG]][GENES[temp_ord_gene],"ts_lfc"] >= 1 & 
                  mean_lfc[[LIG]][GENES[temp_ord_gene],"ts_fdr"] <= 0.1)
  }
  if (!is.null(mean_lfc_CT[[LIG]])){
    DE <- rbind(DE,
                t(sapply(mean_lfc_CT[[LIG]],function(X) 
                  X[GENES[temp_ord_gene],"lfc"] >= 1 & 
                    X[GENES[temp_ord_gene],"fdr"] <= 0.1)))
  }
  if (!is.null(mean_lfc_DS[[LIG]])){
    DE <- rbind(DE,
                t(sapply(mean_lfc_DS[[LIG]],function(X) 
                  X[GENES[temp_ord_gene],"lfc"] >= 1 & 
                    X[GENES[temp_ord_gene],"fdr"] <= 0.1)))
  }
}


par(mar=c(8,7,1,1),mgp=2:0)
plot(x=NULL,y=NULL,xlim=c(0.5,nrow(LFC) + .5),ylim=c(0.5,ncol(LFC) + .5),
     xaxs="i",yaxs="i",xaxt="n",yaxt="n",xlab=NA,ylab=NA,bty="n",asp=1)
# abline(v=1:nrow(LFC),col="grey90")
# abline(h=1:ncol(LFC),col="grey90")

temp_FG <- temp_BG <- diverging_hcl(1000,palette="Blue-Red")[as.vector(LFC)]
temp_FG[as.vector(DE)] <- "mediumseagreen"
symbols(x=rep(1:nrow(FDR),ncol(FDR)),
        y=as.vector(sapply(1:ncol(FDR),function(X) rep(X,nrow(FDR)))),
        circles=as.vector(FDR)/2,inches=F,add=T,xpd=NA,
        fg=temp_FG,bg=temp_BG)

mtext(colnames(FDR),side=2,at=seq_along(colnames(FDR)),las=2,line=0.1,cex=0.8)
text(x=seq_along(rownames(FDR)),y=rep(par("usr")[3] - 0.5,nrow(FDR)),
     labels=rownames(FDR),col=temp_labelcol,xpd=NA,adj=1,srt=45,cex=0.8)
text(x=0,y=par("usr")[3] - 0.5,labels="* Time-series dataset",
     xpd=NA,adj=1,srt=45,cex=0.8)

temp_yh <- par("usr")[4] - par("usr")[3]
segments(x0=rep(line2user(5,2),1000),x1=rep(line2user(4,2),1000),
         y0=seq(0.6 * temp_yh + par("usr")[3],
                0.9 * temp_yh + par("usr")[3],
                length.out=1000),
         y1=seq(0.6 * temp_yh + par("usr")[3],
                0.9 * temp_yh + par("usr")[3],
                length.out=1000),
         xpd=NA,col=diverging_hcl(1000,palette="Blue-Red"))
rect(xleft=line2user(5,2),xright=line2user(4,2),
     ytop=0.9 * temp_yh + par("usr")[3],
     ybottom=(((0.9 - 0.75) * (1 / max(ogLFC))) + 0.75) * temp_yh + par("usr")[3],
     xpd=NA,border="mediumseagreen")
mtext(c(signif(max(abs(ogLFC)),2),"logFC",-1 * signif(max(abs(ogLFC)),2)),
      side=2,line=c(4,5.1,4),adj=c(-0.1,0.5,1.1),
      at=c(0.9,0.75,0.6) * temp_yh + par("usr")[3])

temp_q <- rev(c(0.5,0.2,0.1,0.05,0.01))
symbols(x=rep(line2user(4.5,2),length(temp_q)),
        y=seq(0.4 * temp_yh + par("usr")[3],
              by=-1.5,length.out=length(temp_q)),
        circles=-log10(temp_q)/4,inches=F,add=T,xpd=NA,
        fg=ifelse(temp_q <= 0.1,"mediumseagreen","black"))
mtext(c(paste0("\u2264",temp_q[1] * 100,"%"),
        paste0(temp_q[2:5] * 100,"%")),
      side=2,line=5,cex=0.8,
      at=seq(0.4 * temp_yh + par("usr")[3],
             by=-1.5,length.out=length(temp_q)))
mtext("FDR",side=2,line=5.8,
      at=0.4 * temp_yh + par("usr")[3] - 3)

mtext(LIG,side=3,line=-0.5,at=line2user(4.5,2),font=2,cex=1.5)

