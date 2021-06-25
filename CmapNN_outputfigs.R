library(cmapR)
library(colorspace)
library(scales)
source("~/Dropbox/GDB/line2user.R")

# Load CMap ----
load("~/Dropbox/GDB_archive/CMapCorr_files/lvl4_inputs_allgenes.RData")
lvl4_data <- lvl4_data_all
rm(lvl4_data_all)


# Load NN ----
nn_db <- readRDS(url("https://zenodo.org/record/3260758/files/expression_settings.rds"))
nn_db <- nn_db[sapply(sapply(nn_db,"[[","from"),length) == 1]
nn_db <- sapply(nn_db,"[[","diffexp",simplify=F)
for (X in names(nn_db)) {
  nn_db[[X]] <- as.data.frame(nn_db[[X]])
  if (!any(duplicated(nn_db[[X]]$gene))) {
    rownames(nn_db[[X]]) <- nn_db[[X]]$gene
    nn_db[[X]] <- nn_db[[X]][,colnames(nn_db[[X]]) != "gene"]
  } else {
    for (GENE in nn_db[[X]]$gene[duplicated(nn_db[[X]]$gene)]) {
      temp <- nn_db[[X]][nn_db[[X]]$gene == GENE,]
      nn_db[[X]] <- rbind(temp[which.min(temp$qval),],
                          nn_db[[X]][nn_db[[X]]$gene != GENE,])
    }
    rownames(nn_db[[X]]) <- nn_db[[X]]$gene
    nn_db[[X]] <- nn_db[[X]][,colnames(nn_db[[X]]) != "gene"]
  }
}

load("~/Dropbox/GDB_archive/CMapCorr_files/NN_all_dat.RData")
rm(list=c("X","GENE",grep("^temp",ls(),value=T)))


for (LIG in sort(nn_ligands)) {
# LIG <- "IFNA1"

  if (all(DSinfo[[LIG]]$time_series)) { next }
  if (!LIG %in% lvl4_data@cdesc$pert_iname) { next }
  
  # CMap genes ----
  CM_ogZS <- sapply(unique(lvl4_data@cdesc$cell_id[lvl4_data@cdesc$pert_iname == LIG]),function(CT) 
    rowMeans(lvl4_data@mat[,lvl4_data@cdesc$pert_iname == LIG & lvl4_data@cdesc$cell_id == CT]))
  CM_ogFDR <- apply(CM_ogZS,2,function(Z) p.adjust(pnorm(-abs(Z)),method="fdr"))
  colnames(CM_ogFDR) <- colnames(CM_ogZS) <- sapply(colnames(CM_ogZS),function(X) names(ct14)[ct14 == X])
  
  CM_GENES <- rownames(CM_ogFDR)[apply(CM_ogFDR,1,function(X) sum(X <= 0.1) >= 2)]
  # CM_GENES <- rownames(CM_ogFDR)[apply(CM_ogFDR,1,function(X) any(X <= 0.1))]
  
  
  # NN genes ----
  temp_genes <- sort(Reduce(union,sapply(nn_db[rownames(DSinfo[[LIG]])],rownames,simplify=F)))
  NN_ogZS <- sapply(nn_db[rownames(DSinfo[[LIG]])],"[",temp_genes,"lfc")
  rownames(NN_ogZS) <- temp_genes
  NN_ogFDR <- sapply(nn_db[rownames(DSinfo[[LIG]])],"[",temp_genes,"qval")
  rownames(NN_ogFDR) <- temp_genes
  
  temp_ts <- DSinfo[[LIG]][colnames(NN_ogZS),"time_series"]
  colnames(NN_ogFDR) <- colnames(NN_ogZS) <- paste(DSinfo[[LIG]][colnames(NN_ogZS),"CtAcc"],"*")
  colnames(NN_ogFDR)[temp_ts] <- colnames(NN_ogZS)[temp_ts] <- paste0(colnames(NN_ogZS)[temp_ts],"*")
  if (length(unique(DSinfo[[LIG]]$cell_type)) > 1) {
    NN_GENES <- rownames(NN_ogFDR)[apply(NN_ogFDR,1,function(X) sum(X <= 0.1,na.rm=T) > 1 & all(X <= 0.1,na.rm=T)) & 
                                     apply(NN_ogZS,1,function(X) sum(abs(X) >= 1,na.rm=T) > 1 & all(abs(X) >= 1,na.rm=T))]
  } else {
    NN_GENES <- NULL
  }
  
  # Merge genes & fix gene names ----
  GENES <- unique(c(NN_GENES,lvl4_data@rdesc[CM_GENES,"pr_gene_symbol"]))
  # Fixing SELENOP (formerly SEPP1)
  # temp_CMgenes[temp_CMgenes == "SELENOP"] <- "SEPP1"

    # Fixing missing genes
  # temp_CMgenes <- temp_CMgenes[temp_CMgenes %in% lvl4_data@rdesc$pr_gene_symbol]
  # GENES <- GENES[GENES %in% lvl4_data@rdesc$pr_gene_symbol]
  # temp_CMgenes <- temp_CMgenes[temp_CMgenes %in% rownames(NN_ogZS)]
  # GENES <- GENES[GENES %in% rownames(NN_ogZS)]
  
  # names(temp_CMgenes)[names(temp_CMgenes) == "SEPP1"] <- "SELENOP"
  
  
  # Filter CMap ----
  CM_ogZS <- CM_ogZS[match(GENES,lvl4_data@rdesc$pr_gene_symbol),,drop=F]
  hist(CM_ogZS)
  CM_ogFDR <- CM_ogFDR[match(GENES,lvl4_data@rdesc$pr_gene_symbol),,drop=F]
  hist(-log10(CM_ogFDR))
  rownames(CM_ogZS) <- rownames(CM_ogFDR) <- GENES
  
  
  # Filter NN ----
  NN_ogZS <- NN_ogZS[match(GENES,rownames(NN_ogZS)),,drop=F]
  hist(NN_ogZS)
  NN_ogFDR <- NN_ogFDR[match(GENES,rownames(NN_ogFDR)),,drop=F]
  hist(-log10(NN_ogFDR))
  rownames(NN_ogZS) <- rownames(NN_ogFDR) <- GENES
  
  # Merge data ----
  ogZS <- cbind(CM_ogZS,NN_ogZS)
  ogFDR <- cbind(CM_ogFDR,NN_ogFDR)
  
  # HClust ----
  if (nrow(ogZS) > 1) {
    temp_mat <- -log10(ogFDR)
    temp_mat[temp_mat == 0] <- 1e-10
    temp_mat <- temp_mat * ogZS
    temp_corS <- cor(temp_mat,method="pearson",use="pairwise.complete.obs")
    temp_corG <- cor(t(temp_mat),method="pearson",use="pairwise.complete.obs")
    temp_corS[is.na(temp_corS)] <- 0; temp_corG[is.na(temp_corG)] <- 0
    hS <- hclust(as.dist((1 - abs(temp_corS)) * sign(temp_corS)),method="ward.D2")
    hG <- hclust(as.dist((1 - abs(temp_corG)) * sign(temp_corG)),method="ward.D2")
    ogZS <- ogZS[hG$order,hS$order]
    ogFDR <- ogFDR[hG$order,hS$order]
  }
  
  # Scale CMap ----
  CM_cut <- 5
  # boxplot(CM_ogZS,ylim=c(-CM_cut,CM_cut))
  CM_ZS <- CM_ogZS
  CM_ZS[which(CM_ZS > CM_cut)] <- CM_cut; CM_ZS[which(CM_ZS < -CM_cut)] <- -CM_cut
  CM_ZS <- matrix(cut(c(-1,1,as.vector(CM_ZS / max(abs(CM_ZS),na.rm=T))),
                      breaks=1000,labels=F)[-(1:2)],nrow=nrow(CM_ogZS))
  rownames(CM_ZS) <- rownames(CM_ogZS); colnames(CM_ZS) <- colnames(CM_ogZS)
  
  CM_FDR <- -log10(CM_ogFDR)
  CM_FDR[which(CM_FDR > 3)] <- 3
  CM_FDR <- CM_FDR / max(CM_FDR,na.rm=T)
  CM_FDR[CM_FDR <= 0] <- NA
  
  
  # Scale NN ----
  NN_cut <- CM_cut # if you want these to be different, change the legend.
  # boxplot(NN_ogZS,ylim=c(-NN_cut,NN_cut))
  NN_ZS <- NN_ogZS
  NN_ZS[which(NN_ZS > NN_cut)] <- NN_cut; NN_ZS[which(NN_ZS < -NN_cut)] <- -NN_cut
  NN_ZS <- matrix(cut(c(-1,1,as.vector(NN_ZS / max(abs(NN_ZS),na.rm=T))),
                      breaks=1000,labels=F)[-(1:2)],nrow=nrow(NN_ZS))
  rownames(NN_ZS) <- rownames(NN_ogZS); colnames(NN_ZS) <- colnames(NN_ogZS)
  
  NN_FDR <- -log10(NN_ogFDR)
  NN_FDR[NN_FDR > 3] <- 3
  NN_FDR <- NN_FDR / max(NN_FDR,na.rm=T)
  NN_FDR[NN_FDR <= 0] <- NA
  
  
  # Merge plotting data ----
  ZS <- t(cbind(CM_ZS,NN_ZS)[rownames(ogZS),colnames(ogZS)])
  FDR <- t(cbind(CM_FDR,NN_FDR)[rownames(ogFDR),colnames(ogFDR)])
  GENES <- colnames(FDR)
  GENES[lvl4_data@rdesc[match(GENES,lvl4_data@rdesc$pr_gene_symbol),"pr_is_lm"] %in% "0"] <- 
    paste(GENES[lvl4_data@rdesc[match(GENES,lvl4_data@rdesc$pr_gene_symbol),"pr_is_lm"] %in% "0"],"^")
  
  # PLOT ----
  png(paste0("docs/output_figs/CMapNN_",LIG,".png"),
      width=nrow(FDR) * 0.2 + 2,height=ncol(FDR) * 0.2 + 1.6,
      units="in",res=120)
  
  par(mar=c(7,9,1,1),mgp=2:0)
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
  points(which(is.na(ZS),arr.ind=T),pch=4,col="red",cex=0.5)
  mtext("Genes",side=2,line=0.3,at=line2user(0.1,3),las=2,font=2)
  text(rep(line2user(0.3,1),ncol(FDR)),1:ncol(FDR),GENES,xpd=NA,adj=1,srt=45,cex=0.8)
  text(line2user(0,4),line2user(1.3,2),"Cell lines",xpd=NA,adj=1,srt=45,font=2)
  
  text(1:nrow(FDR),rep(line2user(0.3,2),nrow(FDR)),
       labels=rownames(FDR),xpd=NA,adj=1,srt=45,cex=0.8)
  
  mtext(LIG,side=3,line=0,at=line2user(5.4,2),font=2,cex=1.2)
  
  temp_yh <- par("usr")[4] - par("usr")[3]
  mtext("Z-score /",side=2,line=6.9,at=0.7 * temp_yh + par("usr")[3])
  mtext("LogFC *",side=2,line=6.1,at=0.7 * temp_yh + par("usr")[3])
  segments(x0=rep(line2user(6,2),1000),x1=rep(line2user(4.8,2),1000),
           y0=seq(0.4 * temp_yh + par("usr")[3],
                  temp_yh + par("usr")[3],
                  length.out=1000),
           y1=seq(0.4 * temp_yh + par("usr")[3],
                  temp_yh + par("usr")[3],
                  length.out=1000),
           xpd=NA,col=diverging_hcl(1000,palette="Blue-Red3"))
  text(rep(line2user(5.4,2),2),c(1,0.4) * temp_yh + par("usr")[3],
       labels=c(paste0("\u2265",CM_cut),paste0("\u2264",-CM_cut)),
       col="grey80",pos=c(1,3),offset=0.1,xpd=NA)
  
  mtext("FDR",side=2,line=7.5,las=0,
        at=-2 + 0.1 * temp_yh + par("usr")[3])
  temp_q <- rev(c(0.2,0.1,0.05,0.01,0.001))
  symbols(x=rep(line2user(5.4,2),length(temp_q)),
          y=seq(0.1 * temp_yh + par("usr")[3],
                by=-1,length.out=length(temp_q)),
          circles=-log10(temp_q) / (max(-log10(temp_q)) *2),
          inches=F,add=T,xpd=NA,fg="black")
  mtext(c(paste0("\u2264",temp_q[1] * 100,"%"),
          paste0(temp_q[2:5] * 100,"%")),
        side=2,line=6,cex=0.8,las=2,
        at=seq(0.1 * temp_yh + par("usr")[3],
               by=-1,length.out=length(temp_q)))
  
  legend(x=line2user(-0.05,4) - strwidth("missing data",cex=0.6),
         y=line2user(1,3),xpd=NA,bty="n",
         cex=0.6,pt.cex=0.5,pch=4,col="red",legend="missing data")
  mtext(paste("^ inferred in CMap",
              "* microarray data",
              "** timeseries data",
              sep="\n"),
        side=1,line=6,cex=0.6,at=line2user(0.8,4),adj=1)
  
  dev.off()
}
