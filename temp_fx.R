temp <- load("~/Dropbox/GDB_archive/CMapCorr_files/lig16_DE_ct_FDR.RData")

FDR_all <- FDR_ct
pDE_all <- pDE_ct
data_name="Ligand / cell line"
aCut=0.1
rm(list=temp)
  
  ALPHA <- c(0.1,0.05,0.01,0.001)
  names(ALPHA) <- colnames(pDE_all)
  
  # trim to genes with <Zp% probability of change by chance in at least one ligand
  if (is.na(aCut)) {
    temp_cut <- ALPHA[1]
  } else {
    temp_cut <- aCut
  }
  FDR <- FDR_all[apply(FDR_all,1,function(X) any(X <= temp_cut)),,drop=F]
  if (ncol(FDR) >= 50) {
    FDR <- FDR[,apply(FDR,2,function(X) any(X <= temp_cut))]
  }
  
  if (nrow(FDR) < 2 | ncol(FDR) < 2) {
    stop(paste0("Less than two samples reached FDR threshold of ",aCut * 100,"%"))
  }
  
  # order ligands
  temp_hROW <- hclust(dist(t(FDR)),method="ward.D2") 
  #### Euclidean distance in -log10 space? ####
  
  # order genes
  temp_hGENE <- hclust(dist(FDR),method="ward.D2")
  
  # order genes by DE prior
  # all(rownames(FDR) %in% DEprior$Gene_EntrezID)
  # FDR <- FDR[DEprior$Gene_EntrezID[DEprior$Gene_EntrezID %in% rownames(FDR)],]
  
  FDR <- FDR[temp_hGENE$order,temp_hROW$order]
  countDE <- t(sapply(ALPHA,function(Z) apply(FDR,2,function(X) sum(X <= Z))))
  FDR <- -log10(FDR)
  pDE <- pDE_all[colnames(FDR),]
  scoreDE <- pDE + 1e-4
  scoreDE[scoreDE > 1] <- 1
  scoreDE <- t(-log10(scoreDE))
  
  
  layout(rbind(1:2),widths=c(6,1))
  if (nrow(FDR) >= 100) {
    par(mar=c(4,6,1,0.5),mgp=2:0)
    label_cex <- 0.9
  } else if (data_name == "Ligands") {
    par(mar=c(4,6,1,0.5),mgp=2:0)
    label_cex <- 0.9
  } else if (data_name == "Cell lines") {
    par(mar=c(4,10,1.5,0.5),mgp=2:0)
    label_cex <- 0.9
  } else {
    stop("data_name must be 'Ligands' or 'Cell lines'")
  }
  image(z=FDR / max(abs(FDR)),
        x=1:nrow(FDR),y=1:ncol(FDR),
        col=sequential_hcl(99,palette="Emrld",rev=T),
        breaks=seq(0,1,length.out=100),
        xaxt="n",yaxt="n",xlab=NA,ylab=NA)
  if (nrow(FDR) < 100) {
    mtext(lvl4_data@rdesc[rownames(FDR),"pr_gene_symbol"],
          side=1,las=2,at=1:nrow(FDR),adj=1.1,cex=0.7)
  } else {
    mtext("Genes",side=1,line=1,font=2,cex=1.5)
  }
  if (ncol(FDR) < 100) {
    mtext(colnames(FDR),side=2,las=2,
          at=1:ncol(FDR),adj=1,line=0.1,cex=label_cex)
  }
  if (data_name == "Ligands") {
    mtext("Ligands",side=2,line=4.5,font=2,cex=1.5,las=0)
  } else if (data_name == "Cell lines") {
    mtext("Cell lines",side=3,line=0,font=2,cex=1.5,las=0,
          at=par("usr")[1],adj=1)
  } 
  
  segments(x0=seq(line2user(5,2),line2user(1,2),length.out=1001),
           x1=seq(line2user(5,2),line2user(1,2),length.out=1001),
           y0=rep(line2user(1,1),1001),y1=rep(line2user(1.7,1),1001),
           xpd=NA,col=sequential_hcl(1001,palette="Emrld",rev=T))
  text(c(line2user(5,2),line2user(1,2)),
       rep(line2user(2.2,1),2),xpd=NA,
       labels=c(min(FDR),round(max(FDR))))
  text(line2user(3,2),line2user(3,1),xpd=NA,labels="-log10(FDR)")
  
  if (!is.na(aCut)) {
    rect(xleft=which(FDR >= -log10(aCut),arr.ind=T)[,1] - 0.5,
         xright=which(FDR >= -log10(aCut),arr.ind=T)[,1] + 0.5,
         ybottom=which(FDR >= -log10(aCut),arr.ind=T)[,2] - 0.5,
         ytop=which(FDR >= -log10(aCut),arr.ind=T)[,2] + 0.5,
         border="dodgerblue")
    rect(xleft=seq(line2user(5,2),line2user(1,2),
                   length.out=1001)[round((-log10(aCut) / max(FDR)) * 1001 / 2)],
         xright=line2user(1,2),ybottom=line2user(1.7,1),ytop=line2user(1,1),
         xpd=NA,col=NA,border="dodgerblue")
    text(line2user(1,2),line2user(0.5,1),adj=1,xpd=NA,
         labels=paste0("FDR = ",aCut * 100,"%"),col="dodgerblue",cex=0.9)  
  }
  
  if (data_name == "Ligands") {
    par(mar=c(4,0.5,1,1),mgp=2:0)
  } else if (data_name == "Cell lines") {
    par(mar=c(4,0.5,1.5,1),mgp=2:0)
  }
  image(z=scoreDE / max(scoreDE),
        x=1:nrow(scoreDE),y=1:ncol(scoreDE),
        col=sequential_hcl(99,palette="PinkYl",rev=T),
        breaks=seq(0,1,length.out=100),
        xaxt="n",yaxt="n",xlab=NA,ylab=NA)
  abline(v=seq(1.5,by=1,length.out=3),col="grey50")
  mtext(c("10.0","5.0","1.0","0.1"),
        side=1,at=1:4,las=2,line=0.1,cex=0.9)
  mtext("FDR %",side=1,at=0.2,las=2,line=0.1)
  if (ncol(FDR) < 100) {
    text(as.vector(sapply(1:ncol(countDE),function(X) 1:nrow(countDE))),
         as.vector(sapply(1:ncol(countDE),function(X) rep(X,nrow(countDE)))),
         labels=as.vector(countDE),cex=0.8,col="grey20")
  }
  segments(x0=seq(1.5,4,length.out=1001),
           x1=seq(1.5,4,length.out=1001),
           y0=rep(line2user(2.1,1),1001),y1=rep(line2user(2.6,1),1001),
           xpd=NA,col=sequential_hcl(1001,palette="PinkYl",rev=T))
  text(x=c(1.2,4.3),y=rep(line2user(2.4,1),2),xpd=NA,
       labels=c(0,max(scoreDE)))
  mtext("-log10(P)",side=1,line=2.7,at=2.75)
  mtext("# DE at FDR thresholds",side=4)
  
  par(.PAR)

