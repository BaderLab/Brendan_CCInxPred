inD <- c("lig16_DE","lig295_DE")
names(inD) <- c("16 ligands","295 ligands")
DSname <- inD
D <- DSname[2]

# par(mfrow=c(1,length(DSname)))
# for (D in DSname) {
  load(paste0("~/Dropbox/GDB_archive/CMapCorr_files/",D,"_lig_FDR.RData"))
  scoreDElig <- pDE_lig + 1e-4
  scoreDElig[scoreDElig > 1] <- 1
  scoreDElig <- -log10(scoreDElig)
  load(paste0("~/Dropbox/GDB_archive/CMapCorr_files/",D,"_ct_FDR.RData"))
  scoreDEct <- pDE_ct + 1e-4
  scoreDEct[scoreDEct > 1] <- 1
  scoreDEct <- -log10(scoreDEct)
  load(paste0("~/Dropbox/GDB_archive/CMapCorr_files/",D,"_ligct_FDR.RData"))
  scoreDEligct <- pDE_ligct + 1e-4
  scoreDEligct[scoreDEligct > 1] <- 1
  scoreDEligct <- -log10(scoreDEligct)
  load(paste0("~/Dropbox/GDB_archive/CMapCorr_files/",D,"_rep_FDR.RData"))
  scoreDErep <- pDE_rep + 1e-4
  scoreDErep[scoreDErep > 1] <- 1
  scoreDErep <- -log10(scoreDErep)
  
  par(mar=c(5,3,2,1),mgp=2:0)
  plot(NA,NA,xlim=c(1,19),ylim=c(0,4),
       xaxt="n",xlab=NA,ylab="-log10(P) # DE at FDR thresholds",
       main=names(DSname)[which(DSname == D)])
  boxplot(scoreDEct,at=seq(1,by=5,length.out=4),
          add=T,xaxt="n",yaxt="n",pch=".",cex=2,
          border=qualitative_hcl(4,palette="dynamic")[1],
          col=qualitative_hcl(4,alpha=0.5,palette="dynamic")[1])
  temp_stat <- sapply(colnames(scoreDElig),function(X)
    wilcox.test(scoreDElig[,X],scoreDEct[,X])$p.value)
  if (any(temp_stat <= 2.2e-16)) {
    mtext("*",side=3,line=-0.4,
          at=seq(1,by=5,length.out=4)[temp_stat <= 2.2e-16])
  }
  boxplot(scoreDElig,at=seq(2,by=5,length.out=4),
          add=T,xaxt="n",yaxt="n",pch=".",cex=2,
          border=qualitative_hcl(4,palette="dynamic")[2],
          col=qualitative_hcl(4,alpha=0.5,palette="dynamic")[2])
  boxplot(scoreDEligct,at=seq(3,by=5,length.out=4),
          add=T,xaxt="n",yaxt="n",pch=".",cex=2,
          border=qualitative_hcl(4,palette="dynamic")[3],
          col=qualitative_hcl(4,alpha=0.5,palette="dynamic")[3])
  temp_stat <- sapply(colnames(scoreDElig),function(X) 
    wilcox.test(scoreDElig[,X],scoreDEligct[,X])$p.value)
  if (any(temp_stat <= 2.2e-16)) {
    mtext("*",side=3,line=-0.4,
          at=seq(3,by=5,length.out=4)[temp_stat <= 2.2e-16])
  }
  boxplot(scoreDErep,at=seq(4,by=5,length.out=4),
          add=T,xaxt="n",yaxt="n",pch=".",cex=2,
          border=qualitative_hcl(4,palette="dynamic")[4],
          col=qualitative_hcl(4,alpha=0.5,palette="dynamic")[4])
  temp_stat <- sapply(colnames(scoreDElig),function(X) 
    wilcox.test(scoreDElig[,X],scoreDErep[,X])$p.value)
  if (any(temp_stat <= 2.2e-16)) {
    mtext("*",side=3,line=-0.4,
          at=seq(4,by=5,length.out=4)[temp_stat <= 2.2e-16])
  }
  # points(seq(1,by=5,length.out=4),colMeans(scoreDEct),
  #        pch="-",cex=3,col=alpha("red",0.5))
  # points(seq(2,by=5,length.out=4),colMeans(scoreDElig),
  #        pch="-",cex=3,col=alpha("red",0.5))
  # points(seq(3,by=5,length.out=4),colMeans(scoreDEligct),
  #        pch="-",cex=3,col=alpha("red",0.5))
  # points(seq(4,by=5,length.out=4),colMeans(scoreDErep),
  #        pch="-",cex=3,col=alpha("red",0.5))
  mtext(rep(c("Cell line","Ligand","Lig / Line","Replicate",NA),4),
        col=c(qualitative_hcl(4,palette="dynamic"),NA),
        side=1,line=0.1,las=2,at=1:19,cex=0.8)
  mtext(paste0("FDR = ",c(10,5,1,0.1),"%"),
        side=1,line=3.5,at=seq(2.5,by=5,length.out=4))

