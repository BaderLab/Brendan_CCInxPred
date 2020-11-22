plot(NA,NA,xlim=range(scoreMAT[,LIG]),ylim=c(0,15),bty="n",yaxt="n",
     xlab="Prediction accuracy per cell type",ylab=NA)
axis(2,line=3.5)
mtext("Receptor expression per cell type",
      side=2,line=5.5,cex=0.7)
title(LIG,line=0)
segments(x0=seq(min(scoreMAT),max(scoreMAT),length.out=1000),
         x1=seq(min(scoreMAT),max(scoreMAT),length.out=1000),
         y0=rep(par("usr")[3],1000),
         y1=rep(par("usr")[3] + (par("usr")[4] - par("usr")[3]) * 0.02,1000),
         lwd=2,col=sequential_hcl(1000,palette="PinkYl",rev=T))

temp_order <- rownames(scoreMAT)[order(scoreMAT[,LIG])]
temp_corr <- apply(meanR[lig16R[[LIG]],temp_order],1,function(X)
  cor(scoreMAT[temp_order,LIG],X,method="spearman"))
temp_col <- diverge_hcl(100,palette="Purple-Green")[cut(c(-1,1,temp_corr),100,labels=F)[c(-1,-2)]]
names(temp_col) <- lig16R[[LIG]]
for (GENE in lig16R[[LIG]][order(temp_corr)]) {
  points(scoreMAT[temp_order,LIG],
         meanR[GENE,temp_order],
         type="b",pch=20,col=temp_col[GENE],lwd=2)
}
temp_top <- names(temp_corr)[order(temp_corr,decreasing=T)][1:3]
text(scoreMAT[temp_order[length(temp_order)],LIG],
     meanR[temp_top,temp_order[length(temp_order)]],
     labels=temp_top,col=temp_col[temp_top],
     pos=4,xpd=NA,offset=0.4,
     font=as.integer(sapply(temp_top,function(X) 
       lvl3_ctl@rdesc[lvl3_ctl@rdesc$pr_gene_symbol == X,"pr_is_lm"])) + 1)
text(scoreMAT[temp_order[1],LIG],
     meanR[temp_top,temp_order[1]],
     labels=temp_top,col=temp_col[temp_top],
     pos=2,xpd=NA,offset=0.4,
     font=as.integer(sapply(temp_top,function(X) 
       lvl3_ctl@rdesc[lvl3_ctl@rdesc$pr_gene_symbol == X,"pr_is_lm"])) + 1)

