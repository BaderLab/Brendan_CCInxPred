

# pDE_heatmap <- function(countDE,pDE,scoreDE) {
pDEvector <- pDE_ligct[,"01"]

scoreDE <- pDEvector + 1e-4
scoreDE[scoreDE > 1] <- 1
scoreDE <- -log10(scoreDE)

scoreMAT <- tapply(scoreDE,as.factor(sub("^.+_","",names(scoreDE))),c)
for (X in names(scoreMAT)) {
  names(scoreMAT[[X]]) <- sub(paste0("_",X),"",names(scoreMAT[[X]]))
  scoreMAT[[X]] <- scoreMAT[[X]][sort(names(scoreMAT[[X]]))]
}
scoreMAT <- t(do.call(rbind,scoreMAT))

if (ncol(scoreMAT) > 100) {
  scoreMAT <- scoreMAT[,apply(scoreMAT,2,function(X) any(X >= -log10(0.05 + 1e-4)))]
}

temp_hLIG <- hclust(dist(t(scoreMAT)),method="ward.D2")
temp_hCT <- hclust(dist(scoreMAT),method="ward.D2")
scoreMAT <- scoreMAT[temp_hCT$order,temp_hLIG$order]

par(mar=c(3,2,7,6),mgp=2:0,las=2)
image(scoreMAT,
      col=sequential_hcl(100,palette="PinkYl",rev=T),
      x=1:nrow(scoreMAT),y=1:ncol(scoreMAT),
      xaxt="n",yaxt="n",xlab=NA,ylab=NA)
if (nrow(scoreMAT) < 100) {
  mtext(colnames(scoreMAT),side=4,line=0.1,
        at=1:ncol(scoreMAT),adj=0,cex=0.9)
  mtext("Ligands",side=4,line=4.5,las=0,font=2,cex=1.5,
        at=par("usr")[1] + (par("usr")[2] - par("usr")[1]) / 2,adj=0.5)
} else {
  mtext("Ligands",side=4,line=2,las=0,font=2,cex=1.5,
        at=par("usr")[1] + (par("usr")[2] - par("usr")[1]) / 2,adj=0.5)
}
text(x=1:nrow(scoreMAT),y=rep(line2user(0.3,3),nrow(scoreMAT)),
     labels=rownames(scoreMAT),xpd=NA,adj=0,srt=45,cex=0.8)
text(x=0,y=line2user(1,3),labels="Cell lines",
     xpd=NA,adj=0,srt=45,font=2,cex=1.5)

rect(xleft=which(scoreMAT >= -log10(0.05 + 1e-4),arr.ind=T)[,1] - 0.5,
     xright=which(scoreMAT >= -log10(0.05 + 1e-4),arr.ind=T)[,1] + 0.5,
     ybottom=which(scoreMAT >= -log10(0.05 + 1e-4),arr.ind=T)[,2] - 0.5,
     ytop=which(scoreMAT >= -log10(0.05 + 1e-4),arr.ind=T)[,2] + 0.5,
     border="dodgerblue")

temp_left <- par("usr")[1] + (par("usr")[2] - par("usr")[1]) * 0.3
temp_right <- par("usr")[1] + (par("usr")[2] - par("usr")[1]) * 0.7

segments(x0=seq(temp_left,temp_right,length.out=1000),
         x1=seq(temp_left,temp_right,length.out=1000),
         y0=rep(line2user(1,1),1000),y1=rep(line2user(1.6,1),1000),
         xpd=NA,col=sequential_hcl(1000,palette="PinkYl",rev=T))
mtext(c("100%","< 0.01%"),side=1,line=0.8,las=0,adj=c(1.1,-0.1),
      at=c(temp_left,temp_right))
mtext("Probability of at least # DE occuring by chance",
      las=0,side=1,line=1.7,adj=0.5,
      at=temp_left + (temp_right - temp_left) / 2)

rect(xleft=temp_left + (-log10(0.05 + 1e-4) / max(scoreMAT)) * (temp_right - temp_left),
     xright=temp_right,ybottom=line2user(1,1),ytop=line2user(1.6,1),
     xpd=NA,col=NA,border="dodgerblue")
mtext("p <= 0.05",side=1,las=0,col="dodgerblue",cex=0.9,line=0,adj=1,at=temp_right)
