library(cmapR)
library(umap)
library(ranger)
library(pbapply)
pboptions(type="timer")

# prepping data ----

# ^ common ligands from breast samples only ----
if (file.exists("~/Dropbox/GDB_archive/CMapCorr_files/200731_lvl4_breastonly.RData")) {
  load("~/Dropbox/GDB_archive/CMapCorr_files/200731_lvl4_breastonly.RData")
} else {
  if (exists("lvl4_data")) {
  } else if (file.exists("~/Dropbox/GDB_archive/CMapCorr_files/lvl4_inputs.RData")) {
    load("~/Dropbox/GDB_archive/CMapCorr_files/lvl4_inputs.RData") 
  } else {
    source("lvl4_inputs.R")
  }
  rm(lvl4_data_ctl)
  ctBR <- ct14[grepl("breast_",names(ct14))]
  rm(lig16,ct14)
  BRcdesc <- lvl4_data@cdesc[lvl4_data@cdesc$cell_id %in% ctBR,]
  lig29 <- Reduce(intersect,tapply(BRcdesc$pert_iname,BRcdesc$cell_id,c))
  BRcdesc <- BRcdesc[BRcdesc$pert_iname %in% lig29,]
  BRmat <- lvl4_data@mat[,rownames(BRcdesc)]
  save(BRmat,BRcdesc,lig29,ctBR,file="~/Dropbox/GDB_archive/CMapCorr_files/200731_lvl4_breastonly.RData")
}

# ^ UMAP ----
if (!file.exists("~/Dropbox/GDB_archive/CMapCorr_files/200731_lvl4_breastonly_umap.RData")) {
  temp_param <- umap.defaults
  temp_param$n_neighbors <- 30
  temp_param$metric <- "cosine"
  temp_param$min_dist <- 0.2
  temp_param$n_epochs <- 500
  BRumap <- umap(t(BRmat),config=temp_param)
  save(BRumap,file="~/Dropbox/GDB_archive/CMapCorr_files/200731_lvl4_breastonly_umap.RData")
}
# library(colorspace)
# load("~/Dropbox/GDB_archive/CMapCorr_files/200731_lvl4_breastonly_umap.RData")
# Z <- sample(nrow(BRumap$layout))
# par(mar=c(1.5,1,1,0.5),mgp=c(0,0,0))
# plot(BRumap$layout[Z,],pch=20,
#      col=qualitative_hcl(6,palette="dark3",alpha=0.5)[as.factor(BRcdesc$cell_id[Z])],
#      xaxt="n",yaxt="n",main="Coloured by cell type",xlab="UMAP1",ylab="UMAP2")
# plot(BRumap$layout[Z,],pch=20,
#      col=c(qualitative_hcl(29,palette="dark3",alpha=0.5),
#            "black","black")[factor(BRcdesc$pert_iname[Z],levels=c(lig29,"H2O","PBS"))],
#      xaxt="n",yaxt="n",main="Coloured by ligand",xlab="UMAP1",ylab="UMAP2")
# temp_plate <- as.factor(BRcdesc$rna_plate)
# plot(BRumap$layout[Z,],pch=20,
#      col=qualitative_hcl(length(levels(temp_plate)),palette="dark3",alpha=0.5)[temp_plate[Z]],
#      xaxt="n",yaxt="n",main="Coloured by cell type",xlab="UMAP1",ylab="UMAP2")


# ^ No cell line balancing ----
scores_nobal <- list()
for (LIG in lig29) {
  print(paste0(which(lig29 == LIG),"/",length(lig29)))
  testIDs_true <- rownames(BRcdesc)[BRcdesc$pert_iname == LIG]
  scores_nobal[[LIG]] <- pbsapply(seq_along(testIDs_true),function(X){
    trainIDs_true <- testIDs_true[-X]
    trainIDs_false <- sample(setdiff(rownames(BRcdesc),testIDs_true),
                             length(trainIDs_true))
    # table(BRcdesc[trainIDs_true,c("pert_iname","cell_id")])
    # addmargins(table(BRcdesc[trainIDs_false,c("pert_iname","cell_id")]))
    # plot(BRumap$layout[Z,],pch=20,
    #      col=qualitative_hcl(6,palette="dark3",alpha=0.5)[as.factor(BRcdesc$cell_id[Z])],
    #      xaxt="n",yaxt="n",main="Coloured by cell type",xlab="UMAP1",ylab="UMAP2")
    # rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col=scales::alpha("white",.7))
    # points(BRumap$layout[trainIDs_true,],pch=".",cex=2,col="dodgerblue")
    # points(BRumap$layout[trainIDs_false,],pch=".",cex=2,col="red")
    trainIDs <- sample(c(trainIDs_false,trainIDs_true))
    rfmodel <- ranger(x=t(BRmat[,trainIDs]),
                      y=BRcdesc[trainIDs,"pert_iname"] == LIG,
                      num.threads=8,num.trees=1e3,probability=T,
                      verbose=F)
    testIDs <- setdiff(colnames(BRmat),trainIDs)
    testResults <- predict(rfmodel,
                           t(BRmat[,testIDs]))
    
    temp_train <- as.data.frame(rfmodel$predictions[,which(colnames(rfmodel$predictions) == "TRUE")])
    rownames(temp_train) <- trainIDs
    colnames(temp_train) <- "score"
    temp_train$source <- "train"
    temp_train$label <- BRcdesc[trainIDs,"pert_iname"] == LIG
    
    temp_test <- as.data.frame(testResults$predictions[,which(colnames(rfmodel$predictions) == "TRUE")])
    rownames(temp_test) <- testIDs
    colnames(temp_test) <- "score"
    temp_test$source <- "test"
    temp_test$label <- BRcdesc[testIDs,"pert_iname"] == LIG
    temp <- rbind(temp_train,temp_test)
    return(temp[rownames(BRcdesc),])
  },simplify=F)
}
save(scores_nobal,file="~/Dropbox/GDB_archive/CMapCorr_files/200731_probs_mixall_breast_lvl4_nobal.RData")


