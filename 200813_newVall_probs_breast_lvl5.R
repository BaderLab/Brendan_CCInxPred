library(cmapR)
library(umap)
library(ranger)
library(mccf1)

# ^ common ligands from breast samples only ----
if (file.exists("~/Dropbox/GDB_archive/CMapCorr_files/200817_lvl5_breastonly.RData")) {
  load("~/Dropbox/GDB_archive/CMapCorr_files/200817_lvl5_breastonly.RData")
} else {
  if (exists("lvl5_data")) {
  } else if (file.exists("~/Dropbox/GDB_archive/CMapCorr_files/lvl5_inputs.RData")) {
    load("~/Dropbox/GDB_archive/CMapCorr_files/lvl5_inputs.RData") 
  } else {
    source("lvl5_inputs.R")
  }
  ctBR <- ct14[grepl("breast_",names(ct14))]
  rm(lig16,ct14)
  BRcdesc <- lvl5_data@cdesc[lvl5_data@cdesc$cell_id %in% ctBR,]
  lig29 <- Reduce(intersect,tapply(BRcdesc$pert_iname,BRcdesc$cell_id,c))
  BRcdesc <- BRcdesc[BRcdesc$pert_iname %in% lig29,]
  BRmat <- lvl5_data@mat[,rownames(BRcdesc)]
  save(BRmat,BRcdesc,lig29,ctBR,file="~/Dropbox/GDB_archive/CMapCorr_files/200817_lvl5_breastonly.RData")
}

# ^ UMAP ----
if (!file.exists("~/Dropbox/GDB_archive/CMapCorr_files/200817_lvl5_breastonly_umap.RData")) {
  temp_param <- umap.defaults
  temp_param$n_neighbors <- 30
  temp_param$metric <- "cosine"
  temp_param$min_dist <- 0.2
  temp_param$n_epochs <- 500
  BRumap <- umap(t(BRmat),config=temp_param)
  save(BRumap,file="~/Dropbox/GDB_archive/CMapCorr_files/200817_lvl5_breastonly_umap.RData")
}
# library(colorspace)
# # load("~/Dropbox/GDB_archive/CMapCorr_files/200817_lvl5_breastonly_umap.RData")
# Z <- sample(nrow(BRumap$layout))
# par(mar=c(1.5,1,1,0.5),mgp=c(0,0,0))
# plot(BRumap$layout[Z,],pch=20,
#      col=qualitative_hcl(6,palette="dark3",alpha=0.5)[as.factor(BRcdesc$cell_id[Z])],
#      xaxt="n",yaxt="n",main="Coloured by cell type",xlab="UMAP1",ylab="UMAP2")
# plot(BRumap$layout[Z,],pch=20,
#      col=c(qualitative_hcl(29,palette="dark3",alpha=0.5),
#            "black","black")[factor(BRcdesc$pert_iname[Z],levels=c(lig29,"H2O","PBS"))],
#      xaxt="n",yaxt="n",main="Coloured by ligand",xlab="UMAP1",ylab="UMAP2")


# Train RF model ----
scores_all <- scores_new <- list()

for (LIG in lig29) {
  print(paste0(which(lig29 == LIG),"/",length(lig29)," --------"))
  scores_all[[LIG]] <- scores_new[[LIG]] <- list()
  
  for (CT in ctBR) {
    print(paste0("----- ",which(ctBR == CT),"/",length(ctBR)))
    
    # generalize to novel cell line ----
    new_trainID <- rownames(BRcdesc)[BRcdesc$cell_id != CT]
    new_trainID_true <- new_trainID[BRcdesc[new_trainID,"pert_iname"] == LIG]
    new_trainID_false <- sample(setdiff(new_trainID,new_trainID_true),
                                length(new_trainID_true))
    new_trainID <- sample(c(new_trainID_true,new_trainID_false))
    
    rfmodel <- ranger(x=t(BRmat[,new_trainID]),
                      y=BRcdesc[new_trainID,"pert_iname"] == LIG,
                      num.threads=8,num.trees=1e3,probability=T,
                      verbose=F)
    testResults <- predict(rfmodel,t(BRmat))
    
    new_out <- as.data.frame(testResults$predictions[,which(colnames(rfmodel$predictions) == "TRUE")])
    rownames(new_out) <- colnames(BRmat)
    colnames(new_out) <- "score"
    new_out$label <- BRcdesc$pert_iname == LIG
    new_out$train <- rownames(new_out) %in% new_trainID
    # boxplot(score~label+train,data=new_out,las=2,ylim=0:1)
    scores_new[[LIG]][[CT]] <- new_out

    rm(rfmodel,testResults)
    
    # train on all cell lines ----
    all_trainID_true <- sample(rownames(BRcdesc)[BRcdesc$pert_iname == LIG],
                               length(new_trainID_true))
    all_trainID_false <- sample(rownames(BRcdesc)[BRcdesc$pert_iname != LIG],
                                length(new_trainID_false))
    all_trainID <- sample(c(all_trainID_true,all_trainID_false))
    
    rfmodel <- ranger(x=t(BRmat[,all_trainID]),
                      y=BRcdesc[all_trainID,"pert_iname"] == LIG,
                      num.threads=8,num.trees=1e3,probability=T,
                      verbose=F)
    testResults <- predict(rfmodel,t(BRmat))
    
    all_out <- as.data.frame(testResults$predictions[,which(colnames(rfmodel$predictions) == "TRUE")])
    rownames(all_out) <- colnames(BRmat)
    colnames(all_out) <- "score"
    all_out$label <- BRcdesc$pert_iname == LIG
    all_out$train <- rownames(all_out) %in% all_trainID
    # boxplot(score~label+train,data=all_out,las=2,ylim=0:1)
    scores_all[[LIG]][[CT]] <- all_out

    rm(rfmodel,testResults)
    rm(list=grep("^all",ls(),value=T))
    rm(list=grep("^new",ls(),value=T))

  }
}


# Evaluation metric ----

mccf1_all <- sapply(scores_all,function(X) 
  sapply(X,function(Y)
    summary(
      mccf1(as.integer(Y$label[!Y$train]),
            Y$score[!Y$train])
    )$mccf1_metric
  )
)

mccf1_all_thresh <- sapply(scores_all,function(X) 
  sapply(X,function(Y)
    summary(
      mccf1(as.integer(Y$label[!Y$train]),
            Y$score[!Y$train])
    )$best_threshold
  )
)


mccf1_new <- sapply(scores_new,function(X) 
  sapply(X,function(Y)
    summary(
      mccf1(as.integer(Y$label[!Y$train]),
            Y$score[!Y$train])
    )$mccf1_metric
  )
)

mccf1_new_thresh <- sapply(scores_new,function(X) 
  sapply(X,function(Y)
    summary(
      mccf1(as.integer(Y$label[!Y$train]),
            Y$score[!Y$train])
    )$best_threshold
  )
)


save(scores_all,scores_new,mccf1_all,mccf1_all_thresh,mccf1_new,mccf1_new_thresh,
     file="~/Dropbox/GDB_archive/CMapCorr_files/200817_newVall_probs_breast_lvl5.RData")
