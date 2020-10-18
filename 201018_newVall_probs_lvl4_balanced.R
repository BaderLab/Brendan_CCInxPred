library(cmapR)
library(ranger)
library(mccf1)
library(smotefamily)

load("~/Dropbox/GDB_archive/CMapCorr_files/lvl4_lig16_inputs.RData")

temp <- table(lvl4_data@cdesc[,c("cell_id","pert_iname")])
# temp
# range(temp)
temp_n <- min(temp[temp > 20])


scores_all <- scores_new <- list()
for (LIG in lig16) {
  print(paste0(which(lig16 == LIG),"/",length(lig16)," --------"))
  scores_all[[LIG]] <- scores_new[[LIG]] <- list()
  
  for (CT in ct14) {
    print(paste0("----- ",which(ct14 == CT),"/",length(ct14)))

    # generalize to novel cell line ----
    new_trainID <- rownames(lvl4_data@cdesc)[lvl4_data@cdesc$cell_id != CT]
    
temp_cellID <- lvl4_data@cdesc[new_trainID,"cell_id"]
    
    
    
    new_trainID_true <- new_trainID[lvl4_data@cdesc[new_trainID,"pert_iname"] == LIG]
    new_trainID_false <- sample(setdiff(new_trainID,new_trainID_true),
                                length(new_trainID_true))
    new_trainID <- sample(c(new_trainID_true,new_trainID_false))
    
    rfmodel <- ranger(x=t(lvl4_data@mat[,new_trainID]),
                          y=lvl4_data@cdesc[new_trainID,"pert_iname"] == LIG,
                          num.threads=8,num.trees=1e3,probability=T,
                          verbose=F)
    testResults <- predict(rfmodel,t(lvl4_data@mat))
    
    new_out <- as.data.frame(testResults$predictions[,which(colnames(rfmodel$predictions) == "TRUE")])
    rownames(new_out) <- colnames(lvl4_data@mat)
    colnames(new_out) <- "score"
    new_out$label <- lvl4_data@cdesc$pert_iname == LIG
    new_out$train <- rownames(new_out) %in% new_trainID
    # boxplot(score~label+train,data=new_out,las=2,ylim=0:1)
    scores_new[[LIG]][[CT]] <- new_out
    
    rm(rfmodel,testResults)
    
    # train on all cell lines ----
    all_trainID_true <- sample(rownames(lvl4_data@cdesc)[lvl4_data@cdesc$pert_iname == LIG],
                               length(new_trainID_true))
    all_trainID_false <- sample(rownames(lvl4_data@cdesc)[lvl4_data@cdesc$pert_iname != LIG],
                                length(new_trainID_false))
    all_trainID <- sample(c(all_trainID_true,all_trainID_false))
    
    rfmodel <- ranger(x=t(lvl4_data@mat[,all_trainID]),
                      y=lvl4_data@cdesc[all_trainID,"pert_iname"] == LIG,
                      num.threads=8,num.trees=1e3,probability=T,
                      verbose=F)
    testResults <- predict(rfmodel,t(lvl4_data@mat))
    
    all_out <- as.data.frame(testResults$predictions[,which(colnames(rfmodel$predictions) == "TRUE")])
    rownames(all_out) <- colnames(lvl4_data@mat)
    colnames(all_out) <- "score"
    all_out$label <- lvl4_data@cdesc$pert_iname == LIG
    all_out$train <- rownames(all_out) %in% all_trainID
    # boxplot(score~label+train,data=all_out,las=2,ylim=0:1)
    # mean(all_out$score[all_out$label & all_out$train]) - mean(all_out$score[!all_out$label & all_out$train])
    # summary(mccf1(as.integer(all_out$label[!all_out$train]),all_out$score[!all_out$train]))
    
    scores_all[[LIG]][[CT]] <- all_out

    rm(rfmodel,testResults)
    rm(list=grep("^all",ls(),value=T))
    rm(list=grep("^new",ls(),value=T))
  }
}


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
     file="~/Dropbox/GDB_archive/CMapCorr_files/201018_newVall_probs_lvl4_balanced.RData")
