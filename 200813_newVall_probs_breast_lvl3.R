library(cmapR)
library(ranger)
library(mccf1)
library(pbapply)
pboptions(type="timer")

load("~/Dropbox/GDB_archive/CMapCorr_files/200730_lvl3_breastonly.RData")

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
     file="~/Dropbox/GDB_archive/CMapCorr_files/200813_newVall_probs_breast_lvl3.RData")
