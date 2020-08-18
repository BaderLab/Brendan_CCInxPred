library(cmapR)
library(ranger)
library(pbapply)
pboptions(type="timer")

# prepping data ---
if (exists("lvl4_data")) {
} else if (file.exists("~/Dropbox/GDB_archive/CMapCorr_files/lvl4_inputs.RData")) {
  load("~/Dropbox/GDB_archive/CMapCorr_files/lvl4_inputs.RData") 
} else {
  source("lvl4_inputs.R")
}
temp <- load("~/Dropbox/GDB_archive/CMapCorr_files/lvl5_inputs.RData") 
rm(list=c("temp",temp[!temp %in% c("ct14","lig16")]))
rm(lvl4_data_ctl)


#### TESTING ####
if (file.exists("~/Dropbox/GDB_archive/CMapCorr_files/200721_lvl4_probs_TESTIDs.RData")) {
  load("~/Dropbox/GDB_archive/CMapCorr_files/200721_lvl4_probs_TESTIDs.RData")
} else {
  testIDs_true_TESTING <- sapply(lig16,function(LIG)
    sample(rownames(lvl4_data@cdesc)[lvl4_data@cdesc$pert_iname == LIG],10),
    simplify=F)
  save(testIDs_true_TESTING,file="~/Dropbox/GDB_archive/CMapCorr_files/200721_lvl4_probs_TESTIDs.RData")
}

# ^ No cell line balancing ----
scores_nobal <- list()
for (LIG in lig16) {
  message()
  message(paste0(which(lig16 == LIG),"/",length(lig16)))

  testIDs_true <- testIDs_true_TESTING[[LIG]]  # TESTING
  # testIDs_true <- rownames(lvl4_data@cdesc)[lvl4_data@cdesc$pert_iname == LIG]
  
  scores_nobal[[LIG]] <- pbsapply(seq_along(testIDs_true),function(X) {
    trainIDs_true <- setdiff(rownames(lvl4_data@cdesc)[lvl4_data@cdesc$pert_iname == LIG],
                             testIDs_true[X])
    trainIDs_false <- sample(
      setdiff(rownames(lvl4_data@cdesc)[lvl4_data@cdesc$pert_iname != LIG],
              trainIDs_true),
      length(trainIDs_true))
    trainIDs <- sample(c(trainIDs_false,trainIDs_true))
    
    rfmodel <- ranger(x=t(lvl4_data@mat[,trainIDs]),
                      y=lvl4_data@cdesc[trainIDs,"pert_iname"] == LIG,
                      num.threads=8,num.trees=1e3,probability=T,
                      verbose=F)
    
    testIDs <- setdiff(colnames(lvl4_data@mat),trainIDs)
    testResults <- predict(rfmodel,
                           t(lvl4_data@mat[,testIDs]))
    
    temp_train <- as.data.frame(rfmodel$predictions[,which(colnames(rfmodel$predictions) == "TRUE")])
    rownames(temp_train) <- trainIDs
    colnames(temp_train) <- "score"
    temp_train$source <- "train"
    temp_train$label <- lvl4_data@cdesc[trainIDs,"pert_iname"] == LIG
    
    temp_test <- as.data.frame(testResults$predictions[,which(colnames(rfmodel$predictions) == "TRUE")])
    rownames(temp_test) <- testIDs
    colnames(temp_test) <- "score"
    temp_test$source <- "test"
    temp_test$label <- lvl4_data@cdesc[testIDs,"pert_iname"] == LIG
    
    temp <- rbind(temp_train,temp_test)
    return(temp[rownames(lvl4_data@cdesc),])
  },simplify=F)
}
# save(scores_nobal,file="~/Dropbox/GDB_archive/CMapCorr_files/200721_lvl4_probs_nobal.RData")
save(scores_nobal,file="~/Dropbox/GDB_archive/CMapCorr_files/200721_lvl4_probs_nobal_TESTING.RData")


# ^ Cell line upsampling ----
scores_upbalCT <- list()
for (LIG in lig16) {
  message()
  message(paste0(which(lig16 == LIG),"/",length(lig16)))

  testIDs_true <- testIDs_true_TESTING[[LIG]]  # TESTING
  # testIDs_true <- rownames(lvl4_data@cdesc)[lvl4_data@cdesc$pert_iname == LIG]
  
  scores_upbalCT[[LIG]] <- pbsapply(seq_along(testIDs_true),function(X) {
    trainIDs_true <- setdiff(rownames(lvl4_data@cdesc)[lvl4_data@cdesc$pert_iname == LIG],
                             testIDs_true[X])
    max_ct <- max(table(lvl4_data@cdesc[trainIDs_true,"cell_id"]))
    trainIDs_true <- unlist(lapply(ct14,function(CT)
      sample(trainIDs_true[lvl4_data@cdesc[trainIDs_true,"cell_id"] == CT],
             size=max_ct,replace=T)))
    trainIDs_false <- sample(
      setdiff(rownames(lvl4_data@cdesc)[lvl4_data@cdesc$pert_iname != LIG],
              trainIDs_true),
      length(trainIDs_true))
    trainIDs <- sample(c(trainIDs_false,trainIDs_true))
    
    rfmodel <- ranger(x=t(lvl4_data@mat[,trainIDs]),
                      y=lvl4_data@cdesc[trainIDs,"pert_iname"] == LIG,
                      num.threads=8,num.trees=1e3,probability=T,
                      verbose=F)
    
    testIDs <- setdiff(colnames(lvl4_data@mat),trainIDs)
    testResults <- predict(rfmodel,
                           t(lvl4_data@mat[,testIDs]))
    temp_maxSD <- max(sapply(unique(trainIDs),function(X) 
      sd(rfmodel$predictions[trainIDs == X,
                             which(colnames(rfmodel$predictions) == "TRUE")])),
      na.rm=T)
    if (temp_maxSD > 0.1) {
      message(paste(
        "Max SD between repeated samples:",
        signif(temp_maxSD,digits=3),
        "   Sample:",
        testIDs_true[X]
      ))
    }
    temp_train <- as.data.frame(
      sapply(unique(trainIDs),function(X) 
        mean(rfmodel$predictions[trainIDs == X,
                                 which(colnames(rfmodel$predictions) == "TRUE")])))
    colnames(temp_train) <- "score"
    rownames(temp_train) <- unique(trainIDs)
    temp_train$source <- "train"
    temp_train$label <- lvl4_data@cdesc[rownames(temp_train),"pert_iname"] == LIG
    
    temp_test <- as.data.frame(testResults$predictions[,which(colnames(rfmodel$predictions) == "TRUE")])
    rownames(temp_test) <- testIDs
    colnames(temp_test) <- "score"
    temp_test$source <- "test"
    temp_test$label <- lvl4_data@cdesc[testIDs,"pert_iname"] == LIG
    
    temp <- rbind(temp_train,temp_test)
    return(temp[rownames(lvl4_data@cdesc),])
  },simplify=F)
}
# save(scores_balCT,file="~/Dropbox/GDB_archive/CMapCorr_files/200721_lvl4_probs_balCT.RData")
save(scores_upbalCT,file="~/Dropbox/GDB_archive/CMapCorr_files/200721_lvl4_probs_upbalCT_TESTING.RData")


# ^ Cell line up/downsampling ----
scores_balCT <- list()
for (LIG in lig16) {
  message()
  message(paste0(which(lig16 == LIG),"/",length(lig16)))
  
  testIDs_true <- testIDs_true_TESTING[[LIG]]  # TESTING
  # testIDs_true <- rownames(lvl4_data@cdesc)[lvl4_data@cdesc$pert_iname == LIG]
  
  scores_balCT[[LIG]] <- pbsapply(seq_along(testIDs_true),function(X) {
    trainIDs_true <- setdiff(rownames(lvl4_data@cdesc)[lvl4_data@cdesc$pert_iname == LIG],
                             testIDs_true[X])
    max_ct <- min(table(lvl4_data@cdesc[trainIDs_true,"cell_id"])) * 3
    trainIDs_true <- unlist(lapply(ct14,function(CT)
      sample(trainIDs_true[lvl4_data@cdesc[trainIDs_true,"cell_id"] == CT],
             size=max_ct,replace=sum(lvl4_data@cdesc[trainIDs_true,"cell_id"] == CT) < max_ct)))
    trainIDs_false <- sample(
      setdiff(rownames(lvl4_data@cdesc)[lvl4_data@cdesc$pert_iname != LIG],
              trainIDs_true),
      length(trainIDs_true))
    trainIDs <- sample(c(trainIDs_false,trainIDs_true))
    
    rfmodel <- ranger(x=t(lvl4_data@mat[,trainIDs]),
                      y=lvl4_data@cdesc[trainIDs,"pert_iname"] == LIG,
                      num.threads=8,num.trees=1e3,probability=T,
                      verbose=F)
    
    testIDs <- setdiff(colnames(lvl4_data@mat),trainIDs)
    testResults <- predict(rfmodel,
                           t(lvl4_data@mat[,testIDs]))
    temp_maxSD <- max(sapply(unique(trainIDs),function(X) 
      sd(rfmodel$predictions[trainIDs == X,
                             which(colnames(rfmodel$predictions) == "TRUE")])),
      na.rm=T)
    if (temp_maxSD > 0.1) {
      message(paste(
        "Max SD between repeated samples:",
        signif(temp_maxSD,digits=3),
        "   Sample:",
        testIDs_true[X]
      ))
    }
    temp_train <- as.data.frame(
      sapply(unique(trainIDs),function(X) 
        mean(rfmodel$predictions[trainIDs == X,
                                 which(colnames(rfmodel$predictions) == "TRUE")])))
    colnames(temp_train) <- "score"
    rownames(temp_train) <- unique(trainIDs)
    temp_train$source <- "train"
    temp_train$label <- lvl4_data@cdesc[rownames(temp_train),"pert_iname"] == LIG
    
    temp_test <- as.data.frame(testResults$predictions[,which(colnames(rfmodel$predictions) == "TRUE")])
    rownames(temp_test) <- testIDs
    colnames(temp_test) <- "score"
    temp_test$source <- "test"
    temp_test$label <- lvl4_data@cdesc[testIDs,"pert_iname"] == LIG
    
    temp <- rbind(temp_train,temp_test)
    return(temp[rownames(lvl4_data@cdesc),])
  },simplify=F)
}
# save(scores_balCT,file="~/Dropbox/GDB_archive/CMapCorr_files/200721_lvl4_probs_balCT.RData")
save(scores_balCT,file="~/Dropbox/GDB_archive/CMapCorr_files/200721_lvl4_probs_balCT_TESTING.RData")
