library(cmapR)
library(ranger)

if (exists("lvl3_data")) {
} else if (file.exists("../CMapCorr_files/lvl3_inputs.RData")) {
  load("../CMapCorr_files/lvl3_inputs.RData") 
} else {
  source("lvl3_inputs.R")
}


train_labels <- train_labelsB <- test_labels <- 
  test_results <- test_resultsW <- test_resultsB <-
  model <- modelW <- modelB <- list()
for (L in lig16) {
  message(paste0("---- Ligand ",which(lig16 == L),"/",length(lig16)," ----"))
  tx <- lvl3_data@mat[,lvl3_data@cdesc$pert_iname == L]
  tx_train <- tx[,sample(colnames(tx),round(ncol(tx) / 2))]
  tx_test <- tx[,!colnames(tx) %in% colnames(tx_train)]
  ctl_train <- lvl3_data_ctl@mat[,sample(colnames(lvl3_data_ctl@mat),
                                         round(ncol(lvl3_data_ctl@mat) / 2))]
  ctl_test <- lvl3_data_ctl@mat[,!colnames(lvl3_data_ctl@mat) %in% colnames(ctl_train)]
  train <- t(cbind(ctl_train,tx_train))
  train_labels[[L]] <- as.factor(c(rep("no",ncol(ctl_train)),
                                   rep("yes",ncol(tx_train))))
  test <- t(cbind(ctl_test,tx_test))
  test_labels[[L]] <- as.factor(c(rep("no",ncol(ctl_test)),
                                  rep("yes",ncol(tx_test))))
  
  model[[L]] <- ranger(x=train,
                       y=train_labels[[L]],
                       num.threads=8,
                       verbose=F)
  modelW[[L]] <- ranger(x=train,
                        y=train_labels[[L]],
                        class.weights=length(train_labels[[L]]) /
                          (length(unique(train_labels[[L]])) * 
                             table(train_labels[[L]])),
                        # sklearn.ensemble.RandomForestClassifier(class_weight="balanced")
                        num.threads=8,
                        verbose=F)
  test_results[[L]] <- predict(model[[L]],test)$predictions
  test_resultsW[[L]] <- predict(modelW[[L]],test)$predictions
  
  
  dwn <- sample(colnames(ctl_train),ncol(tx_train))
  train <- t(cbind(ctl_train[,dwn],tx_train))
  train_labelsB[[L]] <- as.factor(c(rep("no",length(dwn)),
                                    rep("yes",ncol(tx_train))))
  modelB[[L]] <- ranger(x=train,
                        y=train_labelsB[[L]],
                        num.threads=8,
                        verbose=F)
  test_resultsB[[L]] <- predict(modelB[[L]],test)$predictions
}

save(train_labels,train_labelsB,test_labels,
     test_results,test_resultsW,test_resultsB,
     model,modelW,modelB,
     file="~/Dropbox/GDB/CMapCorr_files/200326_mixall.RData")
