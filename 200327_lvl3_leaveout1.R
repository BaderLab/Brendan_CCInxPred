library(cmapR)
library(ranger)

if (exists("lvl3_data")) {
} else if (file.exists("../CMapCorr_files/lvl3_inputs.RData")) {
  load("../CMapCorr_files/lvl3_inputs.RData") 
} else {
  source("lvl3_inputs.R")
}


train_labels <- test_labels <- 
  test_results <- model <- list()
for (CT in names(PM)) {
  message(paste0("---- Cell type ",which(names(PM) == CT),"/",length(PM)," ----"))
  train_labels[[CT]] <- test_labels[[CT]] <- 
    test_results[[CT]] <- model[[CT]] <- list()
  for (L in lig15) {
    message(paste0("  -- Ligand ",which(lig15 == L),"/",length(lig15)))
    temp_bal <- sapply(PM[names(PM) != CT],function(X) 
      list(tx=sample(X[[L]]$tx,min(sapply(X[[L]],length))),
           ctl=sample(X[[L]]$ctl,min(sapply(X[[L]],length)))),
      simplify=F)
    tx_train <- lvl3_data@mat[,unlist(sapply(temp_bal,function(X) X$tx),use.names=F)]
    ctl_train <- lvl3_data_ctl@mat[,unlist(sapply(temp_bal,function(X) X$ctl),use.names=F)]
    train <- t(cbind(ctl_train,tx_train))
    train_labels[[CT]][[L]] <- as.factor(c(rep("no",ncol(ctl_train)),
                                           rep("yes",ncol(tx_train))))
    tx_test <- lvl3_data@mat[,sample(PM[[CT]][[L]]$tx,min(sapply(PM[[CT]][[L]],length)))]
    ctl_test <- lvl3_data_ctl@mat[,sample(PM[[CT]][[L]]$ctl,min(sapply(PM[[CT]][[L]],length)))]
    test <- t(cbind(ctl_test,tx_test))
    test_labels[[CT]][[L]] <- as.factor(c(rep("no",ncol(ctl_test)),
                                          rep("yes",ncol(tx_test))))
    
    model[[CT]][[L]] <- ranger(x=train,
                               y=train_labels[[CT]][[L]],
                               num.threads=8,
                               verbose=F)
    test_results[[CT]][[L]] <- predict(model[[CT]][[L]],test)$predictions
    
  }
}
save(train_labels,test_labels,test_results,
     file="~/Dropbox/GDB/CMapCorr_files/200327_leaveout1.RData")

