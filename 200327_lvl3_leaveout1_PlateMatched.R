library(cmapR)
library(ranger)

if (exists("lvl3_data")) {
} else if (file.exists("../CMapCorr_files/lvl3_inputs.RData")) {
  load("../CMapCorr_files/lvl3_inputs.RData") 
} else {
  source("lvl3_inputs.R")
}

PM <- sapply(
  unique(lvl3_data@cdesc[lvl3_data@cdesc$pert_iname %in% lig16,"cell_id"]),
  function(CT) {
    sapply(lig16,function(L) {
      temp_tx <- rownames(lvl3_data@cdesc)[lvl3_data@cdesc$cell_id == CT & 
                                             lvl3_data@cdesc$pert_iname == L]
      temp_pl <- unique(lvl3_data@cdesc[temp_tx,"rna_plate"])
      temp_ctl <- rownames(lvl3_data_ctl@cdesc)[lvl3_data_ctl@cdesc$rna_plate %in% temp_pl &
                                                  lvl3_data_ctl@cdesc$cell_id == CT]
      # message(paste(length(temp_tx),length(temp_ctl)))
      return(list(tx=temp_tx,
                  ctl=temp_ctl))
    },simplify=F)
  },simplify=F)

train_labels <- test_labels <- 
  test_results <- model <- list()
for (CT in names(PM)) {
  message(paste0("---- Cell type ",which(names(PM) == CT),"/",length(PM)," ----"))
  train_labels[[CT]] <- test_labels[[CT]] <- 
    test_results[[CT]] <- model[[CT]] <- list()
  for (L in lig16) {
    message(paste0("  -- Ligand ",which(lig16 == L),"/",length(lig16)))
    temp_bal <- sapply(PM[names(PM) != CT],function(X) 
      list(tx=sample(X[[L]]$tx,min(sapply(X[[L]],length))),
           ctl=sample(X[[L]]$ctl,min(sapply(X[[L]],length)))),
      simplify=F)
    tx_train <- lvl3_data@mat[,unlist(sapply(temp_bal,function(X) X$tx),use.names=F)]
    ctl_train <- lvl3_data_ctl@mat[,unlist(sapply(temp_bal,function(X) X$ctl),use.names=F)]
    train <- t(cbind(ctl_train,tx_train))
    train_labels[[CT]][[L]] <- as.factor(c(rep("no",ncol(ctl_train)),
                                           rep("yes",ncol(tx_train))))
    tx_test <- lvl3_data@mat[,PM[[CT]][[L]]$tx]
    ctl_test <- lvl3_data_ctl@mat[,PM[[CT]][[L]]$ctl]
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
     file="~/Dropbox/GDB/CMapCorr_files/200327_leaveout1_PM.RData")

