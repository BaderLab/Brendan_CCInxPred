library(cmapR)
library(ranger)

# prepping data ---
if (exists("lvl5_data")) {
} else if (file.exists("../CMapCorr_files/lvl5_inputs.RData")) {
  load("../CMapCorr_files/lvl5_inputs.RData") 
} else {
  source("lvl5_inputs.R")
}

temp_id <- rownames(lvl5_data@cdesc)[lvl5_data@cdesc$pert_iname %in% lig15]
trainIDs <- sample(temp_id,round(length(temp_id) / 2))
testIDs <- setdiff(temp_id,trainIDs)

# training ----
rfmodel <- ranger(x=t(lvl5_data@mat[,trainIDs]),
                  y=as.factor(lvl5_data@cdesc[trainIDs,"pert_iname"]),
                  num.threads=8,
                  verbose=F)
# testing ----
rfresults <- predict(rfmodel,t(lvl5_data@mat[,testIDs]))

save(rfmodel,rfresults,trainIDs,testIDs,
     file="../CMapCorr_files/200427_lvl5_mixall.RData")


# balancing data ----
temp_lig_id <- sapply(lig15,function(L) 
  rownames(lvl5_data@cdesc)[lvl5_data@cdesc$pert_iname == L],
  simplify=F)
temp_size <- round(min(sapply(temp_lig_id,length)) / 2)
trainIDs <- sapply(temp_lig_id,function(X) sample(X,temp_size),simplify=F)
testIDs <- mapply(function(all,train) setdiff(all,train),
                  all=temp_lig_id,train=trainIDs)
## balancing test data:
# temp_size <- min(temp_size,sapply(testIDs,length))
# testIDs <-sapply(testIDs,function(X) sample(X,temp_size),simplify=F)
##
trainIDs <- unlist(trainIDs,use.names=F)
testIDs <- unlist(testIDs,use.names=F)

# training ----
rfmodel <- ranger(x=t(lvl5_data@mat[,trainIDs]),
                  y=as.factor(lvl5_data@cdesc[trainIDs,"pert_iname"]),
                  num.threads=8,
                  verbose=F)
# testing ----
rfresults <- predict(rfmodel,t(lvl5_data@mat[,testIDs]))

save(rfmodel,rfresults,trainIDs,testIDs,
     file="../CMapCorr_files/200427_lvl5_mixall_balanced.RData")


