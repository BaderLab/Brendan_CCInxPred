library(cmapR)
library(ranger)
library(pbapply)

# prepping data ---
if (exists("lvl4new_data")) {
} else if (file.exists("../CMapCorr_files/lvl4new.RData")) {
  load("../CMapCorr_files/lvl4new.RData") 
} else {
  source("200706_ZscoreFromAssayed.R")
}
temp <- load("../CMapCorr_files/lvl5_inputs.RData") 
rm(list=c("temp",temp[!temp %in% c("ct14","lig16")]))


# Data saturation test ----
temp_lig_id <- sapply(lig16,function(L) 
  rownames(lvl4new_data@cdesc)[lvl4new_data@cdesc$pert_iname == L],
  simplify=F)
trainIDs <- sapply(
  seq(1,min(sapply(temp_lig_id,length)) - 1,1),
  function(N)
    sapply(temp_lig_id,function(X) sample(X,N),simplify=F),
  simplify=F)
testIDs <- sapply(trainIDs,function(X)
  mapply(function(all,train) setdiff(all,train),
         all=temp_lig_id,train=X),
  simplify=F)

trainIDs <- sapply(trainIDs,unlist,use.names=F)
testIDs <- sapply(testIDs,unlist,use.names=F)

# ^ training ----
rfmodel <- pbsapply(seq_along(trainIDs),function(N) 
  ranger(x=t(lvl4new_data@mat[,trainIDs[[N]]]),
         y=as.factor(lvl4new_data@cdesc[trainIDs[[N]],"pert_iname"]),
         num.threads=8,
         verbose=F),
  simplify=F)

# ^ testing ----
rfresults <- pbsapply(seq_along(rfmodel),function(N)
  predict(rfmodel[[N]],t(lvl4new_data@mat[,testIDs[[N]]])),
  simplify=F)

save(rfmodel,rfresults,trainIDs,testIDs,
     file="../CMapCorr_files/200707_lvl4new_mixall_balanced_saturated.RData")
