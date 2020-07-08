library(cmapR)
library(ranger)
library(pbapply)

# prepping data ---
if (exists("lvl5_data")) {
} else if (file.exists("../CMapCorr_files/lvl5_inputs.RData")) {
  load("../CMapCorr_files/lvl5_inputs.RData") 
} else {
  source("lvl5_inputs.R")
}

# add metadata to matrix
temp_data_frame <- cbind(lvl5_data@cdesc[,c("cell_id","pert_dose","pert_time")],
                         as.data.frame(t(lvl5_data@mat)))


# Data saturation test ----
temp_lig_id <- sapply(lig16,function(L) 
  rownames(lvl5_data@cdesc)[lvl5_data@cdesc$pert_iname == L],
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
pboptions(type="timer")
rfmodel <- pbsapply(seq_along(trainIDs),function(N) 
  ranger(x=temp_data_frame[trainIDs[[N]],],
         y=as.factor(lvl5_data@cdesc[trainIDs[[N]],"pert_iname"]),
         num.threads=8,
         verbose=F),
  simplify=F)

# ^ testing ----
rfresults <- pbsapply(seq_along(rfmodel),function(N)
  predict(rfmodel[[N]],temp_data_frame[testIDs[[N]],]),
  simplify=F)

save(rfmodel,rfresults,trainIDs,testIDs,
     file="../CMapCorr_files/200703_lvl5_mixall_balanced_saturated_metadata.RData")
