library(cmapR)
library(ranger)

# prepping data ---
if (exists("lvl5_data")) {
} else if (file.exists("~/Dropbox/GDB_archive/CMapCorr_files/lvl5_inputs.RData")) {
  load("~/Dropbox/GDB_archive/CMapCorr_files/lvl5_inputs.RData") 
} else {
  source("lvl5_inputs.R")
}

# Leave one out ----
temp_lig_id <- sapply(lig16,function(L) 
  rownames(lvl5_data@cdesc)[lvl5_data@cdesc$pert_iname == L],
  simplify=F)

testIDs <- sapply(1:10,function(Z)
  sapply(temp_lig_id,function(X) 
    c(X,sample(X,size=max(sapply(temp_lig_id,length)) - length(X)))),
  simplify=F)
testIDs <- do.call(rbind,testIDs)
testIDs <- apply(testIDs,2,sample)

trainIDs <- apply(testIDs,1,function(X) {
  mapply(function(TL,AL) 
    sample(setdiff(AL,TL),size=min(sapply(temp_lig_id,length)) - 1),
    TL=X,AL=temp_lig_id)
})

# ^ training ----
rfmodel <- pbapply::pbapply(trainIDs,2,function(X) 
  ranger(x=t(lvl5_data@mat[,X]),
         y=as.factor(lvl5_data@cdesc[X,"pert_iname"]),
         num.threads=8,
         verbose=F))

# ^ testing ----
rfresults <- pbapply::pbsapply(seq_along(rfmodel),function(N)
  predict(rfmodel[[N]],t(lvl5_data@mat[,testIDs[N,]])),
  simplify=F)

save(rfresults,testIDs,
     file="~/Dropbox/GDB_archive/CMapCorr_files/200524_lvl5_mixall_balanced_leave1out_results.RData")
save(rfmodel,trainIDs,
     file="~/Dropbox/GDB_archive/CMapCorr_files/200524_lvl5_mixall_balanced_leave1out_model.RData")
