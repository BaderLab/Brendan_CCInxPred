library(cmapR)
library(ranger)
library(pbapply)

# prepping data ---
if (exists("lvl4_data")) {
} else if (file.exists("~/Dropbox/GDB_archive/CMapCorr_files/lvl4_inputs.RData")) {
  load("~/Dropbox/GDB_archive/CMapCorr_files/lvl4_inputs.RData") 
} else {
  source("lvl4_inputs.R")
}
temp <- load("~/Dropbox/GDB_archive/CMapCorr_files/lvl5_inputs.RData") 
rm(list=c("temp",temp[!temp %in% c("ct14","lig16")]))

# Leave one out ----
temp_lig_id <- sapply(lig16,function(L) 
  rownames(lvl4_data@cdesc)[lvl4_data@cdesc$pert_iname == L],
  simplify=F)

testIDs <- sapply(1:2,function(Z)
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

# add metadata to matrix
temp_data_frame <- cbind(lvl4_data@cdesc[,c("cell_id","pert_dose","pert_time")],
                         as.data.frame(t(lvl4_data@mat)))


# ^ training ----
pboptions(type="timer")
rfmodel <- pbapply::pbapply(trainIDs,2,function(X) 
  ranger(x=temp_data_frame[X,],
         y=as.factor(lvl4_data@cdesc[X,"pert_iname"]),
         num.threads=12,
         verbose=F))

save(rfmodel,trainIDs,
     file="~/Dropbox/GDB_archive/CMapCorr_files/200623_lvl4_mixall_balanced_metadata_leave1out_model.RData")

# ^ testing ----
rfresults <- pbapply::pbsapply(seq_along(rfmodel),function(N)
  predict(rfmodel[[N]],temp_data_frame[testIDs[N,],]),
  simplify=F)

save(rfresults,testIDs,
     file="~/Dropbox/GDB_archive/CMapCorr_files/200623_lvl4_mixall_balanced_metadata_leave1out_results.RData")
