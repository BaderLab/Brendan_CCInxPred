library(cmapR)
library(ranger)
library(pbapply)

# prepping data ---
if (exists("lvl4_data")) {
} else if (file.exists("../CMapCorr_files/lvl4_inputs.RData")) {
  load("../CMapCorr_files/lvl4_inputs.RData") 
} else {
  source("lvl4_inputs.R")
}
temp <- load("../CMapCorr_files/lvl5_inputs.RData") 
rm(list=c("temp",temp[!temp %in% c("ct14","lig16")]))

# t(sapply(ct14,function(CT) sapply(lig16,function(LIG) 
#   sum(lvl4_data@cdesc$cell_id == CT & lvl4_data@cdesc$pert_iname == LIG))))

# Train a model per left-out cell type ----

trainIDs <- testIDs <- list()
for (CT in ct14) {
  temp_lig_id <- sapply(lig16,function(L) 
    rownames(lvl4_data@cdesc)[lvl4_data@cdesc$pert_iname == L &
                                lvl4_data@cdesc$cell_id != CT],
    simplify=F)
  trainIDs[[CT]] <- sapply(temp_lig_id,function(X) sample(X,min(sapply(temp_lig_id,length))),simplify=F)
  trainIDs[[CT]] <- unlist(trainIDs[[CT]],use.names=F)
  testIDs[[CT]] <- rownames(lvl4_data@cdesc)[lvl4_data@cdesc$pert_iname %in% lig16 &
                                               lvl4_data@cdesc$cell_id == CT]
}

# ^ training ----
pboptions(type="timer")
rfmodel <- pbsapply(names(trainIDs),function(N) 
  ranger(x=t(lvl4_data@mat[,trainIDs[[N]]]),
         y=factor(lvl4_data@cdesc[trainIDs[[N]],"pert_iname"],levels=lig16),
         num.threads=8,
         verbose=F),
  simplify=F)

# ^ testing ----
rfresults <- pbsapply(names(trainIDs),function(N)
  predict(rfmodel[[N]],t(lvl4_data@mat[,testIDs[[N]]])),
  simplify=F)

save(rfmodel,rfresults,trainIDs,testIDs,
     file="../CMapCorr_files/200708_lvl4_withhold1CT.RData")
