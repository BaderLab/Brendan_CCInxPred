library(cmapR)
library(ranger)
library(pbapply)

# prepping data ---
if (exists("lvl5_data")) {
} else if (file.exists("~/Dropbox/GDB_archive/CMapCorr_files/lvl5_inputs.RData")) {
  load("~/Dropbox/GDB_archive/CMapCorr_files/lvl5_inputs.RData") 
} else {
  source("lvl5_inputs.R")
}

t(sapply(ct14,function(CT) sapply(lig16,function(LIG)
  sum(lvl5_data@cdesc$cell_id == CT & lvl5_data@cdesc$pert_iname == LIG))))

# Train a model per left-out cell type ----

trainIDs <- testIDs <- list()
for (CT in ct14) {
  temp_lig_id <- sapply(lig16,function(L) 
    rownames(lvl5_data@cdesc)[lvl5_data@cdesc$pert_iname == L &
                                lvl5_data@cdesc$cell_id != CT],
    simplify=F)
  trainIDs[[CT]] <- sapply(temp_lig_id,function(X) sample(X,min(sapply(temp_lig_id,length))),simplify=F)
  trainIDs[[CT]] <- unlist(trainIDs[[CT]],use.names=F)
  testIDs[[CT]] <- rownames(lvl5_data@cdesc)[lvl5_data@cdesc$pert_iname %in% lig16 &
                                               lvl5_data@cdesc$cell_id == CT]
}

# ^ training ----
pboptions(type="timer")
rfmodel <- pbsapply(names(trainIDs),function(N) 
  ranger(x=t(lvl5_data@mat[,trainIDs[[N]]]),
         y=factor(lvl5_data@cdesc[trainIDs[[N]],"pert_iname"],levels=lig16),
         num.threads=8,
         verbose=F),
  simplify=F)

# ^ testing ----
rfresults <- pbsapply(names(trainIDs),function(N)
  predict(rfmodel[[N]],t(lvl5_data@mat[,testIDs[[N]]])),
  simplify=F)

save(rfmodel,rfresults,trainIDs,testIDs,
     file="~/Dropbox/GDB_archive/CMapCorr_files/200708_lvl5_withhold1CT.RData")
