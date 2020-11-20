library(pbapply)
pboptions(type="timer")


pwCOR <- function(MAT,IDs) {
  pbsapply(IDs,function(X) {
    temp <- cor(MAT[,X],method="spearman")
    return(temp[upper.tri(temp)])
  },cl=8,simplify=F)
}


load("~/Dropbox/GDB_archive/CMapCorr_files/lvl4_inputs.RData") 
rm(lvl4_data_ctl)
load("~/Dropbox/GDB_archive/CMapCorr_files/lvl4_inputs_allgenes.RData")

lig16 <- sort(lig16)
IDlig16 <- rownames(lvl4_data@cdesc)[lvl4_data@cdesc$pert_iname %in% lig16 & 
                                       lvl4_data@cdesc$cell_id %in% ct14]

temp_ligcountsperct <- sapply(unique(lvl4_data@cdesc$cell_id),function(CT)
  length(unique(lvl4_data@cdesc[lvl4_data@cdesc$cell_id == CT,"pert_iname"])))
ct9 <- names(temp_ligcountsperct)[temp_ligcountsperct > 100]
names(ct9) <- sapply(ct9,function(X) names(ct14)[ct14 == X])

temp_ctcountsperlig <- sapply(unique(lvl4_data@cdesc$pert_iname),function(LIG)
  length(unique(lvl4_data@cdesc[lvl4_data@cdesc$pert_iname == LIG &
                                  lvl4_data@cdesc$cell_id %in% ct9,"cell_id"])))
lig295 <- names(temp_ctcountsperlig)[temp_ctcountsperlig == 9]
lig295 <- sort(lig295)

IDlig295 <- rownames(lvl4_data@cdesc)[lvl4_data@cdesc$pert_iname %in% lig295 & 
                                        lvl4_data@cdesc$cell_id %in% ct9]

rm(list=grep("^temp",ls(),value=T))



# LIG ----
for (LIG in c("lig16","lig295")) {
  for (DAT in c("lvl4_data","lvl4_data_all")) {
    
    temp <- sapply(get(LIG),function(X) 
      get(paste0("ID",LIG))[get(DAT)@cdesc[get(paste0("ID",LIG)),"pert_iname"] == X],
      simplify=F)
    CORS <- pwCOR(get(DAT)@mat,temp)
    temp <- switch(grepl("all",DAT) + 1,
                   paste0(LIG,"_corr_lig.RData"),
                   paste0(LIG,"_corr_allgenes_lig.RData"))
    save(CORS,file=paste0("~/Dropbox/GDB_archive/CMapCorr_files/",temp))
  }
}


# CT ----
for (LIG in c("lig16","lig295")) {
  for (DAT in c("lvl4_data","lvl4_data_all")) {
    CT <- switch(LIG,
                 lig16=ct14,
                 lig295=ct9)
    temp <- sapply(CT,function(X) 
      get(paste0("ID",LIG))[get(DAT)@cdesc[get(paste0("ID",LIG)),"cell_id"] == X],
      simplify=F)
    CORS <- pwCOR(get(DAT)@mat,temp)
    temp <- switch(grepl("all",DAT) + 1,
                   paste0(LIG,"_corr_ct.RData"),
                   paste0(LIG,"_corr_allgenes_ct.RData"))
    save(CORS,file=paste0("~/Dropbox/GDB_archive/CMapCorr_files/",temp))
  }
}

# LigCT ----
for (LIG in c("lig16","lig295")) {
  for (DAT in c("lvl4_data","lvl4_data_all")) {
    CT <- switch(LIG,
                 lig16=ct14,
                 lig295=ct9)
    temp <- sapply(get(LIG),function(L) 
      sapply(CT,function(C)
        get(paste0("ID",LIG))[get(DAT)@cdesc[get(paste0("ID",LIG)),"cell_id"] == C &
                                get(DAT)@cdesc[get(paste0("ID",LIG)),"pert_iname"] == L]
      ),simplify=F)
    temp <- unlist(temp,recursive=F)
    CORS <- pwCOR(get(DAT)@mat,temp)
    temp <- switch(grepl("all",DAT) + 1,
                   paste0(LIG,"_corr_ligct.RData"),
                   paste0(LIG,"_corr_allgenes_ligct.RData"))
    save(CORS,file=paste0("~/Dropbox/GDB_archive/CMapCorr_files/",temp))
  }
}


# REP ----
for (LIG in c("lig16","lig295")) {
  for (DAT in c("lvl4_data","lvl4_data_all")) {
    temp_tx <- unique(get(DAT)@cdesc[get(paste0("ID",LIG)),c("pert_iname","cell_id","pert_dose","pert_time")])
    rownames(temp_tx) <- apply(temp_tx,1,paste,collapse="_")
    temp <- apply(temp_tx,1,function(X)
      rownames(get(DAT)@cdesc)[get(DAT)@cdesc$pert_iname == X[1] &
                                 get(DAT)@cdesc$cell_id == X[2] &
                                 get(DAT)@cdesc$pert_dose == X[3] &
                                 get(DAT)@cdesc$pert_time == X[4]])
    temp <- temp[sort(names(temp))]
    temp <- temp[sapply(temp,length) > 1]

    CORS <- pwCOR(get(DAT)@mat,temp)
    temp <- switch(grepl("all",DAT) + 1,
                   paste0(LIG,"_corr_rep.RData"),
                   paste0(LIG,"_corr_allgenes_rep.RData"))
    save(CORS,file=paste0("~/Dropbox/GDB_archive/CMapCorr_files/",temp))
  }
}

