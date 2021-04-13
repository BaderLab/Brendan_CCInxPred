library(pbapply)
pboptions(type="timer")


pwCOR <- function(MAT,IDs) {
  pbsapply(IDs,function(X) {
    temp <- cor(MAT[,X],method="spearman")
    return(temp[upper.tri(temp)])
  },cl=8,simplify=F)
}

load("~/Dropbox/GDB_archive/CMapCorr_files/lvl4_inputs.RData")
temp <- load("~/Dropbox/GDB_archive/CMapCorr_files/GSEA_lvl4.RData")
for (X in temp) {
  assign(sub("nes","lvl4",X),get(X))
}
rm(list=temp)

load("~/Dropbox/GDB_archive/CMapCorr_files/lvl5_inputs.RData")
temp <- load("~/Dropbox/GDB_archive/CMapCorr_files/GSEA_lvl5.RData")
for (X in temp) {
  assign(sub("nes","lvl5",X),get(X))
}
rm(list=temp)
rm(X)

temp_ligcountsperct <- sapply(unique(lvl5_data@cdesc$cell_id),function(CT)
  length(unique(lvl5_data@cdesc[lvl5_data@cdesc$cell_id == CT,"pert_iname"])))
ct9 <- names(temp_ligcountsperct)[temp_ligcountsperct > 100]
names(ct9) <- sapply(ct9,function(X) names(ct14)[ct14 == X])

temp_ctcountsperlig <- sapply(unique(lvl5_data@cdesc$pert_iname),function(LIG)
  length(unique(lvl5_data@cdesc[lvl5_data@cdesc$pert_iname == LIG &
                                  lvl5_data@cdesc$cell_id %in% ct9,"cell_id"])))
lig295 <- names(temp_ctcountsperlig)[temp_ctcountsperlig == 9]
lig295 <- sort(lig295)

rm(list=grep("^temp",ls(),value=T))



# LIG ----
print("LIG ------------------")
for (DAT in c("lvl5_data","lvl4_data")) {
  temp <- sapply(lig295,function(LIG) 
    rownames(get(DAT)@cdesc)[get(DAT)@cdesc$pert_iname == LIG & 
                               get(DAT)@cdesc$cell_id %in% ct9],
    simplify=F)
  for (GSEA in grep(paste0("gsea_",sub("_data","",DAT)),ls(),value=T)) {
    print(paste0("-",GSEA," ----"))
    CORS <- pwCOR(get(GSEA),temp)
    temp_filename <- paste0("corrGSEA_",sub("_data","",DAT),"_",sub("gsea_lvl[45]_","",GSEA),"_lig.RData")
    save(CORS,file=paste0("~/Dropbox/GDB_archive/CMapCorr_files/",temp_filename))
  }
  rm(list=grep("^temp",ls(),value=T))
}


# CT ----
print("CT ------------------")
for (DAT in c("lvl5_data","lvl4_data")) {
  temp <- sapply(ct9,function(CT) 
    rownames(get(DAT)@cdesc)[get(DAT)@cdesc$pert_iname %in% lig295 & 
                               get(DAT)@cdesc$cell_id == CT],
    simplify=F)
  for (GSEA in grep(paste0("gsea_",sub("_data","",DAT)),ls(),value=T)) {
    print(paste0("-",GSEA," ----"))
    CORS <- pwCOR(get(GSEA),temp)
    temp_filename <- paste0("corrGSEA_",sub("_data","",DAT),"_",sub("gsea_lvl[45]_","",GSEA),"_ct.RData")
    save(CORS,file=paste0("~/Dropbox/GDB_archive/CMapCorr_files/",temp_filename))
  }
  rm(list=grep("^temp",ls(),value=T))
}


# LigCT ----
print("LIGCT ------------------")
for (DAT in c("lvl5_data","lvl4_data")) {
  temp <- sapply(lig295,function(LIG) {
    sapply(ct9,function(CT) {
      rownames(get(DAT)@cdesc)[get(DAT)@cdesc$pert_iname == LIG & 
                                 get(DAT)@cdesc$cell_id == CT]
    })
  },simplify=F)
  temp_names <- unlist(sapply(names(temp),function(X) paste(X,names(temp[[X]]),sep="_")),use.names=F)
  temp <- unlist(temp,recursive=F,use.names=F)
  names(temp) <- temp_names
  temp <- temp[sapply(temp,length) > 1]
  for (GSEA in grep(paste0("gsea_",sub("_data","",DAT)),ls(),value=T)) {
    print(paste0("-",GSEA," ----"))
    CORS <- pwCOR(get(GSEA),temp)
    temp_filename <- paste0("corrGSEA_",sub("_data","",DAT),"_",sub("gsea_lvl[45]_","",GSEA),"_ligct.RData")
    save(CORS,file=paste0("~/Dropbox/GDB_archive/CMapCorr_files/",temp_filename))
  }
  rm(list=grep("^temp",ls(),value=T))
}


# REP ----
print("REP ------------------")
for (DAT in c("lvl4_data")) {
  temp_tx <- unique(get(DAT)@cdesc[get(DAT)@cdesc$pert_iname %in% lig295 & 
                                     get(DAT)@cdesc$cell_id %in% ct9,
                                   c("pert_iname","cell_id","pert_dose","pert_time")])
  rownames(temp_tx) <- apply(temp_tx,1,paste,collapse="_")
  temp <- apply(temp_tx,1,function(X)
    rownames(get(DAT)@cdesc)[get(DAT)@cdesc$pert_iname == X[1] &
                               get(DAT)@cdesc$cell_id == X[2] &
                               get(DAT)@cdesc$pert_dose == X[3] &
                               get(DAT)@cdesc$pert_time == X[4]])
  temp <- temp[sort(names(temp))]
  temp <- temp[sapply(temp,length) > 1]
  for (GSEA in grep(paste0("gsea_",sub("_data","",DAT)),ls(),value=T)) {
    print(paste0("-",GSEA," ----"))
    CORS <- pwCOR(get(GSEA),temp)
    temp_filename <- paste0("corrGSEA_",sub("_data","",DAT),"_",sub("gsea_lvl[45]_","",GSEA),"_rep.RData")
    save(CORS,file=paste0("~/Dropbox/GDB_archive/CMapCorr_files/",temp_filename))
  }
  rm(list=grep("^temp",ls(),value=T))
}


