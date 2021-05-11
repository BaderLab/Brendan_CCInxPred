library(pbapply)
pboptions(type="timer")


pwCOR <- function(MAT,IDs) {
  pbsapply(IDs,function(X) {
    temp <- cor(MAT[,X],method="spearman",use="complete.obs")
    return(temp[upper.tri(temp)])
  },cl=8,simplify=F)
}

load("~/Dropbox/GDB_archive/CMapCorr_files/lvl4_inputs.RData")
temp <- list.files("~/Dropbox/GDB_archive/CMapCorr_files/","^GSEA_lvl4_")
for (X in temp) {
  temp2 <- load(paste0("~/Dropbox/GDB_archive/CMapCorr_files/",X))
  assign(sub(".RData","",X,fixed=T),sapply(get(temp2)[-1],"[[","NES"))
  rm(list=temp2)
}

temp_ligcountsperct <- sapply(unique(lvl4_data@cdesc$cell_id),function(CT)
  length(unique(lvl4_data@cdesc[lvl4_data@cdesc$cell_id == CT,"pert_iname"])))
ct9 <- names(temp_ligcountsperct)[temp_ligcountsperct > 100]
names(ct9) <- sapply(ct9,function(X) names(ct14)[ct14 == X])

temp_ctcountsperlig <- sapply(unique(lvl4_data@cdesc$pert_iname),function(LIG)
  length(unique(lvl4_data@cdesc[lvl4_data@cdesc$pert_iname == LIG &
                                  lvl4_data@cdesc$cell_id %in% ct9,"cell_id"])))
lig295 <- names(temp_ctcountsperlig)[temp_ctcountsperlig == 9]
lig295 <- sort(lig295)

rm(list=c("X",grep("^temp",ls(),value=T)))



# LIG ----
print("LIG ------------------")
DAT <- "lvl4_data"
temp <- sapply(lig295,function(LIG) 
  rownames(get(DAT)@cdesc)[get(DAT)@cdesc$pert_iname == LIG & 
                             get(DAT)@cdesc$cell_id %in% ct9],
  simplify=F)
for (GSEA in grep(paste0("GSEA_",sub("_data","",DAT)),ls(),value=T)) {
  print(paste0("-",GSEA," ----"))
  CORS <- pwCOR(get(GSEA),temp)
  temp_filename <- paste0("corrGSEA_",sub("GSEA_lvl[45]_","",GSEA),"_lig.RData")
  save(CORS,file=paste0("~/Dropbox/GDB_archive/CMapCorr_files/",temp_filename))
}
rm(list=grep("^temp",ls(),value=T))


# CT ----
print("CT ------------------")
DAT <- "lvl4_data"
temp <- sapply(ct9,function(CT) 
  rownames(get(DAT)@cdesc)[get(DAT)@cdesc$pert_iname %in% lig295 & 
                             get(DAT)@cdesc$cell_id == CT],
  simplify=F)
for (GSEA in grep(paste0("GSEA_",sub("_data","",DAT)),ls(),value=T)) {
  print(paste0("-",GSEA," ----"))
  CORS <- pwCOR(get(GSEA),temp)
  temp_filename <- paste0("corrGSEA_",sub("GSEA_lvl[45]_","",GSEA),"_ct.RData")
  save(CORS,file=paste0("~/Dropbox/GDB_archive/CMapCorr_files/",temp_filename))
}
rm(list=grep("^temp",ls(),value=T))


# LigCT ----
print("LIGCT ------------------")
DAT <- "lvl4_data"
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
for (GSEA in grep(paste0("GSEA_",sub("_data","",DAT)),ls(),value=T)) {
  print(paste0("-",GSEA," ----"))
  CORS <- pwCOR(get(GSEA),temp)
  temp_filename <- paste0("corrGSEA_",sub("GSEA_lvl[45]_","",GSEA),"_ligct.RData")
  save(CORS,file=paste0("~/Dropbox/GDB_archive/CMapCorr_files/",temp_filename))
}
rm(list=grep("^temp",ls(),value=T))



# REP ----
print("REP ------------------")
DAT <- "lvl4_data"
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
for (GSEA in grep(paste0("GSEA_",sub("_data","",DAT)),ls(),value=T)) {
  print(paste0("-",GSEA," ----"))
  CORS <- pwCOR(get(GSEA),temp)
  temp_filename <- paste0("corrGSEA_",sub("GSEA_lvl[45]_","",GSEA),"_rep.RData")
  save(CORS,file=paste0("~/Dropbox/GDB_archive/CMapCorr_files/",temp_filename))
}
rm(list=grep("^temp",ls(),value=T))



