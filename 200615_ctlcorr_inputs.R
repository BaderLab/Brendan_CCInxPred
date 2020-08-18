library(cmapR)
setwd("~/Dropbox/GDB/CMapCorr/")

# lvl3 ----
if (exists("lvl3_data")) {
} else if (file.exists("~/Dropbox/GDB_archive/CMapCorr_files/lvl3_inputs.RData")) {
  load("~/Dropbox/GDB_archive/CMapCorr_files/lvl3_inputs.RData") 
} else {
  source("lvl3_inputs.R")
}

lvl3_ctl_corr <- list()
for (CT in ct14) {
  temp_cond <- unique(lvl3_data_ctl@cdesc[lvl3_data_ctl@cdesc$cell_id == CT,c("pert_iname","pert_time")])
  temp_cond_str <- apply(temp_cond,1,function(X) paste0(paste(X,collapse="_"),"hr"))
  temp_cond_id <- apply(temp_cond,1,
                        function(X) rownames(lvl3_data_ctl@cdesc)[
                          lvl3_data_ctl@cdesc$cell_id == CT &
                            lvl3_data_ctl@cdesc$pert_iname == X[1] &
                            lvl3_data_ctl@cdesc$pert_time == X[2]
                          ])
  names(temp_cond_id) <- temp_cond_str
  temp_cond_id[[paste0(CT,"_all")]] <- unlist(temp_cond_id)
  temp_cond_id <- temp_cond_id[c(paste0(CT,"_all"),
                                 names(temp_cond_id)[names(temp_cond_id) != paste0(CT,"_all")])]
  
  print(paste0(CT,":"))
  lvl3_ctl_corr[[CT]] <- pbapply::pbsapply(temp_cond_id,function(X) {
    temp <- cor(lvl3_data_ctl@mat[,X])
    return(temp[lower.tri(temp)])
  },simplify=F)
}


lvl3_ctl_corr_plate <- apply(unique(lvl3_data_ctl@cdesc[,c("rna_plate","cell_id")]),1,function(X) {
  temp <- cor(lvl3_data_ctl@mat[,lvl3_data_ctl@cdesc$rna_plate == X[1] &
                                  lvl3_data_ctl@cdesc$cell_id == X[2]])
  return(temp[lower.tri(temp)])
})


save(lvl3_ctl_corr,lvl3_ctl_corr_plate,
     file="~/Dropbox/GDB/CMapCorr_files/200615_lvl3ctlcorr.RData")
rm(list=c("CT",grep("^temp",ls(),value=T),"lvl3_ctl_corr"))



lvl3_lig_corr <- pbapply::pbsapply(unique(lvl3_data@cdesc$pert_iname),function(LIG) {
  temp <- cor(lvl3_data@mat[,lvl3_data@cdesc$pert_iname == LIG])
  return(temp[lower.tri(temp)])
})


lvl3_ct_corr <- pbapply::pbsapply(unique(lvl3_data@cdesc$cell_id),function(CT) {
  temp <- cor(lvl3_data@mat[,lvl3_data@cdesc$cell_id == CT])
  return(temp[lower.tri(temp)])
})


lvl3_lig_corr_ct <- list()
for (LIG in unique(lvl3_data@cdesc$pert_iname)) {
  temp_cond <- unique(lvl3_data@cdesc[lvl3_data@cdesc$pert_iname == LIG,"cell_id"])
  temp_cond_id <- sapply(temp_cond,function(CT) 
    rownames(lvl3_data@cdesc)[
      lvl3_data@cdesc$pert_iname == LIG &
        lvl3_data@cdesc$cell_id == CT 
      ],simplify=F)
  
  lvl3_lig_corr_ct[[LIG]] <- sapply(temp_cond_id[sapply(temp_cond_id,length) > 1],
                                    function(X) {
                                      temp <- cor(lvl3_data@mat[,X])
                                      return(temp[lower.tri(temp)])
                                    },simplify=F)
}
lvl3_lig_corr_ct <- lvl3_lig_corr_ct[sapply(lvl3_lig_corr_ct,length) > 0]
lvl3_lig_corr_ct <- sapply(lvl3_lig_corr_ct,unlist,simplify=F)


lvl3_lig_corr_tx <- list()
for (LIG in unique(lvl3_data@cdesc$pert_iname)) {
  temp_cond <- unique(lvl3_data@cdesc[lvl3_data@cdesc$pert_iname == LIG,c("cell_id","pert_dose","pert_time")])
  temp_cond_str <- apply(temp_cond,1,function(X) paste0(X[1],"_",X[2],"ng/mL_",X[3],"hr"))
  temp_cond_id <- sapply(1:nrow(temp_cond),function(X) 
    rownames(lvl3_data@cdesc)[
      lvl3_data@cdesc$pert_iname == LIG &
        lvl3_data@cdesc$cell_id == temp_cond[X,1] &
        lvl3_data@cdesc$pert_dose == temp_cond[X,2] &
        lvl3_data@cdesc$pert_time == temp_cond[X,3]
      ],simplify=F)
  names(temp_cond_id) <- temp_cond_str
  
  lvl3_lig_corr_tx[[LIG]] <- sapply(temp_cond_id[sapply(temp_cond_id,length) > 1],function(X) {
    temp <- cor(lvl3_data@mat[,X])
    return(temp[lower.tri(temp)])
  },simplify=F)
}
lvl3_lig_corr_tx <- lvl3_lig_corr_tx[sapply(lvl3_lig_corr_tx,length) > 0]
lvl3_lig_corr_tx <- sapply(lvl3_lig_corr_tx,unlist,simplify=F)


save(lvl3_ct_corr,lvl3_lig_corr,lvl3_lig_corr_ct,lvl3_lig_corr_tx,
     file="~/Dropbox/GDB/CMapCorr_files/200617_lvl3corr.RData")
rm(list=c("LIG",grep("^temp",ls(),value=T),grep("^lvl3_(lig|ct)",ls(),value=T)))


# lvl4 ----
if (exists("lvl4_data")) {
} else if (file.exists("~/Dropbox/GDB_archive/CMapCorr_files/lvl4_inputs.RData")) {
  load("~/Dropbox/GDB_archive/CMapCorr_files/lvl4_inputs.RData") 
} else {
  source("lvl4_inputs.R")
}

lvl4_ctl_corr_plate <- apply(unique(lvl4_data_ctl@cdesc[,c("rna_plate","cell_id")]),1,function(X) {
  temp <- cor(lvl4_data_ctl@mat[,lvl4_data_ctl@cdesc$rna_plate == X[1] &
                                  lvl4_data_ctl@cdesc$cell_id == X[2]],
              method="spearman")
  return(temp[lower.tri(temp)])
})


lvl4_lig_corr <- pbapply::pbsapply(unique(lvl4_data@cdesc$pert_iname),function(LIG) {
  temp <- cor(lvl4_data@mat[,lvl4_data@cdesc$pert_iname == LIG],
              method="spearman")
  return(temp[lower.tri(temp)])
})


lvl4_ct_corr <- pbapply::pbsapply(unique(lvl4_data@cdesc$cell_id),function(CT) {
  temp <- cor(lvl4_data@mat[,lvl4_data@cdesc$cell_id == CT],
              method="spearman")
  return(temp[lower.tri(temp)])
})


lvl4_lig_corr_ct <- list()
for (LIG in unique(lvl4_data@cdesc$pert_iname)) {
  temp_cond <- unique(lvl4_data@cdesc[lvl4_data@cdesc$pert_iname == LIG,"cell_id"])
  temp_cond_id <- sapply(temp_cond,function(CT) 
    rownames(lvl4_data@cdesc)[
      lvl4_data@cdesc$pert_iname == LIG &
        lvl4_data@cdesc$cell_id == CT 
      ],simplify=F)
  
  lvl4_lig_corr_ct[[LIG]] <- sapply(temp_cond_id[sapply(temp_cond_id,length) > 1],
                                    function(X) {
                                      temp <- cor(lvl4_data@mat[,X],
                                                  method="spearman")
                                      return(temp[lower.tri(temp)])
                                    },simplify=F)
}
lvl4_lig_corr_ct <- lvl4_lig_corr_ct[sapply(lvl4_lig_corr_ct,length) > 0]
lvl4_lig_corr_ct <- sapply(lvl4_lig_corr_ct,unlist,simplify=F)


lvl4_lig_corr_tx <- list()
for (LIG in unique(lvl4_data@cdesc$pert_iname)) {
  temp_cond <- unique(lvl4_data@cdesc[lvl4_data@cdesc$pert_iname == LIG,c("cell_id","pert_dose","pert_time")])
  temp_cond_str <- apply(temp_cond,1,function(X) paste0(X[1],"_",X[2],"ng/mL_",X[3],"hr"))
  temp_cond_id <- sapply(1:nrow(temp_cond),function(X) 
    rownames(lvl4_data@cdesc)[
      lvl4_data@cdesc$pert_iname == LIG &
        lvl4_data@cdesc$cell_id == temp_cond[X,1] &
        lvl4_data@cdesc$pert_dose == temp_cond[X,2] &
        lvl4_data@cdesc$pert_time == temp_cond[X,3]
      ],simplify=F)
  names(temp_cond_id) <- temp_cond_str
  
  lvl4_lig_corr_tx[[LIG]] <- sapply(temp_cond_id[sapply(temp_cond_id,length) > 1],function(X) {
    temp <- cor(lvl4_data@mat[,X],
                method="spearman")
    return(temp[lower.tri(temp)])
  },simplify=F)
}
lvl4_lig_corr_tx <- lvl4_lig_corr_tx[sapply(lvl4_lig_corr_tx,length) > 0]
lvl4_lig_corr_tx <- sapply(lvl4_lig_corr_tx,unlist,simplify=F)


save(list=ls()[grepl("lvl4",ls()) & grepl("corr",ls())],
     file="~/Dropbox/GDB_archive/CMapCorr_files/200625_lvl4corr.RData")
rm(list=c("LIG",grep("^temp",ls(),value=T),grep("^lvl4_(lig|ct)",ls(),value=T)))


# lvl5 ----
if (exists("lvl5_data")) {
} else if (file.exists("~/Dropbox/GDB_archive/CMapCorr_files/lvl5_inputs.RData")) {
  load("~/Dropbox/GDB_archive/CMapCorr_files/lvl5_inputs.RData") 
} else {
  source("lvl5_inputs.R")
}


lvl5_ct_corr <- pbapply::pbsapply(unique(lvl5_data@cdesc$cell_id),function(CT) {
  temp <- cor(lvl5_data@mat[,lvl5_data@cdesc$cell_id == CT],
              method="spearman")
  return(temp[lower.tri(temp)])
})


lvl5_lig_corr <- pbapply::pbsapply(unique(lvl5_data@cdesc$pert_iname),function(LIG) {
  temp <- cor(lvl5_data@mat[,lvl5_data@cdesc$pert_iname == LIG],method="spearman")
  return(temp[lower.tri(temp)])
})


lvl5_lig_corr_ct <- list()
for (LIG in unique(lvl5_data@cdesc$pert_iname)) {
  temp_cond <- unique(lvl5_data@cdesc[lvl5_data@cdesc$pert_iname == LIG,"cell_id"])
  temp_cond_id <- sapply(temp_cond,function(CT) 
    rownames(lvl5_data@cdesc)[
      lvl5_data@cdesc$pert_iname == LIG &
        lvl5_data@cdesc$cell_id == CT 
      ],simplify=F)
  
  lvl5_lig_corr_ct[[LIG]] <- sapply(temp_cond_id[sapply(temp_cond_id,length) > 1],function(X) {
    temp <- cor(lvl5_data@mat[,X],method="spearman")
    return(temp[lower.tri(temp)])
  },simplify=F)
}
lvl5_lig_corr_ct <- lvl5_lig_corr_ct[sapply(lvl5_lig_corr_ct,length) > 0]
lvl5_lig_corr_ct <- sapply(lvl5_lig_corr_ct,unlist,simplify=F)


lvl5_lig_corr_tx <- list()
for (LIG in unique(lvl5_data@cdesc$pert_iname)) {
  temp_cond <- unique(lvl5_data@cdesc[lvl5_data@cdesc$pert_iname == LIG,c("cell_id","pert_dose","pert_time")])
  temp_cond_str <- apply(temp_cond,1,function(X) paste0(X[1],"_",X[2],"ng/mL_",X[3],"hr"))
  temp_cond_id <- sapply(1:nrow(temp_cond),function(X) 
    rownames(lvl5_data@cdesc)[
      lvl5_data@cdesc$pert_iname == LIG &
        lvl5_data@cdesc$cell_id == temp_cond[X,1] &
        lvl5_data@cdesc$pert_dose == temp_cond[X,2] &
        lvl5_data@cdesc$pert_time == temp_cond[X,3]
      ],simplify=F)
  names(temp_cond_id) <- temp_cond_str
  
  lvl5_lig_corr_tx[[LIG]] <- sapply(temp_cond_id[sapply(temp_cond_id,length) > 1],function(X) {
    temp <- cor(lvl5_data@mat[,X],method="spearman")
    return(temp[lower.tri(temp)])
  },simplify=F)
}
lvl5_lig_corr_tx <- lvl5_lig_corr_tx[sapply(lvl5_lig_corr_tx,length) > 0]
lvl5_lig_corr_tx <- sapply(lvl5_lig_corr_tx,unlist,simplify=F)


save(list=ls()[grepl("lvl5",ls()) & grepl("corr",ls())],
     file="~/Dropbox/GDB/CMapCorr_files/200617_lvl5corr.RData")


# lvl4new ----
if (exists("lvl4new_data")) {
} else if (file.exists("~/Dropbox/GDB_archive/CMapCorr_files/lvl4new.RData")) {
  load("~/Dropbox/GDB_archive/CMapCorr_files/lvl4new.RData") 
} else {
  source("200706_ZscoreFromAssayed.R")
}


lvl4new_lig_corr <- pbapply::pbsapply(unique(lvl4new_data@cdesc$pert_iname),function(LIG) {
  temp <- cor(lvl4new_data@mat[,lvl4new_data@cdesc$pert_iname == LIG],
              method="spearman")
  return(temp[lower.tri(temp)])
})


lvl4new_ct_corr <- pbapply::pbsapply(unique(lvl4new_data@cdesc$cell_id),function(CT) {
  temp <- cor(lvl4new_data@mat[,lvl4new_data@cdesc$cell_id == CT],
              method="spearman")
  return(temp[lower.tri(temp)])
})


lvl4new_lig_corr_ct <- list()
for (LIG in unique(lvl4new_data@cdesc$pert_iname)) {
  temp_cond <- unique(lvl4new_data@cdesc[lvl4new_data@cdesc$pert_iname == LIG,"cell_id"])
  temp_cond_id <- sapply(temp_cond,function(CT) 
    rownames(lvl4new_data@cdesc)[
      lvl4new_data@cdesc$pert_iname == LIG &
        lvl4new_data@cdesc$cell_id == CT 
      ],simplify=F)
  
  lvl4new_lig_corr_ct[[LIG]] <- sapply(temp_cond_id[sapply(temp_cond_id,length) > 1],
                                    function(X) {
                                      temp <- cor(lvl4new_data@mat[,X],
                                                  method="spearman")
                                      return(temp[lower.tri(temp)])
                                    },simplify=F)
}
lvl4new_lig_corr_ct <- lvl4new_lig_corr_ct[sapply(lvl4new_lig_corr_ct,length) > 0]
lvl4new_lig_corr_ct <- sapply(lvl4new_lig_corr_ct,unlist,simplify=F)


lvl4new_lig_corr_tx <- list()
for (LIG in unique(lvl4new_data@cdesc$pert_iname)) {
  temp_cond <- unique(lvl4new_data@cdesc[lvl4new_data@cdesc$pert_iname == LIG,c("cell_id","pert_dose","pert_time")])
  temp_cond_str <- apply(temp_cond,1,function(X) paste0(X[1],"_",X[2],"ng/mL_",X[3],"hr"))
  temp_cond_id <- sapply(1:nrow(temp_cond),function(X) 
    rownames(lvl4new_data@cdesc)[
      lvl4new_data@cdesc$pert_iname == LIG &
        lvl4new_data@cdesc$cell_id == temp_cond[X,1] &
        lvl4new_data@cdesc$pert_dose == temp_cond[X,2] &
        lvl4new_data@cdesc$pert_time == temp_cond[X,3]
      ],simplify=F)
  names(temp_cond_id) <- temp_cond_str
  
  lvl4new_lig_corr_tx[[LIG]] <- sapply(temp_cond_id[sapply(temp_cond_id,length) > 1],function(X) {
    temp <- cor(lvl4new_data@mat[,X],
                method="spearman")
    return(temp[lower.tri(temp)])
  },simplify=F)
}
lvl4new_lig_corr_tx <- lvl4new_lig_corr_tx[sapply(lvl4new_lig_corr_tx,length) > 0]
lvl4new_lig_corr_tx <- sapply(lvl4new_lig_corr_tx,unlist,simplify=F)


save(list=ls()[grepl("lvl4new",ls()) & grepl("corr",ls())],
     file="~/Dropbox/GDB_archive/CMapCorr_files/200706_lvl4newcorr.RData")
rm(list=c("LIG",grep("^temp",ls(),value=T),grep("^lvl4new_(lig|ct)",ls(),value=T)))


