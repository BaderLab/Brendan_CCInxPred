library(cmapR)
library(pbapply)
setwd("~/Dropbox/GDB/CMapCorr/")

jaccard_similarity_fx <- function(A,B) {
  length(intersect(A,B)) / length(union(A,B))
}

# lvl4 ----
if (exists("lvl4_data")) {
} else if (file.exists("../CMapCorr_files/lvl4_inputs.RData")) {
  load("../CMapCorr_files/lvl4_inputs.RData") 
} else {
  source("lvl4_inputs.R")
}


lvl4_de <- apply(lvl4_data@mat,2,function(X) names(X)[X > 1.645])


lvl4_jacc_ct <- pbsapply(unique(lvl4_data@cdesc$cell_id),function(CT) {
  temp <- combn(rownames(lvl4_data@cdesc)[lvl4_data@cdesc$cell_id == CT],2)
  mapply(jaccard_similarity_fx,A=lvl4_de[temp[1,]],B=lvl4_de[temp[2,]])
},simplify=F)


lvl4_jacc_lig <- pbsapply(unique(lvl4_data@cdesc$pert_iname),function(LIG) {
  temp <- combn(rownames(lvl4_data@cdesc)[lvl4_data@cdesc$pert_iname == LIG],2)
  mapply(jaccard_similarity_fx,A=lvl4_de[temp[1,]],B=lvl4_de[temp[2,]])
},simplify=F)


lvl4_jacc_lig_ct <- list()
for (LIG in unique(lvl4_data@cdesc$pert_iname)) {
  temp_cond <- unique(lvl4_data@cdesc[lvl4_data@cdesc$pert_iname == LIG,"cell_id"])
  temp_cond_id <- sapply(temp_cond,function(CT) 
    rownames(lvl4_data@cdesc)[
      lvl4_data@cdesc$pert_iname == LIG &
        lvl4_data@cdesc$cell_id == CT 
      ],simplify=F)
  lvl4_jacc_lig_ct[[LIG]] <- sapply(
    temp_cond_id[sapply(temp_cond_id,length) > 1],function(X) {
      temp <- combn(X,2)
      mapply(jaccard_similarity_fx,A=lvl4_de[temp[1,]],B=lvl4_de[temp[2,]])
    },simplify=F)
}  
lvl4_jacc_lig_ct <- lvl4_jacc_lig_ct[sapply(lvl4_jacc_lig_ct,length) > 0]
lvl4_jacc_lig_ct <- sapply(lvl4_jacc_lig_ct,unlist,simplify=F)


lvl4_jacc_lig_tx <- list()
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
  
  lvl4_jacc_lig_tx[[LIG]] <- sapply(
    temp_cond_id[sapply(temp_cond_id,length) > 1],function(X) {
      temp <- combn(X,2)
      mapply(jaccard_similarity_fx,A=lvl4_de[temp[1,]],B=lvl4_de[temp[2,]])
    },simplify=F)
}
lvl4_jacc_lig_tx <- lvl4_jacc_lig_tx[sapply(lvl4_jacc_lig_tx,length) > 0]
lvl4_jacc_lig_tx <- sapply(lvl4_jacc_lig_tx,unlist,simplify=F)

save(list=c("lvl4_de",ls()[grepl("lvl4",ls()) & grepl("jacc",ls())]),
     file="~/Dropbox/GDB/CMapCorr_files/200624_lvl4jacc.RData")
rm(list=c("LIG","CT",grep("^temp",ls(),value=T),grep("^lvl4_",ls(),value=T)))


# lvl5 ----
if (exists("lvl5_data")) {
} else if (file.exists("../CMapCorr_files/lvl5_inputs.RData")) {
  load("../CMapCorr_files/lvl5_inputs.RData") 
} else {
  source("lvl5_inputs.R")
}


lvl5_de <- apply(lvl5_data@mat,2,function(X) names(X)[X > 1.645])


lvl5_jacc_ct <- pbsapply(unique(lvl5_data@cdesc$cell_id),function(CT) {
  temp <- combn(rownames(lvl5_data@cdesc)[lvl5_data@cdesc$cell_id == CT],2)
  mapply(jaccard_similarity_fx,A=lvl5_de[temp[1,]],B=lvl5_de[temp[2,]])
},simplify=F)


lvl5_jacc_lig <- pbsapply(unique(lvl5_data@cdesc$pert_iname),function(LIG) {
  temp <- combn(rownames(lvl5_data@cdesc)[lvl5_data@cdesc$pert_iname == LIG],2)
  mapply(jaccard_similarity_fx,A=lvl5_de[temp[1,]],B=lvl5_de[temp[2,]])
},simplify=F)


lvl5_jacc_lig_ct <- list()
for (LIG in unique(lvl5_data@cdesc$pert_iname)) {
  temp_cond <- unique(lvl5_data@cdesc[lvl5_data@cdesc$pert_iname == LIG,"cell_id"])
  temp_cond_id <- sapply(temp_cond,function(CT) 
    rownames(lvl5_data@cdesc)[
      lvl5_data@cdesc$pert_iname == LIG &
        lvl5_data@cdesc$cell_id == CT 
      ],simplify=F)
  lvl5_jacc_lig_ct[[LIG]] <- sapply(
    temp_cond_id[sapply(temp_cond_id,length) > 1],function(X) {
      temp <- combn(X,2)
      mapply(jaccard_similarity_fx,A=lvl5_de[temp[1,]],B=lvl5_de[temp[2,]])
    },simplify=F)
}  
lvl5_jacc_lig_ct <- lvl5_jacc_lig_ct[sapply(lvl5_jacc_lig_ct,length) > 0]
lvl5_jacc_lig_ct <- sapply(lvl5_jacc_lig_ct,unlist,simplify=F)


lvl5_jacc_lig_tx <- list()
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
  
  lvl5_jacc_lig_tx[[LIG]] <- sapply(
    temp_cond_id[sapply(temp_cond_id,length) > 1],function(X) {
      temp <- combn(X,2)
      mapply(jaccard_similarity_fx,A=lvl5_de[temp[1,]],B=lvl5_de[temp[2,]])
    },simplify=F)
}
lvl5_jacc_lig_tx <- lvl5_jacc_lig_tx[sapply(lvl5_jacc_lig_tx,length) > 0]
lvl5_jacc_lig_tx <- sapply(lvl5_jacc_lig_tx,unlist,simplify=F)

save(list=c("lvl5_de",ls()[grepl("lvl5",ls()) & grepl("jacc",ls())]),
     file="~/Dropbox/GDB/CMapCorr_files/200624_lvl5jacc.RData")
rm(list=c("LIG","CT",grep("^temp",ls(),value=T),grep("^lvl5_",ls(),value=T)))



# lvl4new ----
if (exists("lvl4new_data")) {
} else if (file.exists("../CMapCorr_files/lvl4new.RData")) {
  load("../CMapCorr_files/lvl4new.RData") 
} else {
  source("200706_ZscoreFromAssayed.R")
}


lvl4new_de <- apply(lvl4new_data@mat,2,function(X) names(X)[X > 1.645])


lvl4new_jacc_ct <- pbsapply(unique(lvl4new_data@cdesc$cell_id),function(CT) {
  temp <- combn(rownames(lvl4new_data@cdesc)[lvl4new_data@cdesc$cell_id == CT],2)
  mapply(jaccard_similarity_fx,A=lvl4new_de[temp[1,]],B=lvl4new_de[temp[2,]])
},simplify=F)


lvl4new_jacc_lig <- pbsapply(unique(lvl4new_data@cdesc$pert_iname),function(LIG) {
  temp <- combn(rownames(lvl4new_data@cdesc)[lvl4new_data@cdesc$pert_iname == LIG],2)
  mapply(jaccard_similarity_fx,A=lvl4new_de[temp[1,]],B=lvl4new_de[temp[2,]])
},simplify=F)


lvl4new_jacc_lig_ct <- list()
for (LIG in unique(lvl4new_data@cdesc$pert_iname)) {
  temp_cond <- unique(lvl4new_data@cdesc[lvl4new_data@cdesc$pert_iname == LIG,"cell_id"])
  temp_cond_id <- sapply(temp_cond,function(CT) 
    rownames(lvl4new_data@cdesc)[
      lvl4new_data@cdesc$pert_iname == LIG &
        lvl4new_data@cdesc$cell_id == CT 
      ],simplify=F)
  lvl4new_jacc_lig_ct[[LIG]] <- sapply(
    temp_cond_id[sapply(temp_cond_id,length) > 1],function(X) {
      temp <- combn(X,2)
      mapply(jaccard_similarity_fx,A=lvl4new_de[temp[1,]],B=lvl4new_de[temp[2,]])
    },simplify=F)
}  
lvl4new_jacc_lig_ct <- lvl4new_jacc_lig_ct[sapply(lvl4new_jacc_lig_ct,length) > 0]
lvl4new_jacc_lig_ct <- sapply(lvl4new_jacc_lig_ct,unlist,simplify=F)


lvl4new_jacc_lig_tx <- list()
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
  
  lvl4new_jacc_lig_tx[[LIG]] <- sapply(
    temp_cond_id[sapply(temp_cond_id,length) > 1],function(X) {
      temp <- combn(X,2)
      mapply(jaccard_similarity_fx,A=lvl4new_de[temp[1,]],B=lvl4new_de[temp[2,]])
    },simplify=F)
}
lvl4new_jacc_lig_tx <- lvl4new_jacc_lig_tx[sapply(lvl4new_jacc_lig_tx,length) > 0]
lvl4new_jacc_lig_tx <- sapply(lvl4new_jacc_lig_tx,unlist,simplify=F)

save(list=c("lvl4new_de",ls()[grepl("lvl4new",ls()) & grepl("jacc",ls())]),
     file="~/Dropbox/GDB/CMapCorr_files/200706_lvl4newjacc.RData")
rm(list=c("LIG",grep("^temp",ls(),value=T),grep("^lvl4new_",ls(),value=T)))



