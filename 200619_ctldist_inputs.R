library(cmapR)
setwd("~/Dropbox/GDB/CMapCorr/")

# lvl3 ----
if (exists("lvl3_data")) {
} else if (file.exists("~/Dropbox/GDB_archive/CMapCorr_files/lvl3_inputs.RData")) {
  load("~/Dropbox/GDB_archive/CMapCorr_files/lvl3_inputs.RData") 
} else {
  source("lvl3_inputs.R")
}

lvl3_ctl_dist <- list()
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
  lvl3_ctl_dist[[CT]] <- pbapply::pbsapply(temp_cond_id,function(X) 
    as.numeric(dist(t(lvl3_data_ctl@mat[,X]))),
    simplify=F)
}

save(lvl3_ctl_dist,file="~/Dropbox/GDB/CMapCorr_files/200615_lvl3ctldist.RData")
rm(list=c("CT",grep("^temp",ls(),value=T),"lvl3_ctl_dist"))


lvl3_lig_dist <- pbapply::pbsapply(unique(lvl3_data@cdesc$pert_iname),function(LIG) 
  as.numeric(dist(t(lvl3_data@mat[,lvl3_data@cdesc$pert_iname == LIG]))))


lvl3_ct_dist <- pbapply::pbsapply(unique(lvl3_data@cdesc$cell_id),function(CT) 
  as.numeric(dist(t(lvl3_data@mat[,lvl3_data@cdesc$cell_id == CT]))))


lvl3_lig_dist_ct <- list()
for (LIG in unique(lvl3_data@cdesc$pert_iname)) {
  temp_cond <- unique(lvl3_data@cdesc[lvl3_data@cdesc$pert_iname == LIG,"cell_id"])
  temp_cond_id <- sapply(temp_cond,function(CT) 
    rownames(lvl3_data@cdesc)[
      lvl3_data@cdesc$pert_iname == LIG &
        lvl3_data@cdesc$cell_id == CT 
      ],simplify=F)
  
  lvl3_lig_dist_ct[[LIG]] <- sapply(temp_cond_id[sapply(temp_cond_id,length) > 1],
                                    function(X) as.numeric(dist(t(lvl3_data@mat[,X]))),
                                    simplify=F)
}
lvl3_lig_dist_ct <- lvl3_lig_dist_ct[sapply(lvl3_lig_dist_ct,length) > 0]
lvl3_lig_dist_ct <- sapply(lvl3_lig_dist_ct,unlist,simplify=F)


lvl3_lig_dist_tx <- list()
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
  
  lvl3_lig_dist_tx[[LIG]] <- sapply(temp_cond_id[sapply(temp_cond_id,length) > 1],
                                    function(X) as.numeric(dist(t(lvl3_data@mat[,X]))),
                                    simplify=F)
}
lvl3_lig_dist_tx <- lvl3_lig_dist_tx[sapply(lvl3_lig_dist_tx,length) > 0]
lvl3_lig_dist_tx <- sapply(lvl3_lig_dist_tx,unlist,simplify=F)


save(lvl3_ct_dist,lvl3_lig_dist,lvl3_lig_dist_ct,lvl3_lig_dist_tx,
     file="~/Dropbox/GDB/CMapCorr_files/200617_lvl3dist.RData")
rm(list=c("LIG",grep("^temp",ls(),value=T),grep("^lvl3_(lig|ct)",ls(),value=T)))


