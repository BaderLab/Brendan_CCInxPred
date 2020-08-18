if (exists("lvl3_data")) {
} else if (file.exists("~/Dropbox/GDB_archive/CMapCorr_files/lvl3_inputs.RData")) {
  load("~/Dropbox/GDB_archive/CMapCorr_files/lvl3_inputs.RData") 
} else {
  source("lvl3_inputs.R")
}

if (exists("lvl4_data")) {
} else if (file.exists("~/Dropbox/GDB_archive/CMapCorr_files/lvl4_inputs.RData")) {
  load("~/Dropbox/GDB_archive/CMapCorr_files/lvl4_inputs.RData") 
} else {
  source("lvl4_inputs.R")
}


zscores_list <- pbapply::pbsapply(unique(lvl3_data@cdesc$rna_plate),function(PLATE) {
  sapply(unique(lvl3_data@cdesc[lvl3_data@cdesc$rna_plate == PLATE,"cell_id"]),function(CELL) {
    ctl_med <- apply(lvl3_data_ctl@mat[,lvl3_data_ctl@cdesc$rna_plate == PLATE & 
                                         lvl3_data_ctl@cdesc$cell_id == CELL],1,
                     median)
    ctl_mad <- apply(lvl3_data_ctl@mat[,lvl3_data_ctl@cdesc$rna_plate == PLATE & 
                                         lvl3_data_ctl@cdesc$cell_id == CELL],1,
                    mad)
    ctl_mad[ctl_mad < 0.1] <- 0.1
    return(t(
      sapply(1:nrow(lvl3_data@mat),function(X)
        (lvl3_data@mat[X,lvl3_data@cdesc$rna_plate == PLATE &
                         lvl3_data@cdesc$cell_id == CELL] -
           ctl_med[X]) / ctl_mad[X]
      )
    ))
  },simplify=F)
},simplify=F)
zscores <- do.call(cbind,unlist(zscores_list,recursive=F))
rownames(zscores) <- rownames(lvl3_data@mat)
zscores <- zscores[,colnames(lvl3_data@mat)]

lvl4new_data <- lvl4_data
lvl4new_data@mat <- zscores

# Lost a bunch of samples with no matching controls on the plate:
lvl4new_data@mat <- lvl4new_data@mat[,!apply(lvl4new_data@mat,2,function(X) any(is.na(X)))]
lvl4new_data@cdesc <- lvl4new_data@cdesc[colnames(lvl4new_data@mat),]

save(lvl4new_data,file="~/Dropbox/GDB_archive/CMapCorr_files/lvl4new.RData")
