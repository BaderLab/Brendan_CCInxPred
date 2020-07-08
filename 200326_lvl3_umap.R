library(cmapR)
library(umap)

if (exists("lvl3_data")) {
} else if (file.exists("../CMapCorr_files/lvl3_inputs.RData")) {
  load("../CMapCorr_files/lvl3_inputs.RData") 
} else {
  source("lvl3_inputs.R")
}

# UMAP ----
if (file.exists("~/Dropbox/GDB/CMapCorr_files/200630_lvl3_umap.RData")) {
  load("~/Dropbox/GDB/CMapCorr_files/200630_lvl3_umap.RData") 
} else {
  temp_param <- umap.defaults
  temp_param$n_neighbors <- 30
  temp_param$metric <- "cosine"
  temp_param$min_dist <- 0.2
  temp_param$n_epochs <- 500
  lvl3_umap <- umap(t(cbind(lvl3_data@mat,lvl3_data_ctl@mat)),config=temp_param)
  save(lvl3_umap,file="~/Dropbox/GDB/CMapCorr_files/200630_lvl3_umap.RData")
}

# UMAP per CT ----
if (file.exists("~/Dropbox/GDB/CMapCorr_files/200630_lvl3_ctUMAP.RData")) {
  load("~/Dropbox/GDB/CMapCorr_files/200630_lvl3_ctUMAP.RData") 
} else {
  lvl3_ctUMAP <- list()
  temp_param <- umap.defaults
  temp_param$n_neighbors <- 30
  temp_param$metric <- "cosine"
  temp_param$min_dist <- 0.2
  temp_param$n_epochs <- 500
  for (CT in ct14) {
    message(paste0(which(CT == ct14),"/",length(ct14)))
    temp_tx <- rownames(lvl3_data@cdesc[lvl3_data@cdesc$cell_id == CT &
                                          lvl3_data@cdesc$pert_iname %in% lig16,])
    temp_plate <- unique(lvl3_data@cdesc[temp_tx,"rna_plate"])
    temp_ctrl <- rownames(lvl3_data_ctl@cdesc[lvl3_data_ctl@cdesc$cell_id == CT &
                                                lvl3_data_ctl@cdesc$rna_plate %in% temp_plate,])
    lvl3_ctUMAP[[CT]] <- umap(t(cbind(lvl3_data@mat[,temp_tx],
                                 lvl3_data_ctl@mat[,temp_ctrl])),
                         config=temp_param)
  }
  save(lvl3_ctUMAP,file="~/Dropbox/GDB/CMapCorr_files/200630_lvl3_ctUMAP.RData")
}


# UMAP per LIG ----
if (file.exists("~/Dropbox/GDB/CMapCorr_files/200630_lvl3_ligUMAP.RData")) {
  load("~/Dropbox/GDB/CMapCorr_files/200630_lvl3_ligUMAP.RData") 
} else {
  lvl3_ligUMAP <- list()
  temp_param <- umap.defaults
  temp_param$n_neighbors <- 30
  temp_param$metric <- "cosine"
  temp_param$min_dist <- 0.2
  temp_param$n_epochs <- 500
  for (LIG in lig16) {
    message(paste0(which(LIG == lig16),"/",length(lig16)))
    temp_tx <- rownames(lvl3_data@cdesc[lvl3_data@cdesc$pert_iname == LIG,])
    temp_plate <- unique(lvl3_data@cdesc[temp_tx,"rna_plate"])
    temp_ctrl <- rownames(lvl3_data_ctl@cdesc[lvl3_data_ctl@cdesc$rna_plate %in% temp_plate,])
    lvl3_ligUMAP[[LIG]] <- umap(t(cbind(lvl3_data@mat[,temp_tx],
                                      lvl3_data_ctl@mat[,temp_ctrl])),
                              config=temp_param)
  }
  save(lvl3_ligUMAP,file="~/Dropbox/GDB/CMapCorr_files/200630_lvl3_ligUMAP.RData")
}
