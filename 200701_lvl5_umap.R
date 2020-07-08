library(cmapR)
library(umap)

setwd("~/Dropbox/GDB/CMapCorr/")
if (exists("lvl5_data")) {
} else if (file.exists("../CMapCorr_files/lvl5_inputs.RData")) {
  load("../CMapCorr_files/lvl5_inputs.RData") 
} else {
  source("lvl5_inputs.R")
}

# Not including controls in UMAP data because they're meaningless

# UMAP ----
if (file.exists("~/Dropbox/GDB/CMapCorr_files/200701_lvl5_umap.RData")) {
  load("~/Dropbox/GDB/CMapCorr_files/200701_lvl5_umap.RData") 
} else {
  temp_param <- umap.defaults
  temp_param$n_neighbors <- 30
  temp_param$metric <- "cosine"
  temp_param$min_dist <- 0.2
  temp_param$n_epochs <- 500
  
  lvl5_umap <- umap(t(lvl5_data@mat),config=temp_param)
  save(lvl5_umap,file="~/Dropbox/GDB/CMapCorr_files/200701_lvl5_umap.RData")
}

if (file.exists("~/Dropbox/GDB/CMapCorr_files/200701_lvl5_umap_lig16.RData")) {
  load("~/Dropbox/GDB/CMapCorr_files/200701_lvl5_umap_lig16.RData") 
} else {
  temp_param <- umap.defaults
  temp_param$n_neighbors <- 30
  temp_param$metric <- "cosine"
  temp_param$min_dist <- 0.2
  temp_param$n_epochs <- 500
  
  lvl5_umap_lig16 <- umap(t(lvl5_data@mat[,lvl5_data@cdesc$pert_iname %in% lig16]),config=temp_param)
  save(lvl5_umap_lig16,file="~/Dropbox/GDB/CMapCorr_files/200701_lvl5_umap_lig16.RData")
}


# UMAP by CT ----
if (file.exists("~/Dropbox/GDB/CMapCorr_files/200701_lvl5_ctUMAP.RData")) {
  load("~/Dropbox/GDB/CMapCorr_files/200701_lvl5_ctUMAP.RData") 
} else {
  temp_param <- umap.defaults
  temp_param$n_neighbors <- 5
  temp_param$metric <- "cosine"
  temp_param$min_dist <- 0.2
  temp_param$n_epochs <- 500
  
  lvl5_ctUMAP <- list()
  for (CT in ct14) {
    message(paste0(which(CT == ct14),"/",length(ct14)))
    temp_tx <- rownames(lvl5_data@cdesc[lvl5_data@cdesc$cell_id == CT &
                                          lvl5_data@cdesc$pert_iname %in% lig16,])
    lvl5_ctUMAP[[CT]] <- umap(t(lvl5_data@mat[,temp_tx]),config=temp_param)
  }
  save(lvl5_ctUMAP,file="~/Dropbox/GDB/CMapCorr_files/200701_lvl5_ctUMAP.RData")
}


# UMAP per LIG ----
if (file.exists("~/Dropbox/GDB/CMapCorr_files/200701_lvl5_ligUMAP.RData")) {
  load("~/Dropbox/GDB/CMapCorr_files/200701_lvl5_ligUMAP.RData") 
} else {
  lvl5_ligUMAP <- list()
  temp_param <- umap.defaults
  temp_param$n_neighbors <- 5
  temp_param$metric <- "cosine"
  temp_param$min_dist <- 0.2
  temp_param$n_epochs <- 500
  for (LIG in lig16) {
    message(paste0(which(LIG == lig16),"/",length(lig16)))
    temp_tx <- rownames(lvl5_data@cdesc[lvl5_data@cdesc$pert_iname == LIG,])
    lvl5_ligUMAP[[LIG]] <- umap(t(lvl5_data@mat[,temp_tx]),config=temp_param)
  }
  save(lvl5_ligUMAP,file="~/Dropbox/GDB/CMapCorr_files/200701_lvl5_ligUMAP.RData")
}
