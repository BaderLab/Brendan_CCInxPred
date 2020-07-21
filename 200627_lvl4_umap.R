library(cmapR)
library(umap)

setwd("~/Dropbox/GDB/CMapCorr/")
if (exists("lvl4_data")) {
} else if (file.exists("../CMapCorr_files/lvl4_inputs.RData")) {
  load("../CMapCorr_files/lvl4_inputs.RData") 
} else {
  source("lvl4_inputs.R")
}

# Not including controls in UMAP data because they're meaningless

# UMAP ----
if (file.exists("~/Dropbox/GDB/CMapCorr_files/200630_lvl4_umap.RData")) {
  load("~/Dropbox/GDB/CMapCorr_files/200630_lvl4_umap.RData") 
} else {
  temp_param <- umap.defaults
  temp_param$n_neighbors <- 30
  temp_param$metric <- "cosine"
  temp_param$min_dist <- 0.2
  temp_param$n_epochs <- 500

  lvl4_umap <- umap(t(lvl4_data@mat),config=temp_param)
  save(lvl4_umap,file="~/Dropbox/GDB/CMapCorr_files/200630_lvl4_umap.RData")
}

if (file.exists("~/Dropbox/GDB/CMapCorr_files/200630_lvl4_umap_lig16.RData")) {
  load("~/Dropbox/GDB/CMapCorr_files/200630_lvl4_umap_lig16.RData") 
} else {
  temp_param <- umap.defaults
  temp_param$n_neighbors <- 30
  temp_param$metric <- "cosine"
  temp_param$min_dist <- 0.2
  temp_param$n_epochs <- 500
  
  lvl4_umap_lig16 <- umap(t(lvl4_data@mat[,lvl4_data@cdesc$pert_iname %in% lig16]),config=temp_param)
  save(lvl4_umap_lig16,file="~/Dropbox/GDB/CMapCorr_files/200630_lvl4_umap_lig16.RData")
}


# UMAP by CT ----
if (file.exists("~/Dropbox/GDB/CMapCorr_files/200630_lvl4_ctUMAP.RData")) {
  load("~/Dropbox/GDB/CMapCorr_files/200630_lvl4_ctUMAP.RData") 
} else {
  temp_param <- umap.defaults
  temp_param$n_neighbors <- 30
  temp_param$metric <- "cosine"
  temp_param$min_dist <- 0.2
  temp_param$n_epochs <- 500
  
  lvl4_ctUMAP <- list()
  for (CT in ct14) {
    message(paste0(which(CT == ct14),"/",length(ct14)))
    temp_tx <- rownames(lvl4_data@cdesc[lvl4_data@cdesc$cell_id == CT &
                                          lvl4_data@cdesc$pert_iname %in% lig16,])
    if (temp_tx < 1000) {
      temp_param$n_neighbors <- 5
    } else {
      temp_param$n_neighbors <- 30
    }
    lvl4_ctUMAP[[CT]] <- umap(t(lvl4_data@mat[,temp_tx]),config=temp_param)
  }
  save(lvl4_ctUMAP,file="~/Dropbox/GDB/CMapCorr_files/200630_lvl4_ctUMAP.RData")
}


# UMAP per LIG ----
if (file.exists("~/Dropbox/GDB/CMapCorr_files/200630_lvl4_ligUMAP.RData")) {
  load("~/Dropbox/GDB/CMapCorr_files/200630_lvl4_ligUMAP.RData") 
} else {
  lvl4_ligUMAP <- list()
  temp_param <- umap.defaults
  temp_param$n_neighbors <- 30
  temp_param$metric <- "cosine"
  temp_param$min_dist <- 0.2
  temp_param$n_epochs <- 500
  for (LIG in lig16) {
    message(paste0(which(LIG == lig16),"/",length(lig16)))
    temp_tx <- rownames(lvl4_data@cdesc[lvl4_data@cdesc$pert_iname == LIG,])
    lvl4_ligUMAP[[LIG]] <- umap(t(lvl4_data@mat[,temp_tx]),config=temp_param)
  }
  save(lvl4_ligUMAP,file="~/Dropbox/GDB/CMapCorr_files/200630_lvl4_ligUMAP.RData")
}

# WITH ALL GENES ----
rm(list=ls())
if (exists("lvl4_data_all")) {
} else if (file.exists("../CMapCorr_files/lvl4_inputs_allgenes.RData")) {
  load("../CMapCorr_files/lvl4_inputs_allgenes.RData") 
} else {
  source("lvl4_inputs_allgenes.R")
}

if (file.exists("~/Dropbox/GDB/CMapCorr_files/200721_lvl4all_umap_lig16.RData")) {
  load("~/Dropbox/GDB/CMapCorr_files/200721_lvl4all_umap_lig16.RData") 
} else {
  temp_param <- umap.defaults
  temp_param$n_neighbors <- 30
  temp_param$metric <- "cosine"
  temp_param$min_dist <- 0.2
  temp_param$n_epochs <- 500
  
  lvl4all_umap_lig16 <- umap(t(lvl4_data_all@mat[,lvl4_data_all@cdesc$pert_iname %in% lig16]),config=temp_param)
  save(lvl4all_umap_lig16,file="~/Dropbox/GDB/CMapCorr_files/200721_lvl4all_umap_lig16.RData")
}


# UMAP by CT ----
if (file.exists("~/Dropbox/GDB/CMapCorr_files/200721_lvl4all_ctUMAP.RData")) {
  load("~/Dropbox/GDB/CMapCorr_files/200721_lvl4all_ctUMAP.RData") 
} else {
  temp_param <- umap.defaults
  temp_param$n_neighbors <- 30
  temp_param$metric <- "cosine"
  temp_param$min_dist <- 0.2
  temp_param$n_epochs <- 500
  
  lvl4_ctUMAP <- list()
  for (CT in ct14) {
    message(paste0(which(CT == ct14),"/",length(ct14)))
    temp_tx <- rownames(lvl4_data_all@cdesc[lvl4_data_all@cdesc$cell_id == CT &
                                          lvl4_data_all@cdesc$pert_iname %in% lig16,])
    if (temp_tx < 1000) {
      temp_param$n_neighbors <- 5
    } else {
      temp_param$n_neighbors <- 30
    }
    lvl4_ctUMAP[[CT]] <- umap(t(lvl4_data_all@mat[,temp_tx]),config=temp_param)
  }
  save(lvl4_ctUMAP,file="~/Dropbox/GDB/CMapCorr_files/200721_lvl4all_ctUMAP.RData")
}


# UMAP per LIG ----
if (file.exists("~/Dropbox/GDB/CMapCorr_files/200721_lvl4all_ligUMAP.RData")) {
  load("~/Dropbox/GDB/CMapCorr_files/200721_lvl4all_ligUMAP.RData") 
} else {
  lvl4_ligUMAP <- list()
  temp_param <- umap.defaults
  temp_param$n_neighbors <- 30
  temp_param$metric <- "cosine"
  temp_param$min_dist <- 0.2
  temp_param$n_epochs <- 500
  for (LIG in lig16) {
    message(paste0(which(LIG == lig16),"/",length(lig16)))
    temp_tx <- rownames(lvl4_data_all@cdesc[lvl4_data_all@cdesc$pert_iname == LIG,])
    lvl4_ligUMAP[[LIG]] <- umap(t(lvl4_data_all@mat[,temp_tx]),config=temp_param)
  }
  save(lvl4_ligUMAP,file="~/Dropbox/GDB/CMapCorr_files/200721_lvl4all_ligUMAP.RData")
}
