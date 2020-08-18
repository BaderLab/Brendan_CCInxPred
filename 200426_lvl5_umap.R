library(cmapR)
library(umap)

setwd("~/Dropbox/GDB/CMapCorr/")
if (exists("lvl5_data")) {
} else if (file.exists("~/Dropbox/GDB_archive/CMapCorr_files/lvl5_inputs.RData")) {
  load("~/Dropbox/GDB_archive/CMapCorr_files/lvl5_inputs.RData") 
} else {
  source("lvl5_inputs.R")
}


# UMAP ----
if (file.exists("~/Dropbox/GDB/CMapCorr_files/200426_lvl5_umap.RData")) {
  load("~/Dropbox/GDB/CMapCorr_files/200426_lvl5_umap.RData") 
} else {
  temp_param <- umap.defaults
  temp_param$n_neighbors <- 30
  temp_param$metric <- "cosine"
  temp_param$min_dist <- 0.2
  temp_param$n_epochs <- 500

  lvl5_umap <- umap(t(lvl5_data@mat),config=temp_param)
  save(lvl5_umap,file="~/Dropbox/GDB/CMapCorr_files/200426_lvl5_umap.RData")
}

# UMAP by cell type ----
if (file.exists("~/Dropbox/GDB/CMapCorr_files/200426_lvl5_ctUMAP.RData")) {
  load("~/Dropbox/GDB/CMapCorr_files/200426_lvl5_ctUMAP.RData") 
} else {
  temp_param <- umap.defaults
  temp_param$n_neighbors <- 30
  temp_param$metric <- "cosine"
  temp_param$min_dist <- 0.2
  temp_param$n_epochs <- 500
  
  ctUMAP <- ctUMAPraw <- list()
  for (CT in ct14) {
    ctUMAPraw[[CT]] <- umap(t(lvl5_data@mat[,lvl5_data@cdesc$cell_id == CT]),config=temp_param)
  }
  save(lvl5_ctUMAP,file="~/Dropbox/GDB/CMapCorr_files/200426_lvl5_ctUMAP.RData")
}


# UMAP by cell type with common ligands only ----
if (file.exists("~/Dropbox/GDB/CMapCorr_files/200426_lvl5_ctligUMAP.RData")) {
  load("~/Dropbox/GDB/CMapCorr_files/200426_lvl5_ctligUMAP.RData") 
} else {
  temp_param <- umap.defaults
  temp_param$n_neighbors <- 5  #because many individual cell types only have ~5 samples per cell type.
  temp_param$metric <- "cosine"
  temp_param$min_dist <- 0.2
  temp_param$n_epochs <- 500
  
  temp_ct <- ct14[sapply(ct14,function(CT) 
    ncol(lvl5_data@mat[,lvl5_data@cdesc$cell_id == CT & 
                         lvl5_data@cdesc$pert_iname %in% lig16])) > 100]
  ctligUMAP <- list()
  for (CT in temp_ct) {
    ctligUMAP[[CT]] <- umap(t(lvl5_data@mat[,lvl5_data@cdesc$cell_id == CT & 
                                              lvl5_data@cdesc$pert_iname %in% lig16]),
                            config=temp_param)
  }
  save(ctligUMAP,file="~/Dropbox/GDB/CMapCorr_files/200426_lvl5_ctligUMAP.RData")
}


# UMAP by ligand ----
if (file.exists("~/Dropbox/GDB/CMapCorr_files/200426_lvl5_ligUMAP.RData")) {
  load("~/Dropbox/GDB/CMapCorr_files/200426_lvl5_ligUMAP.RData") 
} else {
  temp_param <- umap.defaults
  temp_param$n_neighbors <- 5  #because many individual cell types only have ~5 samples.
  temp_param$metric <- "cosine"
  temp_param$min_dist <- 0.2
  temp_param$n_epochs <- 500
  
  ligUMAP <- list()
  for (L in lig16) {
    # temp_pca <- prcomp(t(lvl5_data@mat[,lvl5_data@cdesc$pert_iname == L]))
    # ligUMAP[[L]] <- umap(temp_pca$x[,1:50],config=temp_param)
    ligUMAP[[L]] <- umap(t(lvl5_data@mat[,lvl5_data@cdesc$pert_iname == L]),
                         config=temp_param)
  }
  save(ligUMAP,file="~/Dropbox/GDB/CMapCorr_files/200426_lvl5_ligUMAP.RData") 
}


# UMAP of LigPred_mixall ----
if (file.exists("~/Dropbox/GDB/CMapCorr_files/200427_lvl5_UMAPmixall.RData")) {
  load("~/Dropbox/GDB/CMapCorr_files/200427_lvl5_UMAPmixall.RData") 
} else {
  temp_param <- umap.defaults
  temp_param$n_neighbors <- 5  #because many individual cell types only have ~5 samples.
  temp_param$metric <- "cosine"
  temp_param$min_dist <- 0.2
  temp_param$n_epochs <- 500
  
  UMAPmixall <- umap(t(lvl5_data@mat[,lvl5_data@cdesc$pert_iname %in% lig16]),
                     config=temp_param)
  
  plot(UMAPmixall$layout,pch=".")
  
  save(UMAPmixall,file="~/Dropbox/GDB/CMapCorr_files/200427_lvl5_UMAPmixall.RData")
}

