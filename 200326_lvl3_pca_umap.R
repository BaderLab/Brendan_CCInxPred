library(cmapR)
library(umap)

if (exists("lvl3_data")) {
} else if (file.exists("../CMapCorr_files/lvl3_inputs.RData")) {
  load("../CMapCorr_files/lvl3_inputs.RData") 
} else {
  source("lvl3_inputs.R")
}


# pca ----
if (file.exists("~/Dropbox/GDB/CMapCorr_files/200326_lvl3pca.RData")) {
  load("~/Dropbox/GDB/CMapCorr_files/200326_lvl3pca.RData") 
} else {
  lvl3pca <- prcomp(t(cbind(lvl3_data@mat,lvl3_data_ctl@mat)))
  save(lvl3pca,file="~/Dropbox/GDB/CMapCorr_files/200326_lvl3pca.RData")
}

# umap ----
if (file.exists("~/Dropbox/GDB/CMapCorr_files/200326_lvl3umap.RData")) {
  load("~/Dropbox/GDB/CMapCorr_files/200326_lvl3umap.RData") 
} else {
  temp_param <- umap.defaults
  temp_param$n_neighbors <- 30
  temp_param$metric <- "cosine"
  temp_param$min_dist <- 0.3
  lvl3umap <- umap(lvl3pca$x[,1:50],config=temp_param)
  save(lvl3umap,file="~/Dropbox/GDB/CMapCorr_files/200326_lvl3umap.RData")
}

if (file.exists("~/Dropbox/GDB/CMapCorr_files/200326_lvl3umap_raw.RData")) {
  load("~/Dropbox/GDB/CMapCorr_files/200326_lvl3umap_raw.RData") 
} else {
  temp_param <- umap.defaults
  temp_param$n_neighbors <- 30
  temp_param$metric <- "cosine"
  temp_param$min_dist <- 0.3
  lvl3umap_raw <- umap(t(cbind(lvl3_data@mat,lvl3_data_ctl@mat)),config=temp_param)
  save(lvl3umap_raw,file="~/Dropbox/GDB/CMapCorr_files/200326_lvl3umap_raw.RData")
}

# subset DR ----
if (file.exists("~/Dropbox/GDB/CMapCorr_files/200326_ctUMAP.RData")) {
  load("~/Dropbox/GDB/CMapCorr_files/200326_ctUMAP.RData") 
} else {
  ctPCA <- ctUMAP <- list()
  temp_param <- umap.defaults
  temp_param$n_neighbors <- 30
  temp_param$metric <- "cosine"
  temp_param$min_dist <- 0.3
  for (CT in c("BT20","VCAP","MCF7")) {
    ctPCA[[CT]] <- prcomp(t(cbind(lvl3_data@mat[,lvl3_data@cdesc$cell_id == CT],
                                  lvl3_data_ctl@mat[,lvl3_data_ctl@cdesc$cell_id == CT])))
    ctUMAP[[CT]] <- umap(ctPCA[[CT]]$x[,1:50],config=temp_param)
  }
  save(ctPCA,ctUMAP,file="~/Dropbox/GDB/CMapCorr_files/200326_ctUMAP.RData")
}


# leaveout1_PM_followup ----

# PC3, VCAP, MDAMB231, A549 samples used in leaveout1PM
# IL17A, GDNF, TNF, IGF1 samples used in leaveout1PM

if (file.exists("~/Dropbox/GDB/CMapCorr_files/200328_ctpmUMAP.RData")) {
  load("~/Dropbox/GDB/CMapCorr_files/200328_ctpmUMAP.RData") 
} else {
  temp_param <- umap.defaults
  temp_param$n_neighbors <- 30
  temp_param$metric <- "cosine"
  temp_param$min_dist <- 0.3
  
  ctUMAPpm <- list()
  for (CT in c("PC3","VCAP","MDAMB231","A549")) {
    temp_pca <- prcomp(t(cbind(
      lvl3_data@mat[,colnames(lvl3_data@mat) %in% unique(unlist(PM[[CT]],use.names=F))],
      lvl3_data_ctl@mat[,colnames(lvl3_data_ctl@mat) %in% unique(unlist(PM[[CT]],use.names=F))]
    )))
    ctUMAPpm[[CT]] <- umap(temp_pca$x[,1:50],config=temp_param)
  }
  
  ligUMAPpm <- list()
  for (L in c("IL17A","GDNF","TNF","IGF1")) {
    temp_pca <- prcomp(t(cbind(
      lvl3_data@mat[,colnames(lvl3_data@mat) %in% 
                      unique(unlist(sapply(PM,function(X) X[[L]],simplify=F),use.names=F))],
      lvl3_data_ctl@mat[,colnames(lvl3_data_ctl@mat) %in% 
                          unique(unlist(sapply(PM,function(X) X[[L]],simplify=F),use.names=F))]
    )))
    ligUMAPpm[[L]] <- umap(temp_pca$x[,1:50],config=temp_param)
  }
  
  save(ctUMAPpm,ligUMAPpm,file="~/Dropbox/GDB/CMapCorr_files/200328_ctpmUMAP.RData") 
}


# Ligand-wise UMAPs ----
if (file.exists("~/Dropbox/GDB/CMapCorr_files/200328_ligUMAP.RData")) {
  load("~/Dropbox/GDB/CMapCorr_files/200328_ligUMAP.RData") 
} else {
  temp_param <- umap.defaults
  temp_param$n_neighbors <- 5  #because many individual cell types only have ~5 samples.
  temp_param$metric <- "cosine"
  temp_param$min_dist <- 0.3
  
  ligUMAP <- list()
  for (L in lig15) {
    # temp_pca <- prcomp(t(lvl3_data@mat[,lvl3_data@cdesc$pert_iname == L]))
    # ligUMAP[[L]] <- umap(temp_pca$x[,1:50],config=temp_param)
    ligUMAP[[L]] <- umap(t(lvl3_data@mat[,lvl3_data@cdesc$pert_iname == L]),
                         config=temp_param)
  }
  save(ligUMAP,file="~/Dropbox/GDB/CMapCorr_files/200328_ligUMAP.RData") 
}

