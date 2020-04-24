library(cmapR)
library(ranger)
library(umap)

if (!exists("lvl3_data")) {
  setwd("~/Dropbox/GDB/CMapCorr/")
  cmap_file_path <- "~/Data_LINCS/phase1/"
  lvl3_coldata <- read.table(file.path(cmap_file_path,"GSE92742_Broad_LINCS_inst_info.txt"),
                             header=T,sep="\t",row.names=1,colClasses="character",quote="\"")
  # always check # of rows matches .txt with `awk -F "\t" 'END {print NR}' filename`
  # and when it doesn't, make sure R isn't quoting out things based on apostrophes.
  temp_geneinfo <- read.table(file.path(cmap_file_path,"GSE92742_Broad_LINCS_gene_info.txt"),
                              header=T,sep="\t",row.names=1,colClasses="character",quote="\"")
  # Only going to use the landmark genes, because those were actually measured.
  # This reduces the number of features, such that features don't massively outnumber samples.
  
  lvl3_data <- parse.gctx(
    file.path(cmap_file_path,"GSE92742_Broad_LINCS_Level3_INF_mlr12k_n1319138x12328.gctx"),
    cid=rownames(lvl3_coldata)[lvl3_coldata$pert_type == "trt_lig"],
    rid=rownames(temp_geneinfo)[temp_geneinfo$pr_is_lm == "1"])
  temp_id <- lvl3_data@cdesc$id
  lvl3_data@cdesc <- lvl3_coldata[temp_id,]
  # fixing floating-point bullshit
  lvl3_data@cdesc$pert_dose[lvl3_data@cdesc$pert_dose == "0.00999999977648"] <- "0.01"
  lvl3_data@cdesc$pert_dose[lvl3_data@cdesc$pert_dose == "0.029999999329400003"] <- "0.03"
  lvl3_data@cdesc$pert_dose[lvl3_data@cdesc$pert_dose == "0.10000000149"] <- "0.1"
  lvl3_data@cdesc$pert_dose[lvl3_data@cdesc$pert_dose == "0.15000000596"] <- "0.15"
  lvl3_data@cdesc$pert_dose[lvl3_data@cdesc$pert_dose == "0.20000000298"] <- "0.2"
  lvl3_data@cdesc$pert_dose[lvl3_data@cdesc$pert_dose == "0.300000011921"] <- "0.3"
  lvl3_data@cdesc$pert_dose[lvl3_data@cdesc$pert_dose == "1.64999997616"] <- "1.65"
  lvl3_data@cdesc$pert_dose <- sub("\\.?0+$","",lvl3_data@cdesc$pert_dose)
  # fixing dose unit discrepancy
  lvl3_data@cdesc$pert_dose[lvl3_data@cdesc$pert_dose == "0.1" & 
                              lvl3_data@cdesc$pert_dose_unit == "ng/ul"] <- "100"
  lvl3_data@cdesc$pert_dose[lvl3_data@cdesc$pert_dose == "1" & 
                              lvl3_data@cdesc$pert_dose_unit == "ng/ul"] <- "1000"
  lvl3_data@cdesc$pert_dose[lvl3_data@cdesc$pert_dose == "3" & 
                              lvl3_data@cdesc$pert_dose_unit == "ng/ul"] <- "3000"
  lvl3_data@cdesc$pert_dose[lvl3_data@cdesc$pert_dose == "10" & 
                              lvl3_data@cdesc$pert_dose_unit == "ng/ul"] <- "10000"
  lvl3_data@cdesc$pert_dose[lvl3_data@cdesc$pert_dose == "100" & 
                              lvl3_data@cdesc$pert_dose_unit == "ng/ul"] <- "100000"
  lvl3_data@cdesc$pert_dose[lvl3_data@cdesc$pert_dose == "300" & 
                              lvl3_data@cdesc$pert_dose_unit == "ng/ul"] <- "300000"
  lvl3_data@cdesc$pert_dose_unit[lvl3_data@cdesc$pert_dose_unit == "ng/ul"] <- "ng/ml"
  lvl3_data@cdesc$pert_dose_unit[lvl3_data@cdesc$pert_dose_unit == "ng/ml"] <- "ng/mL"
  
  lvl3_data_ctl <- parse.gctx(
    file.path(cmap_file_path,"GSE92742_Broad_LINCS_Level3_INF_mlr12k_n1319138x12328.gctx"),
    cid=rownames(lvl3_coldata)[lvl3_coldata$pert_type == "ctl_vehicle" &
                                 lvl3_coldata$cell_id %in% 
                                 unique(lvl3_coldata[rownames(lvl3_data@cdesc),"cell_id"])],
    rid=rownames(temp_geneinfo)[temp_geneinfo$pr_is_lm == "1"])
  temp_id <- lvl3_data_ctl@cdesc$id
  lvl3_data_ctl@cdesc <- lvl3_coldata[temp_id,]
  # fixing floating-point bullshit
  lvl3_data_ctl@cdesc$pert_dose[lvl3_data_ctl@cdesc$pert_dose == "0.10000000149"] <- "0.1"
  lvl3_data_ctl@cdesc$pert_dose <- sub("\\.?0+$","",lvl3_data_ctl@cdesc$pert_dose)
  
  temp_ct <- sapply(unique(lvl3_data@cdesc$cell_id),function(X) 
    unique(lvl3_data@cdesc$pert_iname[lvl3_data@cdesc$cell_id == X]),
    simplify=F)
  lig15 <- Reduce(intersect,temp_ct)
  
  PM <- sapply(
    unique(lvl3_data@cdesc[lvl3_data@cdesc$pert_iname %in% lig15,"cell_id"]),
    function(CT) {
      sapply(lig15,function(L) {
        temp_tx <- rownames(lvl3_data@cdesc)[lvl3_data@cdesc$cell_id == CT & 
                                               lvl3_data@cdesc$pert_iname == L]
        temp_pl <- unique(lvl3_data@cdesc[temp_tx,"rna_plate"])
        temp_ctl <- rownames(lvl3_data_ctl@cdesc)[lvl3_data_ctl@cdesc$rna_plate %in% temp_pl &
                                                    lvl3_data_ctl@cdesc$cell_id == CT]
        # message(paste(length(temp_tx),length(temp_ctl)))
        return(list(tx=temp_tx,
                    ctl=temp_ctl))
      },simplify=F)
    },simplify=F)
  
  rm(list=grep("^temp",ls(),value=T))
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

