# Landmark genes only ----
load("~/Dropbox/GDB_archive/CMapCorr_files/lvl4_inputs.RData")

hist(lvl4_data_ctl@mat,ylim=c(0,100))
hist(lvl4_data@mat,ylim=c(0,100))

boxplot(sapply(ct14,function(X) rowMeans(lvl4_data_ctl@mat[,lvl4_data_ctl@cdesc$cell_id == X])))


rep_tx <- unique(lvl4_data@cdesc[lvl4_data@cdesc$pert_iname %in% lig16 &
                                   lvl4_data@cdesc$cell_id %in% ct14,
                                 c("pert_iname","cell_id","pert_dose","pert_time")])
rownames(rep_tx) <- apply(rep_tx,1,paste,collapse="_")

X <- unlist(rep_tx[1,])
temp_repL <- lvl4_data@cdesc$pert_iname == X[1] &
  lvl4_data@cdesc$cell_id == X[2] &
  lvl4_data@cdesc$pert_dose == X[3] &
  lvl4_data@cdesc$pert_time == X[4]

temp_P <- pnorm(-abs(rowMeans(lvl4_data@mat[,temp_repL])))
temp_Q <- p.adjust(temp_P,"fdr")                
sum(temp_Q1 <= 0.01)



# All genes (including inferred) ----
load("~/Dropbox/GDB_archive/CMapCorr_files/lvl4_inputs_allgenes_CTL.RData")
