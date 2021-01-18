library(cmapR)
library(colorspace)
library(scales)
library(pbapply)

temp <- load("~/Dropbox/GDB_archive/CMapCorr_files/lvl4_inputs.RData") 

temp_ligcountsperct <- sapply(unique(lvl4_data@cdesc$cell_id),function(CT)
  length(unique(lvl4_data@cdesc[lvl4_data@cdesc$cell_id == CT,"pert_iname"])))
ct9 <- names(temp_ligcountsperct)[temp_ligcountsperct > 100]
names(ct9) <- sapply(ct9,function(X) names(ct14)[ct14 == X])

lvl4_data@cdesc <- lvl4_data@cdesc[lvl4_data@cdesc$cell_id %in% ct9,]
lvl4_data@mat <- lvl4_data@mat[,rownames(lvl4_data@cdesc)]

temp_ctcountsperlig <- sapply(unique(lvl4_data@cdesc$pert_iname),function(LIG)
  length(unique(lvl4_data@cdesc[lvl4_data@cdesc$pert_iname == LIG,"cell_id"])))

lig295 <- names(temp_ctcountsperlig)[temp_ctcountsperlig == 9]

lvl4_data@cdesc <- lvl4_data@cdesc[lvl4_data@cdesc$pert_iname %in% lig295,]
lvl4_data@mat <- lvl4_data@mat[,rownames(lvl4_data@cdesc)]

rm(list=temp[temp != "lvl4_data"])
rm(list=grep("^temp",ls(),value=T))



meanZlig <- sapply(lig295,function(LIG)
  rowMeans(lvl4_data@mat[,lvl4_data@cdesc$pert_iname == LIG]))

qqnorm(meanZlig[,"ACHE"],pch=".",cex=2)
qqline(meanZlig[,"ACHE"])
qqnorm(as.vector(meanZlig),pch=".",cex=2)
qqline(as.vector(meanZlig))

