library(cmapR)

rm(list=ls())
setwd("~/Dropbox/GDB/CMapCorr/")
temp_cmap_path <- "~/Data_LINCS/phase1"

# sample metadata ----
temp_coldata <- read.gctx.meta(
  file.path(temp_cmap_path,"annotated_GSE92742_Broad_LINCS_Level5_COMPZ_n473647x12328.gctx"),
  dimension="column")

gene_info <- read.table(file.path(temp_cmap_path,"GSE92742_Broad_LINCS_gene_info.txt"),
                        header=T,sep="\t",row.names=1,colClasses="character",quote="\"")
# Only going to use the landmark genes, because those were actually measured.
# This reduces the number of features, such that features don't massively outnumber samples.

cell_info <- read.table("~/Data_LINCS/phase1/GSE92742_Broad_LINCS_cell_info.txt",
                        header=T,sep="\t",row.names=1,colClasses="character",quote="\"")

lvl5_data <- parse.gctx(
  file.path(temp_cmap_path,"annotated_GSE92742_Broad_LINCS_Level5_COMPZ_n473647x12328.gctx"),
  cid=temp_coldata$id[temp_coldata$pert_type == "trt_lig"],
  rid=rownames(gene_info)[gene_info$pr_is_lm == "1"])
if (all(lvl5_data@cdesc$id == colnames(lvl5_data@mat))) {
  rownames(lvl5_data@cdesc) <- lvl5_data@cdesc$id
  lvl5_data@cdesc <- lvl5_data@cdesc[,colnames(lvl5_data@cdesc) != "id"]
}
if (all(lvl5_data@rdesc$id == rownames(lvl5_data@mat))) {
  rownames(lvl5_data@rdesc) <- lvl5_data@rdesc$id
  lvl5_data@rdesc <- lvl5_data@rdesc[,colnames(lvl5_data@rdesc) != "id"]
}
# fixing dose unit discrepancy
lvl5_data@cdesc$pert_dose_unit[lvl5_data@cdesc$pert_dose_unit == "ng/<fd><fd>L"] <- "ng/uL"
lvl5_data@cdesc$pert_dose[lvl5_data@cdesc$pert_dose == "0.1" & 
                            lvl5_data@cdesc$pert_dose_unit == "ng/uL"] <- "100.0"
lvl5_data@cdesc$pert_dose[lvl5_data@cdesc$pert_dose == "1.0" & 
                            lvl5_data@cdesc$pert_dose_unit == "ng/uL"] <- "1000.0"
lvl5_data@cdesc$pert_dose[lvl5_data@cdesc$pert_dose == "3.0" & 
                            lvl5_data@cdesc$pert_dose_unit == "ng/uL"] <- "3000.0"
lvl5_data@cdesc$pert_dose[lvl5_data@cdesc$pert_dose == "10.0" & 
                            lvl5_data@cdesc$pert_dose_unit == "ng/uL"] <- "10000.0"
lvl5_data@cdesc$pert_dose[lvl5_data@cdesc$pert_dose == "100.0" & 
                            lvl5_data@cdesc$pert_dose_unit == "ng/uL"] <- "100000.0"
lvl5_data@cdesc$pert_dose[lvl5_data@cdesc$pert_dose == "300.0" & 
                            lvl5_data@cdesc$pert_dose_unit == "ng/uL"] <- "300000.0"
lvl5_data@cdesc$pert_dose_unit[lvl5_data@cdesc$pert_dose_unit == "ng/uL"] <- "ng/mL"

# lvl5_data_ctl <- parse.gctx(
#   file.path(temp_cmap_path,"annotated_GSE92742_Broad_LINCS_Level5_COMPZ_n473647x12328.gctx"),
#   cid=temp_coldata$id[temp_coldata$pert_type == "ctl_vehicle" & 
#                         temp_coldata$cell_id %in% unique(lvl5_data@cdesc$cell_id)],
#   rid=rownames(gene_info)[gene_info$pr_is_lm == "1"])

temp_ct <- sapply(unique(lvl5_data@cdesc$cell_id),function(X) 
  unique(lvl5_data@cdesc$pert_iname[lvl5_data@cdesc$cell_id == X]),
  simplify=F)
lig15 <- Reduce(intersect,temp_ct)

temp_ctype <- apply(cell_info[unique(lvl5_data@cdesc$cell_id),c("primary_site","sample_type")],1,
                    function(X) paste(X,collapse=" "))
temp_md <- data.frame(sets=names(temp_ctype),source.tissue=temp_ctype)
ct14 <- names(sort(table(lvl5_data@cdesc$cell_id[lvl5_data@cdesc$pert_iname %in% lig15]),decreasing=T))
names(ct14) <- paste(temp_md[ct14,"source.tissue"],ct14)
names(ct14) <- gsub(" ","_",names(ct14))

rm(list=grep("^temp",ls(),value=T))

save(list=ls(),file="~/Dropbox/GDB/CMapCorr_files/lvl5_inputs.RData")
