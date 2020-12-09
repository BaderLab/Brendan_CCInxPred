library(cmapR)

temp <- load("~/Dropbox/GDB_archive/CMapCorr_files/lvl3_inputs.RData")
rm(list=ls()[!ls() %in% c("ct14","lig16")])

temp_cmap_path <- "~/Data_LINCS/phase1"

# sample metadata ----
temp_coldata <- read.gctx.meta(
  file.path(temp_cmap_path,"annotated_GSE92742_Broad_LINCS_Level5_COMPZ_n473647x12328.gctx"),
  dim="column")

lvl5_data <- parse.gctx(
  file.path(temp_cmap_path,"annotated_GSE92742_Broad_LINCS_Level5_COMPZ_n473647x12328.gctx"),
  cid=temp_coldata$id[temp_coldata$pert_type == "trt_lig"])

if (all(lvl5_data@cdesc$id == colnames(lvl5_data@mat))) {
  rownames(lvl5_data@cdesc) <- lvl5_data@cdesc$id
  lvl5_data@cdesc <- lvl5_data@cdesc[,colnames(lvl5_data@cdesc) != "id"]
}
if (all(lvl5_data@rdesc$id == rownames(lvl5_data@mat))) {
  rownames(lvl5_data@rdesc) <- lvl5_data@rdesc$id
  lvl5_data@rdesc <- lvl5_data@rdesc[,colnames(lvl5_data@rdesc) != "id"]
}
#fixing fucked-up ligand names
lvl5_data@cdesc$pert_iname <- toupper(lvl5_data@cdesc$pert_iname)
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

rm(list=grep("^temp",ls(),value=T))

save(list=ls(),file="~/Dropbox/GDB_archive/CMapCorr_files/lvl5_inputs_allgenes.RData")
