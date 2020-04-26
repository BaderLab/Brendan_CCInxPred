library(cmapR)

setwd("~/Dropbox/GDB/CMapCorr/")
temp_cmap_path <- "~/Data_LINCS/phase1"

# sample metadata ----
temp_coldata <- read.gctx.meta(
  file.path(temp_cmap_path,"annotated_GSE92742_Broad_LINCS_Level5_COMPZ_n473647x12328.gctx"),
  dimension="column")
lvl5_data <- parse.gctx(
  file.path(temp_cmap_path,"annotated_GSE92742_Broad_LINCS_Level5_COMPZ_n473647x12328.gctx"),
  cid=temp_coldata$id[temp_coldata$pert_type == "trt_lig"])
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

lvl5_data_ctl <- parse.gctx(
  file.path(temp_cmap_path,"annotated_GSE92742_Broad_LINCS_Level5_COMPZ_n473647x12328.gctx"),
  cid=temp_coldata$id[temp_coldata$pert_type == "ctl_vehicle" & 
                        temp_coldata$cell_id %in% unique(lvl5_data@cdesc$cell_id)])
