library(cmapR)

rm(list=ls())
setwd("~/Dropbox/GDB/CMapCorr/")
temp_cmap_path <- "~/Data_LINCS/phase1"
temp_coldata <- read.table(file.path(temp_cmap_path,"GSE92742_Broad_LINCS_inst_info.txt"),
                           header=T,sep="\t",row.names=1,colClasses="character",quote="\"")
# always check # of rows matches .txt with `awk -F "\t" 'END {print NR}' filename`
# and when it doesn't, make sure R isn't quoting out things based on apostrophes.
temp_geneinfo <- read.table(file.path(temp_cmap_path,"GSE92742_Broad_LINCS_gene_info.txt"),
                            header=T,sep="\t",row.names=1,colClasses="character",quote="\"")
# Only going to use the landmark genes, because those were actually measured.
# This reduces the number of features, such that features don't massively outnumber samples.

temp_cell_info <- read.table("~/Data_LINCS/phase1/GSE92742_Broad_LINCS_cell_info.txt",
                        header=T,sep="\t",row.names=1,colClasses="character",quote="\"")

lvl3_data_all <- parse.gctx(
  file.path(temp_cmap_path,"GSE92742_Broad_LINCS_Level3_INF_mlr12k_n1319138x12328.gctx"),
  cid=rownames(temp_coldata)[temp_coldata$pert_type == "trt_lig"])
temp_id <- lvl3_data_all@cdesc$id
lvl3_data_all@cdesc <- temp_coldata[temp_id,]
#fixing fucked-up ligand names (except NRG/ALPHA/BETA)
lvl3_data_all@cdesc$pert_iname <- toupper(lvl3_data_all@cdesc$pert_iname)
# fixing floating-point bullshit
lvl3_data_all@cdesc$pert_dose[lvl3_data_all@cdesc$pert_dose == "0.00999999977648"] <- "0.01"
lvl3_data_all@cdesc$pert_dose[lvl3_data_all@cdesc$pert_dose == "0.029999999329400003"] <- "0.03"
lvl3_data_all@cdesc$pert_dose[lvl3_data_all@cdesc$pert_dose == "0.10000000149"] <- "0.1"
lvl3_data_all@cdesc$pert_dose[lvl3_data_all@cdesc$pert_dose == "0.15000000596"] <- "0.15"
lvl3_data_all@cdesc$pert_dose[lvl3_data_all@cdesc$pert_dose == "0.20000000298"] <- "0.2"
lvl3_data_all@cdesc$pert_dose[lvl3_data_all@cdesc$pert_dose == "0.300000011921"] <- "0.3"
lvl3_data_all@cdesc$pert_dose[lvl3_data_all@cdesc$pert_dose == "1.64999997616"] <- "1.65"
lvl3_data_all@cdesc$pert_dose <- sub("\\.?0+$","",lvl3_data_all@cdesc$pert_dose)
# fixing dose unit discrepancy
lvl3_data_all@cdesc$pert_dose[lvl3_data_all@cdesc$pert_dose == "0.1" & 
                            lvl3_data_all@cdesc$pert_dose_unit == "ng/ul"] <- "100"
lvl3_data_all@cdesc$pert_dose[lvl3_data_all@cdesc$pert_dose == "1" & 
                            lvl3_data_all@cdesc$pert_dose_unit == "ng/ul"] <- "1000"
lvl3_data_all@cdesc$pert_dose[lvl3_data_all@cdesc$pert_dose == "3" & 
                            lvl3_data_all@cdesc$pert_dose_unit == "ng/ul"] <- "3000"
lvl3_data_all@cdesc$pert_dose[lvl3_data_all@cdesc$pert_dose == "10" & 
                            lvl3_data_all@cdesc$pert_dose_unit == "ng/ul"] <- "10000"
lvl3_data_all@cdesc$pert_dose[lvl3_data_all@cdesc$pert_dose == "100" & 
                            lvl3_data_all@cdesc$pert_dose_unit == "ng/ul"] <- "100000"
lvl3_data_all@cdesc$pert_dose[lvl3_data_all@cdesc$pert_dose == "300" & 
                            lvl3_data_all@cdesc$pert_dose_unit == "ng/ul"] <- "300000"
lvl3_data_all@cdesc$pert_dose_unit[lvl3_data_all@cdesc$pert_dose_unit == "ng/ul"] <- "ng/ml"
lvl3_data_all@cdesc$pert_dose_unit[lvl3_data_all@cdesc$pert_dose_unit == "ng/ml"] <- "ng/mL"

lvl3_data_all_ctl <- parse.gctx(
  file.path(temp_cmap_path,"GSE92742_Broad_LINCS_Level3_INF_mlr12k_n1319138x12328.gctx"),
  cid=rownames(temp_coldata)[temp_coldata$pert_type == "ctl_vehicle" &
                               temp_coldata$rna_plate %in% 
                               unique(temp_coldata[rownames(lvl3_data_all@cdesc),"rna_plate"])])
temp_id <- lvl3_data_all_ctl@cdesc$id
lvl3_data_all_ctl@cdesc <- temp_coldata[temp_id,]
#fixing fucked-up ligand names
lvl3_data_all_ctl@cdesc$pert_iname <- toupper(lvl3_data_all_ctl@cdesc$pert_iname)
# fixing floating-point bullshit
lvl3_data_all_ctl@cdesc$pert_dose[lvl3_data_all_ctl@cdesc$pert_dose == "0.10000000149"] <- "0.1"
lvl3_data_all_ctl@cdesc$pert_dose <- sub("\\.?0+$","",lvl3_data_all_ctl@cdesc$pert_dose)


rm(list=grep("^temp",ls(),value=T))

save(list=ls(),file="~/Dropbox/GDB/CMapCorr_files/lvl3_inputs_allgenes.RData")
