library(cmapR)

temp <- load("~/Dropbox/GDB_archive/CMapCorr_files/lvl3_inputs.RData")
rm(list=ls()[!ls() %in% c("ct14","lig16")])

temp_cmap_path <- "~/Data_LINCS/phase1"
temp_coldata <- read.table(file.path(temp_cmap_path,"GSE92742_Broad_LINCS_inst_info.txt"),
                           header=T,sep="\t",row.names=1,colClasses="character",quote="\"")
# always check # of rows matches .txt with `awk -F "\t" 'END {print NR}' filename`
# and when it doesn't, make sure R isn't quoting out things based on apostrophes.
temp_geneinfo <- read.table(file.path(temp_cmap_path,"GSE92742_Broad_LINCS_gene_info.txt"),
                            header=T,sep="\t",row.names=1,colClasses="character",quote="\"")
# Only going to use the landmark genes, because those were actually measured.
# This reduces the number of features, such that features don't massively outnumber samples.

lvl4_data <- parse.gctx(
  file.path(temp_cmap_path,"GSE92742_Broad_LINCS_Level4_ZSPCINF_mlr12k_n1319138x12328.gctx"),
  cid=rownames(temp_coldata)[temp_coldata$pert_type == "trt_lig"],
  rid=rownames(temp_geneinfo)[temp_geneinfo$pr_is_lm == "1"])
lvl4_data@rdesc <- temp_geneinfo[lvl4_data@rdesc$id,]

temp_id <- lvl4_data@cdesc$id
lvl4_data@cdesc <- temp_coldata[temp_id,]
#fixing fucked-up ligand names (except NRG/ALPHA/BETA)
lvl4_data@cdesc$pert_iname <- toupper(lvl4_data@cdesc$pert_iname)
# fixing floating-point bullshit
lvl4_data@cdesc$pert_dose[lvl4_data@cdesc$pert_dose == "0.00999999977648"] <- "0.01"
lvl4_data@cdesc$pert_dose[lvl4_data@cdesc$pert_dose == "0.029999999329400003"] <- "0.03"
lvl4_data@cdesc$pert_dose[lvl4_data@cdesc$pert_dose == "0.10000000149"] <- "0.1"
lvl4_data@cdesc$pert_dose[lvl4_data@cdesc$pert_dose == "0.15000000596"] <- "0.15"
lvl4_data@cdesc$pert_dose[lvl4_data@cdesc$pert_dose == "0.20000000298"] <- "0.2"
lvl4_data@cdesc$pert_dose[lvl4_data@cdesc$pert_dose == "0.300000011921"] <- "0.3"
lvl4_data@cdesc$pert_dose[lvl4_data@cdesc$pert_dose == "1.64999997616"] <- "1.65"
lvl4_data@cdesc$pert_dose <- sub("\\.?0+$","",lvl4_data@cdesc$pert_dose)
# fixing dose unit discrepancy
lvl4_data@cdesc$pert_dose[lvl4_data@cdesc$pert_dose == "0.1" & 
                            lvl4_data@cdesc$pert_dose_unit == "ng/ul"] <- "100"
lvl4_data@cdesc$pert_dose[lvl4_data@cdesc$pert_dose == "1" & 
                            lvl4_data@cdesc$pert_dose_unit == "ng/ul"] <- "1000"
lvl4_data@cdesc$pert_dose[lvl4_data@cdesc$pert_dose == "3" & 
                            lvl4_data@cdesc$pert_dose_unit == "ng/ul"] <- "3000"
lvl4_data@cdesc$pert_dose[lvl4_data@cdesc$pert_dose == "10" & 
                            lvl4_data@cdesc$pert_dose_unit == "ng/ul"] <- "10000"
lvl4_data@cdesc$pert_dose[lvl4_data@cdesc$pert_dose == "100" & 
                            lvl4_data@cdesc$pert_dose_unit == "ng/ul"] <- "100000"
lvl4_data@cdesc$pert_dose[lvl4_data@cdesc$pert_dose == "300" & 
                            lvl4_data@cdesc$pert_dose_unit == "ng/ul"] <- "300000"
lvl4_data@cdesc$pert_dose_unit[lvl4_data@cdesc$pert_dose_unit == "ng/ul"] <- "ng/ml"
lvl4_data@cdesc$pert_dose_unit[lvl4_data@cdesc$pert_dose_unit == "ng/ml"] <- "ng/mL"

# Subsetting to lig16 data only
lvl4_data@mat <- lvl4_data@mat[,lvl4_data@cdesc$pert_iname %in% lig16]
lvl4_data@cdesc <- lvl4_data@cdesc[lvl4_data@cdesc$pert_iname %in% lig16,]

rm(list=grep("^temp",ls(),value=T))

save(list=ls(),file="~/Dropbox/GDB_archive/CMapCorr_files/lvl4_lig16_inputs.RData")
