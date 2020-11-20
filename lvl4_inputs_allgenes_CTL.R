library(cmapR)

rm(list=ls())

temp <- load("~/Dropbox/GDB_archive/CMapCorr_files/lvl4_inputs_allgenes.RData")

temp_cmap_path <- "~/Data_LINCS/phase1"
temp_coldata <- read.table(file.path(temp_cmap_path,"GSE92742_Broad_LINCS_inst_info.txt"),
                           header=T,sep="\t",row.names=1,colClasses="character",quote="\"")
# always check # of rows matches .txt with `awk -F "\t" 'END {print NR}' filename`
# and when it doesn't, make sure R isn't quoting out things based on apostrophes.
temp_geneinfo <- read.table(file.path(temp_cmap_path,"GSE92742_Broad_LINCS_gene_info.txt"),
                            header=T,sep="\t",row.names=1,colClasses="character",quote="\"")
# need to label landmark genes

lvl4_data_all_ctl <- parse.gctx(
  file.path(temp_cmap_path,"GSE92742_Broad_LINCS_Level4_ZSPCINF_mlr12k_n1319138x12328.gctx"),
  cid=rownames(temp_coldata)[temp_coldata$pert_type == "ctl_vehicle" &
                               temp_coldata$cell_id %in% 
                               unique(temp_coldata[rownames(lvl4_data_all@cdesc),"cell_id"])])
lvl4_data_all_ctl@rdesc <- temp_geneinfo[lvl4_data_all_ctl@rdesc$id,]

temp_id <- lvl4_data_all_ctl@cdesc$id
lvl4_data_all_ctl@cdesc <- temp_coldata[temp_id,]
#fixing fucked-up ligand names (except NRG/ALPHA/BETA)
lvl4_data_all_ctl@cdesc$pert_iname <- toupper(lvl4_data_all_ctl@cdesc$pert_iname)
# fixing floating-point bullshit
lvl4_data_all_ctl@cdesc$pert_dose[lvl4_data_all_ctl@cdesc$pert_dose == "0.10000000149"] <- "0.1"
lvl4_data_all_ctl@cdesc$pert_dose <- sub("\\.?0+$","",lvl4_data_all_ctl@cdesc$pert_dose)


save(lvl4_data_all_ctl,file="~/Dropbox/GDB_archive/CMapCorr_files/lvl4_inputs_allgenes_CTL.RData")
