library(cmapR)

rm(list=ls())
setwd("~/Dropbox/GDB/CMapCorr/")
temp <- load("~/Dropbox/GDB/CMapCorr_files/lvl5_inputs.RData")
rm(list=c("temp",temp[!temp %in% c("lig15","ct14")]))

temp_cmap_path <- "~/Data_LINCS/phase1"
temp_coldata <- read.table(file.path(temp_cmap_path,"GSE92742_Broad_LINCS_inst_info.txt"),
                           header=T,sep="\t",row.names=1,colClasses="character",quote="\"")
# always check # of rows matches .txt with `awk -F "\t" 'END {print NR}' filename`
# and when it doesn't, make sure R isn't quoting out things based on apostrophes.
temp_geneinfo <- read.table(file.path(temp_cmap_path,"GSE92742_Broad_LINCS_gene_info.txt"),
                            header=T,sep="\t",row.names=1,colClasses="character",quote="\"")

lvl3_ctl <- parse.gctx(
  file.path(temp_cmap_path,"GSE92742_Broad_LINCS_Level3_INF_mlr12k_n1319138x12328.gctx"),
  cid=rownames(temp_coldata)[temp_coldata$pert_type == "ctl_vehicle" &
                               temp_coldata$cell_id %in% ct14])
temp_id <- lvl3_ctl@cdesc$id
lvl3_ctl@cdesc <- temp_coldata[temp_id,]
# fixing floating-point bullshit
lvl3_ctl@cdesc$pert_dose[lvl3_ctl@cdesc$pert_dose == "0.10000000149"] <- "0.1"
lvl3_ctl@cdesc$pert_dose <- sub("\\.?0+$","",lvl3_ctl@cdesc$pert_dose)
lvl3_ctl@rdesc <- cbind(lvl3_ctl@rdesc,temp_geneinfo[lvl3_ctl@rdesc$id,])

rm(list=grep("^temp",ls(),value=T))

save(list=ls(),file="~/Dropbox/GDB/CMapCorr_files/lvl3_allgenesctrl.RData")
