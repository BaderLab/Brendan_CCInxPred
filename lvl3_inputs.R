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

cell_info <- read.table("~/Data_LINCS/phase1/GSE92742_Broad_LINCS_cell_info.txt",
                        header=T,sep="\t",row.names=1,colClasses="character",quote="\"")

lvl3_data <- parse.gctx(
  file.path(temp_cmap_path,"GSE92742_Broad_LINCS_Level3_INF_mlr12k_n1319138x12328.gctx"),
  cid=rownames(temp_coldata)[temp_coldata$pert_type == "trt_lig"],
  rid=rownames(temp_geneinfo)[temp_geneinfo$pr_is_lm == "1"])
temp_id <- lvl3_data@cdesc$id
lvl3_data@cdesc <- temp_coldata[temp_id,]
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
  file.path(temp_cmap_path,"GSE92742_Broad_LINCS_Level3_INF_mlr12k_n1319138x12328.gctx"),
  cid=rownames(temp_coldata)[temp_coldata$pert_type == "ctl_vehicle" &
                               temp_coldata$cell_id %in% 
                               unique(temp_coldata[rownames(lvl3_data@cdesc),"cell_id"])],
  rid=rownames(temp_geneinfo)[temp_geneinfo$pr_is_lm == "1"])
temp_id <- lvl3_data_ctl@cdesc$id
lvl3_data_ctl@cdesc <- temp_coldata[temp_id,]
# fixing floating-point bullshit
lvl3_data_ctl@cdesc$pert_dose[lvl3_data_ctl@cdesc$pert_dose == "0.10000000149"] <- "0.1"
lvl3_data_ctl@cdesc$pert_dose <- sub("\\.?0+$","",lvl3_data_ctl@cdesc$pert_dose)

temp_ct <- sapply(unique(lvl3_data@cdesc$cell_id),function(X) 
  unique(lvl3_data@cdesc$pert_iname[lvl3_data@cdesc$cell_id == X]),
  simplify=F)
lig15 <- Reduce(intersect,temp_ct)

temp_ctype <- apply(cell_info[unique(lvl3_data@cdesc$cell_id),c("primary_site","sample_type")],1,
                    function(X) paste(X,collapse=" "))
temp_md <- data.frame(sets=names(temp_ctype),source.tissue=temp_ctype)
ct14 <- names(sort(table(lvl3_data@cdesc$cell_id[lvl3_data@cdesc$pert_iname %in% lig15]),decreasing=T))
names(ct14) <- paste(temp_md[ct14,"source.tissue"],ct14)
names(ct14) <- gsub(" ","_",names(ct14))

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

save(list=ls(),file="~/Dropbox/GDB/CMapCorr_files/lvl3_inputs.RData")
