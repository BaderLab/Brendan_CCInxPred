load("~/Dropbox/GDB_archive/CMapCorr_files/lvl4_inputs_allgenes.RData")
lvl4_data <- lvl4_data_all
rm(lvl4_data_all)

temp_ligcountsperct <- sapply(unique(lvl4_data@cdesc$cell_id),function(CT)
  length(unique(lvl4_data@cdesc[lvl4_data@cdesc$cell_id == CT,"pert_iname"])))
ct9 <- names(temp_ligcountsperct)[temp_ligcountsperct > 100]
names(ct9) <- sapply(ct9,function(X) names(ct14)[ct14 == X])

temp_ctcountsperlig <- sapply(unique(lvl4_data@cdesc$pert_iname),function(LIG)
  length(unique(lvl4_data@cdesc[lvl4_data@cdesc$pert_iname == LIG,"cell_id"])))

lig295 <- names(temp_ctcountsperlig)[temp_ctcountsperlig >= 9]

rm(list=grep("^temp",ls(),value=T))


# Unique marker genes per ligand ----
load("~/Dropbox/GDB_archive/CMapCorr_files/lig295_DE_lig_FDR.RData")
temp_logica_lig <- FDR_lig <= 0.05
temp_logica_lig
temp_genes <- apply(temp_logica_lig,2,which)
temp_genes
temp_genes <- temp_genes[sapply(temp_genes,length) > 0]
temp_genes
sapply(temp_genes,names)
temp_genes <- sapply(temp_genes,names)

unique_lig <- sapply(names(temp_genes),function(X) 
  setdiff(temp_genes[[X]],unlist(temp_genes[names(temp_genes) != X])))
unique_lig <- unique_lig[sapply(unique_lig,length) > 0]

unique_lig <- sapply(unique_lig,function(X) lvl4_data@rdesc[X,"pr_gene_symbol"])

rm(list=grep("^temp",ls(),value=T))


# LigCT ----
load("~/Dropbox/GDB_archive/CMapCorr_files/lig295_DE_ligct_FDR.RData")
temp_lig <- sapply(strsplit(colnames(FDR_ligct),"_"),function(X) X[length(X)])
table(sapply(lig295,function(LIG) sum(FDR_ligct[,temp_lig == LIG] <= 0.05)) > 0)
table(sapply(lig295,function(LIG) sum(pDE_ligct[temp_lig == LIG,] <= 0.05)) > 0)


# Corr CMap vs NN ----
load("~/Dropbox/GDB_archive/CMapCorr_files/NN_ALLvCT_corr.RData")
load("~/Dropbox/GDB_archive/CMapCorr_files/lig295_corr_lig.RData")
CORS_lig <- CORS; rm(CORS)
load("~/Dropbox/GDB_archive/CMapCorr_files/lig295_corr_ligct.RData")
CORS_ligct <- CORS; rm(CORS)

median(unlist(CORS_lig))
median(unlist(corr_all))
wilcox.test(unlist(CORS_lig),unlist(corr_all))
t.test(unlist(CORS_lig),unlist(corr_all))

median(unlist(CORS_ligct))
median(unlist(corr_ct))
wilcox.test(unlist(CORS_ligct),unlist(corr_ct))
t.test(unlist(CORS_ligct),unlist(corr_ct))


# TAS ----
temp_sig <- read.table("~/Data_LINCS/phase1/GSE92742_Broad_LINCS_sig_metrics.txt",
                       header=T,sep="\t",row.names=1,colClasses="character",quote="\"")

TASdf <- data.frame(CC=as.numeric(temp_sig[lvl5_data@cid,"distil_cc_q75"]),
                    SS=as.numeric(temp_sig[lvl5_data@cid,"distil_ss"]),
                    TAS=as.numeric(temp_sig[lvl5_data@cid,"tas"]),
                    cell_id=lvl5_data@cdesc$cell_id,
                    pert_iname=lvl5_data@cdesc$pert_iname,
                    pert_dose=lvl5_data@cdesc$pert_dose,
                    pert_time=lvl5_data@cdesc$pert_time,
                    row.names=lvl5_data@cid)
TASdf$CC[TASdf$CC < -1] <- NA
save(TASdf,file="~/Dropbox/GDB_archive/CMapCorr_files/lvl5_TAS.RData")
