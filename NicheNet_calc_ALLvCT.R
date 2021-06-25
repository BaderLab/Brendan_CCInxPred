library(pbapply)

# load data ----
# nn_model <- readRDS(url("https://zenodo.org/record/3260758/files/ligand_target_matrix.rds"))
nn_db <- readRDS(url("https://zenodo.org/record/3260758/files/expression_settings.rds"))

nn_db <- nn_db[sapply(sapply(nn_db,"[[","from"),length) == 1]

nn_ligands <- unique(sapply(nn_db,"[[","from"))
nn_lig_dsNames <- sapply(nn_ligands,function(LIG) names(nn_db)[sapply(nn_db,function(X) LIG %in% X$from)])
nn_ligands <- nn_ligands[sapply(nn_lig_dsNames,length) > 1]
nn_lig_dsNames <- nn_lig_dsNames[nn_ligands]


# map accessions ----
nn_DEinfo <- read.table("~/Dropbox/GDB/CMapCorr/NicheNet_DBinfo_nDE.txt",
                        header=T,sep="\t",as.is=T)
rownames(nn_DEinfo) <- nn_DEinfo$setting_id
nn_DEinfo <- nn_DEinfo[,-which(colnames(nn_DEinfo) == "setting_id")]

DSinfo <- list()
for (LIG in nn_ligands) {
  temp_DE <- sapply(nn_lig_dsNames[[LIG]],function(X) 
    nn_db[[X]]$diffexp[nn_db[[X]]$diffexp$lfc >= 1 & nn_db[[X]]$diffexp$qval <= 0.1,"gene"])
  temp_DSname <- sapply(sapply(temp_DE,length),function(N) {
    temp <- rownames(nn_DEinfo)[nn_DEinfo$ligand == LIG][nn_DEinfo[nn_DEinfo$ligand == LIG,"nr_upgenes"] %in% N]
    if (length(temp) == 0) {
      temp <- rownames(nn_DEinfo)[nn_DEinfo$ligand == LIG][nn_DEinfo[nn_DEinfo$ligand == LIG,"nr_degenes"] %in% N]
    }
    return(temp)
  })
  if (!is.character(temp_DSname)) { stop(LIG) }
  DSinfo[[LIG]] <- nn_DEinfo[temp_DSname,c("time_series","cell_type","ligand")]
  rownames(DSinfo[[LIG]]) <- names(temp_DSname)
  rm(list=grep("^temp",ls(),value=T))
}
# DSinfo$WNT1$cell_type <- paste(DSinfo$WNT1$cell_type,
#                                sapply(strsplit(rownames(DSinfo$WNT1),"_"),"[[",3),
#                                sep="_")
for (LIG in names(DSinfo)) {
  DSinfo[[LIG]]$CtAcc <- paste(DSinfo[[LIG]]$accession,DSinfo[[LIG]]$cell_type)
}


# collate replicate data ----
nn_lig_rep <- pbsapply(nn_ligands,function(LIG) {
  temp_gene <- Reduce(intersect,sapply(nn_lig_dsNames[[LIG]],
                                       function(X) nn_db[[X]]$diffexp$gene,
                                       simplify=F))
  
  temp_diffexp <- sapply(nn_lig_dsNames[[LIG]],function(X) {
    temp <- nn_db[[X]]$diffexp
    temp <- temp[temp$gene %in% temp_gene,]
    if (any(duplicated(temp$gene))) {
      temp_dup <- unique(temp$gene[duplicated(temp$gene)])
      temp_rm <- unlist(lapply(temp_dup,function(DUP) {
        temp_id <- rownames(temp[temp$gene == DUP,])
        return(temp_id[-which.min(temp[temp_id,"qval"])])
      }))
      temp <- temp[!rownames(temp) %in% temp_rm,]
    }
    rownames(temp) <- temp$gene
    return(list(lfc=temp[temp_gene,"lfc"] * 
                  ifelse(DSinfo[[LIG]][X,"time_series"],1,-1), ## BECAUSE NN FLIPPED logFC ##
                qval=temp[temp_gene,"qval"]))
  },simplify=F)
  
  temp_lfc <- sapply(temp_diffexp,"[[","lfc")
  temp_qval <- sapply(temp_diffexp,"[[","qval")
  rownames(temp_lfc) <- rownames(temp_qval) <- temp_gene
  
  return(list(lfc=temp_lfc,
              qval=temp_qval))
},simplify=F)


save(nn_ligands,DSinfo,nn_lig_rep,
     file="~/Dropbox/GDB_archive/CMapCorr_files/NN_ALLvCT_dat.RData")



# correlation calc ----
corr_all <- corr_ct <- list()
for (LIG in nn_ligands) {
  temp_ts <- rownames(DSinfo[[LIG]])[DSinfo[[LIG]]$time_series]
  temp_nts <- rownames(DSinfo[[LIG]])[!DSinfo[[LIG]]$time_series]
  
  if (length(temp_ts) > 1) {
    temp_ct <- sapply(unique(DSinfo[[LIG]][temp_ts,"cell_type"]),function(X) 
      temp_ts[DSinfo[[LIG]][temp_ts,"cell_type"] == X],
      simplify=F)
    
    temp_cor2 <- cor(nn_lig_rep[[LIG]]$lfc[,temp_ts],method="spearman")
    temp_cor_names <- sapply(colnames(temp_cor2),function(X) 
      sapply(rownames(temp_cor2), function(Y) paste(X,Y,sep=".")))
    temp_cor <- temp_cor2[lower.tri(temp_cor2)]
    names(temp_cor) <- temp_cor_names[lower.tri(temp_cor_names)]
    rm(temp_cor2,temp_cor_names)
    
    if (length(temp_cor) > 0) {
      for (Z in names(temp_ct)[sapply(temp_ct,length) > 1]) {
        temp_hits <- sapply(strsplit(names(temp_cor),".",fixed=T),function(X) all(X %in% temp_ct[[Z]]))
        corr_ct[[LIG]] <- append(corr_ct[[LIG]],temp_cor[temp_hits])
        temp_cor <- temp_cor[!temp_hits]
        rm(temp_hits)
      }
      rm(Z)
    }
    if (length(temp_cor) > 0) {
      corr_all[[LIG]] <- append(corr_all[[LIG]],temp_cor)
    }
    rm(temp_cor,temp_ct)
  }
  
  if (length(temp_nts) > 1) {
    temp_ct <- sapply(unique(DSinfo[[LIG]][temp_nts,"cell_type"]),function(X) 
      temp_nts[DSinfo[[LIG]][temp_nts,"cell_type"] == X],
      simplify=F)
    
    temp_cor2 <- cor(nn_lig_rep[[LIG]]$lfc[,temp_nts],method="spearman")
    temp_cor_names <- sapply(colnames(temp_cor2),function(X) 
      sapply(rownames(temp_cor2), function(Y) paste(X,Y,sep=".")))
    temp_cor <- temp_cor2[lower.tri(temp_cor2)]
    names(temp_cor) <- temp_cor_names[lower.tri(temp_cor_names)]
    rm(temp_cor2,temp_cor_names)
    
    if (length(temp_cor) > 0) {
      for (Z in names(temp_ct)[sapply(temp_ct,length) > 1]) {
        temp_hits <- sapply(strsplit(names(temp_cor),".",fixed=T),function(X) all(X %in% temp_ct[[Z]]))
        corr_ct[[LIG]] <- append(corr_ct[[LIG]],temp_cor[temp_hits])
        temp_cor <- temp_cor[!temp_hits]
        rm(temp_hits)
      }
      rm(Z)
    }
    if (length(temp_cor) > 0) {
      corr_all[[LIG]] <- append(corr_all[[LIG]],temp_cor)
    }
    rm(temp_cor,temp_ct)
  }
  rm(temp_nts,temp_ts)
}


save(corr_all,corr_ct,
     file="~/Dropbox/GDB_archive/CMapCorr_files/NN_ALLvCT_corr.RData")



# DE overlap background ----
nn_DE <- list()
temp_commongenes <- Reduce(intersect,sapply(nn_lig_rep,function(X) rownames(X$lfc)))
for (LIG in nn_ligands) {
    nn_DE[[LIG]] <- mapply(intersect,
                           # apply(nn_lig_rep[[LIG]]$lfc,2,function(X) names(which(X >= 1))),
                           apply(nn_lig_rep[[LIG]]$lfc,2,function(X) names(which(X >= 1 | X <= -1))),
                           apply(nn_lig_rep[[LIG]]$qval,2,function(X) names(which(X <= 0.1))))
    # nn_DE[[LIG]] <- sapply(nn_DE[[LIG]],function(X) X[X %in% temp_commongenes])
}

par(mar=c(3,3,1,1),mgp=2:0)
hist(sapply(unlist(nn_DE,recursive=F),length),
     xlab="# DE per dataset (|logFC| >= 1, FDR <= 0.1)",main=NA)

temp_DE <- unlist(nn_DE,recursive=F)
DEbkgd <- pbsapply(min(sapply(DSinfo,nrow)):max(sapply(DSinfo,nrow)),function(X)
  sapply(1:1e6,function(Y)
    length(Reduce(intersect,sample(temp_DE,X)))
  ),cl=4)
colnames(DEbkgd) <- min(sapply(DSinfo,nrow)):max(sapply(DSinfo,nrow))
rm(temp_DE)


# DE overlap calc ----
de_all <- de_ct <- list()
P_all <- P_ct <- list()
for (LIG in nn_ligands) {
  temp_ct <- sapply(unique(DSinfo[[LIG]]$cell_type),function(X) 
    rownames(DSinfo[[LIG]])[DSinfo[[LIG]]$cell_type == X],
    simplify=F)
  for (X in names(temp_ct)[sapply(temp_ct,length) > 1]) {
    Y <- paste(temp_ct[[X]],collapse=".")
    de_ct[[LIG]][[Y]] <- Reduce(intersect,nn_DE[[LIG]][temp_ct[[X]]])
    P_ct[[LIG]][[Y]] <- sum(DEbkgd[,as.character(length(temp_ct[[X]]))] >= 
                              length(de_ct[[LIG]][[Y]])) / nrow(DEbkgd)
    P_ct[[LIG]][[Y]][P_ct[[LIG]][[Y]] == 0] <- 0.1 / nrow(DEbkgd)
  }
  Y <- paste(rownames(DSinfo[[LIG]]),collapse=".")
  if (!Y %in% names(de_ct[[LIG]])) {
    de_all[[LIG]][[Y]] <- Reduce(intersect,nn_DE[[LIG]])
    P_all[[LIG]][[Y]] <- sum(DEbkgd[,as.character(length(nn_DE[[LIG]]))] >= 
                               length(de_all[[LIG]][[Y]])) / nrow(DEbkgd)
    P_all[[LIG]][[Y]][P_all[[LIG]][[Y]] == 0] <- 0.1 / nrow(DEbkgd)
  }
}


rm(list=c("X","Y","LIG",grep("^temp",ls(),value=T)))
save(nn_DE,de_all,de_ct,P_all,P_ct,
     file="~/Dropbox/GDB_archive/CMapCorr_files/NN_ALLvCT_DEovlp.RData")



## DE by mean LFC background ----
# temp_nts <- unlist(sapply(DSinfo,function(X) rownames(X)[!X$time_series]))
# temp_nts_genes <- Reduce(intersect,sapply(nn_db[temp_nts],function(X) X$diffexp$gene))
# mean_lfc_nts_bkgd <- pbsapply(
#   2:max(sapply(DSinfo,function(X) sum(!X$time_series))),
#   function(X)
#     sapply(1:1e3,function(Y)
#       rowMeans(sapply(nn_db[sample(temp_nts,X)],
#                       function(X) {
#                         temp <- X$diffexp$lfc
#                         names(temp) <- X$diffexp$gene
#                         return(temp[temp_nts_genes])
#                       }))
#     ),cl=3)
# colnames(mean_lfc_nts_bkgd) <- 2:max(sapply(DSinfo,function(X) sum(!X$time_series)))
# 
# temp_ts <- unlist(sapply(DSinfo,function(X) rownames(X)[X$time_series]))
# temp_ts_genes <- Reduce(intersect,sapply(nn_db[temp_ts],function(X) X$diffexp$gene))
# mean_lfc_ts_bkgd <- pbsapply(
#   2:max(sapply(DSinfo,function(X) sum(X$time_series))),
#   function(X)
#     sapply(1:1e3,function(Y)
#       rowMeans(sapply(nn_db[sample(temp_ts,X)],
#                       function(X) {
#                         temp <- X$diffexp$lfc
#                         names(temp) <- X$diffexp$gene
#                         return(temp[temp_ts_genes])
#                       }))
#     ),cl=3)
# colnames(mean_lfc_ts_bkgd) <- 2:max(sapply(DSinfo,function(X) sum(X$time_series)))
# 
# 
## DE by mean LFC ----
# mean_lfc <- mean_lfc_CT <- list()
# for (LIG in nn_ligands) {
#   print(paste(which(nn_ligands == LIG),"/",length(nn_ligands)))
#   temp_ts <- DSinfo[[LIG]][colnames(nn_lig_rep[[LIG]]$lfc),"time_series"]
#   
#   if (sum(!temp_ts) > 1) {
#     mean_lfc[[LIG]] <- data.frame(nts_lfc=rowMeans(nn_lig_rep[[LIG]]$lfc[,!temp_ts]))
#     mean_lfc[[LIG]]$nts_pval <- pbsapply(mean_lfc[[LIG]]$nts_lfc,function(X) 
#       2 * min(
#         c(sum(mean_lfc_nts_bkgd[,as.character(sum(!temp_ts))] >= X) / nrow(mean_lfc_nts_bkgd),
#           sum(mean_lfc_nts_bkgd[,as.character(sum(!temp_ts))] <= X) / nrow(mean_lfc_nts_bkgd))
#       ),cl=8)
#     mean_lfc[[LIG]]$nts_pval[mean_lfc[[LIG]]$nts_pval == 0] <- 2 / nrow(mean_lfc_nts_bkgd)
#     mean_lfc[[LIG]]$nts_fdr <- p.adjust(mean_lfc[[LIG]]$nts_pval,"fdr")
#     
#     temp_CT2 <- unique(DSinfo[[LIG]][!temp_ts,"cell_type"])
#     temp_CT <- sapply(temp_CT2,function(X) 
#       rownames(DSinfo[[LIG]])[DSinfo[[LIG]]$cell_type == X & !temp_ts],
#       simplify=F)
#     temp_CT <- temp_CT[sapply(temp_CT,length) > 1]
#     for (CT in names(temp_CT)) {
#       mean_lfc_CT[[LIG]][[CT]] <- data.frame(lfc=rowMeans(nn_lig_rep[[LIG]]$lfc[,temp_CT[[CT]]]))
#       mean_lfc_CT[[LIG]][[CT]]$pval <- pbsapply(mean_lfc_CT[[LIG]][[CT]]$lfc,function(X) 
#         2 * min(
#           c(sum(mean_lfc_nts_bkgd[,as.character(length(temp_CT[[CT]]))] >= X) / nrow(mean_lfc_nts_bkgd),
#             sum(mean_lfc_nts_bkgd[,as.character(length(temp_CT[[CT]]))] <= X) / nrow(mean_lfc_nts_bkgd))
#         ),cl=8)
#       mean_lfc_CT[[LIG]][[CT]]$pval[mean_lfc_CT[[LIG]][[CT]]$pval == 0] <- 2 / nrow(mean_lfc_nts_bkgd)
#       mean_lfc_CT[[LIG]][[CT]]$fdr <- p.adjust(mean_lfc_CT[[LIG]][[CT]]$pval,"fdr")
#     }
#   }
#   if (sum(temp_ts) > 1) {
#     if (is.null(mean_lfc[[LIG]])) {
#       mean_lfc[[LIG]] <- data.frame(ts_lfc=rowMeans(nn_lig_rep[[LIG]]$lfc[,temp_ts]))
#     } else {
#       mean_lfc[[LIG]]$ts_lfc <- rowMeans(nn_lig_rep[[LIG]]$lfc[,temp_ts])
#     }
#     mean_lfc[[LIG]]$ts_pval <- pbsapply(mean_lfc[[LIG]]$ts_lfc,function(X) 
#       sum(mean_lfc_ts_bkgd[,as.character(sum(temp_ts))] >= X) / nrow(mean_lfc_ts_bkgd),
#       cl=8)
#     mean_lfc[[LIG]]$ts_pval[mean_lfc[[LIG]]$ts_pval == 0] <- 1 / nrow(mean_lfc_ts_bkgd)
#     mean_lfc[[LIG]]$ts_fdr <- p.adjust(mean_lfc[[LIG]]$ts_pval,"fdr")
#   }
# }
# 
# 
# save(mean_lfc,mean_lfc_CT,mean_lfc_DS,
#      file="~/Dropbox/GDB_archive/CMapCorr_files/NN_ALLvCT_DEmeanlfc.RData")

