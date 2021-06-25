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
  if (is.list(temp_DSname)) { stop(LIG) }
  DSinfo[[LIG]] <- nn_DEinfo[temp_DSname,c("time_series","cell_type","ligand")]
  DSinfo[[LIG]]$accession <- sub(paste0("-",LIG,".*$"),"",rownames(DSinfo[[LIG]]))
  rownames(DSinfo[[LIG]]) <- names(temp_DSname)
  rm(list=grep("^temp",ls(),value=T))
}
# DSinfo$WNT1$cell_type <- paste(DSinfo$WNT1$cell_type,
#                                sapply(strsplit(rownames(DSinfo$WNT1),"_"),"[[",3),
#                                sep="_")
for (LIG in names(DSinfo)) {
  DSinfo[[LIG]]$CtAcc <- paste(DSinfo[[LIG]]$accession,DSinfo[[LIG]]$cell_type)
}


DSinfo <- sapply(DSinfo,function(X) X[!X$time_series,],simplify=F)
DSinfo <- DSinfo[sapply(DSinfo,nrow) > 1]
nn_ligands <- names(DSinfo)

nn_lig_dsNames <- sapply(DSinfo,rownames)


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
     file="~/Dropbox/GDB_archive/CMapCorr_files/NN_noTS_dat.RData")



# correlation calc ----
corr_all <- corr_ct <- corr_ds <- list()
for (LIG in nn_ligands) {
  temp_ts <- rownames(DSinfo[[LIG]])[DSinfo[[LIG]]$time_series]
  temp_nts <- rownames(DSinfo[[LIG]])[!DSinfo[[LIG]]$time_series]
  
  if (length(temp_ts) > 1) {
    temp_acc <- sapply(unique(DSinfo[[LIG]][temp_ts,"CtAcc"]),function(X) 
      temp_ts[DSinfo[[LIG]][temp_ts,"CtAcc"] == X],
      simplify=F)
    temp_ct <- sapply(unique(DSinfo[[LIG]][temp_ts,"cell_type"]),function(X) 
      temp_ts[DSinfo[[LIG]][temp_ts,"cell_type"] == X],
      simplify=F)
    
    temp_cor2 <- cor(nn_lig_rep[[LIG]]$lfc[,temp_ts],method="spearman")
    temp_cor_names <- sapply(colnames(temp_cor2),function(X) 
      sapply(rownames(temp_cor2), function(Y) paste(X,Y,sep=".")))
    temp_cor <- temp_cor2[lower.tri(temp_cor2)]
    names(temp_cor) <- temp_cor_names[lower.tri(temp_cor_names)]
    rm(temp_cor2,temp_cor_names)
    
    for (Y in names(temp_acc)[sapply(temp_acc,length) > 1]) {
      temp_hits <- sapply(strsplit(names(temp_cor),".",fixed=T),function(X) all(X %in% temp_acc[[Y]]))
      corr_ds[[LIG]] <- append(corr_ds[[LIG]],temp_cor[temp_hits])
      temp_cor <- temp_cor[!temp_hits]
      rm(temp_hits)
    }
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
    rm(Y,temp_cor,temp_acc,temp_ct)
  }
  
  if (length(temp_nts) > 1) {
    temp_acc <- sapply(unique(DSinfo[[LIG]][temp_nts,"CtAcc"]),function(X) 
      temp_nts[DSinfo[[LIG]][temp_nts,"CtAcc"] == X],
      simplify=F)
    temp_ct <- sapply(unique(DSinfo[[LIG]][temp_nts,"cell_type"]),function(X) 
      temp_nts[DSinfo[[LIG]][temp_nts,"cell_type"] == X],
      simplify=F)
    
    temp_cor2 <- cor(nn_lig_rep[[LIG]]$lfc[,temp_nts],method="spearman")
    temp_cor_names <- sapply(colnames(temp_cor2),function(X) 
      sapply(rownames(temp_cor2), function(Y) paste(X,Y,sep=".")))
    temp_cor <- temp_cor2[lower.tri(temp_cor2)]
    names(temp_cor) <- temp_cor_names[lower.tri(temp_cor_names)]
    rm(temp_cor2,temp_cor_names)
    
    for (Y in names(temp_acc)[sapply(temp_acc,length) > 1]) {
      temp_hits <- sapply(strsplit(names(temp_cor),".",fixed=T),function(X) all(X %in% temp_acc[[Y]]))
      corr_ds[[LIG]] <- append(corr_ds[[LIG]],temp_cor[temp_hits])
      temp_cor <- temp_cor[!temp_hits]
      rm(temp_hits)
    }
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
    rm(Y,temp_cor,temp_acc,temp_ct)
  }
  rm(temp_nts,temp_ts)
}
rm(LIG)

save(corr_all,corr_ct,corr_ds,
     file="~/Dropbox/GDB_archive/CMapCorr_files/NN_noTS_corr.RData")



# DE overlap background ----
nn_DE <- list()
temp_commongenes <- Reduce(intersect,sapply(nn_lig_rep,function(X) rownames(X$lfc)))
for (LIG in nn_ligands) {
  nn_DE[[LIG]] <- mapply(intersect,
                         apply(nn_lig_rep[[LIG]]$lfc,2,function(X) names(which(X >= 1))),
                         apply(nn_lig_rep[[LIG]]$qval,2,function(X) names(which(X <= 0.1))))
  nn_DE[[LIG]] <- sapply(nn_DE[[LIG]],function(X) X[X %in% temp_commongenes])
}

par(mar=c(3,3,1,1),mgp=2:0)
hist(sapply(unlist(nn_DE,recursive=F),length),
     xlab="# DE per dataset (logFC >= 1, FDR <= 0.1)",main=NA)

temp_DE <- unlist(nn_DE,recursive=F)
DEbkgd <- pbsapply(min(sapply(DSinfo,nrow)):max(sapply(DSinfo,nrow)),function(X)
  sapply(1:1e6,function(Y)
    length(Reduce(intersect,sample(temp_DE,X)))
  ),cl=4)
colnames(DEbkgd) <- min(sapply(DSinfo,nrow)):max(sapply(DSinfo,nrow))
rm(temp_DE)


# DE overlap calc ----
de_all <- de_ct <- de_ds <- list()
P_all <- P_ct <- P_ds <- list()
for (LIG in nn_ligands) {
  temp_acc <- sapply(unique(DSinfo[[LIG]]$CtAcc),function(X) 
    rownames(DSinfo[[LIG]])[DSinfo[[LIG]]$CtAcc == X],
    simplify=F)
  temp_ct <- sapply(unique(DSinfo[[LIG]]$cell_type),function(X) 
    rownames(DSinfo[[LIG]])[DSinfo[[LIG]]$cell_type == X],
    simplify=F)
  temp_ct <- temp_ct[!sapply(temp_ct,function(CT) 
    any(sapply(temp_acc,function(ACC) all(CT == ACC))))]
  
  for (X in names(temp_acc)[sapply(temp_acc,length) > 1]) {
    Y <- paste(temp_acc[[X]],collapse=".")
    de_ds[[LIG]][[Y]] <- Reduce(intersect,nn_DE[[LIG]][temp_acc[[X]]])
    P_ds[[LIG]][[Y]] <- sum(DEbkgd[,as.character(length(temp_acc[[X]]))] >= 
                              length(de_ds[[LIG]][[Y]])) / nrow(DEbkgd)
    P_ds[[LIG]][[Y]][P_ds[[LIG]][[Y]] == 0] <- 0.1 / nrow(DEbkgd)
  }
  for (X in names(temp_ct)[sapply(temp_ct,length) > 1]) {
    Y <- paste(temp_ct[[X]],collapse=".")
    de_ct[[LIG]][[Y]] <- Reduce(intersect,nn_DE[[LIG]][temp_ct[[X]]])
    P_ct[[LIG]][[Y]] <- sum(DEbkgd[,as.character(length(temp_ct[[X]]))] >= 
                              length(de_ct[[LIG]][[Y]])) / nrow(DEbkgd)
    P_ct[[LIG]][[Y]][P_ct[[LIG]][[Y]] == 0] <- 0.1 / nrow(DEbkgd)
  }
  Y <- paste(rownames(DSinfo[[LIG]]),collapse=".")
  if (!Y %in% c(names(de_ct[[LIG]]),names(de_ds[[LIG]]))) {
    de_all[[LIG]][[Y]] <- Reduce(intersect,nn_DE[[LIG]])
    P_all[[LIG]][[Y]] <- sum(DEbkgd[,as.character(length(nn_DE[[LIG]]))] >= 
                               length(de_all[[LIG]][[Y]])) / nrow(DEbkgd)
    P_all[[LIG]][[Y]][P_all[[LIG]][[Y]] == 0] <- 0.1 / nrow(DEbkgd)
  }
}
rm(list=c("X","Y","LIG",grep("^temp",ls(),value=T)))

save(nn_DE,de_all,de_ct,de_ds,P_all,P_ct,P_ds,
     file="~/Dropbox/GDB_archive/CMapCorr_files/NN_noTS_DEovlp.RData")


# # DE by mean LFC background ----
# temp_nts <- unlist(sapply(DSinfo,rownames))
# temp_nts_genes <- Reduce(intersect,sapply(nn_db[temp_nts],function(X) X$diffexp$gene))
# mean_lfc_nts_bkgd <- pbsapply(
#   2:max(sapply(DSinfo,function(X) sum(!X$time_series))),
#   function(X)
#     sapply(1:1e2,function(Y)
#       rowMeans(sapply(nn_db[sample(temp_nts,X)],
#                       function(X) {
#                         temp <- X$diffexp$lfc
#                         names(temp) <- X$diffexp$gene
#                         return(temp[temp_nts_genes])
#                       }))
#     ),cl=3)
# colnames(mean_lfc_nts_bkgd) <- 2:max(sapply(DSinfo,function(X) sum(!X$time_series)))
# 
# 
# # DE by mean LFC ----
# mean_lfc <- list()
# for (LIG in nn_ligands) {
#   print(paste(which(nn_ligands == LIG),"/",length(nn_ligands)))
#   
#   mean_lfc[[LIG]] <- data.frame(nts_lfc=rowMeans(nn_lig_rep[[LIG]]$lfc))
#   mean_lfc[[LIG]]$nts_pval <- pbsapply(mean_lfc[[LIG]]$nts_lfc,function(X) 
#     sum(mean_lfc_nts_bkgd[,as.character(nrow(DSinfo[[LIG]]))] >= X) / nrow(mean_lfc_nts_bkgd),
#     cl=12)
#   mean_lfc[[LIG]]$nts_fdr <- p.adjust(mean_lfc[[LIG]]$nts_pval,"fdr")
# }
# 
# 
# save(mean_lfc,
#      file="~/Dropbox/GDB_archive/CMapCorr_files/NN_noTS_DEmeanlfc.RData")
