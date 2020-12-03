library(colorspace)
library(pbapply)
pboptions(type="timer")
.PAR <- par(no.readonly=T)


# Re-analysis of NicheNet dataset ----

nn_model <- readRDS(url("https://zenodo.org/record/3260758/files/ligand_target_matrix.rds"))
nn_db <- readRDS(url("https://zenodo.org/record/3260758/files/expression_settings.rds"))

nn_ligands <- unique(unlist(sapply(nn_db,"[[","from")))
nn_lig_dsNames <- sapply(nn_ligands,function(LIG) names(nn_db)[sapply(nn_db,function(X) LIG %in% X$from)])
nn_ligands <- nn_ligands[sapply(nn_lig_dsNames,length) > 1]
nn_lig_dsNames <- nn_lig_dsNames[nn_ligands]



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
    return(list(lfc=temp[temp_gene,"lfc"],
                qval=temp[temp_gene,"qval"]))
  },simplify=F)
  
  temp_lfc <- sapply(temp_diffexp,"[[","lfc")
  temp_qval <- sapply(temp_diffexp,"[[","qval")
  rownames(temp_lfc) <- rownames(temp_qval) <- temp_gene
  
  return(list(lfc=temp_lfc,
              qval=temp_qval))
},simplify=F)

nn_cor <- pbsapply(nn_lig_rep,function(X)
  cor(X$lfc,method="spearman")[upper.tri(matrix(nrow=ncol(X$lfc),ncol=ncol(X$lfc)))],
  simplify=F)
boxplot(nn_cor,las=3,ylab="Spearman's CC")

# DE per NicheNet paper: lfc >= 1 & qval <= 0.1





# Prob DE overlap from Supplementary Table ----
DEcounts <- read.table("~/Dropbox/GDB/CMapCorr/NicheNet_LigDEoverlap_count.txt",
                       sep="\t",header=T,row.names=1,as.is=T)
DEratio <- read.table("~/Dropbox/GDB/CMapCorr/NicheNet_LigDEoverlap_ratio.txt",
                      sep="\t",header=T,row.names=1,as.is=T)
numDE <- rowSums(DEcounts,na.rm=T)
nGene <- c(length(Reduce(intersect,lapply(nn_db,function(X) X$diffexp$gene))),
           range(sapply(nn_db,function(X) nrow(X$diffexp))))



# Mapping datasets to GEO accession ----

nn_DEinfo <- read.table("~/Dropbox/GDB/CMapCorr/NicheNet_DBinfo_nDE.txt",
                        header=T,sep="\t",as.is=T)

DEgenes <- list()
for (LIG in names(nn_lig_dsNames)) {
  DEgenes[[LIG]] <- sapply(nn_lig_dsNames[[LIG]],function(X) 
    nn_db[[X]]$diffexp[nn_db[[X]]$diffexp$lfc >= 1 & nn_db[[X]]$diffexp$qval <= 0.1,"gene"])
  sapply(sapply(DEgenes[[LIG]],length),function(N) 
    c(nn_DEinfo[nn_DEinfo$ligand == LIG,"setting_id"][nn_DEinfo[nn_DEinfo$ligand == LIG,"nr_upgenes"] %in% N],
      nn_DEinfo[nn_DEinfo$ligand == LIG,"setting_id"][nn_DEinfo[nn_DEinfo$ligand == LIG,"nr_degenes"] %in% N]))
}

nn_DEinfo[nn_DEinfo$ligand == LIG,"nr_upgenes"]








