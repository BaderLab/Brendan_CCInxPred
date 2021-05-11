library(pbapply)
library(fgsea)
library(fastmatch)
library(colorspace)
.PAR <- par(no.readonly=T)
pboptions(type="timer")

# Create pathways lists ----
# http://www.baderlab.org/GeneSets
# ^ ALL ----
temp_file <- "~/Dropbox/GDB_archive/Human_GOBP_AllPathways_no_GO_iea_March_01_2021_entrezgene.gmt" 
temp_lines <- strsplit(readLines(temp_file), "\t")
pathways_ALL <- lapply(temp_lines,tail,-2)
pathways_ALL_names <- sapply(temp_lines,function(X) X[2])
names(pathways_ALL_names) <- names(pathways_ALL) <- sapply(temp_lines,head,1)

# ^ GO biological process ----
temp_file <- "~/Dropbox/GDB_archive/Human_GO_bp_no_GO_iea_entrezgene.gmt" # http://www.baderlab.org/GeneSets
temp_lines <- strsplit(readLines(temp_file), "\t")
pathways_GOBP <- lapply(temp_lines,tail,-2)
pathways_GOBP_names <- sapply(temp_lines,function(X) X[2])
names(pathways_GOBP_names) <- names(pathways_GOBP) <- sapply(temp_lines,head,1)

# ^ Reactome ----
temp_file <- "~/Dropbox/GDB_archive/Human_Reactome_March_01_2021_Entrezgene.gmt" # http://www.baderlab.org/GeneSets
temp_lines <- strsplit(readLines(temp_file), "\t")
pathways_REACT <- lapply(temp_lines,tail,-2)
pathways_REACT_names <- sapply(temp_lines,function(X) X[2])
names(pathways_REACT_names) <- names(pathways_REACT) <- sapply(temp_lines,head,1)

rm(list=grep("^temp",ls(),value=T))


# load lvl5_data ----
load("~/Dropbox/GDB_archive/CMapCorr_files/lvl5_inputs_allgenes.RData")
temp_ligcountsperct <- sapply(unique(lvl5_data@cdesc$cell_id),function(CT)
  length(unique(lvl5_data@cdesc[lvl5_data@cdesc$cell_id == CT,"pert_iname"])))
ct9 <- names(temp_ligcountsperct)[temp_ligcountsperct > 100]
names(ct9) <- sapply(ct9,function(X) names(ct14)[ct14 == X])

temp_ctcountsperlig <- sapply(unique(lvl5_data@cdesc$pert_iname),function(LIG)
  length(unique(lvl5_data@cdesc[lvl5_data@cdesc$pert_iname == LIG &
                                  lvl5_data@cdesc$cell_id %in% ct9,"cell_id"])))
lig295 <- names(temp_ctcountsperlig)[temp_ctcountsperlig == 9]
lig295 <- sort(lig295)

temp_id <- rownames(lvl5_data@cdesc)[lvl5_data@cdesc$pert_iname %in% lig295 & 
                                        lvl5_data@cdesc$cell_id %in% ct9]
lvl5_data@mat <- lvl5_data@mat[,temp_id]
lvl5_data@cdesc <- lvl5_data@cdesc[temp_id,]
lvl5_data@cid <- temp_id
rm(list=c("lig16","ct14",grep("^temp",ls(),value=T)))


# filter pathways ----
pathways_REACT <- sapply(pathways_REACT,function(X) 
  X[X %in% lvl5_data@rid],simplify=F)
pathways_REACT <- pathways_REACT[sapply(pathways_REACT,length) >= 10 & 
                                   sapply(pathways_REACT,length) <= 200]
pathways_REACT_names <- pathways_REACT_names[names(pathways_REACT)]

pathways_GOBP <- sapply(pathways_GOBP,function(X) 
  X[X %in% lvl5_data@rid],simplify=F)
pathways_GOBP <- pathways_GOBP[sapply(pathways_GOBP,length) >= 10 & 
                                 sapply(pathways_GOBP,length) <= 200]
pathways_GOBP_names <- pathways_GOBP_names[names(pathways_GOBP)]

pathways_ALL <- sapply(pathways_ALL,function(X) 
  X[X %in% lvl5_data@rid],simplify=F)
pathways_ALL <- pathways_ALL[sapply(pathways_ALL,length) >= 10 & 
                               sapply(pathways_ALL,length) <= 200]
pathways_ALL_names <- pathways_ALL_names[names(pathways_ALL)]





# gene set permutation background calc ----
TEST <- lvl5_data@cid[lvl5_data@cdesc$cell_id == ct9[1] | 
                        lvl5_data@cdesc$pert_iname == lig295[1]]


ES <- pbapply(lvl5_data@mat[,TEST],2,function(DATA) {
  temp_DATA <- sort(DATA,decreasing=T)
  sapply(pathways_GOBP,function(X) 
    calcGseaStat(stats=temp_DATA,
                 selectedStats=na.omit(fmatch(X,names(temp_DATA))),
                 scoreType="std"))
})



  

# testing gene set permutation ----
# Do we need to do permutations per sample, or 
# can all samples serve as a background distribution for any sample?

temp_genes <- lvl5_data@rid[lvl5_data@rid %in% unlist(pathways_GOBP)]
temp_setsize <- sapply(pathways_GOBP,length)

temp_sets <- sapply(1:100,function(I)
  sapply(temp_setsize,function(X) sample(temp_genes,X)),
  simplify=F)

temp_es <- pbsapply(temp_sets,function(REP) {
  apply(lvl5_data@mat[,TEST[1:5]],2,function(DATA) {
    temp_DATA <- sort(DATA,decreasing=T)
    sapply(REP,function(SET) 
      calcGseaStat(stats=temp_DATA,
                   selectedStats=na.omit(fmatch(SET,names(temp_DATA))),
                   scoreType="std")
    )
  })
},simplify=F)

sample_bkgd <- sapply(1:ncol(temp_es[[1]]),function(ID) {
  temp <- sapply(temp_es,function(X) X[,ID])
  data.frame(MEAN=rowMeans(temp),
             SD=apply(temp,1,sd),
             COL=rep(
               qualitative_hcl(ncol(temp_es[[1]]),
                               palette="dark3",
                               alpha=0.5),
               nrow(temp)))
},simplify=F)
sample_bkgd <- do.call(rbind,sample_bkgd)
sample_bkgd <- sample_bkgd[sample(nrow(sample_bkgd)),]

par(mar=c(3,3,1,1),mgp=2:0)
plot(SD~MEAN,data=sample_bkgd,
     pch=".",cex=2,col=sample_bkgd$COL,
     xlab="MEAN",ylab="SD")
# doesn't appear to be any sample-related structure?


