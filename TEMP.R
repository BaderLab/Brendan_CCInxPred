library(pbapply)
pboptions(type="timer")

allCuts <- seq(4,0,by=-0.01)

for (GSEA in c("REACT","GOBP")) {
  load(paste0("~/Dropbox/GDB_archive/CMapCorr_files/pfdrGSEA_",GSEA,"_DS.RData"))
  print(GSEA)
  
  for (X in grep("Qscore",ls(),value=T)) {
    assign(paste0("Qall_",strsplit(X,"_")[[1]][2]),get(X))
  }
  rm(list=grep("Qscore",ls(),value=T))
  
  for (X in grep("NQhits",ls(),value=T)) {
    assign(paste0("NQall_",strsplit(X,"_")[[1]][2]),get(X))
  }
  rm(list=grep("NQhits",ls(),value=T))
  
  
  for (X in grep("Pscore",ls(),value=T)) {
    tempQ <- pbsapply(get(X),function(Y)
      apply(Y,2,p.adjust,method="fdr"),
      simplify=F)
    assign(paste0("Qper_",strsplit(X,"_")[[1]][2]),tempQ)
    
    tempN <- pbsapply(10^-allCuts,function(CUT) 
      sapply(tempQ,function(X) mean(X <= CUT)))
    assign(paste0("NQper_",strsplit(X,"_")[[1]][2]),tempN)
    rm(tempQ,tempN)
  }
  rm(X)
  
  save(list=ls()[!ls() %in% c("allCuts","GSEA")],
       file=paste0("~/Dropbox/GDB_archive/CMapCorr_files/pfdrGSEA_",GSEA,"_DS.RData"))
}