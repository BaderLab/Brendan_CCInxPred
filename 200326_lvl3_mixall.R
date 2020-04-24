library(cmapR)
library(ranger)

if (!exists("lvl3_data")) {
  setwd("~/Dropbox/GDB/CMapCorr/")
  cmap_file_path <- "~/Data_LINCS/phase1/"
  lvl3_coldata <- read.table(file.path(cmap_file_path,"GSE92742_Broad_LINCS_inst_info.txt"),
                             header=T,sep="\t",row.names=1,colClasses="character",quote="\"")
  # always check # of rows matches .txt with `awk -F "\t" 'END {print NR}' filename`
  # and when it doesn't, make sure R isn't quoting out things based on apostrophes.
  temp_geneinfo <- read.table(file.path(cmap_file_path,"GSE92742_Broad_LINCS_gene_info.txt"),
                              header=T,sep="\t",row.names=1,colClasses="character",quote="\"")
  # Only going to use the landmark genes, because those were actually measured.
  # This reduces the number of features, such that features don't massively outnumber samples.
  
  lvl3_data <- parse.gctx(
    file.path(cmap_file_path,"GSE92742_Broad_LINCS_Level3_INF_mlr12k_n1319138x12328.gctx"),
    cid=rownames(lvl3_coldata)[lvl3_coldata$pert_type == "trt_lig"],
    rid=rownames(temp_geneinfo)[temp_geneinfo$pr_is_lm == "1"])
  temp_id <- lvl3_data@cdesc$id
  lvl3_data@cdesc <- lvl3_coldata[temp_id,]
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
    file.path(cmap_file_path,"GSE92742_Broad_LINCS_Level3_INF_mlr12k_n1319138x12328.gctx"),
    cid=rownames(lvl3_coldata)[lvl3_coldata$pert_type == "ctl_vehicle" &
                                 lvl3_coldata$cell_id %in% 
                                 unique(lvl3_coldata[rownames(lvl3_data@cdesc),"cell_id"])],
    rid=rownames(temp_geneinfo)[temp_geneinfo$pr_is_lm == "1"])
  temp_id <- lvl3_data_ctl@cdesc$id
  lvl3_data_ctl@cdesc <- lvl3_coldata[temp_id,]
  # fixing floating-point bullshit
  lvl3_data_ctl@cdesc$pert_dose[lvl3_data_ctl@cdesc$pert_dose == "0.10000000149"] <- "0.1"
  lvl3_data_ctl@cdesc$pert_dose <- sub("\\.?0+$","",lvl3_data_ctl@cdesc$pert_dose)
  
  temp_ct <- sapply(unique(lvl3_data@cdesc$cell_id),function(X) 
    unique(lvl3_data@cdesc$pert_iname[lvl3_data@cdesc$cell_id == X]),
    simplify=F)
  lig15 <- Reduce(intersect,temp_ct)

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
} 


train_labels <- train_labelsB <- test_labels <- 
  test_results <- test_resultsW <- test_resultsB <-
  model <- modelW <- modelB <- list()
for (L in lig15) {
  message(paste0("---- Ligand ",which(lig15 == L),"/",length(lig15)," ----"))
  tx <- lvl3_data@mat[,lvl3_data@cdesc$pert_iname == L]
  tx_train <- tx[,sample(colnames(tx),round(ncol(tx) / 2))]
  tx_test <- tx[,!colnames(tx) %in% colnames(tx_train)]
  ctl_train <- lvl3_data_ctl@mat[,sample(colnames(lvl3_data_ctl@mat),
                                         round(ncol(lvl3_data_ctl@mat) / 2))]
  ctl_test <- lvl3_data_ctl@mat[,!colnames(lvl3_data_ctl@mat) %in% colnames(ctl_train)]
  train <- t(cbind(ctl_train,tx_train))
  train_labels[[L]] <- as.factor(c(rep("no",ncol(ctl_train)),
                                   rep("yes",ncol(tx_train))))
  test <- t(cbind(ctl_test,tx_test))
  test_labels[[L]] <- as.factor(c(rep("no",ncol(ctl_test)),
                                  rep("yes",ncol(tx_test))))
  
  model[[L]] <- ranger(x=train,
                       y=train_labels[[L]],
                       num.threads=8,
                       verbose=F)
  modelW[[L]] <- ranger(x=train,
                        y=train_labels[[L]],
                        class.weights=length(train_labels[[L]]) /
                          (length(unique(train_labels[[L]])) * 
                             table(train_labels[[L]])),
                        # sklearn.ensemble.RandomForestClassifier(class_weight="balanced")
                        num.threads=8,
                        verbose=F)
  test_results[[L]] <- predict(model[[L]],test)$predictions
  test_resultsW[[L]] <- predict(modelW[[L]],test)$predictions
  
  
  dwn <- sample(colnames(ctl_train),ncol(tx_train))
  train <- t(cbind(ctl_train[,dwn],tx_train))
  train_labelsB[[L]] <- as.factor(c(rep("no",length(dwn)),
                                    rep("yes",ncol(tx_train))))
  modelB[[L]] <- ranger(x=train,
                        y=train_labelsB[[L]],
                        num.threads=8,
                        verbose=F)
  test_resultsB[[L]] <- predict(modelB[[L]],test)$predictions
}

save(train_labels,train_labelsB,test_labels,
     test_results,test_resultsW,test_resultsB,
     model,modelW,modelB,
     file="~/Dropbox/GDB/CMapCorr_files/200326_mixall.RData")
