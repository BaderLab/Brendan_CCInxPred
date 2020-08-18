if (exists("lvl3_data")) {
} else if (file.exists("~/Dropbox/GDB_archive/CMapCorr_files/lvl3_inputs.RData")) {
  load("~/Dropbox/GDB_archive/CMapCorr_files/lvl3_inputs.RData") 
} else {
  source("lvl3_inputs.R")
}



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



par(mfrow=c(7,2),mar=c(3,3,2,1),mgp=2:0)
for (X in PM) {
  plot(sapply(sapply(X,"[[","ctl",simplify=F),length),
       sapply(sapply(X,"[[","tx",simplify=F),length),
       xlab="CTRL",ylab="TX",pch=20)
}
