library(cmapR)
library(pbapply)
setwd("~/Dropbox/GDB/CMapCorr/")

# LVL4 ----

if (exists("lvl4_data")) {
} else if (file.exists("../CMapCorr_files/lvl4_inputs.RData")) {
  load("../CMapCorr_files/lvl4_inputs.RData") 
} else {
  source("lvl4_inputs.R")
}


lvl4_de <- apply(lvl4_data@mat,2,function(X) names(X)[X > 1.645]) # 2.326 for 99th %ile

# ^ CT ----
if (!file.exists("~/Dropbox/GDB/CMapCorr_files/200702_lvl4_OVLP_ct.RData")) {
  lvl4_OVLPact_ct <- lvl4_OVLPexp_ct <- lvl4_OVLPprob_ct <- list()
  for (CT in ct14) {
    message(paste(which(ct14 == CT),"/",length(ct14)))
    temp_ids <- rownames(lvl4_data@cdesc)[lvl4_data@cdesc$cell_id == CT]
    temp_expected <- temp_actual <- temp_prob <- list()
    temp_threestrikes <- 0
    for (NSET in 2:length(temp_ids)) {
      
      if (choose(length(temp_ids),NSET) > 3e6) {
        temp_combos <- list()
        for (Z in 1:1000) {
          repeat {
            temp <- sample(temp_ids,NSET)
            if (!any(sapply(sapply(temp_combos,setdiff,y=temp),length) == 0)) {
              break
            }
          }
          temp_combos[[Z]] <- temp
        }
        temp_combos <- as.data.frame(temp_combos,stringsAsFactors=F)
      } else {
        temp_combos <- combn(temp_ids,NSET)
        if (ncol(temp_combos) > 1000) {
          temp_combos <- temp_combos[,sample(ncol(temp_combos),1000)]
        }
      }
      
      temp_expected[[as.character(NSET)]] <- apply(temp_combos,2,function(IDS)
        do.call(prod,as.list(sapply(lvl4_de[IDS],length) / nrow(lvl4_data@mat))) * nrow(lvl4_data@mat))
      temp_actual[[as.character(NSET)]] <- apply(temp_combos,2,function(IDS) 
        length(Reduce(intersect,lvl4_de[IDS])))
      temp_prob[[as.character(NSET)]] <-   apply(temp_combos,2,function(IDS) {
        temp <- length(Reduce(intersect,lvl4_de[IDS]))
        if (temp > 0) {
          SuperExactTest::cpsets(x=temp - 1,
                                 L=sapply(lvl4_de[IDS],length),
                                 n=nrow(lvl4_data@mat),
                                 lower.tail=F)
        } else {
          return(1)
        }
      })
      
      if (all(temp_actual[[as.character(NSET)]] == 0)) {
        temp_threestrikes <- temp_threestrikes + 1
      }
      if (temp_threestrikes >= 2) { break }
    }
    lvl4_OVLPexp_ct[[CT]] <- as.data.frame(temp_expected)
    lvl4_OVLPact_ct[[CT]] <- as.data.frame(temp_actual)
    lvl4_OVLPprob_ct[[CT]] <- as.data.frame(temp_prob)
  }
  save(lvl4_OVLPact_ct,lvl4_OVLPexp_ct,lvl4_OVLPprob_ct,
       file="~/Dropbox/GDB/CMapCorr_files/200702_lvl4_OVLP_ct.RData")
}

# ^ LIG ----
if (!file.exists("~/Dropbox/GDB/CMapCorr_files/200702_lvl4_OVLP_lig.RData")) {
  lvl4_OVLPact_lig <- lvl4_OVLPexp_lig <- lvl4_OVLPprob_lig <- list()
  for (LIG in lig16) {
    message(paste(which(lig16 == LIG),"/",length(lig16)))
    temp_ids <- rownames(lvl4_data@cdesc)[lvl4_data@cdesc$pert_iname == LIG]
    temp_expected <- temp_actual <- temp_prob <- list()
    temp_threestrikes <- 0
    for (NSET in 2:length(temp_ids)) {
      
      if (choose(length(temp_ids),NSET) > 3e6) {
        temp_combos <- list()
        for (Z in 1:1000) {
          repeat {
            temp <- sample(temp_ids,NSET)
            if (!any(sapply(sapply(temp_combos,setdiff,y=temp),length) == 0)) {
              break
            }
          }
          temp_combos[[Z]] <- temp
        }
        temp_combos <- as.data.frame(temp_combos,stringsAsFactors=F)
      } else {
        temp_combos <- combn(temp_ids,NSET)
        if (ncol(temp_combos) > 1000) {
          temp_combos <- temp_combos[,sample(ncol(temp_combos),1000)]
        }
      }
      
      temp_expected[[as.character(NSET)]] <- apply(temp_combos,2,function(IDS)
        do.call(prod,as.list(sapply(lvl4_de[IDS],length) / nrow(lvl4_data@mat))) * nrow(lvl4_data@mat))
      temp_actual[[as.character(NSET)]] <- apply(temp_combos,2,function(IDS) 
        length(Reduce(intersect,lvl4_de[IDS])))
      temp_prob[[as.character(NSET)]] <-   apply(temp_combos,2,function(IDS) {
        temp <- length(Reduce(intersect,lvl4_de[IDS]))
        if (temp > 0) {
          SuperExactTest::cpsets(x=temp - 1,
                                 L=sapply(lvl4_de[IDS],length),
                                 n=nrow(lvl4_data@mat),
                                 lower.tail=F)
        } else {
          return(1)
        }
      })
      
      if (all(temp_actual[[as.character(NSET)]] == 0)) {
        temp_threestrikes <- temp_threestrikes + 1
      }
      if (temp_threestrikes >= 2) { break }
    }
    lvl4_OVLPexp_lig[[LIG]] <- as.data.frame(temp_expected)
    lvl4_OVLPact_lig[[LIG]] <- as.data.frame(temp_actual)
    lvl4_OVLPprob_lig[[LIG]] <- as.data.frame(temp_prob)
  }
  save(lvl4_OVLPact_lig,lvl4_OVLPexp_lig,lvl4_OVLPprob_lig,
       file="~/Dropbox/GDB/CMapCorr_files/200702_lvl4_OVLP_lig.RData")
}

# ^ LIG by CT ----
if (!file.exists("~/Dropbox/GDB/CMapCorr_files/200702_lvl4_OVLP_lig_ct.RData")) {
  lvl4_OVLPact_lig_ct <- lvl4_OVLPexp_lig_ct <- lvl4_OVLPprob_lig_ct <- list()
  for (LIG in lig16) {
    message(paste(which(lig16 == LIG),"/",length(lig16)))
    temp_expected <- temp_actual <- temp_prob <- list()
    for (CT in ct14) {
      temp_expected[[CT]] <- temp_actual[[CT]] <- temp_prob[[CT]] <- list()
      temp_ids <- rownames(lvl4_data@cdesc)[lvl4_data@cdesc$pert_iname == LIG &
                                              lvl4_data@cdesc$cell_id == CT]
      temp_threestrikes <- 0
      for (NSET in 2:length(temp_ids)) {
        if (choose(length(temp_ids),NSET) > 3e6) {
          temp_combos <- list()
          for (Z in 1:1000) {
            repeat {
              temp <- sample(temp_ids,NSET)
              if (!any(sapply(sapply(temp_combos,setdiff,y=temp),length) == 0)) {
                break
              }
            }
            temp_combos[[Z]] <- temp
          }
          temp_combos <- as.data.frame(temp_combos,stringsAsFactors=F)
        } else {
          temp_combos <- combn(temp_ids,NSET)
          if (ncol(temp_combos) > 1000) {
            temp_combos <- temp_combos[,sample(ncol(temp_combos),1000)]
          }
        }
        
        temp_expected[[CT]][[as.character(NSET)]] <- apply(temp_combos,2,function(IDS)
          do.call(prod,as.list(sapply(lvl4_de[IDS],length) / nrow(lvl4_data@mat))) * nrow(lvl4_data@mat))
        temp_actual[[CT]][[as.character(NSET)]] <- apply(temp_combos,2,function(IDS) 
          length(Reduce(intersect,lvl4_de[IDS])))
        temp_prob[[CT]][[as.character(NSET)]] <-   apply(temp_combos,2,function(IDS) {
          temp <- length(Reduce(intersect,lvl4_de[IDS]))
          if (temp > 0) {
            SuperExactTest::cpsets(x=temp - 1,
                                   L=sapply(lvl4_de[IDS],length),
                                   n=nrow(lvl4_data@mat),
                                   lower.tail=F)
          } else {
            return(1)
          }
        })
        
        if (all(temp_actual[[CT]][[as.character(NSET)]] == 0)) {
          temp_threestrikes <- temp_threestrikes + 1
        }
        if (temp_threestrikes >= 2) { break }
      }
    }
    lvl4_OVLPexp_lig_ct[[LIG]] <- temp_expected
    lvl4_OVLPact_lig_ct[[LIG]] <- temp_actual
    lvl4_OVLPprob_lig_ct[[LIG]] <- temp_prob
  }
  save(lvl4_OVLPact_lig_ct,lvl4_OVLPexp_lig_ct,lvl4_OVLPprob_lig_ct,
       file="~/Dropbox/GDB/CMapCorr_files/200702_lvl4_OVLP_lig_ct.RData")
}


# LVL5 ----

if (exists("lvl5_data")) {
} else if (file.exists("../CMapCorr_files/lvl5_inputs.RData")) {
  load("../CMapCorr_files/lvl5_inputs.RData") 
} else {
  source("lvl5_inputs.R")
}


lvl5_de <- apply(lvl5_data@mat,2,function(X) names(X)[X > 1.645]) # 2.326 for 99th %ile

# ^ CT ----
if (!file.exists("~/Dropbox/GDB/CMapCorr_files/200702_lvl5_OVLP_ct.RData")) {
  lvl5_OVLPact_ct <- lvl5_OVLPexp_ct <- lvl5_OVLPprob_ct <- list()
  for (CT in ct14) {
    message(paste(which(ct14 == CT),"/",length(ct14)))
    temp_ids <- rownames(lvl5_data@cdesc)[lvl5_data@cdesc$cell_id == CT]
    temp_expected <- temp_actual <- temp_prob <- list()
    temp_threestrikes <- 0
    for (NSET in 2:length(temp_ids)) {
      
      if (choose(length(temp_ids),NSET) > 3e6) {
        temp_combos <- list()
        for (Z in 1:1000) {
          repeat {
            temp <- sample(temp_ids,NSET)
            if (!any(sapply(sapply(temp_combos,setdiff,y=temp),length) == 0)) {
              break
            }
          }
          temp_combos[[Z]] <- temp
        }
        temp_combos <- as.data.frame(temp_combos,stringsAsFactors=F)
      } else {
        temp_combos <- combn(temp_ids,NSET)
        if (ncol(temp_combos) > 1000) {
          temp_combos <- temp_combos[,sample(ncol(temp_combos),1000)]
        }
      }
      
      temp_expected[[as.character(NSET)]] <- apply(temp_combos,2,function(IDS)
        do.call(prod,as.list(sapply(lvl5_de[IDS],length) / nrow(lvl5_data@mat))) * nrow(lvl5_data@mat))
      temp_actual[[as.character(NSET)]] <- apply(temp_combos,2,function(IDS) 
        length(Reduce(intersect,lvl5_de[IDS])))
      temp_prob[[as.character(NSET)]] <-   apply(temp_combos,2,function(IDS) {
        temp <- length(Reduce(intersect,lvl5_de[IDS]))
        if (temp > 0) {
          SuperExactTest::cpsets(x=temp - 1,
                                 L=sapply(lvl5_de[IDS],length),
                                 n=nrow(lvl5_data@mat),
                                 lower.tail=F)
        } else {
          return(1)
        }
      })
      
      if (all(temp_actual[[as.character(NSET)]] == 0)) {
        temp_threestrikes <- temp_threestrikes + 1
      }
      if (temp_threestrikes >= 2) { break }
    }
    lvl5_OVLPexp_ct[[CT]] <- as.data.frame(temp_expected)
    lvl5_OVLPact_ct[[CT]] <- as.data.frame(temp_actual)
    lvl5_OVLPprob_ct[[CT]] <- as.data.frame(temp_prob)
  }
  save(lvl5_OVLPact_ct,lvl5_OVLPexp_ct,lvl5_OVLPprob_ct,
       file="~/Dropbox/GDB/CMapCorr_files/200702_lvl5_OVLP_ct.RData")
}

# ^ LIG ----
if (!file.exists("~/Dropbox/GDB/CMapCorr_files/200702_lvl5_OVLP_lig.RData")) {
  lvl5_OVLPact_lig <- lvl5_OVLPexp_lig <- lvl5_OVLPprob_lig <- list()
  for (LIG in lig16) {
    message(paste(which(lig16 == LIG),"/",length(lig16)))
    temp_ids <- rownames(lvl5_data@cdesc)[lvl5_data@cdesc$pert_iname == LIG]
    temp_expected <- temp_actual <- temp_prob <- list()
    temp_threestrikes <- 0
    for (NSET in 2:length(temp_ids)) {
      
      if (choose(length(temp_ids),NSET) > 3e6) {
        temp_combos <- list()
        for (Z in 1:1000) {
          repeat {
            temp <- sample(temp_ids,NSET)
            if (!any(sapply(sapply(temp_combos,setdiff,y=temp),length) == 0)) {
              break
            }
          }
          temp_combos[[Z]] <- temp
        }
        temp_combos <- as.data.frame(temp_combos,stringsAsFactors=F)
      } else {
        temp_combos <- combn(temp_ids,NSET)
        if (ncol(temp_combos) > 1000) {
          temp_combos <- temp_combos[,sample(ncol(temp_combos),1000)]
        }
      }
      
      temp_expected[[as.character(NSET)]] <- apply(temp_combos,2,function(IDS)
        do.call(prod,as.list(sapply(lvl5_de[IDS],length) / nrow(lvl5_data@mat))) * nrow(lvl5_data@mat))
      temp_actual[[as.character(NSET)]] <- apply(temp_combos,2,function(IDS) 
        length(Reduce(intersect,lvl5_de[IDS])))
      temp_prob[[as.character(NSET)]] <-   apply(temp_combos,2,function(IDS) {
        temp <- length(Reduce(intersect,lvl5_de[IDS]))
        if (temp > 0) {
          SuperExactTest::cpsets(x=temp - 1,
                                 L=sapply(lvl5_de[IDS],length),
                                 n=nrow(lvl5_data@mat),
                                 lower.tail=F)
        } else {
          return(1)
        }
      })
      
      if (all(temp_actual[[as.character(NSET)]] == 0)) {
        temp_threestrikes <- temp_threestrikes + 1
      }
      if (temp_threestrikes >= 2) { break }
    }
    lvl5_OVLPexp_lig[[LIG]] <- as.data.frame(temp_expected)
    lvl5_OVLPact_lig[[LIG]] <- as.data.frame(temp_actual)
    lvl5_OVLPprob_lig[[LIG]] <- as.data.frame(temp_prob)
  }
  save(lvl5_OVLPact_lig,lvl5_OVLPexp_lig,lvl5_OVLPprob_lig,
       file="~/Dropbox/GDB/CMapCorr_files/200702_lvl5_OVLP_lig.RData")
}

# ^ LIG by CT ----
if (!file.exists("~/Dropbox/GDB/CMapCorr_files/200702_lvl5_OVLP_lig_ct.RData")) {
  lvl5_OVLPact_lig_ct <- lvl5_OVLPexp_lig_ct <- lvl5_OVLPprob_lig_ct <- list()
  for (LIG in lig16) {
    message(paste(which(lig16 == LIG),"/",length(lig16)))
    temp_expected <- temp_actual <- temp_prob <- list()
    for (CT in ct14) {
      temp_expected[[CT]] <- temp_actual[[CT]] <- temp_prob[[CT]] <- list()
      temp_ids <- rownames(lvl5_data@cdesc)[lvl5_data@cdesc$pert_iname == LIG &
                                              lvl5_data@cdesc$cell_id == CT]
      temp_threestrikes <- 0
      if (length(temp_ids) < 2) { next }
      for (NSET in 2:length(temp_ids)) {
        if (choose(length(temp_ids),NSET) > 3e6) {
          temp_combos <- list()
          for (Z in 1:1000) {
            repeat {
              temp <- sample(temp_ids,NSET)
              if (!any(sapply(sapply(temp_combos,setdiff,y=temp),length) == 0)) {
                break
              }
            }
            temp_combos[[Z]] <- temp
          }
          temp_combos <- as.data.frame(temp_combos,stringsAsFactors=F)
        } else {
          temp_combos <- combn(temp_ids,NSET)
          if (ncol(temp_combos) > 1000) {
            temp_combos <- temp_combos[,sample(ncol(temp_combos),1000)]
          }
        }
        
        temp_expected[[CT]][[as.character(NSET)]] <- apply(temp_combos,2,function(IDS)
          do.call(prod,as.list(sapply(lvl5_de[IDS],length) / nrow(lvl5_data@mat))) * nrow(lvl5_data@mat))
        temp_actual[[CT]][[as.character(NSET)]] <- apply(temp_combos,2,function(IDS) 
          length(Reduce(intersect,lvl5_de[IDS])))
        temp_prob[[CT]][[as.character(NSET)]] <-   apply(temp_combos,2,function(IDS) {
          temp <- length(Reduce(intersect,lvl5_de[IDS]))
          if (temp > 0) {
            SuperExactTest::cpsets(x=temp - 1,
                                   L=sapply(lvl5_de[IDS],length),
                                   n=nrow(lvl5_data@mat),
                                   lower.tail=F)
          } else {
            return(1)
          }
        })
        
        if (all(temp_actual[[CT]][[as.character(NSET)]] == 0)) {
          temp_threestrikes <- temp_threestrikes + 1
        }
        if (temp_threestrikes >= 2) { break }
      }
    }
    lvl5_OVLPexp_lig_ct[[LIG]] <- temp_expected
    lvl5_OVLPact_lig_ct[[LIG]] <- temp_actual
    lvl5_OVLPprob_lig_ct[[LIG]] <- temp_prob
  }
  save(lvl5_OVLPact_lig_ct,lvl5_OVLPexp_lig_ct,lvl5_OVLPprob_lig_ct,
       file="~/Dropbox/GDB/CMapCorr_files/200702_lvl5_OVLP_lig_ct.RData")
}



# lvl4new ----

if (exists("lvl4new_data")) {
} else if (file.exists("../CMapCorr_files/lvl4new.RData")) {
  load("../CMapCorr_files/lvl4new.RData") 
} else {
  source("200706_ZscoreFromAssayed.R")
}


lvl4new_de <- apply(lvl4new_data@mat,2,function(X) names(X)[X > 1.645]) # 2.326 for 99th %ile

# ^ CT ----
if (!file.exists("~/Dropbox/GDB/CMapCorr_files/200707_lvl4new_OVLP_ct.RData")) {
  lvl4new_OVLPact_ct <- lvl4new_OVLPexp_ct <- lvl4new_OVLPprob_ct <- list()
  for (CT in ct14) {
    message(paste(which(ct14 == CT),"/",length(ct14)))
    temp_ids <- rownames(lvl4new_data@cdesc)[lvl4new_data@cdesc$cell_id == CT]
    temp_expected <- temp_actual <- temp_prob <- list()
    temp_threestrikes <- 0
    for (NSET in 2:length(temp_ids)) {
      
      if (choose(length(temp_ids),NSET) > 3e6) {
        temp_combos <- list()
        for (Z in 1:1000) {
          repeat {
            temp <- sample(temp_ids,NSET)
            if (!any(sapply(sapply(temp_combos,setdiff,y=temp),length) == 0)) {
              break
            }
          }
          temp_combos[[Z]] <- temp
        }
        temp_combos <- as.data.frame(temp_combos,stringsAsFactors=F)
      } else {
        temp_combos <- combn(temp_ids,NSET)
        if (ncol(temp_combos) > 1000) {
          temp_combos <- temp_combos[,sample(ncol(temp_combos),1000)]
        }
      }
      
      temp_expected[[as.character(NSET)]] <- apply(temp_combos,2,function(IDS)
        do.call(prod,as.list(sapply(lvl4new_de[IDS],length) / nrow(lvl4new_data@mat))) * nrow(lvl4new_data@mat))
      temp_actual[[as.character(NSET)]] <- apply(temp_combos,2,function(IDS) 
        length(Reduce(intersect,lvl4new_de[IDS])))
      temp_prob[[as.character(NSET)]] <-   apply(temp_combos,2,function(IDS) {
        temp <- length(Reduce(intersect,lvl4new_de[IDS]))
        if (temp > 0) {
          SuperExactTest::cpsets(x=temp - 1,
                                 L=sapply(lvl4new_de[IDS],length),
                                 n=nrow(lvl4new_data@mat),
                                 lower.tail=F)
        } else {
          return(1)
        }
      })
      
      if (all(temp_actual[[as.character(NSET)]] == 0)) {
        temp_threestrikes <- temp_threestrikes + 1
      }
      if (temp_threestrikes >= 2) { break }
    }
    lvl4new_OVLPexp_ct[[CT]] <- as.data.frame(temp_expected)
    lvl4new_OVLPact_ct[[CT]] <- as.data.frame(temp_actual)
    lvl4new_OVLPprob_ct[[CT]] <- as.data.frame(temp_prob)
  }
  save(lvl4new_OVLPact_ct,lvl4new_OVLPexp_ct,lvl4new_OVLPprob_ct,
       file="~/Dropbox/GDB/CMapCorr_files/200707_lvl4new_OVLP_ct.RData")
}

# ^ LIG ----
if (!file.exists("~/Dropbox/GDB/CMapCorr_files/200707_lvl4new_OVLP_lig.RData")) {
  lvl4new_OVLPact_lig <- lvl4new_OVLPexp_lig <- lvl4new_OVLPprob_lig <- list()
  for (LIG in lig16) {
    message(paste(which(lig16 == LIG),"/",length(lig16)))
    temp_ids <- rownames(lvl4new_data@cdesc)[lvl4new_data@cdesc$pert_iname == LIG]
    temp_expected <- temp_actual <- temp_prob <- list()
    temp_threestrikes <- 0
    for (NSET in 2:length(temp_ids)) {
      
      if (choose(length(temp_ids),NSET) > 3e6) {
        temp_combos <- list()
        for (Z in 1:1000) {
          repeat {
            temp <- sample(temp_ids,NSET)
            if (!any(sapply(sapply(temp_combos,setdiff,y=temp),length) == 0)) {
              break
            }
          }
          temp_combos[[Z]] <- temp
        }
        temp_combos <- as.data.frame(temp_combos,stringsAsFactors=F)
      } else {
        temp_combos <- combn(temp_ids,NSET)
        if (ncol(temp_combos) > 1000) {
          temp_combos <- temp_combos[,sample(ncol(temp_combos),1000)]
        }
      }
      
      temp_expected[[as.character(NSET)]] <- apply(temp_combos,2,function(IDS)
        do.call(prod,as.list(sapply(lvl4new_de[IDS],length) / nrow(lvl4new_data@mat))) * nrow(lvl4new_data@mat))
      temp_actual[[as.character(NSET)]] <- apply(temp_combos,2,function(IDS) 
        length(Reduce(intersect,lvl4new_de[IDS])))
      temp_prob[[as.character(NSET)]] <-   apply(temp_combos,2,function(IDS) {
        temp <- length(Reduce(intersect,lvl4new_de[IDS]))
        if (temp > 0) {
          SuperExactTest::cpsets(x=temp - 1,
                                 L=sapply(lvl4new_de[IDS],length),
                                 n=nrow(lvl4new_data@mat),
                                 lower.tail=F)
        } else {
          return(1)
        }
      })
      
      if (all(temp_actual[[as.character(NSET)]] == 0)) {
        temp_threestrikes <- temp_threestrikes + 1
      }
      if (temp_threestrikes >= 2) { break }
    }
    lvl4new_OVLPexp_lig[[LIG]] <- as.data.frame(temp_expected)
    lvl4new_OVLPact_lig[[LIG]] <- as.data.frame(temp_actual)
    lvl4new_OVLPprob_lig[[LIG]] <- as.data.frame(temp_prob)
  }
  save(lvl4new_OVLPact_lig,lvl4new_OVLPexp_lig,lvl4new_OVLPprob_lig,
       file="~/Dropbox/GDB/CMapCorr_files/200707_lvl4new_OVLP_lig.RData")
}

# ^ LIG by CT ----
if (!file.exists("~/Dropbox/GDB/CMapCorr_files/200707_lvl4new_OVLP_lig_ct.RData")) {
  lvl4new_OVLPact_lig_ct <- lvl4new_OVLPexp_lig_ct <- lvl4new_OVLPprob_lig_ct <- list()
  for (LIG in lig16) {
    message(paste(which(lig16 == LIG),"/",length(lig16)))
    temp_expected <- temp_actual <- temp_prob <- list()
    for (CT in ct14) {
      temp_expected[[CT]] <- temp_actual[[CT]] <- temp_prob[[CT]] <- list()
      temp_ids <- rownames(lvl4new_data@cdesc)[lvl4new_data@cdesc$pert_iname == LIG &
                                              lvl4new_data@cdesc$cell_id == CT]
      temp_threestrikes <- 0
      for (NSET in 2:length(temp_ids)) {
        if (choose(length(temp_ids),NSET) > 3e6) {
          temp_combos <- list()
          for (Z in 1:1000) {
            repeat {
              temp <- sample(temp_ids,NSET)
              if (!any(sapply(sapply(temp_combos,setdiff,y=temp),length) == 0)) {
                break
              }
            }
            temp_combos[[Z]] <- temp
          }
          temp_combos <- as.data.frame(temp_combos,stringsAsFactors=F)
        } else {
          temp_combos <- combn(temp_ids,NSET)
          if (ncol(temp_combos) > 1000) {
            temp_combos <- temp_combos[,sample(ncol(temp_combos),1000)]
          }
        }
        
        temp_expected[[CT]][[as.character(NSET)]] <- apply(temp_combos,2,function(IDS)
          do.call(prod,as.list(sapply(lvl4new_de[IDS],length) / nrow(lvl4new_data@mat))) * nrow(lvl4new_data@mat))
        temp_actual[[CT]][[as.character(NSET)]] <- apply(temp_combos,2,function(IDS) 
          length(Reduce(intersect,lvl4new_de[IDS])))
        temp_prob[[CT]][[as.character(NSET)]] <-   apply(temp_combos,2,function(IDS) {
          temp <- length(Reduce(intersect,lvl4new_de[IDS]))
          if (temp > 0) {
            SuperExactTest::cpsets(x=temp - 1,
                                   L=sapply(lvl4new_de[IDS],length),
                                   n=nrow(lvl4new_data@mat),
                                   lower.tail=F)
          } else {
            return(1)
          }
        })
        
        if (all(temp_actual[[CT]][[as.character(NSET)]] == 0)) {
          temp_threestrikes <- temp_threestrikes + 1
        }
        if (temp_threestrikes >= 2) { break }
      }
    }
    lvl4new_OVLPexp_lig_ct[[LIG]] <- temp_expected
    lvl4new_OVLPact_lig_ct[[LIG]] <- temp_actual
    lvl4new_OVLPprob_lig_ct[[LIG]] <- temp_prob
  }
  save(lvl4new_OVLPact_lig_ct,lvl4new_OVLPexp_lig_ct,lvl4new_OVLPprob_lig_ct,
       file="~/Dropbox/GDB/CMapCorr_files/200707_lvl4new_OVLP_lig_ct.RData")
}
