library(cmapR)
library(ranger)
library(pbapply)

# prepping data ---
if (exists("lvl4_data")) {
} else if (file.exists("../CMapCorr_files/lvl4_inputs.RData")) {
  load("../CMapCorr_files/lvl4_inputs.RData") 
} else {
  source("lvl4_inputs.R")
}
temp <- load("../CMapCorr_files/lvl5_inputs.RData") 
rm(list=c("temp",temp[!temp %in% c("ct14","lig16")]))

