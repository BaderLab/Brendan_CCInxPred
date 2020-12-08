library(pbapply)
setwd("~/Dropbox/GDB_archive/CMapCorr_files/GSE47542_RAW/")

temp_TNF_p <- c("GSM1151945_Agilent_251485053931_S01_GE2_107_Sep09_1_4.txt.gz",
                "GSM1151946_Agilent_251485053932_S01_GE2_107_Sep09_1_4.txt.gz",
                "GSM1151950_Agilent_251485054074_S01_GE2_107_Sep09_1_4.txt.gz")
temp_TNF_m <- c("GSM1151933_Agilent_251485053931_S01_GE2_107_Sep09_1_2.txt.gz",
                "GSM1151934_Agilent_251485053932_S01_GE2_107_Sep09_1_2.txt.gz",
                "GSM1151938_Agilent_251485054074_S01_GE2_107_Sep09_1_2.txt.gz")
temp_CTRL_p <- c("GSM1151939_Agilent_251485053931_S01_GE2_107_Sep09_1_3.txt.gz",
                 "GSM1151940_Agilent_251485053932_S01_GE2_107_Sep09_1_3.txt.gz",
                 "GSM1151947_Agilent_251485054010_S01_GE2_107_Sep09_1_4.txt.gz")
temp_CTRL_m <- c("GSM1151928_Agilent_251485053932_S01_GE2_107_Sep09_1_1.txt.gz",
                 "GSM1151935_Agilent_251485054010_S01_GE2_107_Sep09_1_2.txt.gz",
                 "GSM1151951_Agilent_251485053931_S01_GE2_107_Sep09_1_1.txt.gz")

TNFp <- sapply(temp_TNF_p,function(X) {
  DS <- read.table(X,sep="\t",header=T,
                   skip=which(grepl("^FEAT",scan(X,"character",sep="\n"))) - 1)
  return(DS[DS$SystematicName == "NM_000201","LogRatio"])
},USE.NAMES=F)

TNFm <- sapply(temp_TNF_m,function(X) {
  DS <- read.table(X,sep="\t",header=T,
                   skip=which(grepl("^FEAT",scan(X,"character",sep="\n"))) - 1)
  return(DS[DS$SystematicName == "NM_000201","LogRatio"])
},USE.NAMES=F)

CTRLp <- sapply(temp_CTRL_p,function(X) {
  DS <- read.table(X,sep="\t",header=T,
                   skip=which(grepl("^FEAT",scan(X,"character",sep="\n"))) - 1)
  return(DS[DS$SystematicName == "NM_000201","LogRatio"])
},USE.NAMES=F)

CTRLm <- sapply(temp_CTRL_m,function(X) {
  DS <- read.table(X,sep="\t",header=T,
                   skip=which(grepl("^FEAT",scan(X,"character",sep="\n"))) - 1)
  return(DS[DS$SystematicName == "NM_000201","LogRatio"])
},USE.NAMES=F)


par(mar=c(3,3,3,1),mgp=2:0)
boxplot(list(TNFp=as.vector(TNFp),
             CTRLp=as.vector(CTRLp),
             TNFm=as.vector(TNFm),
             CTRLm=as.vector(CTRLm)),
        main="ICAM1 expression",
        ylab="Log-ratio vs reference RNA")
