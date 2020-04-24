temp1 <- read.table("~/Data_LINCS/phase1/GSE92742_Broad_LINCS_inst_info.txt",
                     header=T,sep="\t",row.names=1,colClasses="character",quote="\"")
temp2 <- read.table("~/Data_LINCS/phase2/GSE70138_Broad_LINCS_inst_info.txt",
                     header=T,sep="\t",row.names=1,colClasses="character",quote="\"")

table(temp1$pert_type)
table(temp2$pert_type)

rm(list=grep("^temp",ls(),value=T))
