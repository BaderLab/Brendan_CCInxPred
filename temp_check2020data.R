library(rhdf5)
library(cmapR)

lvl4_data <- parse_gctx("~/Data_LINCS/2020beta/level4_beta_trt_misc_n26428x12328.gctx")

# always check # of rows matches .txt with `awk -F "\t" 'END {print NR}' filename`
# and when it doesn't, make sure R isn't quoting out things based on apostrophes.
info_cell <- read.table("~/Data_LINCS/2020beta/cellinfo_beta.txt",
                        header=T,sep="\t",row.names=1,colClasses="character",
                        quote="\"",comment.char="")
info_gene <- read.table("~/Data_LINCS/2020beta/geneinfo_beta.txt",
                        header=T,sep="\t",row.names=1,colClasses="character",
                        quote="\"",comment.char="")
info_trt <- read.table("~/Data_LINCS/2020beta/compoundinfo_beta.txt",
                       header=T,sep="\t",row.names=1,colClasses="character",
                       quote="\"",comment.char="")
info_coldata <- read.table("~/Data_LINCS/2020beta/siginfo_beta.txt",
                           header=T,sep="\t",colClasses="character",
                           quote="\"",comment.char="")
