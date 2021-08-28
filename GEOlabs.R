temp <- read.table("GEOlabs.txt",sep="\t",row.names=1,quote="")
GEO <- temp[,1]
names(GEO) <- rownames(temp)
GEO <- sapply(GEO,strsplit,split=",")
GEO <- sapply(GEO,function(X) gsub("^ | $","",X))
GEO <- sapply(GEO,function(X) gsub("[0-9]","",X))
GEO <- sapply(GEO,function(X) gsub("^ | $","",X))

dups <- unlist(GEO)[which(duplicated(unlist(GEO)))]

sapply(dups,function(X) which(sapply(GEO,function(Y) X %in% Y)))
