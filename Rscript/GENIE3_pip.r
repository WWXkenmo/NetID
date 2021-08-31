#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)
# loading Parameters
if(length(args)!=2) {
    stop("Rscript GENIE3_pip.R ExpressionProfile, out_dir")
} else if (length(args == 2)) {
    cat(paste("Input: ",args[1]," \n", "outdir: ", args[2], " \n",sep=""))
}
input <- args[1];
output <- args[2];

###GENIE3 processing
suppressPackageStartupMessages(library(GENIE3))
suppressPackageStartupMessages(library(Matrix))

ks.final <- read.csv(input,header=TRUE,row.names=1)
ks.final <- as.matrix(ks.final)
system.time({
g3 <- GENIE3(ks.final,nCores=12,verbose=TRUE,nTrees=500,regulators = rownames(ks.final), targets = rownames(ks.final) )
}
)

g3_pair <- as(as.matrix(g3), "dgCMatrix")
g3_pair <- summary(g3_pair)
g3_pair <- as.data.frame(g3_pair)
g3_pair$i<- rownames(g3)[g3_pair$i]
g3_pair$j <- colnames(g3)[g3_pair$j]
colnames(g3_pair) <- c("Gene1","Gene2","EdgeWeight")
g3_pair <- g3_pair[order(g3_pair[,3],decreasing=TRUE),]
rownames(g3_pair) <- NULL

dir.create(paste(output,"GENIE3",sep=""))
setwd(paste(output,"GENIE3",sep=""))
write.table(g3_pair,"rankedEdges.csv",sep="\t",quote=FALSE,row.names=FALSE)