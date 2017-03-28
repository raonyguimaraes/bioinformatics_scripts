## try http:// if https:// URLs are not supported
source("https://bioconductor.org/biocLite.R")
biocLite("trackViewer")

library(Gviz)
library(rtracklayer)
library(trackViewer)

SNP <- c(10, 100, 105, 108, 400, 410, 420, 600, 700, 805, 840, 1400, 1402)
x <- sample.int(100, length(SNP))
sample.gr <- GRanges("chr1", IRanges(SNP, width=1, names=paste0("snp", SNP)))
features <- GRanges("chr1", IRanges(c(1, 501, 1001), 
                                    width=c(120, 400, 405),
                                    names=paste0("block", 1:3)))
lolliplot(sample.gr, features)