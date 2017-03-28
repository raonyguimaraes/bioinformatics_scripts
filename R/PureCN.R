
#RTFM
#http://bioconductor.org/packages/devel/bioc/vignettes/PureCN/inst/doc/PureCN.pdf

## try http:// if https:// URLs are not supported
source("https://bioconductor.org/biocLite.R")
biocLite("PureCN")
library("PureCN")

browseVignettes("PureCN")

bam.file <- system.file("extdata", "ex1.bam", package="PureCN", mustWork = TRUE)
interval.file <- system.file("extdata", "ex1_intervals.txt", package="PureCN", mustWork = TRUE)
calculateBamCoverageByInterval(bam.file=bam.file, interval.file=interval.file, output.file="ex1_coverage.txt")
