## ----style-knitr, eval=TRUE, echo=FALSE, results="asis"--------------------
BiocStyle::latex2()

## ----load-purecn, echo=FALSE, message=FALSE--------------------------------
library(PureCN)
set.seed(1234)

## ----examplecoverage-------------------------------------------------------
bam.file <- system.file("extdata", "ex1.bam", package="PureCN", 
    mustWork = TRUE)
interval.file <- system.file("extdata", "ex1_intervals.txt", 
    package="PureCN", mustWork = TRUE)

calculateBamCoverageByInterval(bam.file=bam.file, 
 interval.file=interval.file, output.file="ex1_coverage.txt")

## ----examplegc-------------------------------------------------------------
interval.file <- system.file("extdata", "ex2_intervals.txt", 
        package = "PureCN", mustWork = TRUE)
reference.file <- system.file("extdata", "ex2_reference.fa", 
    package = "PureCN", mustWork = TRUE)
calculateGCContentByInterval(interval.file, reference.file, 
    output.file = "ex2_gc_file.txt")

## ----examplegc2------------------------------------------------------------
bed.file <- system.file("extdata", "ex2_intervals.bed", 
        package = "PureCN", mustWork = TRUE)

intervals <- import(bed.file)

calculateGCContentByInterval(intervals, reference.file, 
    output.file = "ex2_gc_file.txt")

## ----example_files, message=FALSE, warning=FALSE, results='hide'-----------
library(PureCN)

normal.coverage.file <- system.file("extdata", "example_normal.txt",
    package="PureCN") 
normal2.coverage.file <- system.file("extdata", "example_normal2.txt", 
    package="PureCN") 
normal.coverage.files <- c(normal.coverage.file, normal2.coverage.file) 
tumor.coverage.file <- system.file("extdata", "example_tumor.txt", 
    package="PureCN") 
seg.file <- system.file("extdata", "example_seg.txt",
    package = "PureCN")
vcf.file <- system.file("extdata", "example_vcf.vcf", package="PureCN")
gc.gene.file <- system.file("extdata", "example_gc.gene.file.txt", 
    package="PureCN")

## ----figuregccorrect, fig.show='hide', fig.width=7, fig.height=4, warning=FALSE----
correctCoverageBias(normal.coverage.file, gc.gene.file, 
    output.file="example_normal_loess.txt", plot.gc.bias=TRUE)

## ----normaldb--------------------------------------------------------------
normalDB <- createNormalDatabase(normal.coverage.files)
# serialize, so that we need to do this only once for each assay
saveRDS(normalDB, file="normalDB.rds")

## ----normaldbpca-----------------------------------------------------------
normalDB <- readRDS("normalDB.rds")
# get the best normal
best.normal.coverage.file <- findBestNormal(tumor.coverage.file, 
    normalDB)

## ----normaldbpcapool-------------------------------------------------------
# get the best 2 normals and average them
pool <- findBestNormal(tumor.coverage.file, normalDB, 
    num.normals=2, pool=TRUE, remove.chrs=c("chrX", "chrY"))

## ----targetweightfile1, message=FALSE--------------------------------------
target.weight.file <- "target_weights.txt"
createTargetWeights(tumor.coverage.file, normal.coverage.files, 
    target.weight.file)

## ----ucsc_segmental--------------------------------------------------------
# Instead of using a pool of normals to find low quality regions,
# we use suitable BED files, for example from the UCSC genome browser.

# We do not download these in this vignette to avoid build failures 
# due to internet connectivity problems.
downloadFromUCSC <- FALSE
if (downloadFromUCSC) {
    library(rtracklayer)
    mySession <- browserSession("UCSC")
    genome(mySession) <- "hg19"
    tbl.segmentalDups <- getTable( ucscTableQuery(mySession, 
        track="Segmental Dups", table="genomicSuperDups"))
    tbl.simpleRepeats <- getTable( ucscTableQuery(mySession, 
        track="Simple Repeats", table="simpleRepeat"))
    export(tbl.segmentalDups, "hg19_segmentalDuplications.bed")
    export(tbl.simpleRepeats, "hg19_simpleRepeats.bed")
}

snp.blacklist <- c("hg19_segmentalDuplications.bed", 
    "hg19_simpleRepeats.bed")

## ----runpurecn-------------------------------------------------------------
ret <-runAbsoluteCN(normal.coverage.file=pool,
# normal.coverage.file=normal.coverage.file, 
    tumor.coverage.file=tumor.coverage.file, vcf.file=vcf.file, 
    genome="hg19", sampleid="Sample1", 
    gc.gene.file=gc.gene.file, normalDB=normalDB,
# args.setMappingBiasVcf=list(normal.panel.vcf.file=normal.panel.vcf.file),
# args.filterVcf=list(snp.blacklist=snp.blacklist, 
# stats.file=mutect.stats.file), 
    args.segmentation=list(target.weight.file=target.weight.file), 
    post.optimize=FALSE, plot.cnv=FALSE, verbose=FALSE)

## ----createoutput----------------------------------------------------------
file.rds <- "Sample1_PureCN.rds"
saveRDS(ret, file=file.rds)
pdf("Sample1_PureCN.pdf", width=10, height=11)
plotAbs(ret, type="all")
dev.off()

## ----figureexample1, fig.show='hide', fig.width=6, fig.height=6------------
plotAbs(ret, type="overview")

## ----figureexample2, fig.show='hide', fig.width=6, fig.height=6------------
plotAbs(ret, 1, type="hist")

## ----figureexample3, fig.show='hide', fig.width=8, fig.height=8------------
plotAbs(ret, 1, type="BAF")

## ----figureexample4, fig.show='hide', fig.width=8, fig.height=8------------
plotAbs(ret, 1, type="AF")

## ----output1---------------------------------------------------------------
names(ret)

## ----output3---------------------------------------------------------------
head(predictSomatic(ret), 3)

## ----output4---------------------------------------------------------------
vcf <- predictSomatic(ret, return.vcf=TRUE)
writeVcf(vcf, file="Sample1_PureCN.vcf") 

## ----calling2--------------------------------------------------------------
gene.calls <- callAlterations(ret)
head(gene.calls)

## ----loh-------------------------------------------------------------------
loh <- callLOH(ret)
head(loh)

## ----curationfile----------------------------------------------------------
createCurationFile(file.rds) 

## ----readcurationfile------------------------------------------------------
ret <- readCurationFile(file.rds)

## ----curationfileshow------------------------------------------------------
read.csv("Sample1_PureCN.csv")

## ----customseg-------------------------------------------------------------
retSegmented <- runAbsoluteCN(seg.file=seg.file, 
    gc.gene.file=gc.gene.file, vcf.file=vcf.file, 
    max.candidate.solutions=1, genome="hg19",
    test.purity=seq(0.3,0.7,by=0.05), verbose=FALSE, 
    plot.cnv=FALSE)

## ----figurecustombaf, fig.show='hide', fig.width=8, fig.height=8-----------
plotAbs(retSegmented, 1, type="BAF")

## ----customlogratio, message=FALSE-----------------------------------------
# We still use the log-ratio exactly as normalized by PureCN for this
# example
log.ratio <- calculateLogRatio(readCoverageGatk(normal.coverage.file),
    readCoverageGatk(tumor.coverage.file))

retLogRatio <- runAbsoluteCN(log.ratio=log.ratio,
    gc.gene.file=gc.gene.file, vcf.file=vcf.file, 
    max.candidate.solutions=1, genome="hg19",
    test.purity=seq(0.3,0.7,by=0.05), verbose=FALSE, 
    normalDB=normalDB, plot.cnv=FALSE)

## ----power1, fig.show='hide', fig.width=6, fig.height=6--------------------
purity <- c(0.1,0.15,0.2,0.25,0.4,0.6,1)
coverage <- seq(5,35,1)
power <- lapply(purity, function(p) sapply(coverage, function(cv) 
    calculatePowerDetectSomatic(coverage=cv, purity=p, ploidy=2, 
    verbose=FALSE)$power))

# Figure S7b in Carter et al.
plot(coverage, power[[1]], col=1, xlab="Sequence coverage", 
    ylab="Detection power", ylim=c(0,1), type="l")

for (i in 2:length(power)) lines(coverage, power[[i]], col=i)
abline(h=0.8, lty=2, col="grey")     
legend("bottomright", legend=paste("Purity", purity),
    fill=seq_along(purity))

## ----power2, fig.show='hide', fig.width=6, fig.height=6--------------------
coverage <- seq(5,350,1)
power <- lapply(purity, function(p) sapply(coverage, function(cv) 
    calculatePowerDetectSomatic(coverage=cv, purity=p, ploidy=2, 
    cell.fraction=0.2, verbose=FALSE)$power))
plot(coverage, power[[1]], col=1, xlab="Sequence coverage", 
    ylab="Detection power", ylim=c(0,1), type="l")

for (i in 2:length(power)) lines(coverage, power[[i]], col=i)
abline(h=0.8, lty=2, col="grey")     
legend("bottomright", legend=paste("Purity", purity),
    fill=seq_along(purity))

## ----annotatesymbols, message=FALSE, warning=FALSE-------------------------
biocLite("TxDb.Hsapiens.UCSC.hg19.knownGene")
biocLite("org.Hs.eg.db")

library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db) 

.annotateIntervals <- function(gc.gene.file, txdb, output.file = NULL) {
    gc <- read.delim(gc.gene.file, as.is=TRUE)
    # misuse this function to convert interval string into data.frame
    gc.data <- readCoverageGatk(gc.gene.file)
    grGC <- GRanges(seqnames=gc.data$chr, 
        IRanges(start=gc.data$probe_start, end=gc.data$probe_end))
    id <- transcriptsByOverlaps(txdb, ranges=grGC, columns = "GENEID")
    id$SYMBOL <-select(org.Hs.eg.db, as.character(id$GENEID), "SYMBOL")[,2]
    gc$Gene <- id$SYMBOL[findOverlaps(grGC, id, select="first")]
    if (!is.null(output.file)) {
        write.table(gc, file = output.file, row.names = FALSE, 
            quote = FALSE, sep = "\t")
    }
    invisible(gc)
}

.annotateIntervals(gc.gene.file, TxDb.Hsapiens.UCSC.hg19.knownGene)

## ----sessioninfo, results='asis', echo=FALSE-------------------------------
toLatex(sessionInfo())

