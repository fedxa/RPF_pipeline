library(GenomicAlignments)
library(GenomicFeatures)
library(rtracklayer)
library(RiboProfiling)

plotLengthDistribution <- function(aln, lens=NULL, base='data', save=TRUE, plot=TRUE) {
    matchLenDistr <- histMatchLength(aln, 0)
    ## Should in principle plot the matchLenDistr ggplot object.... But we are lazy:
    reads <- matchLenDistr[[1]]$counts
    readslen <- matchLenDistr[[1]]$matchSize
    col <- rep('gray', length(reads))
    if (is.null(lens)) {
        mlen <- getPlotLengthMaximum(matchLenDistr[[1]])
        lens <- c(mlen-1, mlen, mlen+1)
    }
    col[readslen %in% lens] <- 'darkgreen'
    if (plot)
        barplot(reads, names.arg=readslen, col=col, xlab="Read length", ylab="Mapped reads")
    if (save) {
        pdf(paste0(base,"_Length.pdf"), width=174/25.4, height=85/25.4)
        barplot(reads, names.arg=readslen, col=col, xlab="Read length", ylab="Mapped reads")
        dev.off()
    }
    matchLenDistr[[1]]
}

getPlotLengthMaximum <- function(pp) {
    as.integer(as.character(pp$matchSize[which.max(pp$counts)]))
}


## Dirty fix for old version of RiboProfiling package
if (R.version$major=="3" && R.version$minor=="2.2")
    readsToStartOrEnd <- function(aln, what) {
        stopifnot(what=="start")
        readsToReadStart(aln)
    }

plotASiteShifts <- function(aln, lens, txdb) {
    if (is(aln, "GAlignments"))
        alnGRanges <- readsToStartOrEnd(aln, what="start")
    else if (is(aln, "GRanges"))
        alnGRanges <- aln
    else
        stop(paste("Do not know what to do with typeof(aln)==", typeof(aln)))
    oneBinRanges <- aroundPromoter(txdb, alnGRanges, percBestExpressed=0.03)
    ##the coverage in the TSS flanking region for the reads with match sizes 29:31
    listPromoterCov <- readStartCov(
        alnGRanges,
        oneBinRanges,
        matchSize=lens,
        fixedInterval=c(-20, 20),
        renameChr="aroundTSS",
        charPerc="perc")
    plot <- plotSummarizedCov(listPromoterCov)
    shifts <- lapply(listPromoterCov, function(iSumCov) {
        maxPeak <- max(iSumCov$values)
        maxPeakPos <- start(iSumCov)[which(iSumCov$values == 
                                           maxPeak)]
        maxPeakPos
        })
    list(-as.integer(shifts[c(-1)]),
         plot)
}

genShiftedRanges <- function(aln, lens, shifts) {
    ## Shifting properly  -- get information from the previous plot!
    asites <- GRanges()
    for (i in seq_along(lens)) {
        l <- lens[[i]]
        alnsel <- aln[qwidth(aln)==l] # qwidth to get width of the read, not the target!
        pl <- granges(qnarrow(alnsel[strand(alnsel)=='+'], start=shifts[[i]]+1, width=1))
        pl$score <- l
        asites <- c(asites, pl)
        ## Negative strand
        pl <- granges(qnarrow(alnsel[strand(alnsel)=='-'], end=l-shifts[[i]], width=1))
        pl$score <- l
        asites <- c(asites, pl)
    }
    asites <- sort(asites)
    asites
}


## Reads WIG pair representing read counts at a given site
##  Note: the resulting GRangeList may have erroneous score field
## readWigsWithCounts <- function(plus,minus) {
##     wp <- import.wig(plus)
##     strand(wp) <- '+'
##     wm <- import.wig(minus)
##     strand(wm) <- '-'
##     c(rep(wp,score(wp)),rep(wm,score(wm)))
## }

genWigs <- function(asites, wigbase) {
    pluscov <- coverage(asites[strand(asites)=='+'])
    export.bedGraph(pluscov, paste0(wigbase,".plus.bedgraph"))
    system(paste("bgzip -f", paste0(wigbase,".plus.bedgraph")))
    ## system(paste("~/Bio/userApps/bin/bedGraphToBigWig",
    ##              paste0(wigbase,".plus.bedgraph"),
    ##              "~/net/frt2/data/hg19/hg19.fa.fai",
    ##              paste0(wigbase,".plus.bw")))

    minuscov <- coverage(asites[strand(asites)=='-'])
    export.bedGraph(minuscov, paste0(wigbase,".minus.bedgraph"))
    system(paste("bgzip -f", paste0(wigbase,".minus.bedgraph")))
    ## system(paste("~/Bio/userApps/bin/bedGraphToBigWig",
    ##              paste0(wigbase,".minus.bedgraph"),
    ##              "~/net/frt2/data/hg19/hg19.fa.fai",
    ##              paste0(wigbase,".minus.bw")))
}



stripEncodeVersion <- function(namelist, stripencode=TRUE) {
    if (stripencode)
        sub('\\.[0-9]*$', '', namelist)
    else
        namelist
}


genResultTable <- function(se=secds, stripencode=TRUE) {
    m <- data.frame(name=stripEncodeVersion(rownames(se),stripencode),
                    fullname=rownames(se),
                    cdslen=sum(width(rowRanges(se))), assay(se))
    m$rpkm <- m$reads / (sum(m$reads)/1e6) / (m$cdslen/1000)
    return(m)
}



doAllSample <- function(bam, base, txdb, lens=NULL, shifts=NULL, plot=TRUE) {
    shitScanBamParam <- ScanBamParam(what=c("flag","strand","pos","qwidth","cigar"),
                                     flag=scanBamFlag(isUnmappedQuery=F),
                                     tagFilter=list(NH=1))
    message("Reading ", bam)
    aln <- readGAlignments(BamFile(bam), param=shitScanBamParam)
    message("Counting lengths")
    pp <- plotLengthDistribution(aln, lens, base, save=TRUE, plot=plot)
    if (is.null(lens)) {
        mlen <- getPlotLengthMaximum(pp)
        lens <- c(mlen-1, mlen, mlen+1)
    }
    message("Searching A-sites")
    # pdf(paste0(base,"_tmp.pdf"))
    rr <- plotASiteShifts(aln, lens, txdb)
    if (plot)
        print(rr[2])
    ggbio::ggsave(paste0(base,"_Asiteshifts.pdf"), rr[[2]])
    if (is.null(shifts))
        shifts <- rr[[1]]
    write.csv(data.frame(bam=bam,
                         lens=paste(lens,collapse=" "),
                         shifts=paste(shifts,collapse=" ")),
              paste0(base,".shifts.csv"), row.names=FALSE)
    message("Shifting")
    asites <- genShiftedRanges(aln, lens, shifts)
    message("Checking shifts")
    rr2 <- plotASiteShifts(asites, lens, txdb)
    if (plot)
        print(rr2[2])
    ggbio::ggsave(paste0(base,"_shiftedCheck.pdf"), rr2[[2]])
    message("Generating wigs")
    genWigs(asites, base)
    cds <- cdsBy(txdb, by="gene")
    message("Counting RPKMS")
    secds <- summarizeOverlaps(features=cds, reads=asites, mode="Union")
    genResultTable(secds)
}


## Run the script from snakemake
if (exists("snakemake")) {
    ## save.image("debug.RData")
    options(device=pdf)
    message <- function(...) { cat(...) } ## Fix partially the crazy rpy2 warnings
    txdb <- makeTxDbFromGFF(snakemake@config[["GENOME_GFF"]])

    bam <- snakemake@input$bam

    lens <- NULL
    shifts <- NULL
    try({
        df <- read.csv(snakemake@config[["MANUAL_SHIFTS_FILE"]], as.is=TRUE);
        lens <- as.integer(strsplit(df[df$bam==bam, "lens"], " ")[[1]]);
        shifts <- as.integer(strsplit(df[df$bam==bam, "shifts"], " ")[[1]]);
    })
    ## lens <- snakemake@config[["RETAINED_READ_LENGTHS"]]
    ## shifts <- snakemake@config[["RETAINED_READ_ASITE_OFFSETS"]]
    
    base <- file.path("wigs", snakemake@wildcards$sample)
    res <- doAllSample(snakemake@input$bam, base, txdb, lens=lens, shifts=shifts, plot=FALSE)
    write.csv(res, paste0(base,".csv"), row.names=FALSE)
}