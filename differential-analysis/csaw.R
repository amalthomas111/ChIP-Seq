#R script for differential analysis without replicates using csaw
#Author: Jenny
#!/usr/bin/R
library(csaw)
library(edgeR)

bam.files <- c("/home/cmb-00/as/amalthom/DSB_laeiti/bams/bowtie1bams/phH2Av_Treat_R1.sorted.nodup.bam", "/home/cmb-00/as/amalthom/DSB_laeiti/bams/bowtie1bams/phH2Av_UNT_R1.sorted.nodup.bam")
bamnames <- c("phH2Av_Treat_R1", "phH2Av_UNT_R1")

win.width <- 500
wide.win.width <- 10000
param <- readParam(minq=0, restrict=c("chr2L", "chr2R", "chr3L", "chr3R", "chr4", "chrX")) # chrY not included

knownlength <- FALSE

if (!knownlength) {
  ## estimate fragment length
  collected.frag.len <- list()
  for (curbam in bam.files) {
    filename <- rev(unlist(strsplit(curbam, "/")))[1]
    x <- correlateReads(curbam, max.dist = 500, param = param)
    frag.len <- maximizeCcf(x)
    collected.frag.len[[curbam]] <- frag.len
    ## correlation plot
    pdf(paste(filename, "_CCF.pdf", sep=""))
    plot(0:500, x, type="l", ylab="CCF", xlab="Delay (bp)")
    abline(v=frag.len)
    text(x= frag.len, y=x[frag.len], labels = maximizeCcf(x), pos=4)
    dev.off()
    write(x= frag.len, file= paste(filename, "_fraglen.txt", sep="") )
  }
} eles {
  ## change to proper fragment lengths
  collected.frag.len <- list(200, 200)
}

#################################
#  Converting reads to counts   #
#################################
## sliding window
# data <- windowCounts(bam.files, spacing=100, ext=unlist(collected.frag.len), width=win.width, param=param)
## bins
data <- windowCounts(bam.files, bin=T, ext=unlist(collected.frag.len), width=win.width, param=param)

gr <- rowRanges(data)
df <- cbind(data.frame(seqnames=seqnames(gr),
  starts=start(gr)-1, ends=end(gr),
  names=c(rep(".", length(gr))),
  scores=c(rep(".", length(gr))),
  strands=strand(gr)),
  phH2Av_Treat = assay(data)[,1],
  phH2Av_UNT = assay(data)[,2]  )
write.table(df, paste("count_win", win.width,".txt",sep=""), quote=F)

#######################################
# Filtering out uninteresting windows #
#######################################
# filtere out windows with low counts
abundances <- aveLogCPM(asDGEList(data))
summary(abundances)
keep.simple <- abundances > -1  # arbitrary cutoff
filtered.data <- data[keep.simple,]
summary(keep.simple)

#####################################
# Calculating normalization factors #
#####################################
binned <- windowCounts(bam.files, bin=TRUE, width=wide.win.width, param=param)
normfacs <- normOffsets(binned)

# Visualizing normalization efforts with MA plots
adj.counts <- cpm(asDGEList(binned), log=TRUE)
pdf("MA.pdf")
for (i in 1:(length(bam.files)-1)) {
  cur.x <- adj.counts[,1]
  cur.y <- adj.counts[,1+i]
  smoothScatter(x=(cur.x+cur.y)/2+6*log2(10), y=cur.x-cur.y,
  xlab="A", ylab="M", main=paste(bamnames[1], "vs", bamnames[2]))
  all.dist <- diff(log2(normfacs[c(i+1, 1)]))
  abline(h=all.dist, col="red")
}
dev.off()


####################################
# Testing for differential binding #
####################################
# design matrix
grouping <- factor(c('Treat', 'UNT'))
design <- model.matrix(~0 + grouping)
colnames(design) <- levels(grouping)
y <- asDGEList(filtered.data, norm.factors=normfacs)

## Estimating the dispersions (needs multiple replicates)
# y <- estimateDisp(y, design)
# fit <- glmQLFit(y, design, robust=TRUE)

## Without replicates (option1)
# EdgeR manual pp.23: pick a reasonable dispersion value
#           0.1 for data on genetically identical model organisms
# bcv <- 0.1
# y <- asDGEList(filtered.data, norm.factors=normfacs, group=1:2)
# results <- exactTest(y, dispersion=bcv^2)

## Without replicates (option2)
## https://www.bioconductor.org/help/course-materials/2015/BioC2015/csaw_lab.html
## Interesting questions
## What do we do when we don't have any replicates (caution required)?
contrast <- makeContrasts(Treat - UNT, levels=design)
norep.fit <- glmFit(y, design, dispersion=0.05)
norep.results <- glmLRT(norep.fit, contrast=contrast)

## write out the raw test results
gr <- rowRanges(filtered.data)
df <- cbind(data.frame(seqnames=seqnames(gr),
  starts=start(gr)-1, ends=end(gr),
  names=c(rep(".", length(gr))),
  scores=c(rep(".", length(gr))),
  strands=strand(gr)),
  norep.results$table )
write.table(df, "glmLRT.result.txt", quote=F, row.names=F, sep="\t")


indup <- df$PValue < 0.01 & df$logFC>0
write.table(df[indup, ], "phH2Av_Treat_vs_UNT.500bp_bins_up_p0.01.txt", quote=F, row.names=F, sep="\t")
inddown <- df$PValue < 0.01 & df$logFC<0
write.table(df[inddown, ], "phH2Av_Treat_vs_UNT.500bp_bins_down_p0.01.txt", quote=F, row.names=F, sep="\t")


