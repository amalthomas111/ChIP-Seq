#Author: Thomas
#Rscript for Differential analysis using csaw (data with Replicates)
#Input files: csv file, outputname
#Change blacklist file, frag.length if reqd.
args = commandArgs(trailingOnly=TRUE)
if (length(args)!=2)
{
  stop("Usage: Rscript differential_csaw.R <inputcsvfile> <name>", call.=FALSE)
}
if (!file.exists(args[1])
{
  stop("Input file not found!\
Usage: Rscript differential_csaw.R <inputcsvfile> <name>", call.=FALSE)
}

suppressPackageStartupMessages({
library(rtracklayer)
library(csaw)
library(edgeR)
library(session)
})

#read input csv file
inputfile = args[1]
outputname = args[2]

####TEST
#inputfile = "data_Brx42Brainmets-vs-42other.csv"
#outputname = "Brx42Brainmets-vs-42other"
#######

######### CHANGE BLACKLIST FILE ##############
#read blacklist regions for hg19
blacklistdata <- read.table("/media/master/hdd/WORK/blacklist/hg19_consensusBlacklist.bed",header=F) 
head(blacklistdata)
colnames(blacklistdata) <- c('chr','start','end','id','score','strand')
head(blacklistdata)
blacklist_granges <- with(blacklistdata, GRanges(chr, IRanges(start+1, end),strand, score, id=id))
head(blacklist_granges)
param <- readParam(discard=blacklist_granges)

data = read.csv(inputfile,header=T)

bam.files = as.vector(data$bamReads)
condition = as.vector(data$Condition)
#exp. data
data.frame(BAM=basename(bam.files), Condition = condition)

#frag.len=
bins <- windowCounts(bam.files, bin=TRUE, width=2000L, param=param)
win.data <- windowCounts(bam.files, param=param, width=150L, ext = NA)#ext=frag.len)

win.data
win.data$totals 

#remove low reads
abundances <- aveLogCPM(asDGEList(win.data))
summary(abundances)
keep <- abundances > aveLogCPM(5, lib.size=mean(win.data$totals))
print("keep1 summary")
summary(keep)

#filter by fold change above global avg
filter.stat <- filterWindows(win.data, background=bins, type="global")
min.fc <- 3
keep1 <- filter.stat$filter > log2(min.fc)
print("keep1 summary")
summary(keep1)

save.session(paste0(outputname,'_keep1keep.Rdata'))

#data filtered
filtered.data <- win.data[keep & keep1,]
print("#######Filtered data over ###############")
nrow(filtered.data)

offsets <- normOffsets(filtered.data, type="loess")
head(offsets) 
#adjc <- log2(assay(filtered.data)+0.5)
#norm.adjc <- adjc - offsets/log(2)

condition <- factor(condition)
design <- model.matrix(~0+condition)
colnames(design) <- levels(condition)
design 

print("######## Design over #########")

y <- asDGEList(filtered.data)
y$offset <- offsets 

y <- estimateDisp(y, design)
summary(y$trended.dispersion) 

#plotBCV(y)

fit <- glmQLFit(y, design, robust=TRUE)
#plotQLDisp(fit)

summary(fit$df.prior)
pdf(file = paste0(outputname,"_mds.pdf"))
plotMDS(norm.adjc, labels=gsub(".sorted.chr.nodup.filt.bam","",basename(bam.files)),
        col=c("red", "blue")[as.integer(condition)]) 

dev.off()

contrast <- makeContrasts(group1-group2, levels=design)
res <- glmQLFTest(fit, contrast=contrast)
head(res$table)

merged <- mergeWindows(rowRanges(filtered.data), tol=100, max.width=5000)

tabcom <- combineTests(merged$id, res$table)
head(tabcom) 

#FDR cut-off
is.sig <- tabcom$FDR <= 0.05
summary(is.sig) 

has.up <- tabcom$logFC.up > 0
has.down <- tabcom$logFC.down > 0
data.frame(Total.DB=sum(is.sig), DB.up=sum(is.sig & has.up & !has.down),
           DB.down=sum(is.sig & !has.up & has.down), DB.both=sum(is.sig & has.up & has.down)) 


createtable = function(obj){
  obj = obj[order(-elementMetadata(obj)$score)]
  table = data.frame(chr=seqnames(obj),
                     start=start(obj),
                     end=end(obj),
                     score=elementMetadata(obj)$score)
  return(table)
}

out.ranges <- merged$region

simplified_up <- out.ranges[is.sig & has.up & !has.down]
simplified_up$score <- -10*log10(tabcom$PValue[is.sig & has.up & !has.down])

simplified_down <- out.ranges[is.sig & !has.up & has.down]
simplified_down$score <- -10*log10(tabcom$PValue[is.sig & !has.up & has.down])

simplified_both <- out.ranges[is.sig & has.up & has.down]
simplified_both$score <- -10*log10(tabcom$PValue[is.sig & has.up & has.down])

if(length(simplified_up)>0){
  up = createtable(simplified_up)
  write.table(up, file=paste0(outputname,"_up_fdr0.05.bed"), quote=F, sep="\t", row.names=F, col.names=T)
}

if(length(simplified_down)>0){
  down = createtable(simplified_down)
  write.table(down,file=paste0(outputname,"_down_fdr0.05.bed"), quote=F, sep="\t", row.names=F, col.names=T)
}
if(length(simplified_both)>0){
  both = createtable(simplified_both)
  write.table(both,file=paste0(outputname,"_both_fdr0.05.bed"), quote=F, sep="\t", row.names=F, col.names=T)
}

