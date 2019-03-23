args = commandArgs(trailingOnly=TRUE)
if (length(args)!=2)
{
  stop("Usage: Rscript chipQC.R <genome> <samplename(without ext)", call.=FALSE)
}
suppressPackageStartupMessages({
library(GenomicAlignments)
library(ChIPQC)
library(GenomicRanges)
})

inputbam=paste0("bams/",args[2],".sorted.bam")

if (!file.exists(inputbam))
{
  stop("Input file not found!\
Usage: Rscript chipQC.R <genome>  <samplename(without ext)", call.=FALSE)
}
mainDir=getwd()
subDir="fragmentlength"
dir.create(file.path(mainDir,subDir), showWarnings = FALSE)

if(args[1]=="hg38"){
cat("\nGenome:hg38\n")
load("/home/rcf-40/amalthom/panases_soft/genomes/blacklist/hg38_blacklist.Robj")
blacklist = hg38_blacklist
chroms = paste0("chr",c(1:22,'X','Y'))
}else if(args[1]=="hg19"){
cat("\nGenome:hg19\n")
data(blacklist_hg19)
blacklist=blacklist.hg19
chroms = paste0("chr",c(1:22,'X','Y'))
}else if(args[1]=="mm10"){
cat("\nGenome:mm10\n")
load("/home/rcf-40/amalthom/panases_soft/genomes/blacklist/mm10_blacklist.Robj")
chroms = paste0("chr",c(1:19,'X','Y'))
blacklist=mm10_blacklist
}else if(args[1]=="dm6"){
blacklist=NULL
chroms = c('chr3R','chr3L','chr2R','chrX','chr2L','chrY','chr4','chrM')
}else{
stop("Unknown genome!Exiting")
}


exampleExp = ChIPQCsample(inputbam,peaks=NULL,blacklist = blacklist, chromosomes = chroms,shift=1:400)
QCmetrics(exampleExp)
name=paste0("fragmentlength/fraglen_",args[2])
file.create(name)
fl=as.data.frame(QCmetrics(exampleExp))["FragL",]
write(fl,file=name)
metrics = t(as.data.frame(QCmetrics(exampleExp)))
metrics
rownames(metrics) = args[2]
write.table(metrics,file="fragmentlength/BamQC_details.txt", sep = "\t",col.names=F, append= T,quote = F)
