#Rscript to find overlap b/w two BED3 files
#Author: Thomas
args = commandArgs(trailingOnly=TRUE)
if (length(args)==0)
{
  stop("Usage: Rscript findoverlap_granges.R <bedfile1> <bedfile2>", call.=FALSE)
}
if (!file.exists(args[1]) || !file.exists(args[2]))
{
  stop("Input file not found!\
Usage: Rscript findoverlap_granges.R <bedfile1> <bedfile2>", call.=FALSE)
}
suppressPackageStartupMessages({
library(ChIPQC)
library(rtracklayer)
})

bed1=args[1]
bed2=args[2]

name1=gsub("_peaks.narrowPeak.blacklistcleared","",bed1)
name1=gsub("_merged.bed","",bed1)
name1=gsub(".bed","",bed1)
name2=gsub("_peaks.narrowPeak.blacklistcleared","",bed2)
name2=gsub("_merged.bed","",bed2)
name2=gsub(".bed","",bed2)
#name1=gsub(".R0_peaks.bed","",basename(bed1))
#name2=gsub(".R0_peaks.bed","",basename(bed2))
outputname=paste0(name1,"-vs-",name2)
outputname
#name1=gsub("R0_peaks.xls","",bed1)
#name2=gsub("R0_peaks.xls","",bed2)

#create Granges from Bed file
gr1 = read.table(bed1,sep="\t",header = F)
head=c("chr","start","end")
names(gr1)=head
firstPeakSet = GRanges(seqnames=gr1$chr,
            ranges=IRanges(start=gr1$start+1,end=gr1$end),
            strand=rep("*",nrow(gr1)))

gr2 = read.table(bed2,sep="\t",header = F)
names(gr2)=head
secondPeakSet = GRanges(seqnames=gr2$chr,
                       ranges=IRanges(start=gr2$start+1,end=gr2$end),
                       strand=rep("*",nrow(gr2)))

#firstPeakSet = toGRanges(bed1,format="BED")
#secondPeakSet = toGRanges(bed2,format="BED")
#firstPeakSet <- ChIPQC:::GetGRanges(bed1, sep="\t", simple=F)
#secondPeakSet <- ChIPQC:::GetGRanges(bed2, sep="\t", simple=F)

OnlyfirstPeakSet <- firstPeakSet[!(firstPeakSet %over% secondPeakSet)]
firstANDsecondPeakSets <- firstPeakSet[firstPeakSet %over% secondPeakSet]
onlysecondPeakset <- secondPeakSet[!(secondPeakSet %over% firstPeakSet)]

len1=length(OnlyfirstPeakSet)
#comm1=length(firstANDsecondPeakSets)
len2=length(onlysecondPeakset)

allPeaks <- c(firstPeakSet,secondPeakSet)
allPeaksReduced <- reduce(allPeaks)
commonPeaks <- allPeaksReduced[allPeaksReduced %over% firstPeakSet 
                               & allPeaksReduced %over% secondPeakSet]
comm1=length(commonPeaks)

write.table(cbind(outputname,len1,comm1,len2),file="overlap_peakstats.tsv",append=T,sep="\t",quote = F ,col.names = F,row.names = F)

if(len1 > 0)
export.bed(OnlyfirstPeakSet, paste0(outputname,"_",name1,"_uniq.bed"))
if(len2 > 0)
export.bed(onlysecondPeakset, paste0(outputname,"_",name2,"_uniq.bed"))
if(comm1 > 0)
export.bed(commonPeaks,paste0(outputname,"_common.bed"))
