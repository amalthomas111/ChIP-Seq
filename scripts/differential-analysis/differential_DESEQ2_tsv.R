#Author: Thomas
#Script to perfrom differential analysis for atac-seq/chip-seq 
#using DESEQ2 (with replicates)
#arguments: bed regions and count-data (tsv files)
#Change treatment factor according to the data
args = commandArgs(trailingOnly=TRUE)
if (length(args)!=3)
{
  stop("Usage: Rscript differential_DSEQ2.R <bedfile> <countable(tsv file)> <outputname>", call.=FALSE)
}
if (!file.exists(args[1]) || !file.exists(args[2]))
{
  stop("Input file not found!\
Usage: Rscript differential_DSEQ2.R <bedfile> <countable(tsv file)> <outputname>", call.=FALSE)
}

suppressPackageStartupMessages({
  library("DESeq2")
  library(ggplot2)
})

#read arguments
bedfile=args[1]
count_table_file=args[2]
outputname=args[3]

####TEST RUN##
# bedfile="Brx68_totalpeaks.merged.bed"
# count_table_file="countables/Brx68_count-table_H-N8.tsv"
# outputname="test"
#############

#function to convert BED3 to Granges
create_granges_frombed3 = function(bed){
  gr1 = read.table(bed,sep="\t",header = F)
  head=c("chr","start","end")
  names(gr1)=head
  gr = GRanges(seqnames=gr1$chr,
                 ranges=IRanges(start=gr1$start+1,end=gr1$end),
                 strand=rep("*",nrow(gr1)))
  return(gr)
}
#granges created for regions
totalwindows=create_granges_frombed3(bedfile)

#read input count-table
df= read.table(count_table_file,sep="\t",row.names=1,header = T)
print("head counttable")
head(df)
#create design table.
###### CHANGE THIS ACCORDINGLY #######
treatment = c("group2","group2","group1","group1")
#,"group2","group2","group2","group2","group2","group2")
########################################
design.matrix = data.frame(row.names = colnames(df),
                            treatment = treatment)
print("DESEQ2 Design Matrix")
design.matrix

dds = DESeqDataSetFromMatrix(countData = df,
                              colData = design.matrix,
                              design = ~ treatment,
                              rowRanges=totalwindows)


dds = DESeq(dds)
dds_results = results(dds, contrast=c("treatment","group1","group2"),
                         format="GRanges")

#differential regions
up_fdr = dds_results[dds_results$padj < 0.05 & !is.na(dds_results$padj) & dds_results$log2FoldChange > 0]

down_fdr = dds_results[dds_results$padj < 0.05 & !is.na(dds_results$padj) & dds_results$log2FoldChange < 0]
cat("#####\n No. of Differential regions identified:\nfdr_up\n")
length(up_fdr)
print("fdr_down")
length(down_fdr)

#function to create data frame for output generation
create_output =function(df,flag){
  if(flag==1) df=df[order(df$padj)]
  if(flag==2) df=df[order(df$pvalue)]
  gr1 = data.frame(chr=seqnames(df),
                       start=start(df)-1,
                       end=end(df),
                       log2Foldchange=elementMetadata(df)$log2FoldChange,
                       pvalue=elementMetadata(df)$pvalue,
                       padj=elementMetadata(df)$padj)  
  
}

# FDR<0.05 output files
if(length(up_fdr)){
  up.regions=create_output(up_fdr,1)
  write.table(up.regions, file=paste0(outputname,"_DESEQ2_up_fdr0.05.bed"), quote=F, sep="\t", row.names=F, col.names=T)
}
if(length(down_fdr))
{
  down.regions=create_output(down_fdr,1)
  write.table(down.regions,file=paste0(outputname,"_DESEQ2_down_fdr0.05.bed"), quote=F, sep="\t", row.names=F, col.names=T)
}

#create other outputfiles

up = dds_results[dds_results$log2FoldChange > 0]
top_up = create_output(up,2)
down = dds_results[dds_results$log2FoldChange < 0]
top_down = create_output(down,2)

#create new directory otheroutput and write the following
#files to that directory
mainDir=getwd()
subDir="otheroutput"
dir.create(file.path(mainDir,subDir), showWarnings = FALSE)

write.table(top_up, file=paste0("otheroutput/",outputname,"_DESEQ2_up_sortpval.bed"), quote=F, sep="\t", row.names=F, col.names=T)
write.table(top_up, file=paste0("otheroutput/",outputname,"_DESEQ2_up_sortpval.bed"), quote=F, sep="\t", row.names=F, col.names=F)
write.table(head(top_up,n=100L), file=paste0("otheroutput/",outputname,"_DESEQ2_up_sortpval_top100.bed"), 
            quote=F, sep="\t", row.names=F, col.names=F)

write.table(top_down, file=paste0("otheroutput/",outputname,"_DESEQ2_down_sortpval.bed"), quote=F, sep="\t", row.names=F, col.names=T)
write.table(top_down, file=paste0("otheroutput/",outputname,"_DESEQ2_down_sortpval.bed"), quote=F, sep="\t", row.names=F, col.names=F)
write.table(head(top_down,n=100L), file=paste0("otheroutput/",outputname,"_DESEQ2_down_sortpval_top100.bed"),
            quote=F, sep="\t", row.names=F, col.names=F)

#create new directory plots
mainDir=getwd()
subDir="plots"
dir.create(file.path(mainDir,subDir), showWarnings = FALSE)
vsd = varianceStabilizingTransformation(dds, blind = FALSE)
write.table(assay(vsd), file = paste0("plots/",outputname,"_count.table.vsd.txt"), sep = '\t', quote = FALSE)
pdf(file=paste0("plots/",outputname,"_pca.pdf"))
plotPCA(vsd, intgroup="treatment") + theme_bw()
dev.off()
