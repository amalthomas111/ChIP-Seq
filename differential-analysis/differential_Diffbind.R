#Author: Thomas
#Script to perform differential analysis using Diffbind (2 replicates min)
#Input: csv file containing all metadata
#Up regions: first group in the contrast (it varies with input files)
#so check which is up and which is down
args = commandArgs(trailingOnly=TRUE)
if (length(args)!=2)
{
  stop("Usage: Rscript differential_diffbind.R <csv file> <outputname> ", call.=FALSE)
}
if (!file.exists(args[1]))
{
  stop("Input csv file not found!\
Usage: Rscript differential_diffbind.R <csv file> <outputname>", call.=FALSE)
}
suppressPackageStartupMessages({
require(DiffBind)
})
#packageVersion("DiffBind")
#samples = read.csv("data_info.csv")
#setwd("/media/master/hdd/WORK/egor/atac/oihana_data/0.newpreprocessing_peaks/diffbind")
inputfile = args[1]
name = args[2]

##### TEST #####
#inputfile= "data_68_H-N8.csv"
#outputname="BrxH-N8"
################

atac = dba(sampleSheet = inputfile)
atac = dba.count(atac, score=DBA_SCORE_READS,fragmentSize = 0)
atac
dba.overlap(atac,mode = DBA_OLAP_RATE)

#count reads in peak
counts <- dba.peakset(atac, bRetrieve=TRUE, DataType=DBA_DATA_FRAME)
df = counts[,c(4,5,6,7)]
rownames(df)=c(paste0(counts$CHR,":",counts$START,"-",counts$END))
head(df)
#write the count-data to file
write.csv(df,file=paste0("countdata_",name,".csv"),quote = F, row.names = T)

#heatmap based on counts
pdf(file = paste0(name,"_headtmap_counts.pdf"))
plot(atac)
dev.off()

#differential analysis
atac = dba.contrast(atac, categories = DBA_CONDITION, minMembers=2)
atac
atac = dba.analyze(atac,method = DBA_DESEQ2)
atac
atac.db = dba.report(atac,bCounts = T)
print("Up FDR 0.05")
sum(atac.db$Fold > 0)
print("Down FDR 0.05")
sum(atac.db$Fold < 0)

if(length(atac.db)!=0) {
  #heat map based on differential regions
  pdf(file = paste0(name,"_heatmap_counts_differential.pdf"))
  plot(atac,contrast=1)
  dev.off()
  
  #create data frame for output writing
  createtable = function(obj){
    obj = obj[order(elementMetadata(obj)$FDR)]
    table = data.frame(chr=seqnames(obj),
                       start=start(obj),
                       end=end(obj),
                       log2Foldchange=elementMetadata(obj)$Fold,
                       fdr=elementMetadata(obj)$FDR)
    return(table)
  }
  
  atac.db_up = atac.db[elementMetadata(atac.db)$Fold>0]
  atac.db_down = atac.db[elementMetadata(atac.db)$Fold<0]
  
  #FDR <=0.05 output
  if(length(atac.db_up)>0){
    up = createtable(atac.db_up)
    write.table(up, file=paste0(name,"_up_fdr0.05.bed"), quote=F, sep="\t", row.names=F, col.names=T)
  }
  #FDR <=0.05 output
  if(length(atac.db_down)>0){
    down = createtable(atac.db_down)
    write.table(down,file=paste0(name,"_down_fdr0.05.bed"), quote=F, sep="\t", row.names=F, col.names=T)
  }
  
  a = atac.db[elementMetadata(atac.db)$Fold > 0 & elementMetadata(atac.db)$FDR <= 0.01]
  print("FDR 0.01 up")
  print(length(a))
  b = atac.db[elementMetadata(atac.db)$Fold < 0 & elementMetadata(atac.db)$FDR <=0.01]
  print("FDR 0.01 down")
  print(length(b))
  
  #FDR <=0.01 output
  if(length(a)>0)
  {
    df.top <- createtable(a)
    write.table(df.top, file=paste0(name,"_up_fdr0.01.bed"), quote=F, sep="\t", row.names=F, col.names=T)
  }
  
  #FDR <=0.01 output
  if(length(b)>0)
  {
    df.down = createtable(b)
    write.table(df.down, file=paste0(name,"_down_fdr0.01.bed"), quote=F, sep="\t", row.names=F, col.names=T)
  }

  #Create plots based on differential regions
  print("pca")
  pdf(file = paste0(name,"_pca.pdf"))
  dba.plotPCA(atac,DBA_CONDITION,label = DBA_CONDITION)
  dev.off()
  print("pca after differential")
  pdf(file = paste0(name,"_pca_differential.pdf"))
  dba.plotPCA(atac,contrast = 1,label = DBA_CONDITION)
  dev.off()
  print("heatmap")
  pdf( file =paste0(name,"_heatmap_bindingaff.pdf"))
  dba.plotHeatmap(atac,contrast = 1,correlations = F,scale="row")
  dev.off()
}else {
  print("No differential region identified!!")
}
# png("maplot.png")
# dba.plotMA(atac)
# dev.off()

#XY plots (with raw and normalized data)
#par(mfrow=c(1,2))
#dba.plotMA(atac,bNormalized=FALSE)
#dba.plotMA(atac,bNormalized=TRUE)

