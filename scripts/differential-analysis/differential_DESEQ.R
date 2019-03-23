#Author: Thomas
#Code for differential count-data analysis using DESEQ (no replicate case)
#Inputfile reqd: count-data file (tsv file),mention outputname
suppressPackageStartupMessages({
require(DESeq)
require(GenomicRanges)
})

#inputfiles 
#bedfile="AsiSi.flank5kb.win.bed"
count_table_file="dsb_pH2Av_UNT-vs-Treat.tsv"
outputname="pH2Av_UNT-vs-Treat"


#read input count-table
df= read.table(count_table_file,sep="\t",row.names=1,header = T)
print("head counttable")
head(df)
df = df[rowSums(df)>0,]
head(df)
condition <- factor( c( "untreated", "treated") )
cds <- newCountDataSet( df, condition )
cds <- estimateSizeFactors( cds )
cds <- estimateDispersions( cds, method="blind", sharingMode="fit-only" )
vsd <- varianceStabilizingTransformation( cds )
res <- nbinomTest( cds,  "untreated" , "treated")

#check MA plot
plotMA(res)

up_fold <- res[res$log2FoldChange > 0,]
down_fold = res[res$log2FoldChange < 0,]

#create dataframe to write output
create_output =function(df){
  df=df[order(df$pval),]
  df=head(df,n=1000L)
  gr1 <- data.frame(region=(df)$id,
                    log2Foldchange=(df)$log2FoldChange,
                    pvalue=(df)$pval)
  return(gr1)
}

up=create_output(up_fold)
down=create_output(down_fold)

#write output
write.table(up, file=paste0(outputname,"_DESEQ_up_sort.bed"), quote=F, sep="\t", row.names=F, col.names=T)
write.table(down, file=paste0(outputname,"_DESEQ_down_sort.bed"), quote=F, sep="\t", row.names=F, col.names=T)
