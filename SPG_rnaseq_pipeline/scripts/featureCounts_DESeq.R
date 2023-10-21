#getting variables from bash
genome<-Sys.getenv("GENOME")
spikein_genome<-Sys.getenv("SPIKEIN_GENOME")
data_info_file<-Sys.getenv("DATA_INFO_FILE")
aligned_bam_dir<-Sys.getenv("ALIGNED_BAM_DIR")
thread_count<-as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK"))
gtf_file<-Sys.getenv("GTF_FILE")
fcounts_deseq_dir<-Sys.getenv("FC_DSQ_DIR")
refr_group_name<-Sys.getenv("REFR_GROUP_NAME")
se_or_pe<-Sys.getenv("SE_or_PE")
seq_stranded<-Sys.getenv("SEQ_STRANDED")
# factor_name<-Sys.getenv("FACTOR_NAME")

#set up sampleSheet
library(data.table)
pre_samplesheet<-fread(data_info_file)
pre_samplesheet
samplesheet<-unique(pre_samplesheet[,.(SampleID,Condition,Replicate)])
samplesheet[,bamReads:=paste(aligned_bam_dir,SampleID,"_2pass_Aligned.sortedByCoord.out.bam",sep=""),by=SampleID]
samplesheet[,fcounts_colname:=sub("_2pass_Aligned.sortedByCoord.out.bam","",basename(bamReads)),by=SampleID]

#samplesheet[,Factor:=factor_name]
#samplesheet[,Spikein:=paste(aligned_bam_dir,SampleID,".",spikein_genome,".bam",sep=""),by=SampleID]

# setting variables for featureCounts counting
if(se_or_pe=="PE"){
	paired_setting=T
} else if (se_or_pe=="SE"){
	paired_setting=F
}

if(seq_stranded=="Yes"){
	strand_setting<-2
} else if(seq_stranded=="No"){
	strand_setting<-0
}

cat("Starting featureCounts\n")
library(Rsubread)
fcounts<-featureCounts(files=samplesheet$bamReads, 
	annot.ext=gtf_file, isPairedEnd=paired_setting, nthreads=thread_count,
	isGTFAnnotationFile=T, strandSpecific=strand_setting)

saveRDS(fcounts,paste0(fcounts_deseq_dir,"fcounts_output.rds"))

# adjust colnames for DESeq2 setup
samplesheet[,fcounts_colname:=basename(bamReads),by=SampleID]
colnames(fcounts$counts)<- sapply(colnames(fcounts$counts),
                        function(x){samplesheet[fcounts_colname==x,SampleID]})

# making colData info for DESeq2
colData<-as.data.frame(samplesheet[,.(Condition)])
rownames(colData)<-samplesheet$SampleID
colData$Condition<-factor(colData$Condition)
colData$Condition<-relevel(colData$Condition,ref=refr_group_name)
colData

library(DESeq2)
dds<-DESeqDataSetFromMatrix(countData=fcounts$counts,colData=colData,design=~Condition)

# In case there is a spike-in, do DESeq on spike-in genome first to get norm factor and supply it here
# # supply pre-existing Size Factor
# norm_factors_dt<-fread("../results/diffbind/S2P_norm_factors.txt")
# norm_factors<-1/norm_factors_dt$NormFactors
# names(norm_factors)<-norm_factors_dt$SampleID
# sizeFactors(dds)<-norm_factors

# run DESeq
dds<-DESeq(dds)
res<-results(dds)

res_dt<-as.data.table(res,keep.rownames="gene_id")
norm_counts<-as.data.table(counts(dds,normalized=T),keep.rownames="gene_id")
res_dt<-merge(norm_counts,res_dt,by="gene_id")

# get gene name and merge with table
library(rtracklayer)
annot_gtf<-import(gtf_file)
annot_genes<-annot_gtf[annot_gtf$type=="gene"]
names_table<-data.table(gene_id=annot_genes$gene_id, gene_name=annot_genes$gene_name)
res_dt<-merge(names_table,res_dt)
#---------------------------------------------------------------------------#
#extract up, down, unchanged and unexpressed genes
padj_thr<-0.1
res_dt[is.na(padj),leg:="0_unexpr"]
res_dt[padj<padj_thr & log2FoldChange > log(1.5,2) ,leg:="2_upreg"]
res_dt[padj<padj_thr & log2FoldChange < (-log(1.5,2)) ,leg:="3_downreg"]
res_dt[is.na(leg), leg:="1_unchanged"]
res_dt[,.N,by=leg]


# write results to a data.table
write.table(res_dt,paste0(fcounts_deseq_dir,"normCounts.txt"),sep="\t",quote=F,row.name=F)
#res_dt<-fread("normCounts.txt")
#---------------------------------------------------------------------------#
#do MA plot
res_dt$leg<-factor(res_dt$leg,levels=c("1_unchanged","0_unexpr","2_upreg","3_downreg"))
pdf(paste0(fcounts_deseq_dir,"MA_plot.pdf"),useDingbats=F,height=3,width=5)
library(ggplot2)
require(egg)
require(ggrastr)
p1<-ggplot(res_dt[order(leg)],aes(x=baseMean,y=log2FoldChange,color=leg))+
  rasterize(geom_point(size=0.8,shape=16),dpi=2400)+
  #geom_point(size=0.8,shape=16)+
  scale_x_continuous("Mean of normalized counts",trans="log10")+
  scale_y_continuous("log2 fold change")+
  scale_color_manual(values=c("grey","darkgrey","salmon","skyblue"))+
  geom_hline(linewidth=0.25,yintercept=0)+
  theme_classic()+
  theme(text=element_text(family="sans"),
        axis.line=element_line(color='black',linewidth=0.25),
        axis.ticks=element_line(colour="black",linewidth = 0.25),
        axis.title=element_text(size=7),
        axis.text.y=element_text(size=5),
        axis.text.x=element_text(size=5))
grid.arrange(grobs=lapply(list(p1),set_panel_size,width=unit(3.5,"cm"),height=unit(2.5,"cm")))
dev.off()
#---------------------------------------------------------------------------#
#doing PC1 vs PC2 analysis. As in DESeq2 vignette: http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
vsd <- vst(dds, blind=FALSE)
#trick to add labels to the plot from: https://support.bioconductor.org/p/90791/
z<-plotPCA(vsd, intgroup=c("Condition"))
pdf(paste0(fcounts_deseq_dir,"PC1_vs_PC2_plot.pdf"),useDingbats=F,height=3,width=5)
z+geom_label(aes(label = name)) 
dev.off()
#---------------------------------------------------------------------------#
