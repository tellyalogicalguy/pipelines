#getting variables from bash
spikein_genome<-Sys.getenv("SPIKEIN_GENOME")
data_info_file<-Sys.getenv("DATA_INFO_FILE")
filt_bam_dir<-Sys.getenv("FILT_BAM_DIR")
peaks_dir<-Sys.getenv("PEAKS_DIR")
db_out_dir<-Sys.getenv("DB_DIR")
factor_name<-Sys.getenv("FACTOR_NAME")
qc_out_dir<-Sys.getenv("QC_DIR")

#set up sampleSheet
library(data.table)
pre_samplesheet<-fread(data_info_file)
pre_samplesheet
samplesheet<-unique(pre_samplesheet[,.(SampleID,Condition,Replicate)])
samplesheet[,Tissue:="CH12"]
samplesheet[,Treatment:="Full-media"]
samplesheet[,Factor:=factor_name]
samplesheet[,bamReads:=paste(filt_bam_dir,SampleID,".bam",sep=""),by=SampleID]
samplesheet[,Peaks:=paste(peaks_dir,SampleID,"_peaks.narrowPeak",sep=""),by=SampleID]
samplesheet[,PeakCaller:="narrow"]
samplesheet[,Spikein:=paste(filt_bam_dir,SampleID,".",spikein_genome,".bam",sep=""),by=SampleID]
db_samplesheet_file<-paste(db_out_dir,factor_name,"_db_samplesheet.csv",sep="")
write.csv(samplesheet,db_samplesheet_file,row.names=F)

#do DiffBind
library(DiffBind)
sample_db<-dba(sampleSheet=samplesheet)
sample_counts<-dba.count(sample_db)

#do only if we have a KO sample
if ("KO"%in%samplesheet[,Condition]){

#if only one WT or KO sample
if (max(samplesheet[Condition=="WT",Replicate])<2 || max(samplesheet[Condition=="KO",Replicate])<2) {

	if (spikein_genome != "none"){
		sample_model<-dba.normalize(sample_counts,spikein=T)
	} else if (spikein_genome == "none"){
		sample_model<-dba.normalize(sample_counts)
	}
}

#do only if we have two samples each
else if (max(samplesheet[Condition=="WT",Replicate])>=2 && max(samplesheet[Condition=="KO",Replicate])>=2) {

	sample_model<-dba.contrast(sample_counts, minMembers=2,reorderMeta=list(Condition="WT"))

	if (spikein_genome == "dm6"){
		sample_model<-dba.normalize(sample_model,normalize="RLE",background = T,spikein = T)
	} else if (spikein_genome == "k12"){
		sample_model<-dba.normalize(sample_model,spikein = T)
	} else if (spikein_genome == "none"){
		sample_model<-dba.normalize(sample_model,normalize="RLE",background = T)
	}

	sample_model<-dba.analyze(sample_model,bBlacklist = F)
	cat("Writing files from R\n")

	#save objects
	rdataobj_name<-paste(db_out_dir,factor_name,"_spikeinNorm_csawBg_DBobjs.rda",sep="")
	save(sample_model,sample_counts,sample_db,file=rdataobj_name)

	#save MA plot
	ma_plot_name=paste(db_out_dir,factor_name,"_KOvsWT_MAplot.pdf",sep="")
	pdf(file = ma_plot_name,width = 3,height = 3)
	dba.plotMA(sample_model, bNormalized = T, sub="Norm: Spikein RLE (csaw bg)",th=0.1)
	dev.off()

	#save DB 
	sample_db_report<-dba.report(sample_model,th=1)
	db_all_GRanges_name<-paste(db_out_dir,factor_name,"_DB_all_GRanges.rds",sep="")
	saveRDS(sample_db_report,file=db_all_GRanges_name)
	#only diff bound
	sample_db_report<-dba.report(sample_model,th=0.1)
	db_only_GRanges_name<-paste(db_out_dir,factor_name,"_DB_only_GRanges.rds",sep="")
	saveRDS(sample_db_report,file=db_only_GRanges_name)

	db_only_txt_name<-paste(db_out_dir,factor_name,"_DB_only.txt",sep="")
	write.table(as.data.table(sample_db_report),db_only_txt_name,quote=F,row.names=F,sep="\t")
}

#save norm factors
norm_factors<-data.table(SampleID=sample_db$samples$SampleID,NormFactors=(1/sample_model$norm$DESeq2$norm.facs))
norm_factors_name<-paste(db_out_dir,factor_name,"_norm_factors.txt",sep="")
write.table(norm_factors,norm_factors_name,quote=F,row.names=F,sep="\t")
} else {
#make dummy file for facilitating following steps. values are not used
norm_factors<-data.table(SampleID=sample_db$samples$SampleID,NormFactors=1)
norm_factors_name<-paste(db_out_dir,factor_name,"_norm_factors.txt",sep="")
write.table(norm_factors,norm_factors_name,quote=F,row.names=F,sep="\t")
}

cat("Writing consensus peaks files from R\n")
#save "all" peaks
sample_peakset<-dba.peakset(sample_db,bRetrieve=T)
peakset_GRangesfile_name<-paste(db_out_dir,factor_name,"_all_peakset_GRanges.rds",sep="")
saveRDS(sample_peakset,file=peakset_GRangesfile_name)
peakset_df<-as.data.frame(sample_peakset)
peakset_bed_name<-paste(db_out_dir,factor_name,"_all_peaks.bed",sep="")
library(rtracklayer)
export(sample_peakset,con=peakset_bed_name)

cat("Finished writing DiffBind files from R\n")


cat("Doing ChIP-seq QC with ChIPQC\n")
library(ChIPQC)
samplesheet<-samplesheet[,!"Spikein"]
data(blacklist_mm10,package="Signac")
chr_list_vector<-seqlevels(blacklist_mm10)
chipObj<-ChIPQC(samplesheet, annotation="mm10",chromosomes=chr_list_vector,
consensus=T,bCount=T,summits=250)
chipqc_report_name<-paste(factor_name,"_ChIPQC_report",sep="")
chipqc_report_folder<-paste(qc_out_dir,"/ChIPQC_report/",sep="")
ChIPQCreport(chipObj,reportName=chipqc_report_name,reportFolder=chipqc_report_folder)
cat("ChIPQC Done.\n")

