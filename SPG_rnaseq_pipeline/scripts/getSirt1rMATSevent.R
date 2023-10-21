data_info_dir<-Sys.getenv("DATA_INFO_DIR")
data_info_file<-Sys.getenv("DATA_INFO_FILE")
pipeline_vars_file<-paste0(data_info_dir,"pipeline_config_vars.sh")
refr_group_name<-Sys.getenv("REFR_GROUP_NAME")
test_group_name<-Sys.getenv("TEST_GROUP_NAME")

rmats_dir<-Sys.getenv("RMATS_DIR")
jc_file<-paste0(rmats_dir,"SE.MATS.JCEC.txt")

library(data.table)
jc_dt<-fread(jc_file)
jc_dt<-jc_dt[FDR<0.1]
jc_sirt1_dt<-jc_dt[grep("Sirt1",geneSymbol,ignore.case=T)]

if (nrow(jc_sirt1_dt)>0) {
	write("export SIRT1_RMATS_EVENT=Yes", file=pipeline_vars_file, append=T)
	write.table(jc_sirt1_dt, paste0(rmats_dir,"Sirt1_rMATS_SE_event.txt"),
		row.names=F, quote=F, sep="\t")

	# make group_info.gf file for rMATS2sashimiplot
	group_info_file<-paste0(rmats_dir,"groupInfoForSashimi.gf")
	data_info_dt<-fread(data_info_file)
	data_info_dt<-unique(data_info_dt[,.(SampleID,Condition,Replicate)])
	
	noOfRefSamples<-nrow(data_info_dt[Condition==refr_group_name])
	refr_indices<-1:noOfRefSamples
	test_indices<-1:nrow(data_info_dt[Condition==test_group_name]) + noOfRefSamples

	write(paste(paste0(refr_group_name,":"), paste(refr_indices,collapse=",")),
		group_info_file)
	write(paste(paste0(test_group_name,":"), paste(test_indices,collapse=",")),
		group_info_file, append=T)

} else if (nrow(jc_sirt1_dt)==0) {
	write("export SIRT1_RMATS_EVENT=No", file=pipeline_vars_file, append=T)
}


