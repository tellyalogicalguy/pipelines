#!/bin/bash
#SBATCH --account=rpp-jfcote11
#SBATCH --mail-user=poorani.subramani@mail.mcgill.ca
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=8000
#SBATCH --time=2:00:0
#SBATCH --tmp=100G
#SBATCH --output="./slurm_outputs/%x_slurm-%j.out"

set -e
set -u
set -o pipefail

#---------------------------------------------------------#
module load gcc/9.3.0
module load r/4.2.1

# do diffbind only if not control
if [[ $CONTROL_OR_NOT == "no" ]]
then

echo $(date -Iseconds)" SPG_ChIPseq_pipeline: 5. Starting DiffBind and ChIPQC.."
Rscript ${SCRIPT_CUR_DIR}/diffBind_ChIPQC.R
echo $(date -Iseconds)" SPG_ChIPseq_pipeline: 5. Done with DiffBind and ChIPQC"

fi
#---------------------------------------------------------#
module load python/3.8.10
virtualenv --no-download $SLURM_TMPDIR/env
set +u
source $SLURM_TMPDIR/env/bin/activate
set -u
pip install --no-index --upgrade pip
pip install --no-index -r ${SCRIPT_CUR_DIR}/deepTools_requirements.txt

NORM_FACTORS_FILE="${DB_DIR}/${FACTOR_NAME}_norm_factors.txt"

echo $(date -Iseconds)" SPG_ChIPseq_pipeline: 6. Making bigwig files.."
# make bw files
SAMPLE_IDS=($(awk 'NR>1 {print $1}' ${DATA_INFO_FILE} | uniq))
for SAMPLE_ID in ${SAMPLE_IDS[@]}
do
	FILTERED_BAM="${FILT_BAM_DIR}/${SAMPLE_ID}.bam"
	echo  $(date -Iseconds)" SPG_ChIPseq_pipeline: 6. Making bigwig files for "$SAMPLE_ID
	
	if [[ $SE_or_PE == "SE" ]]
	then
		FRAG_LENGTH=$(grep -Po '(?<=(predicted fragment length is ))[0-9]*(?= bps)' ${PEAKS_DIR}/${SAMPLE_ID}_macs2.log)
		echo "MACS2 predicted fragment length is "$FRAG_LENGTH
		EXTEND_OPT="--extendReads ${FRAG_LENGTH} --centerReads"
	elif [[ $SE_or_PE == "PE" ]]
	then
        	EXTEND_OPT="--extendReads"
	fi

	if [[ $SPIKEIN_GENOME != "none" ]]
	then
		NORM_FACTOR=$(awk -v S_ID=$SAMPLE_ID '$1==S_ID {print $2}' $NORM_FACTORS_FILE)
		echo "with scalingFactor: "$NORM_FACTOR
		bamCoverage -b $FILTERED_BAM -o ${BW_DIR}/${SAMPLE_ID}.bw \
		-of "bigwig" --scaleFactor ${NORM_FACTOR} --binSize 20 \
		-p $SLURM_CPUS_PER_TASK --normalizeUsing "None" \
		${EXTEND_OPT}
	elif [[ $SPIKEIN_GENOME == "none" ]]
	then
		echo "with lib size normalization (RPGC)"
		bamCoverage -b $FILTERED_BAM -o ${BW_DIR}/${SAMPLE_ID}.bw \
		-of "bigwig" --binSize 20 -p $SLURM_CPUS_PER_TASK \
		--normalizeUsing "RPGC" --effectiveGenomeSize 2652783500 \
		${EXTEND_OPT}
	fi
done
#---------------------------------------------------------#
# setting up for bigwig merging
module load gcc/9.3.0
module load htslib/1.16
module load gsl/2.6
module load kentutils/401
module load python/3.8.10
# adding paths to libBwigWig
export CPATH=$CPATH:/home/subrampg/local_installs/usr_local_substitute/libBigWig/include/
export LIBRARY_PATH=$LIBRARY_PATH:/home/subrampg/local_installs/usr_local_substitute/libBigWig/lib/
export LD_LIBRARY_PATH=/home/subrampg/local_installs/usr_local_substitute/libBigWig/lib/
# adding wiggletools to $PATH
export PATH=$PATH:/home/subrampg/local_installs/WiggleTools/bin
CHROM_SIZE_FILE=/lustre07/scratch/subrampg/binf_analyses_data/customCH12genome/reformed_igk/GRCm38.primary_assembly.genome_reformed_reformed.chrom.sizes
#wiggletools mean ${CTCF_FILES[@]} | wigToBigWig stdin $CHROM_SIZE_FILE CTCF_WT_merged.bw
#--------------------------------#
# do bw merging only if not control
if [[ $CONTROL_OR_NOT == "yes" ]]
then
#if more than one control (1 WT and 1 KO IgG for e.g.), merge bam files
echo $(date -Iseconds)" SPG_ChIPseq_pipeline: 6. Merging control BAM files.."
module load samtools/1.16.1
SLRM_LOCAL_TMP_DIR=$SLURM_TMPDIR/$RANDOM/${FACTOR_NAME}/
mkdir -p ${SLRM_LOCAL_TMP_DIR}

INPUT_SAMPLES=($(awk 'NR>1 {print $1}' $DATA_INFO_FILE | uniq))
for i in ${!INPUT_SAMPLES[@]}; do
	INPUT_BAMS[$i]="${FILT_BAM_DIR}/${INPUT_SAMPLES[$i]}.bam"
	INPUT_BWS[$i]="${BW_DIR}/${INPUT_SAMPLES[$i]}.bw"
done

INPUT_MERGED_BAM="${FILT_BAM_DIR}/${FACTOR_NAME}_merged.bam"
INPUT_MERGED_BAM_LOCAL="${SLRM_LOCAL_TMP_DIR}/${FACTOR_NAME}_merged.bam"
samtools merge -o ${INPUT_MERGED_BAM_LOCAL} ${INPUT_BAMS[@]} \
-@ $SLURM_CPUS_PER_TASK

cp ${INPUT_MERGED_BAM_LOCAL} ${FILT_BAM_DIR}/

samtools index ${INPUT_MERGED_BAM} ${INPUT_MERGED_BAM}.bai

#merging input bw files
if [ ${#INPUT_SAMPLES[@]} -eq 1 ]
then
        echo "There is only 1 input sample."
	INPUT_MERGED_BW="${BW_DIR}/${INPUT_SAMPLES[0]}.bw"
	#exit 0;
else
	echo $(date -Iseconds)" SPG_ChIPseq_pipeline: 6. Merging ${#INPUT_SAMPLES[@]} control bigwig files.."
	INPUT_MERGED_BW="${BW_DIR}/${FACTOR_NAME}_merged.bw"
	wiggletools mean ${INPUT_BWS[@]} | wigToBigWig stdin $CHROM_SIZE_FILE ${INPUT_MERGED_BW}
fi
#---------------------------------------------------------#
#CONTROL_OR_NOT check
elif [[ $CONTROL_OR_NOT == "no" ]]
then
#Merge the same genotypes
#Will only work for the first 2 replicates of WT and KO
#Need to be extended for more replicates if necessary!!

#merging WT bw files
WT_SAMPLES=($(awk '$4=="WT" {print $1}' $DATA_INFO_FILE | uniq))

if [ ${#WT_SAMPLES[@]} -eq 1 ]
then
	echo "There is only 1 WT sample."
	WT_MERGED_BW="${BW_DIR}/${WT_SAMPLES[0]}.bw"
else
	for i in ${!WT_SAMPLES[@]}; do 
		WT_BWS[$i]="${BW_DIR}/${WT_SAMPLES[$i]}.bw"
	done

	echo $(date -Iseconds)" SPG_ChIPseq_pipeline: 6. Merging ${#WT_SAMPLES[@]} WT bigwig files.."
	WT_MERGED_BW="${BW_DIR}/${FACTOR_NAME}_WT_merged.bw"
	wiggletools mean ${WT_BWS[@]} | wigToBigWig stdin $CHROM_SIZE_FILE ${WT_MERGED_BW}
fi

# merging KO bw files
KO_SAMPLES=($(awk '$4=="KO" {print $1}' $DATA_INFO_FILE | uniq))

if [ ${#KO_SAMPLES[@]} -eq 0 ]
then
	echo "No KO samples detected. Prcoessed only WT samples."
elif [ ${#KO_SAMPLES[@]} -eq 1 ]
then
	echo "There is only 1 KO sample."
	KO_MERGED_BW="${BW_DIR}/${KO_SAMPLES[0]}.bw"
else

	for i in ${!KO_SAMPLES[@]}; do 
		KO_BWS[$i]="${BW_DIR}/${KO_SAMPLES[$i]}.bw"
	done

	echo $(date -Iseconds)" SPG_ChIPseq_pipeline: 6. Merging ${#KO_SAMPLES[@]} KO bigwig files.."
	KO_MERGED_BW="${BW_DIR}/${FACTOR_NAME}_KO_merged.bw"
	wiggletools mean ${KO_BWS[@]} | wigToBigWig stdin $CHROM_SIZE_FILE ${KO_MERGED_BW}
fi

#---------------------------------------------------------#
# do KO-WT comparison only if KO sample is present
if [ ${#KO_SAMPLES[@]} -ne 0 ]
then
	echo $(date -Iseconds)" SPG_ChIPseq_pipeline: 6. Making KO-WT bigwig files.."

	if [ ${#WT_SAMPLES[@]} -ne ${#KO_SAMPLES[@]} ]
	then
		echo "--------------------------------------------------------------"
		echo "Warning: The number of WT and KO replicates are not the same!!"
		echo "--------------------------------------------------------------"
	fi

	#do diff
	bigwigCompare -b1 ${KO_MERGED_BW} -b2 ${WT_MERGED_BW} \
	-o "${BW_DIR}/${FACTOR_NAME}_KO-WT.bw" -of "bigwig" \
	--scaleFactors 1:1 --operation "subtract" --binSize 20 \
	-p $SLURM_CPUS_PER_TASK --verbose
fi
#---------------------------------------------------------#

fi #end of CONTROL_OR_NOT check

#---------------------------------------------------------#

deactivate
rm -rf $SLURM_TMPDIR/env
#---------------------------------------------------------#
echo $(date -Iseconds)" SPG_ChIPseq_pipeline: 7. Compiling QC metrics.."
module load python/3.8.10
virtualenv --no-download $SLURM_TMPDIR/env
set +u
source $SLURM_TMPDIR/env/bin/activate
set -u
pip install --no-index --upgrade pip
pip install --no-index -r ${SCRIPT_CUR_DIR}/multiqc_requirements.txt

multiqc ${BW_DIR}/../../ -o ${QC_DIR}/ -n "${QC_DIR}/${FACTOR_NAME}_multiQC_report"

#making symlinks
for SYML_DIR in DB_SYML_DIR BW_SYML_DIR QC_SYML_DIR FILT_BAM_SYML_DIR
do
	SCRATCH_DIR_VARNAME="$(sed "s:_SYML::" <(echo ${SYML_DIR}))"
	eval SCRATCH_DIR='$'"${SCRATCH_DIR_VARNAME}"
	echo $SCRATCH_DIR_VARNAME" = "$SCRATCH_DIR
	for file in $(ls ${SCRATCH_DIR})
	do
	if [ ! -L "${!SYML_DIR}${file}" ]; then
		ln -s ${SCRATCH_DIR}${file} ${!SYML_DIR}${file}
	fi
	done

find ${!SYML_DIR} -xtype l -delete

done


