#!/bin/bash
set -e
set -u
set -o pipefail

#Setting up directories
PROJECT="/home/subrampg/"
PROJECTPARENT_DIR="${PROJECT}binf_analyses/"
SCRATCHPARENT_DIR="/lustre07/scratch/subrampg/binf_analyses_data/"
SCRIPT_CUR_DIR="$(cd -P "$(dirname "${BASH_SOURCE[0]}")" 2>&1 >/dev/null && pwd)/"

DATA_INFO_REL_DIR="../data/"
DATA_INFO_DIR="$(cd "$DATA_INFO_REL_DIR" 2>&1 >/dev/null && pwd)/"
DATA_INFO_FILE="${DATA_INFO_DIR}/data_info.txt"

#Setting up pipeline_config_vars file that stores pipeline variables to be sourced by other scritps
export PIPELINE_VARS_FILE=${DATA_INFO_DIR}/pipeline_config_vars.sh
echo '#!/bin/bash' > ${PIPELINE_VARS_FILE}
echo 'set -e' >> ${PIPELINE_VARS_FILE}
echo 'set -u' >> ${PIPELINE_VARS_FILE}
echo 'set -o pipefail' >> ${PIPELINE_VARS_FILE}
echo "export SCRIPT_CUR_DIR="$SCRIPT_CUR_DIR >> ${PIPELINE_VARS_FILE}
echo "export DATA_INFO_DIR="$DATA_INFO_DIR >> ${PIPELINE_VARS_FILE}
echo "export DATA_INFO_FILE="$DATA_INFO_FILE >> ${PIPELINE_VARS_FILE}

#getting samples and sequencing info
#-----------------------------------
GRP_NAMES=($(awk 'NR>1 {print $4}' $DATA_INFO_FILE | sort | uniq))
if [[ ${#GRP_NAMES[@]} -ne 2 ]]
then
	echo "This number of sample groups cannot be handled by this pipeline: "${GRP_NAMES[@]}
	echo "This pipeline currently can only handle exactly two sample groups."
	echo "Please modify the data_info.txt accordingly to proceed further."
	exit 1
fi

echo "What is the reference group name? ( ${GRP_NAMES[1]} / ${GRP_NAMES[0]} )"
read REFR_GROUP_NAME
if [[ $REFR_GROUP_NAME == ${GRP_NAMES[0]} ]]
then
	TEST_GROUP_NAME=${GRP_NAMES[1]}
elif [[ $REFR_GROUP_NAME == ${GRP_NAMES[1]} ]]
then
	TEST_GROUP_NAME=${GRP_NAMES[0]}
else
	echo "Reference group name entered incorrectly. Not proceeding with alignments."
	exit 1
fi
echo "Reference group name is: "$REFR_GROUP_NAME"; and test group name is: "$TEST_GROUP_NAME
echo "export REFR_GROUP_NAME="$REFR_GROUP_NAME >> ${PIPELINE_VARS_FILE}
echo "export TEST_GROUP_NAME="$TEST_GROUP_NAME >> ${PIPELINE_VARS_FILE}

# SE_or_PE
echo "Is this dataset single-end or paired-end? (SE/PE)"
read SE_or_PE
case $SE_or_PE in
	SE | PE )
		echo "This data set is: "$SE_or_PE
		echo "export SE_or_PE="$SE_or_PE >> ${PIPELINE_VARS_FILE}
		;;
	*)
		echo "Enter valid option. Not proceeding with alignments."
		exit 1
		;;
esac

# stranded or not
echo "Is the sequencing in this dataset stranded? ( Yes / No )"
read SEQ_STRANDED
case $SEQ_STRANDED in
	Yes | No )
		echo "This data set is: "$SEQ_STRANDED
		echo "export SEQ_STRANDED="$SEQ_STRANDED >> ${PIPELINE_VARS_FILE}
		;;
	*)
		echo "Enter valid option. Not proceeding with alignments."
		exit 1
		;;
esac

# Genome
echo "Enter genome to align to (custom_ch12/mm10/hg38):"
read GENOME
case $GENOME in
	mm10 )
		GTF_FILE="${SCRATCHPARENT_DIR}common/mm10/gencode.vM23.primary_assembly.annotation.gtf"
		GFF3_FILE="${SCRATCHPARENT_DIR}common/mm10/gencode.vM23.primary_assembly.annotation.gff3"
		CHROM_SIZE_FILE="${SCRATCHPARENT_DIR}common/mm10/GRCm38.primary_assembly.genome.fa.fai"
		;;
	hg38 )
		GTF_FILE="${SCRATCHPARENT_DIR}common/hg38/gencode.v43.primary_assembly.annotation.gtf"
		GFF3_FILE="${SCRATCHPARENT_DIR}common/hg38/gencode.v43.primary_assembly.annotation.gff3"
		CHROM_SIZE_FILE="${SCRATCHPARENT_DIR}common/hg38/GRCh38.primary_assembly.genome.fa.fai"
		;;
	custom_ch12 )
		GTF_FILE="${SCRATCHPARENT_DIR}customCH12genome/reformed_igk/gencode.vM23.primary_assembly.customCH12.annotation.9cols.gtf"
		GFF3_FILE="${SCRATCHPARENT_DIR}customCH12genome/reformed_igk/gencode.vM23.primary_assembly.annotation_reformed_reformed.gff3"
		CHROM_SIZE_FILE="${SCRATCHPARENT_DIR}customCH12genome/reformed_igk/GRCm38.primary_assembly.genome_reformed_reformed.fa.fai"
		;;
	*)
		echo "Enter valid genome option. Not proceeding with alignments."
		exit 1
		;;
esac
echo "Genome to align to is: "$GENOME

# sjdbOverhang
echo "Enter sjdbOverhang value to use for STAR alignment to $GENOME genome (49 / 75 / 99 / 149):"
read SJDB_OVHG_VAL
case $SJDB_OVHG_VAL in
	149 | 99 | 75 | 49 )
		echo "The chosen sjdbOverhang value is: "$SJDB_OVHG_VAL
		;;
	*)
		echo "Enter valid option. Not proceeding with alignments."
		exit 1
		;;
esac

# spike-in genome
echo "Enter spike-in genome (dm6/k12/none):"
read SPIKEIN_GENOME
case $SPIKEIN_GENOME in
	dm6 | k12 | none )
		echo "Spike-in genome is: "$SPIKEIN_GENOME
		echo "export SPIKEIN_GENOME="$SPIKEIN_GENOME >> ${PIPELINE_VARS_FILE}
		# add chromosome file size and gtf file
		#export SPIKEIN_GENOME
		;;
	*)
		echo "Enter valid spike-in genome option. Not proceeding with alignments."
		exit 1
		;;
esac

# adding genome info to pipeline_config_vars.sh 
echo "export GENOME="$GENOME >> ${PIPELINE_VARS_FILE}
echo "export SJDB_OVHG_VAL="$SJDB_OVHG_VAL >> ${PIPELINE_VARS_FILE}
GENOME_DIR="${SCRATCHPARENT_DIR}/common/star_indexes/${GENOME}/sjdbOH_${SJDB_OVHG_VAL}/"
echo "export GENOME_DIR="$GENOME_DIR >> ${PIPELINE_VARS_FILE}
echo "export GTF_FILE="$GTF_FILE >> ${PIPELINE_VARS_FILE}
echo "export GFF3_FILE="$GFF3_FILE >> ${PIPELINE_VARS_FILE}
echo "export CHROM_SIZE_FILE="$CHROM_SIZE_FILE >> ${PIPELINE_VARS_FILE}


#raw fastq input directory
RAW_FASTQ_REL_DIR="$SCRIPT_CUR_DIR/../data/raw_fastq/"
RAW_FASTQ_ABS_DIR="$(cd "$RAW_FASTQ_REL_DIR" 2>&1 >/dev/null && pwd)/"
RAW_FASTQ_DIR="$(sed "s:$PROJECTPARENT_DIR:$SCRATCHPARENT_DIR:" <(echo $RAW_FASTQ_ABS_DIR))"
echo "export RAW_FASTQ_DIR="$RAW_FASTQ_DIR >> ${PIPELINE_VARS_FILE}

#make the desired output directory structures here
FILT_FASTQ_REL_DIR="${SCRIPT_CUR_DIR}/../data/filteredFastq/"
ALIGNED_BAM_REL_DIR="${SCRIPT_CUR_DIR}/../data/alignedBAM/"
BW_REL_DIR="../results/bw_files/"

for DIR in FILT_FASTQ_REL_DIR ALIGNED_BAM_REL_DIR BW_REL_DIR
do
	if [ ! -d ${!DIR} ]
	then
		mkdir -p "${!DIR}" || exit 1;
	fi

	SYML_DIR_ABS_PATH="$(cd "${!DIR}" 2>&1 >/dev/null && pwd)/"
	OUT_SUBDIRS="$(sed "s:$PROJECTPARENT_DIR::" <(echo $SYML_DIR_ABS_PATH))"
	SCRATCH_DIR_ABS_PATH="${SCRATCHPARENT_DIR}${OUT_SUBDIRS}"
	if [ ! -d ${SCRATCH_DIR_ABS_PATH} ]
	then
		mkdir -p ${SCRATCH_DIR_ABS_PATH}
	fi
	SYML_DIR_VARNAME="$(sed "s:_REL:_SYML:" <(echo $DIR))"
	SCRATCH_DIR_VARNAME="$(sed "s:_REL::" <(echo $DIR))"
	eval export ${SYML_DIR_VARNAME}="${SYML_DIR_ABS_PATH}"
	eval export ${SCRATCH_DIR_VARNAME}="${SCRATCH_DIR_ABS_PATH}"
done
echo "local data folder structure:"
#tree $DATA_INFO_DIR
SCRATCH_DATA_DIR="$(sed "s:$PROJECTPARENT_DIR:$SCRATCHPARENT_DIR:" <(echo $DATA_INFO_DIR))"
echo "scratch data folder structure:"
#tree $SCRATCH_DATA_DIR

echo "export FILT_FASTQ_DIR="$FILT_FASTQ_DIR >> ${PIPELINE_VARS_FILE}
echo "export FILT_FASTQ_SYML_DIR="$FILT_FASTQ_SYML_DIR >> ${PIPELINE_VARS_FILE}
echo "export ALIGNED_BAM_DIR="$ALIGNED_BAM_DIR >> ${PIPELINE_VARS_FILE}
echo "export ALIGNED_BAM_SYML_DIR="$ALIGNED_BAM_SYML_DIR >> ${PIPELINE_VARS_FILE}
echo "export BW_DIR="$BW_DIR >> ${PIPELINE_VARS_FILE}
echo "export BW_SYML_DIR="$BW_SYML_DIR >> ${PIPELINE_VARS_FILE}

BAM_PROG_FILE="${SCRIPT_CUR_DIR}/slurm_outputs/BAM_progress_file.txt"
BAM_PROG_2P_FILE="${SCRIPT_CUR_DIR}/slurm_outputs/BAM_progress_2Pass_file.txt"
echo "export BAM_PROG_FILE="$BAM_PROG_FILE >> ${PIPELINE_VARS_FILE}
echo "export BAM_PROG_2P_FILE="$BAM_PROG_2P_FILE >> ${PIPELINE_VARS_FILE}
if [ -e ${BAM_PROG_FILE} ]
then
	rm ${BAM_PROG_FILE}
	rm ${BAM_PROG_2P_FILE}
fi

NORM_FACTOR_FILE="${SCRIPT_CUR_DIR}../data/normFactor_file.txt"
echo "export NORM_FACTOR_FILE="$NORM_FACTOR_FILE >> ${PIPELINE_VARS_FILE}
if [ -e ${NORM_FACTOR_FILE} ]
then
	rm ${NORM_FACTOR_FILE}
fi
echo -e "SampleID\tNorm" >> ${NORM_FACTOR_FILE}

SAMPLE_IDS=($(awk 'NR>1 {print $1}' ${DATA_INFO_FILE} | uniq))
#export FACTOR_NAME=${SAMPLE_IDS[0]/_*/}
for SAMPLE_ID in "${SAMPLE_IDS[@]}"
do
	echo $SAMPLE_ID >> ${BAM_PROG_FILE}
	echo $SAMPLE_ID >> ${BAM_PROG_2P_FILE}
	echo -e "$SAMPLE_ID\t1" >> ${NORM_FACTOR_FILE}
done

source ${DATA_INFO_DIR}/pipeline_config_vars.sh
cat ${DATA_INFO_DIR}/pipeline_config_vars.sh

# Write code to accomodate spike-in. Here, if spike-in, do alignment and DESeq first for spikein genome and then do normFactor calculation to overwrite the normFactor file..
# Write files to a differenct folder like ../data/spikein_aligned/ folder

#for each read pair, send the variables to sbatch
for SAMPLE_ID in "${SAMPLE_IDS[@]}"
do
	echo "Submitting: "$SAMPLE_ID
	export SAMPLE_ID
	SLRM_NAME="12_trimAlign_${SAMPLE_ID}"
	if [[ $SE_or_PE == "PE" ]]
	then
		echo "Starting paired-end alignment procedure..."
		sbatch --job-name "$SLRM_NAME" 12_trimAlign_sbatch.sh
	elif [[ $SE_or_PE == "SE" ]]
	then
		echo "Starting single-end alignment procedure..."
		sbatch --job-name "$SLRM_NAME" 12_trimAlign_sbatch_SE.sh
	fi
	sleep 0.2
done

