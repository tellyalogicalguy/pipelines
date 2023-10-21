#!/bin/bash
set -e
set -u
set -o pipefail

#Setting up directories
PROJECT="/home/subrampg/"
PROJECTPARENT_DIR="${PROJECT}binf_analyses/"
SCRATCHPARENT_DIR="/lustre07/scratch/subrampg/binf_analyses_data/"
CUR_DIR="$(cd -P "$(dirname "${BASH_SOURCE[0]}")" 2>&1 >/dev/null && pwd)/"
export SCRIPT_CUR_DIR=$CUR_DIR

DATA_INFO_REL_DIR="../data/"
DATA_INFO_DIR="$(cd "$DATA_INFO_REL_DIR" 2>&1 >/dev/null && pwd)/"
export DATA_INFO_FILE="${DATA_INFO_DIR}/data_info.txt"

echo "Is this dataset single-end or paired-end? (SE/PE)"
read SE_or_PE
case $SE_or_PE in
	SE | PE )
		echo "This data set is: "$SE_or_PE
		export SE_or_PE
		;;
	*)
		echo "Enter valid option. Not proceeding with alignments."
		exit 1
		;;
esac

echo "Enter spike-in genome (dm6/k12/none):"
read SPIKEIN_GENOME
case $SPIKEIN_GENOME in
	dm6 | k12 | none)
		echo "Spike-in genome is: "$SPIKEIN_GENOME
		export SPIKEIN_GENOME
		;;
	*)
		echo "Enter valid spike-in genome option. Not proceeding with alignments."
		exit 1
		;;
esac

echo "Is this a IgG or Input control? (yes/no)"
read CONTROL_OR_NOT
case $CONTROL_OR_NOT in
	yes )
		echo "These are control samples."
		;;
	no )
		echo "There are not control samples."
		echo "Is this ChIPseq or CUTandRUN? (ChIPseq/CUTandRUN)"
		read TECHNIQUE
		case $TECHNIQUE in
			ChIPseq )
				echo "These are "$TECHNIQUE" samples."
				#Input sample BAM for peak calling
				CHIP_INPUT_REL_DIR="../../input/data/filteredBAM/"
				CHIP_INPUT_DIR="$(cd "$CHIP_INPUT_REL_DIR" 2>&1 >/dev/null && pwd)/"
				export CHIP_INPUT_BAM="${CHIP_INPUT_DIR}/Input.bam"
				;;
			CUTandRUN )
				echo "These are "$TECHNIQUE" samples."
				#Input sample BAM for peak calling
				CHIP_INPUT_REL_DIR="../../igg_candr/data/filteredBAM/"
				CHIP_INPUT_DIR="$(cd "$CHIP_INPUT_REL_DIR" 2>&1 >/dev/null && pwd)/"
				export CHIP_INPUT_BAM="${CHIP_INPUT_DIR}/IgG_merged.bam"
				;;
			* )
				echo "Enter valid technique"
				exit 1
				;;
		esac
		;;
	* )
		echo "Enter 'yes' or 'no'"
		exit 1
		;;
esac
export CONTROL_OR_NOT


export PIPELINE_VARS_FILE=${DATA_INFO_DIR}/pipeline_config_vars.sh
echo '#!/bin/bash' > ${PIPELINE_VARS_FILE}
echo 'set -e' >> ${PIPELINE_VARS_FILE}
echo 'set -u' >> ${PIPELINE_VARS_FILE}
echo 'set -o pipefail' >> ${PIPELINE_VARS_FILE}

echo "export SE_or_PE="$SE_or_PE >> ${PIPELINE_VARS_FILE}
echo "export TECHNIQUE="$TECHNIQUE >> ${PIPELINE_VARS_FILE}
echo "export SPIKEIN_GENOME="$SPIKEIN_GENOME >> ${PIPELINE_VARS_FILE}
echo "export CONTROL_OR_NOT="$CONTROL_OR_NOT >> ${PIPELINE_VARS_FILE}


export GENOME_DIR="${SCRATCHPARENT_DIR}/common/bt2_indexes/customCH12genome/" #bt2 index loc

#raw fastq input directory
RAW_FASTQ_REL_DIR="$CUR_DIR/../data/raw_fastq/"
RAW_FASTQ_ABS_DIR="$(cd "$RAW_FASTQ_REL_DIR" 2>&1 >/dev/null && pwd)/"
export RAW_FASTQ_DIR="$(sed "s:$PROJECTPARENT_DIR:$SCRATCHPARENT_DIR:" <(echo $RAW_FASTQ_ABS_DIR))"

#make the desired output directory structures here
FILT_FASTQ_REL_DIR="${CUR_DIR}/../data/filteredFastq/"
UNFILT_BAM_REL_DIR="${CUR_DIR}/../data/alignedUnfiltQNAMEsortedBAM/"
INTERMED_FILES_REL_DIR="${CUR_DIR}/../data/intermediateFiles/"
FILT_BAM_REL_DIR="${CUR_DIR}/../data/filteredBAM/"
PEAKS_REL_DIR="../results/peaks/"

for DIR in FILT_FASTQ_REL_DIR UNFILT_BAM_REL_DIR INTERMED_FILES_REL_DIR FILT_BAM_REL_DIR PEAKS_REL_DIR
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
tree $DATA_INFO_DIR
SCRATCH_DATA_DIR="$(sed "s:$PROJECTPARENT_DIR:$SCRATCHPARENT_DIR:" <(echo $DATA_INFO_DIR))"
echo "scratch data folder structure:"
tree $SCRATCH_DATA_DIR

echo "GENOME_DIR= "$GENOME_DIR
echo "DATA_INFO_FILE= "$DATA_INFO_FILE
echo "FILT_FASTQ_DIR = "$FILT_FASTQ_DIR
echo "FILT_FASTQ_SYML_DIR = "$FILT_FASTQ_SYML_DIR
echo "UNFILT_BAM_DIR = "$UNFILT_BAM_DIR
echo "UNFILT_BAM_SYML_DIR = "$UNFILT_BAM_SYML_DIR
echo "INTERMED_FILES_DIR = "$INTERMED_FILES_DIR
echo "INTERMED_FILES_SYML_DIR = "$INTERMED_FILES_SYML_DIR
echo "FILT_BAM_DIR = "$FILT_BAM_DIR
echo "FILT_BAM_SYML_DIR = "$FILT_BAM_SYML_DIR
echo "PEAKS_DIR = "$PEAKS_DIR
echo "PEAKS_SYML_DIR = "$PEAKS_SYML_DIR
#echo "ChIP Input BAM file = "$CHIP_INPUT_BAM

export BAM_PROG_FILE="${SCRIPT_CUR_DIR}/slurm_outputs/BAM_progress_file.txt"
if [ -e ${BAM_PROG_FILE} ]
then
	rm ${BAM_PROG_FILE}
fi

SAMPLE_IDS=($(awk 'NR>1 {print $1}' ${DATA_INFO_FILE} | uniq))
export FACTOR_NAME=${SAMPLE_IDS[0]/_*/}
for SAMPLE_ID in "${SAMPLE_IDS[@]}"
do
	echo $SAMPLE_ID >> ${BAM_PROG_FILE}
done


#for each read pair, send the variables to sbatch
for SAMPLE_ID in "${SAMPLE_IDS[@]}"
do
	echo "Submitting: "$SAMPLE_ID
	export SAMPLE_ID
	SLRM_NAME="1234_trimAlignFilterMACS_${SAMPLE_ID}"
	if [[ $SE_or_PE == "PE" ]]
	then
		echo "Starting paired-end alignment procedure..."
		sbatch --job-name "$SLRM_NAME" 1234_trimAlignFilterMACS_sbatch.sh
	elif [[ $SE_or_PE == "SE" ]]
	then
		echo "Starting single-end alignment procedure..."
		sbatch --job-name "$SLRM_NAME" 1234_trimAlignFilterMACS_SE_sbatch.sh
	fi
	sleep 0.2
done
