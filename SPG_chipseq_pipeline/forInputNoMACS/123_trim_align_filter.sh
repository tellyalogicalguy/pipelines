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

export GENOME_DIR="${SCRATCHPARENT_DIR}/common/bt2_indexes/customCH12genome/" #bt2 index loc

DATA_INFO_REL_DIR="../data/"
DATA_INFO_DIR="$(cd "$DATA_INFO_REL_DIR" 2>&1 >/dev/null && pwd)/"
export DATA_INFO_FILE="${DATA_INFO_DIR}/data_info.txt"

#input directory
RAW_FASTQ_REL_DIR="$CUR_DIR/../data/raw_fastq/"
RAW_FASTQ_ABS_DIR="$(cd "$RAW_FASTQ_REL_DIR" 2>&1 >/dev/null && pwd)/"
#optional. changes path to original files in scratch instead of the symlinks in home
export RAW_FASTQ_DIR="$(sed "s:$PROJECTPARENT_DIR:$SCRATCHPARENT_DIR:" <(echo $RAW_FASTQ_ABS_DIR))"

#make the desired output directory structures here
FILT_FASTQ_REL_DIR="${CUR_DIR}/../data/filteredFastq/"
UNFILT_BAM_REL_DIR="${CUR_DIR}/../data/alignedUnfiltQNAMEsortedBAM/"
INTERMED_FILES_REL_DIR="${CUR_DIR}/../data/intermediateFiles/"
FILT_BAM_REL_DIR="${CUR_DIR}/../data/filteredBAM/"

for DIR in FILT_FASTQ_REL_DIR UNFILT_BAM_REL_DIR INTERMED_FILES_REL_DIR FILT_BAM_REL_DIR
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

#for each read pair, send the variables to sbatch
SAMPLE_IDS=($(awk 'NR>1 {print $1}' ${DATA_INFO_FILE} | uniq))
for SAMPLE_ID in "${SAMPLE_IDS[@]}"
do
	echo "Submitting: "$SAMPLE_ID
	export SAMPLE_ID
	SLRM_NAME="123_trim_align_filter_${SAMPLE_ID}"
	sbatch --job-name "$SLRM_NAME" 123_trim_align_filter_sbatch.sh
	sleep 0.2
done
