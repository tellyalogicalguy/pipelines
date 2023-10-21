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

DATA_INFO_REL_DIR="${CUR_DIR}/../data/"
DATA_INFO_DIR="$(cd "$DATA_INFO_REL_DIR" 2>&1 >/dev/null && pwd)/"
export DATA_INFO_FILE="${DATA_INFO_DIR}/data_info.txt"
source ${DATA_INFO_DIR}/pipeline_config_vars.sh
SAMPLE_IDS=($(awk 'NR>1 {print $1}' ${DATA_INFO_FILE} | uniq))
export FACTOR_NAME=${SAMPLE_IDS[0]/_*/}

FILT_BAM_REL_DIR="../data/filteredBAM/"
PEAKS_REL_DIR="../results/peaks/"
#FILT_BAM_SYML_DIR="$(cd "$FILT_BAM_REL_DIR" 2>&1 >/dev/null && pwd)/"
#export FILT_BAM_DIR="$(sed "s:$PROJECTPARENT_DIR:$SCRATCHPARENT_DIR:" <(echo $FILT_BAM_SYML_DIR))"

DB_REL_DIR="../results/diffbind/"
BW_REL_DIR="../results/bw_files/"
QC_REL_DIR="../results/qc_stats/"

for DIR in DB_REL_DIR BW_REL_DIR QC_REL_DIR FILT_BAM_REL_DIR PEAKS_REL_DIR
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

echo "DB_DIR = "$DB_DIR
echo "DB_SYML_DIR = "$DB_SYML_DIR
echo "BW_DIR = "$BW_DIR
echo "BW_SYML_DIR = "$BW_SYML_DIR
echo "QC_DIR = "$QC_DIR
echo "QC_SYML_DIR = "$QC_SYML_DIR
echo "FACTOR_NAME = "$FACTOR_NAME
echo "FILT_BAM_DIR = "$FILT_BAM_DIR
echo "PEAKS_DIR = "$PEAKS_DIR
echo "SPIKEIN_GENOME = "$SPIKEIN_GENOME

if [[ $SE_or_PE == "PE" ]]
then
	echo "Submitting for DiffBind and making bigwig tracks.. (PE)"
	sbatch --job-name "567_DBbwQC" 567_DBbwQC_sbatch.sh
elif [[ $SE_or_PE == "SE" ]]
then
	echo "Submitting for DiffBind and making bigwig tracks.. (SE)"
	sbatch --job-name "567_DBbwQC" 567_DBbwQC_sbatch.sh
fi

