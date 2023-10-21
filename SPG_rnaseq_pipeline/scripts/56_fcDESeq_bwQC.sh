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
source ${DATA_INFO_DIR}/pipeline_config_vars.sh
SAMPLE_IDS=($(awk 'NR>1 {print $1}' ${DATA_INFO_FILE} | uniq))
#export FACTOR_NAME=${SAMPLE_IDS[0]/_*/}

ALIGNED_BAM_REL_DIR="../data/alignedBAM/"
FC_DSQ_REL_DIR="../results/fcounts_deseq/"
RMATS_REL_DIR="../results/rMATS/"
BW_REL_DIR="../results/bw_files/"
QC_REL_DIR="../results/qc_stats/"

for DIR in FC_DSQ_REL_DIR BW_REL_DIR QC_REL_DIR ALIGNED_BAM_REL_DIR RMATS_REL_DIR
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

echo "ALIGNED_BAM_DIR="$ALIGNED_BAM_DIR
echo "export FC_DSQ_DIR="$FC_DSQ_DIR >> ${DATA_INFO_DIR}/pipeline_config_vars.sh
echo "export FC_DSQ_SYML_DIR="$FC_DSQ_SYML_DIR >> ${DATA_INFO_DIR}/pipeline_config_vars.sh
echo "export RMATS_DIR="$RMATS_DIR >> ${DATA_INFO_DIR}/pipeline_config_vars.sh
echo "export RMATS_SYML_DIR="$RMATS_SYML_DIR >> ${DATA_INFO_DIR}/pipeline_config_vars.sh
echo "export BW_DIR="$BW_DIR >> ${DATA_INFO_DIR}/pipeline_config_vars.sh
echo "export BW_SYML_DIR="$BW_SYML_DIR >> ${DATA_INFO_DIR}/pipeline_config_vars.sh
echo "export QC_DIR="$QC_DIR >> ${DATA_INFO_DIR}/pipeline_config_vars.sh
echo "export QC_SYML_DIR="$QC_SYML_DIR >> ${DATA_INFO_DIR}/pipeline_config_vars.sh

source ${DATA_INFO_DIR}/pipeline_config_vars.sh
export SCRIPT_CUR_DIR
cat ${DATA_INFO_DIR}/pipeline_config_vars.sh

#echo "FACTOR_NAME = "$FACTOR_NAME
#echo "SPIKEIN_GENOME = "$SPIKEIN_GENOME

if [[ $SEQ_STRANDED == "Yes" ]]
then
	echo "Submitting for DiffBind and making bigwig tracks.. (PE)"
	sbatch --job-name "56_fcDESeq_bwQC" 56_fcDESeq_bwQC_sbatch.sh
elif [[ $SEQ_STRANDED == "No" ]]
then
	echo "Submitting for DiffBind and making bigwig tracks.. (PE)"
	sbatch --job-name "56_fcDESeq_bwQC_unstranded" 56_fcDESeq_bwQC_sbatch_unstranded.sh
fi

