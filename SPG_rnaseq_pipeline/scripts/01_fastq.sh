#!/bin/bash
set -e
set -u
set -o pipefail

# this script will accept one folder or multiple fastq.gz files as inputs

#check arguments
if [[ $# -eq 1 ]]
then
        if [[ -d "$1" ]]
	then
		echo "FASTQ will be performed for files in the directory: "$1
		FILES=($(find $1 -maxdepth 1 -name "*fastq.gz" | sort))
		if [[ ${#FILES[@]} -eq 0 ]]
		then
			echo "No files fastq.gz files present in the directory."
			exit 1
		fi
	elif [[ -f "$1" ]] && [[ $1 =~ 'fastq.gz' ]]
	then
		FILES=($@)
	else
		echo "Invalid input"
	fi
elif [[ $# -gt 1 ]]
then
	FILES=($@)
	for FILE in ${FILES[@]}
	do
		if [[ ! -f $FILE ]] || [[ ! $FILE =~ 'fastq.gz' ]]
		then
			echo $FILE" is not a valid entry."
			echo "Supply a single directory or one or more files."
			exit 1
		fi
	done
elif [[ $# -eq 0 ]]
then
        echo "Supply a single directory or one or more files."
        exit 1
fi

echo "FASTQ will be performed for the following files:"
printf '%s\n' "${FILES[@]}"

#---------------------------------------------------------#
echo "Enter name of folder to save FASTQC analysis to:"
read FAST_MULTIQC_REPORT_NAME
#---------------------------------------------------------#
# setting up directories
SCRIPT_CUR_DIR="$(cd -P "$(dirname "${BASH_SOURCE[0]}")" 2>&1 >/dev/null && pwd)/"

OUT_REL_DIR="../results/fastqc/${FAST_MULTIQC_REPORT_NAME}/"
if [[ ! -d "${OUT_REL_DIR}" ]]
then
	mkdir -p "${OUT_REL_DIR}"
fi
OUT_DIR="$(cd "$OUT_REL_DIR" 2>&1 >/dev/null && pwd)/"
echo ${OUT_DIR}
LOG_FILE="${OUT_DIR}/Logfile.log"

#---------------------------------------------------------#
module load fastqc/0.11.9

echo $(date -Iseconds)" SPG_RNAseq_pipeline: 0. Starting FASTQC.." | tee "${LOG_FILE}" 
fastqc "${FILES[@]}" --outdir "${OUT_DIR}" --threads 8 2>&1 | tee -a "${LOG_FILE}"
echo $(date -Iseconds)" SPG_RNAseq_pipeline: 0. Done FASTQC." | tee -a "${LOG_FILE}"

#---------------------------------------------------------#
echo $(date -Iseconds)" SPG_RNAseq_pipeline: 0. Compiling FASTQC results.."
module load python/3.8.10
ENVDIR="$TMPDIR/$(whoami)/$RANDOM/env"
virtualenv --no-download $ENVDIR
set +u
source $ENVDIR/bin/activate
set -u
pip install --no-index --upgrade pip
pip install --no-index -r ${SCRIPT_CUR_DIR}/multiqc_requirements.txt

multiqc ${OUT_DIR} -o ${OUT_DIR}/ -n "${OUT_DIR}/multiQC_report"

deactivate
rm $ENVDIR -r

