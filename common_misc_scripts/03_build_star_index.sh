#!/bin/bash
set -e
set -u
set -o pipefail

echo "Edit paths for fasta, gtf and output directory location before beginning!!"
exit 0

#get overhang information
if [[ $# -eq 1 ]] && [[ "$1" =~ ^[0-9]+$ ]] && [[ "$1" -gt 38 ]] && [[ "$1" -lt 250 ]]
then
	export SJDB_OVHG_VAL="$1"
	echo "Making STAR index for sjdbOverhang value of "$SJDB_OVHG_VAL
else
	echo "Invalid input. Enter the sjdbOverhang value for STAR index between 39 and 249bp."
	exit 1
fi

#Setting up directories
PROJECT="/home/subrampg/"
PROJECTPARENT_DIR="${PROJECT}binf_analyses/"
SCRATCHPARENT_DIR="/lustre07/scratch/subrampg/binf_analyses_data/"
CUR_DIR="$(cd -P "$(dirname "${BASH_SOURCE[0]}")" 2>&1 >/dev/null && pwd)"

# make slurm_outputs folder
if [[ ! -d $CUR_DIR/slurm_outputs ]]
then
	mkdir -p $CUR_DIR/slurm_outputs
fi

#genome dir - make sure this is correct!
export FASTA_LOC="$SCRATCHPARENT_DIR/common/hg38/GRCh38.primary_assembly.genome.fa"
export GTF_LOC="$SCRATCHPARENT_DIR/common/hg38/gencode.v43.primary_assembly.annotation.gtf"

#output directory - make sure this is correct!
OUT_DIR_REL="${PROJECTPARENT_DIR}/common/star_indexes/hg38/sjdbOH_${SJDB_OVHG_VAL}/" #make the desired directory structure here
echo "OUT_DIR_REL = "$OUT_DIR_REL

if [ ! -d $OUT_DIR_REL ]; then
        mkdir -p "$OUT_DIR_REL"
fi
OUT_DIR_ABS="$(cd "$OUT_DIR_REL" 2>&1 >/dev/null && pwd)/"
OUT_SUBDIRS="$(sed "s:$PROJECTPARENT_DIR::" <(echo $OUT_DIR_ABS))/"
export OUT_DIR="${SCRATCHPARENT_DIR}${OUT_SUBDIRS}"
export SYMLINK_DIR="${PROJECTPARENT_DIR}${OUT_SUBDIRS}"

echo "OUT_DIR = "$OUT_DIR
if [ ! -d $OUT_DIR ]; then
        mkdir -p ${OUT_DIR}
fi


SLRM_NAME="03_build_star_index_overhang_"${SJDB_OVHG_VAL}
sbatch --job-name "$SLRM_NAME" 03_build_star_index_sbatch.sh
sleep 0.2

