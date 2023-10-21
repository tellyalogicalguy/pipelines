#!/bin/bash
set -e
set -u
set -o pipefail

#Setting up directories
PROJECT="/home/subrampg/"
PROJECTPARENT_DIR="${PROJECT}binf_analyses/"
SCRATCHPARENT_DIR="/lustre07/scratch/subrampg/binf_analyses_data/"
CUR_DIR="$(cd -P "$(dirname "${BASH_SOURCE[0]}")" 2>&1 >/dev/null && pwd)"

#genome dir (STAR indexed directory)
export GENOME_LOC="$SCRATCHPARENT_DIR/customCH12genome/reformed_igk/GRCm38.primary_assembly.genome_reformed_reformed.fa"

#output directory
OUT_DIR_REL="${PROJECTPARENT_DIR}/common/bt2_indexes/customCH12genome/" #make the desired directory structure here
if [ ! -d $OUT_DIR_REL ]; then
        mkdir -p "$OUT_DIR_REL"
fi
OUT_DIR_ABS="$(cd "$OUT_DIR_REL" 2>&1 >/dev/null && pwd)/"
OUT_SUBDIRS="$(sed "s:$PROJECTPARENT_DIR::" <(echo $OUT_DIR_ABS))/"
export OUT_DIR="${SCRATCHPARENT_DIR}${OUT_SUBDIRS}"
export SYMLINK_DIR="${PROJECTPARENT_DIR}${OUT_SUBDIRS}"

echo "out-dir: "$OUT_DIR
if [ ! -d $OUT_DIR ]; then
        mkdir -p ${OUT_DIR}
fi


SLRM_NAME="02_build_bt2_index"
sbatch --job-name "$SLRM_NAME" 02_build_bt2_index_sbatch.sh
sleep 0.2

