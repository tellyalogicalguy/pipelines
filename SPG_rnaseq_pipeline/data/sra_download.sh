#!/bin/bash
set -e
set -u
set -o pipefail

#Setting up directories
PROJECT="/home/subrampg/"
PROJECTPARENT_DIR="${PROJECT}binf_analyses/"
SCRATCHPARENT_DIR="/lustre07/scratch/subrampg/binf_analyses_data/"
CUR_DIR="$(cd -P "$(dirname "${BASH_SOURCE[0]}")" 2>&1 >/dev/null && pwd)/"

OUT_REL_DIR="./raw_fastq/"
OUT_SYML_DIR="$(cd "$OUT_REL_DIR" 2>&1 >/dev/null && pwd)/"
OUT_DIR="$(sed "s:$PROJECTPARENT_DIR:$SCRATCHPARENT_DIR:" <(echo $OUT_SYML_DIR))"

SRA_FILE="./sra_list.txt"
SRA_LIST=($(cat $SRA_FILE))

cd ${OUT_DIR}
echo "In folder: "$(pwd)
#exit 0

module load gcc/9.3.0
module load sra-toolkit/3.0.0

echo "Downloading fastq files to :"$OUT_DIR

for SRA in ${SRA_LIST[@]}
do
	echo "Downloading: "$SRA
	fasterq-dump -p -t ${OUT_DIR} --outdir ${OUT_DIR} $SRA
done

# compress fastq files for STAR consistency
echo "gzipping .fastq files.."
find ${OUT_DIR} -name "*.fastq" | xargs -I "{}" -n 1 -P 8 gzip {}

#---------------------------------------------------------#
#making symlinks
for SYML_DIR in OUT_SYML_DIR
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
