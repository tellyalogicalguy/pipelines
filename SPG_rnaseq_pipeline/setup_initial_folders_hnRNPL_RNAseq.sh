#!/bin/bash
set -e
set -u
set -o pipefail

#check number of arguments
if [[ $# -eq 1 ]]
then
	CELLTYPE_NAME=($(basename "$1"))
	# PREFIX="$(sed "s:$CELLTYPE_NAME::" <(echo $1))"
else
	echo "Incorrect number of arguments. Supply CELLTYPE_NAME"
	exit 1
fi


#Setting up directories
# in local
PROJECT="/home/subrampg/"
PROJECTPARENT_DIR="${PROJECT}binf_analyses/"
SCRATCHPARENT_DIR="/lustre07/scratch/subrampg/binf_analyses_data/"
CUR_DIR="$(pwd)/"

echo "Setting up folders in "${CUR_DIR}" (local)"

mkdir -p ${CUR_DIR}/$CELLTYPE_NAME/rna_seq/hnrnpl/data/raw_fastq/
mkdir -p ${CUR_DIR}/$CELLTYPE_NAME/rna_seq/hnrnpl/scripts/slurm_outputs/
mkdir -p ${CUR_DIR}/$CELLTYPE_NAME/rna_seq/hnrnpl/scripts/seff_outputs/
cp ~/script_templates/SPG_rnaseq_pipeline/scripts/* ${CUR_DIR}/$CELLTYPE_NAME/rna_seq/hnrnpl/scripts/
cp ~/script_templates/SPG_rnaseq_pipeline/data/* ${CUR_DIR}/$CELLTYPE_NAME/rna_seq/hnrnpl/data/



tree ${CUR_DIR}${CELLTYPE_NAME}
# in scratch
SCRATCH_EQV_OF_CUR_DIR="$(sed "s:$PROJECTPARENT_DIR:$SCRATCHPARENT_DIR:" <(echo $CUR_DIR))"
echo "Setting up folders in "${SCRATCH_EQV_OF_CUR_DIR}" (scratch)"
mkdir -p ${SCRATCH_EQV_OF_CUR_DIR}/$CELLTYPE_NAME/rna_seq/hnrnpl/data/raw_fastq/

tree ${SCRATCH_EQV_OF_CUR_DIR}${CELLTYPE_NAME}

echo "Transfer raw FASTQ files with rsync from local computer to:"
echo "subrampg@narval.computecanada.ca:${SCRATCH_EQV_OF_CUR_DIR}/$CELLTYPE_NAME/rna_seq/hnrnpl/data/raw_fastq/"
echo ""
echo "Then do commands in ./$CELLTYPE_NAME/rna_seq/hnrnpl/data/initial_raw_fastq_setup_commands.txt, i.e.:"
echo ""
cat ./$CELLTYPE_NAME/rna_seq/hnrnpl/data/initial_raw_fastq_setup_commands.txt

