#!/bin/bash
#SBATCH --account=rpp-jfcote11
#SBATCH --mail-user=poorani.subramani@mail.mcgill.ca
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=4000
#SBATCH --tmp=100G
#SBATCH --time=1:00:0
#SBATCH --output="./slurm_outputs/%x_slurm-%j.out"

set -e
set -u
set -o pipefail


module load star/2.7.9a

STAR --runMode genomeGenerate \
--runThreadN $SLURM_CPUS_PER_TASK \
--genomeDir "$OUT_DIR" \
--outTmpDir "${OUT_DIR}/tmp/" \
--genomeFastaFiles "$FASTA_LOC" \
--sjdbGTFfile "$GTF_LOC" \
--sjdbOverhang "$SJDB_OVHG_VAL"


#create symlinks to the indexed files in the scratch-based folder; create folder with the same subdirectory structure if the folder doesn't already exist.

if [ ! -d $SYMLINK_DIR ]; then
    mkdir -p "$SYMLINK_DIR"
fi

ln -s ${OUT_DIR}* "$SYMLINK_DIR"

