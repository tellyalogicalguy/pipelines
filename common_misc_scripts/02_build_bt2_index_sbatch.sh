#!/bin/bash
#SBATCH --account=rpp-jfcote11
#SBATCH --mail-user=poorani.subramani@mail.mcgill.ca
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=8000
#SBATCH --time=0:50:0
#SBATCH --output="%x_slurm-%j.out"

set -e
set -u
set -o pipefail


module load bowtie2/2.4.4

bowtie2-build --threads 16 ${GENOME_LOC} ${OUT_DIR}/customCH12


#create symlinks to the indexed files in the scratch-based folder; create folder with the same subdirectory structure if the folder doesn't already exist.

if [ ! -d $SYMLINK_DIR ]; then
    mkdir -p "$SYMLINK_DIR"
fi

ln -s ${OUT_DIR}* "$SYMLINK_DIR"

