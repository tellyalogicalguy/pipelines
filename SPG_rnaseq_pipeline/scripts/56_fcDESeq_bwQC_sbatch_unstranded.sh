#!/bin/bash
#SBATCH --account=rpp-jfcote11
#SBATCH --mail-user=poorani.subramani@mail.mcgill.ca
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=4000
#SBATCH --time=3:00:0
#SBATCH --tmp=100G
#SBATCH --output="./slurm_outputs/%x_slurm-%j.out"

set -e
set -u
set -o pipefail

SLRM_LOCAL_TMP_DIR=$SLURM_TMPDIR/$RANDOM/
mkdir -p ${SLRM_LOCAL_TMP_DIR}

#---------------------------------------------------------#
module load gcc/9.3.0
module load r/4.2.1

echo $(date -Iseconds)" SPG_RNAseq_pipeline: 5. Starting featureCounts and DESeq2.."
Rscript ${SCRIPT_CUR_DIR}/featureCounts_DESeq.R
echo $(date -Iseconds)" SPG_RNAseq_pipeline: 5. Done with featureCounts and DESeq2."
#---------------------------------------------------------#
# setting up for bigwig merging
module load gcc/9.3.0
module load htslib/1.16
module load gsl/2.6
module load kentutils/401
module load python/3.8.10
# adding paths to libBwigWig
export CPATH=$CPATH:/home/subrampg/local_installs/usr_local_substitute/libBigWig/include/
export LIBRARY_PATH=$LIBRARY_PATH:/home/subrampg/local_installs/usr_local_substitute/libBigWig/lib/
export LD_LIBRARY_PATH=/home/subrampg/local_installs/usr_local_substitute/libBigWig/lib/
# adding wiggletools to $PATH
export PATH=$PATH:/home/subrampg/local_installs/WiggleTools/bin
#CHROM_SIZE_FILE=/lustre07/scratch/subrampg/binf_analyses_data/common/mm10/GRCm38.primary_assembly.genome.fa.fai
#--------------------------------#
#Merge the same genotypes
#merging WT bw files
#REFR_SAMPLES=($(awk '$4=="WT" {print $1}' $DATA_INFO_FILE | uniq))
REFR_SAMPLES=($(awk -v groupName="$REFR_GROUP_NAME" '$4==groupName {print $1}' $DATA_INFO_FILE | sort -u ))

if [ ${#REFR_SAMPLES[@]} -eq 1 ]
then
	echo "There is only 1 WT sample."
	REFR_MERGED_BW="${BW_DIR}/${REFR_SAMPLES[0]}.bw"
else
	for i in ${!REFR_SAMPLES[@]}; do 
		REFR_BWS[$i]="${BW_DIR}/${REFR_SAMPLES[$i]}.bw"
	done

	echo $(date -Iseconds)" SPG_RNAseq_pipeline: 6. Merging ${#REFR_SAMPLES[@]} WT bigwig files.."
	REFR_MERGED_BW="${BW_DIR}/${REFR_GROUP_NAME}_merged_unstranded.bw"
	wiggletools mean ${REFR_BWS[@]} | wigToBigWig stdin $CHROM_SIZE_FILE ${REFR_MERGED_BW}
fi

# merging KO bw files
#TEST_SAMPLES=($(awk '$4=="KO" {print $1}' $DATA_INFO_FILE | uniq))
TEST_SAMPLES=($(awk -v groupName="$TEST_GROUP_NAME" '$4==groupName {print $1}' $DATA_INFO_FILE | sort -u ))

if [ ${#TEST_SAMPLES[@]} -eq 0 ]
then
	echo "No KO samples detected. Prcoessed only WT samples."
elif [ ${#TEST_SAMPLES[@]} -eq 1 ]
then
	echo "There is only 1 KO sample."
	TEST_MERGED_BW="${BW_DIR}/${TEST_SAMPLES[0]}.bw"
else

	for i in ${!TEST_SAMPLES[@]}; do 
		TEST_BWS[$i]="${BW_DIR}/${TEST_SAMPLES[$i]}.bw"
	done

	echo $(date -Iseconds)" SPG_RNAseq_pipeline: 6. Merging ${#TEST_SAMPLES[@]} KO bigwig files.."
	TEST_MERGED_BW="${BW_DIR}/${TEST_GROUP_NAME}_merged_unstranded.bw"
	wiggletools mean ${TEST_BWS[@]} | wigToBigWig stdin $CHROM_SIZE_FILE ${TEST_MERGED_BW}
fi

##---------------------------------------------------------#
# do KO-WT comparison only if KO sample is present - dont think this is wise in the case of RNAseq, unless if the strand infomration is lost. Do with wiggletools instead of bigwigCompare.

if [ ${#TEST_SAMPLES[@]} -ne 0 ]
then
	echo $(date -Iseconds)" SPG_RNAseq_pipeline: 6. Making KO-WT bigwig files with wiggletools.."

	if [ ${#REFR_SAMPLES[@]} -ne ${#TEST_SAMPLES[@]} ]
	then
		echo "--------------------------------------------------------------"
		echo "Warning: The number of WT and KO replicates are not the same!!"
		echo "--------------------------------------------------------------"
	fi

	#do diff
	wiggletools diff ${TEST_MERGED_BW} ${REFR_MERGED_BW} | \
	wigToBigWig stdin $CHROM_SIZE_FILE "${BW_DIR}/${TEST_GROUP_NAME}-${REFR_GROUP_NAME}_unstranded.bw"
fi

##---------------------------------------------------------#
## do KO-WT comparison only if KO sample is present - dont think this is wise in the case of RNAseq, unless if the strand infomration is lost
#module load python/3.8.10
#virtualenv --no-download $SLURM_TMPDIR/env
#set +u
#source $SLURM_TMPDIR/env/bin/activate
#set -u
#pip install --no-index --upgrade pip
#pip install --no-index -r ${SCRIPT_CUR_DIR}/deepTools_requirements.txt
#
#if [ ${#TEST_SAMPLES[@]} -ne 0 ]
#then
#	echo $(date -Iseconds)" SPG_RNAseq_pipeline: 6. Making KO-WT bigwig files.."
#
#	if [ ${#REFR_SAMPLES[@]} -ne ${#TEST_SAMPLES[@]} ]
#	then
#		echo "--------------------------------------------------------------"
#		echo "Warning: The number of WT and KO replicates are not the same!!"
#		echo "--------------------------------------------------------------"
#	fi
#
#	#do diff
#	bigwigCompare -b1 ${TEST_MERGED_BW} -b2 ${REFR_MERGED_BW} \
#	-o "${BW_DIR}/${TEST_GROUP_NAME}-${REFR_GROUP_NAME}_unstranded.bw" -of "bigwig" \
#	--scaleFactors 1:1 --operation "subtract" --binSize 1 \
#	-p $SLURM_CPUS_PER_TASK --verbose
#fi
#
#deactivate
#rm -rf $SLURM_TMPDIR/env
##---------------------------------------------------------#
echo $(date -Iseconds)" SPG_RNAseq_pipeline: 7. Compiling QC metrics.."
module load python/3.8.10
virtualenv --no-download $SLURM_TMPDIR/env
set +u
source $SLURM_TMPDIR/env/bin/activate
set -u
pip install --no-index --upgrade pip
pip install --no-index -r ${SCRIPT_CUR_DIR}/multiqc_requirements.txt

multiqc ${BW_DIR}/../../ -o ${QC_DIR}/ -n "${QC_DIR}/multiQC_report"

#making symlinks
for SYML_DIR in FC_DSQ_SYML_DIR BW_SYML_DIR QC_SYML_DIR
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

#---------------------------------------------------------#
# start rMATS
sbatch --job-name="7_rMATS" ${SCRIPT_CUR_DIR}/7_rMATS_sbatch.sh

