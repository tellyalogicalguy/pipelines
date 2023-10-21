#!/bin/bash
#SBATCH --account=rpp-jfcote11
#SBATCH --mail-user=poorani.subramani@mail.mcgill.ca
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=8000
#SBATCH --time=2:00:0
#SBATCH --tmp=100G
#SBATCH --output="./slurm_outputs/%x_slurm-%j.out"

set -e
set -u
set -o pipefail

SLRM_LOCAL_TMP_DIR=$SLURM_TMPDIR/$RANDOM/${SAMPLE_ID}/
mkdir -p ${SLRM_LOCAL_TMP_DIR}

#---------------------------------------------------------#
FILT_FASTQ_R1=${FILT_FASTQ_DIR}/${SAMPLE_ID}_R1.filteredFastq.gz
FILT_FASTQ_R2=${FILT_FASTQ_DIR}/${SAMPLE_ID}_R2.filteredFastq.gz

echo $(date -Iseconds)" SPG_RNAseq_pipeline: 3. Starting second pass with STAR"

module load star/2.7.9a

SAMPLE_IDS=($(awk 'NR>1 {print $1}' ${DATA_INFO_FILE} | uniq))
SJDB_FILENAMES=(${SAMPLE_IDS[@]/%/_SJ.out.tab})
SJDB_FILES_LOC=${SJDB_FILENAMES[@]/#/$ALIGNED_BAM_DIR}

STAR --runThreadN "$SLURM_CPUS_PER_TASK" \
	--genomeDir "$GENOME_DIR" \
	--readFilesIn $FILT_FASTQ_R1 $FILT_FASTQ_R2 \
	--readFilesCommand gunzip -c \
	--outSAMtype BAM SortedByCoordinate \
	--outFileNamePrefix "${ALIGNED_BAM_DIR}/${SAMPLE_ID}_2pass_" \
	--sjdbFileChrStartEnd ${SJDB_FILES_LOC}

ALIGNED_BAM="${ALIGNED_BAM_DIR}/${SAMPLE_ID}_2pass_Aligned.sortedByCoord.out.bam"
ALIGNED_BAM_BAI="${ALIGNED_BAM_DIR}/${SAMPLE_ID}_2pass_Aligned.sortedByCoord.out.bam.bai"

echo $(date -Iseconds)" SPG_RNAseq_pipeline: 3. Indexing BAM"
module load samtools/1.16.1
samtools index -@ "$SLURM_CPUS_PER_TASK" -o $ALIGNED_BAM_BAI $ALIGNED_BAM 
echo $(date -Iseconds)" SPG_RNAseq_pipeline: 3. Done STAR alignment and indexing."
#---------------------------------------------------------#
module load python/3.8.10
virtualenv --no-download $SLURM_TMPDIR/env
set +u
source $SLURM_TMPDIR/env/bin/activate
set -u
pip install --no-index --upgrade pip
pip install --no-index -r ${SCRIPT_CUR_DIR}/deepTools_requirements.txt

NORM_FACTORS_FILE="${DATA_INFO_DIR}/normFactor_file.txt"

echo $(date -Iseconds)" SPG_RNAseq_pipeline: 4. Making individual bigwig files.."

if [[ $SEQ_STRANDED == "Yes" ]]
then

	if [[ $SPIKEIN_GENOME != "none" ]]
	then
		NORM_FACTOR=$(awk -v S_ID=$SAMPLE_ID '$1==S_ID {print $2}' $NORM_FACTORS_FILE)
		echo "with scalingFactor: "$NORM_FACTOR
		bamCoverage -b $ALIGNED_BAM \
		--filterRNAstrand forward -o ${BW_DIR}/${SAMPLE_ID}_fwd.bw \
		-of "bigwig" --scaleFactor ${NORM_FACTOR} --binSize 1 \
		-p "max" --normalizeUsing "None" --exactScaling

		bamCoverage -b $ALIGNED_BAM \
		--filterRNAstrand reverse -o ${BW_DIR}/${SAMPLE_ID}_rev.bw \
		-of "bigwig" --scaleFactor ${NORM_FACTOR} --binSize 1 \
		-p "max" --normalizeUsing "None" --exactScaling

	elif [[ $SPIKEIN_GENOME == "none" ]]
	then
		echo "with lib size normalization (RPGC)"
		bamCoverage -b $ALIGNED_BAM \
		--filterRNAstrand forward -o ${BW_DIR}/${SAMPLE_ID}_fwd.bw \
		-of "bigwig" --binSize 1 -p "max" --exactScaling \
		--normalizeUsing "RPGC" --effectiveGenomeSize 2652783500

		bamCoverage -b $ALIGNED_BAM \
		--filterRNAstrand reverse -o ${BW_DIR}/${SAMPLE_ID}_rev.bw \
		-of "bigwig" --binSize 1 -p "max" --exactScaling \
		--normalizeUsing "RPGC" --effectiveGenomeSize 2652783500
	fi

elif [[ $SEQ_STRANDED == "No" ]]
then

	if [[ $SPIKEIN_GENOME != "none" ]]
	then
		NORM_FACTOR=$(awk -v S_ID=$SAMPLE_ID '$1==S_ID {print $2}' $NORM_FACTORS_FILE)
		echo "with scalingFactor: "$NORM_FACTOR
		bamCoverage -b $ALIGNED_BAM \
		-o ${BW_DIR}/${SAMPLE_ID}.bw \
		-of "bigwig" --scaleFactor ${NORM_FACTOR} --binSize 1 \
		-p "max" --normalizeUsing "None" --exactScaling

	elif [[ $SPIKEIN_GENOME == "none" ]]
	then
		echo "with lib size normalization (RPGC)"
		bamCoverage -b $ALIGNED_BAM \
		-o ${BW_DIR}/${SAMPLE_ID}.bw \
		-of "bigwig" --binSize 1 -p "max" --exactScaling \
		--normalizeUsing "RPGC" --effectiveGenomeSize 2652783500
	fi

fi

echo $(date -Iseconds)" SPG_RNAseq_pipeline: 4. Done making individual bigwig files.."
#---------------------------------------------------------#
#making symlinks
for SYML_DIR in FILT_FASTQ_SYML_DIR ALIGNED_BAM_SYML_DIR BW_SYML_DIR
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
# update progress file and start the next step in the pipeline if appropriate
echo $(date -Iseconds)" SPG_RNAseq_pipeline: Updating BAM progress file"
gawk -i inplace -v S_ID="${SAMPLE_ID}" 'BEGIN{OFS="\t"} $1==S_ID {$2="Done"}'1 ${BAM_PROG_2P_FILE}

cat ${BAM_PROG_2P_FILE}

while IFS=$'\t' read -r s_id status; do
        if [[ $status != "Done" ]]
        then
                echo "Waiting for "$s_id
                exit 0
        fi
done < ${BAM_PROG_2P_FILE}


echo $(date -Iseconds)" SPG_RNAseq_pipeline: Starting featureCounts, DESeq2, rMATS and bamCoverage.."
${SCRIPT_CUR_DIR}/56_fcDESeq_bwQC.sh


