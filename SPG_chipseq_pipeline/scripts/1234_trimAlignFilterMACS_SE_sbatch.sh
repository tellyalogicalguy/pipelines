#!/bin/bash
#SBATCH --account=rpp-jfcote11
#SBATCH --mail-user=poorani.subramani@mail.mcgill.ca
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=8000
#SBATCH --time=3:00:0
#SBATCH --tmp=100G
#SBATCH --output="./slurm_outputs/%x_slurm-%j.out"

set -e
set -u
set -o pipefail

SLRM_LOCAL_TMP_DIR=$SLURM_TMPDIR/$RANDOM/${SAMPLE_ID}/
mkdir -p ${SLRM_LOCAL_TMP_DIR}

#---------------------------------------------------------#
RAW_FASTQ_FILES=($(awk -v S_ID="${SAMPLE_ID}" '$1==S_ID {print $3}' ${DATA_INFO_FILE}))
declare -a RAW_FASTQ
RAW_FASTQ[0]=$(find ${RAW_FASTQ_DIR} -name "${RAW_FASTQ_FILES[0]}")

echo $(date -Iseconds)" SPG_ChIPseq_pipeline: 1. Starting FASTQ filtering with fastp"
module load fastp/0.23.1

FILT_FASTQ_R1=${SLRM_LOCAL_TMP_DIR}/${SAMPLE_ID}_R1.filteredFastq.gz

fastp --in1 ${RAW_FASTQ[0]} --out1 $FILT_FASTQ_R1 \
--html ${FILT_FASTQ_DIR}/${SAMPLE_ID}.html --json /dev/null

cp $FILT_FASTQ_R1 ${FILT_FASTQ_DIR}/

echo $(date -Iseconds)" SPG_ChIPseq_pipeline: 1. Done FASTQ filtering."
#---------------------------------------------------------#

echo $(date -Iseconds)" SPG_ChIPseq_pipeline: 2. Starting aligning with bowtie2"
module load bowtie2/2.4.4
module load samtools/1.16.1

UNFILT_BAM=${SLRM_LOCAL_TMP_DIR}${SAMPLE_ID}.aligned.bam
UNFILT_SORT_BAM=${SLRM_LOCAL_TMP_DIR}${SAMPLE_ID}.unfilt.bam
BT2_LOG=${SLRM_LOCAL_TMP_DIR}${SAMPLE_ID}_bowtie2.log

bowtie2 --mm -x $GENOME_DIR/customCH12 --threads $SLURM_CPUS_PER_TASK \
--rg-id "${SAMPLE_ID}" \
--rg "ID:${SAMPLE_ID}" --rg "PL:ILLUMINA" --rg "LB:CH12" --rg "PU:unit1" --rg "SM:${SAMPLE_ID}" \
-U <(zcat -f ${FILT_FASTQ_R1}) 2> $BT2_LOG | samtools view -h -b > $UNFILT_BAM

cp $BT2_LOG ${UNFILT_BAM_DIR}/
rm $FILT_FASTQ_R1
rm $BT2_LOG

echo $(date -Iseconds)" SPG_ChIPseq_pipeline: 2. Sorting aligned BAM.."

samtools sort -T ${SLRM_LOCAL_TMP_DIR} -o $UNFILT_SORT_BAM \
-@ $SLURM_CPUS_PER_TASK -m 2G -n $UNFILT_BAM

rm $UNFILT_BAM

echo $(date -Iseconds)" SPG_ChIPseq_pipeline: 2. Copying BAM to scratch.."
cp $UNFILT_SORT_BAM ${UNFILT_BAM_DIR}/

echo $(date -Iseconds)" SPG_ChIPseq_pipeline: 2. Done aligning. Output is unfiltered name-sorted BAM."
#---------------------------------------------------------#
echo $(date -Iseconds)" SPG_ChIPseq_pipeline: 3. Starting MarkduplicatesSpark.."
module load gatk/4.2.4.0
module load java/1.8.0_192

MARKDUP_BAM=${SLRM_LOCAL_TMP_DIR}${SAMPLE_ID}.markdup.bam
MARKDUP_METRICS=${SLRM_LOCAL_TMP_DIR}${SAMPLE_ID}.dup.stats.qc

gatk --java-options '-Xmx30G' MarkDuplicatesSpark -I $UNFILT_SORT_BAM \
-O $MARKDUP_BAM -M $MARKDUP_METRICS --read-validation-stringency LENIENT \
-- --spark-runner LOCAL --spark-master local[$SLURM_CPUS_PER_TASK] \
--conf $(echo "'spark.local.dir=${SLRM_LOCAL_TMP_DIR}'")

rm $UNFILT_SORT_BAM
cp $MARKDUP_BAM ${INTERMED_FILES_DIR}/
cp $MARKDUP_METRICS ${INTERMED_FILES_DIR}/
rm $MARKDUP_METRICS

echo $(date -Iseconds)" SPG_ChIPseq_pipeline: 3. Done Markduplicates."

echo $(date -Iseconds)" SPG_ChIPseq_pipeline: 3. Starting quality filtering and removing duplicates.."
module load samtools/1.16.1

FILTERED_BAM=${SLRM_LOCAL_TMP_DIR}${SAMPLE_ID}.bam

samtools view -h -F 1804 -q 1 -b $MARKDUP_BAM -@ $SLURM_CPUS_PER_TASK | \
samtools sort - -T ${SLRM_LOCAL_TMP_DIR} -o $FILTERED_BAM -@ $SLURM_CPUS_PER_TASK -m 2G

cp $FILTERED_BAM ${FILT_BAM_DIR}/
echo $(date -Iseconds)" SPG_ChIPseq_pipeline: 3. Done filtering. Output is filtered coord-sorted BAM."

echo $(date -Iseconds)" SPG_ChIPseq_pipeline: 3. Starting indexing.."
samtools index ${FILT_BAM_DIR}${SAMPLE_ID}.bam ${FILT_BAM_DIR}${SAMPLE_ID}.bam.bai -@ $SLURM_CPUS_PER_TASK
echo $(date -Iseconds)" SPG_ChIPseq_pipeline: 3. Done indexing."

echo $(date -Iseconds)" SPG_ChIPseq_pipeline: 3. Starting stats."
samtools stats $MARKDUP_BAM > ${INTERMED_FILES_DIR}${SAMPLE_ID}.unfilt.dupmark.stats.txt \
-@ $SLURM_CPUS_PER_TASK

samtools stats $FILTERED_BAM > ${FILT_BAM_DIR}${SAMPLE_ID}.filtBAMstats.txt \
-@ $SLURM_CPUS_PER_TASK
echo $(date -Iseconds)" SPG_ChIPseq_pipeline: 3. Done with stats."
#---------------------------------------------------------#
# call peaks only if not control
if [[ $CONTROL_OR_NOT == "no" ]]
then

echo $(date -Iseconds)" SPG_ChIPseq_pipeline: 4. Starting MACS2 peak calling.."
module load python/3.8.10
virtualenv --no-download $SLURM_TMPDIR/env
set +u
source $SLURM_TMPDIR/env/bin/activate
set -u
pip install --no-index --upgrade pip
pip install --no-index -r ${SCRIPT_CUR_DIR}/macs2_requirements.txt
chmod +x ${SLURM_TMPDIR}/env/bin/macs2

GENOMESIZE='mm'

macs2 callpeak -t ${FILTERED_BAM} -c ${CHIP_INPUT_BAM} \
-f BAM --outdir ${PEAKS_DIR} -n ${PEAKS_DIR}/${SAMPLE_ID} \
-g ${GENOMESIZE} 2> ${PEAKS_DIR}/${SAMPLE_ID}_macs2.log

FRAG_LENGTH=$(grep -Po '(?<=(predicted fragment length is ))[0-9]*(?= bps)' ${PEAKS_DIR}/${SAMPLE_ID}_macs2.log)
echo "MACS2 predicted fragment length is "$FRAG_LENGTH

echo $(date -Iseconds)" SPG_ChIPseq_pipeline: 4. Done peak calling."

fi
#---------------------------------------------------------#
# align to spikein genome if appropriate
export SLURM_CPUS_PER_TASK
if [[ "$SPIKEIN_GENOME" = "dm6" || "$SPIKEIN_GENOME" = "k12" ]]
then
echo $(date -Iseconds)" SPG_ChIPseq_pipeline: Detected "$SPIKEIN_GENOME" as spikein genome."
echo $(date -Iseconds)" SPG_ChIPseq_pipeline: Aligning to $SPIKEIN_GENOME genome.."
${SCRIPT_CUR_DIR}/23_spikein_align_filter.sh
fi

#---------------------------------------------------------#
#making symlinks
for SYML_DIR in FILT_FASTQ_SYML_DIR UNFILT_BAM_SYML_DIR INTERMED_FILES_SYML_DIR FILT_BAM_SYML_DIR PEAKS_SYML_DIR
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
echo $(date -Iseconds)" SPG_ChIPseq_pipeline: Updating BAM progress file"
gawk -i inplace -v S_ID="${SAMPLE_ID}" 'BEGIN{OFS="\t"} $1==S_ID {$2="Done"}'1 ${BAM_PROG_FILE}

cat ${BAM_PROG_FILE}

while IFS=$'\t' read -r s_id status; do 
	if [[ $status != "Done" ]]
	then
		echo "Waiting for "$s_id
		exit 0
	fi 
done < ${BAM_PROG_FILE}


echo $(date -Iseconds)" SPG_ChIPseq_pipeline: Starting DiffBind and bamCoverage.."
${SCRIPT_CUR_DIR}/567_DBbwQC.sh


