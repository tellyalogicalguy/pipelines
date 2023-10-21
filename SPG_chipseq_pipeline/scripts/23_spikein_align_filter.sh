#!/bin/bash
set -e
set -u
set -o pipefail

#Setting up directories
PROJECT="/home/subrampg/"
PROJECTPARENT_DIR="${PROJECT}binf_analyses/"
SCRATCHPARENT_DIR="/lustre07/scratch/subrampg/binf_analyses_data/"
CUR_DIR="$(cd -P "$(dirname "${BASH_SOURCE[0]}")" 2>&1 >/dev/null && pwd)/"

GENOME_DIR="${SCRATCHPARENT_DIR}/common/bt2_indexes/${SPIKEIN_GENOME}/" #bt2 index loc

FILT_FASTQ_R1="${FILT_FASTQ_DIR}/${SAMPLE_ID}_R1.filteredFastq.gz"
FILT_FASTQ_R2="${FILT_FASTQ_DIR}/${SAMPLE_ID}_R2.filteredFastq.gz"

SLRM_LOCAL_TMP_DIR=$SLURM_TMPDIR/$RANDOM/${SAMPLE_ID}/
mkdir -p ${SLRM_LOCAL_TMP_DIR}

#---------------------------------------------------------#
echo $(date -Iseconds)" SPG_ChIPseq_pipeline(${SPIKEIN_GENOME}): 2. Starting aligning with bowtie2"
module load bowtie2/2.4.4
module load samtools/1.16.1

UNFILT_BAM=${SLRM_LOCAL_TMP_DIR}${SAMPLE_ID}.aligned.${SPIKEIN_GENOME}.bam
UNFILT_SORT_BAM=${SLRM_LOCAL_TMP_DIR}${SAMPLE_ID}.unfilt.${SPIKEIN_GENOME}.bam
BT2_LOG=${SLRM_LOCAL_TMP_DIR}${SAMPLE_ID}_bowtie2.${SPIKEIN_GENOME}.log

bowtie2 --mm -x $GENOME_DIR/${SPIKEIN_GENOME} --threads $SLURM_CPUS_PER_TASK \
--rg-id "${SAMPLE_ID}" \
--rg "ID:${SAMPLE_ID}" --rg "PL:ILLUMINA" --rg "LB:CH12" --rg "PU:unit1" --rg "SM:${SAMPLE_ID}" \
-1 <(zcat -f ${FILT_FASTQ_R1}) -2 <(zcat -f ${FILT_FASTQ_R2}) \
2> $BT2_LOG | samtools view -h -b -q 1 > $UNFILT_BAM

cp $BT2_LOG ${FILT_BAM_DIR}${SAMPLE_ID}_bowtie2.${SPIKEIN_GENOME}.log
echo $(date -Iseconds)" SPG_ChIPseq_pipeline(${SPIKEIN_GENOME}): 2. Sorting aligned BAM.."

#samtools view -h -q 1 -b -@ $SLURM_CPUS_PER_TASK $UNFILT_BAM | \
samtools sort -T ${SLRM_LOCAL_TMP_DIR} -o $UNFILT_SORT_BAM \
-@ $SLURM_CPUS_PER_TASK -m 2G -n $UNFILT_BAM

rm $UNFILT_BAM
echo $(date -Iseconds)" SPG_ChIPseq_pipeline(${SPIKEIN_GENOME}): 2. Done aligning."
#---------------------------------------------------------#
echo $(date -Iseconds)" SPG_ChIPseq_pipeline(${SPIKEIN_GENOME}): 3. Starting MarkduplicatesSpark.."
module load gatk/4.2.4.0
module load java/1.8.0_192

MARKDUP_BAM=${SLRM_LOCAL_TMP_DIR}${SAMPLE_ID}.markdup.${SPIKEIN_GENOME}.bam
MARKDUP_METRICS=${SLRM_LOCAL_TMP_DIR}${SAMPLE_ID}.dup.${SPIKEIN_GENOME}.stats.qc

gatk --java-options '-Xmx30G' MarkDuplicatesSpark -I $UNFILT_SORT_BAM \
-O $MARKDUP_BAM -M $MARKDUP_METRICS --read-validation-stringency LENIENT \
-- --spark-runner LOCAL --spark-master local[$SLURM_CPUS_PER_TASK] \
--conf $(echo "'spark.local.dir=${SLRM_LOCAL_TMP_DIR}'")

rm $UNFILT_SORT_BAM
rm $MARKDUP_METRICS
echo $(date -Iseconds)" SPG_ChIPseq_pipeline(${SPIKEIN_GENOME}): 3. Done Markduplicates."

echo $(date -Iseconds)" SPG_ChIPseq_pipeline(${SPIKEIN_GENOME}): 3. Starting quality filtering and removing duplicates.."
module load samtools/1.16.1

FILTERED_BAM=${SLRM_LOCAL_TMP_DIR}${SAMPLE_ID}.${SPIKEIN_GENOME}.bam

samtools view -h -F 1804 -q 1 -b $MARKDUP_BAM -@ $SLURM_CPUS_PER_TASK | \
samtools sort - -T ${SLRM_LOCAL_TMP_DIR} -o $FILTERED_BAM -@ $SLURM_CPUS_PER_TASK -m 2G

cp $FILTERED_BAM ${FILT_BAM_DIR}/
echo $(date -Iseconds)" SPG_ChIPseq_pipeline(${SPIKEIN_GENOME}): 3. Done filtering. Output is filtered coord-sorted BAM."

echo $(date -Iseconds)" SPG_ChIPseq_pipeline(${SPIKEIN_GENOME}): 3. Starting indexing.."
samtools index ${FILT_BAM_DIR}${SAMPLE_ID}.${SPIKEIN_GENOME}.bam \
${FILT_BAM_DIR}${SAMPLE_ID}.${SPIKEIN_GENOME}.bam.bai -@ $SLURM_CPUS_PER_TASK
echo $(date -Iseconds)" SPG_ChIPseq_pipeline(${SPIKEIN_GENOME}): 3. Done indexing."

echo $(date -Iseconds)" SPG_ChIPseq_pipeline(${SPIKEIN_GENOME}): 3. Starting stats."
samtools stats $FILTERED_BAM > ${FILT_BAM_DIR}${SAMPLE_ID}.${SPIKEIN_GENOME}.filtBAMstats.txt \
-@ $SLURM_CPUS_PER_TASK
echo $(date -Iseconds)" SPG_ChIPseq_pipeline(${SPIKEIN_GENOME}): 3. Done with stats."
#---------------------------------------------------------#

