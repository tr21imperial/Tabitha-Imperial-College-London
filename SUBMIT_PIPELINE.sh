#SUBMIT_PIPELINE.sh


#!/bin/bash

set -euo pipefail

# ============================================================
# Ribo-seq Two-step pipeline master submission script
# ============================================================
# Stage 1: upstream processing
#   - FastQC
#   - Trim Galore paired-end trimming
#   - Bowtie2 contaminant removal
#   - STAR genome alignment
#   - Samtools stats
#
# Stage 2: Ribo-seq analysis
#   - Ribo-seQC
#   - RiboTaper
#   - BEDTools / MetamORF counting
#   - QC summary
#   - MultiQC
#
# Stage 2 jobs start only if corresponding Stage 1 job finishes
# successfully using PBS dependency: -W depend=afterok
# ============================================================

# ============================================================
# SECTION 1: SAMPLE INFORMATION
# ============================================================

SAMPLE_1_NAME="MT1155"
SAMPLE_2_NAME="MT1156"
SAMPLE_3_NAME="MT1157"

SAMPLE_1_R1="MT1155_R1"
SAMPLE_2_R1="MT1156_R1"
SAMPLE_3_R1="MT1157_R1"

SAMPLE_1_R2="MT1155_R2"
SAMPLE_2_R2="MT1156_R2"
SAMPLE_3_R2="MT1157_R2"

# ============================================================
# SECTION 2: PROJECT PATHS
# ============================================================

PROJECT=/rds/general/project/bbsrc_sva/live/Tabitha
USER_DIR=/rds/general/project/bbsrc_sva/live
FASTQ_DIR=/rds/general/project/bbsrc_sva/live/fastq
LOG_DIR=$PROJECT/logfiles

STAGE1_SCRIPT=$PROJECT/scripts/STAGE1_UPSTREAM_PROCESSING.sh
STAGE2_SCRIPT=$PROJECT/scripts/STAGE2_RIBOSEQ_ANALYSIS.sh

GENOME_DIR=/rds/general/project/bbsrc_sva/live/hg38
GENOME_FASTA=/rds/general/project/bbsrc_sva/live/hg38/hg38.fa
ANNOTATION_GTF=/rds/general/project/bbsrc_sva/live/annotations/gencode.v47.annotation.gtf
CONTAMINANT_IDX=$PROJECT/indexes/contaminants/contaminants
METAMORF_BED=$PROJECT/MetamORF/MetamORF_library.bed

# ============================================================
# SECTION 3: AUTOMATIC SETUP
# ============================================================

#mkdir -p "$LOG_DIR"
#mkdir -p "$PROJECT/aln"
#mkdir -p "$PROJECT/stats"
#mkdir -p "$PROJECT/QC/raw"
#mkdir -p "$PROJECT/QC/trimmed"
#mkdir -p "$PROJECT/QC/aligned"
#mkdir -p "$PROJECT/QC/summaries"
#mkdir -p "$PROJECT/QC/multiqc"
#mkdir -p "$PROJECT/riboseqc"
#mkdir -p "$PROJECT/psites"
#mkdir -p "$PROJECT/ribotaper"
#mkdir -p "$PROJECT/ribotaper_annot"
#mkdir -p "$PROJECT/bedtools"
#mkdir -p "$PROJECT/MetamORF"
#mkdir -p "$PROJECT/indexes/contaminants"
#mkdir -p "$PROJECT/scripts"


ERRORS=0
check_file() {
    if [ ! -f "$1" ]; then
        echo "ERROR: File not found: $1"
        ERRORS=1
    else
        echo "OK: $1"
    fi
}

check_dir() {
    if [ ! -d "$1" ]; then
        echo "ERROR: Directory not found: $1"
        ERRORS=1
    else
        echo "OK: $1"
    fi
}

# Scripts
check_file "$STAGE1_SCRIPT"
check_file "$STAGE2_SCRIPT"

# FASTQ files
check_file "$FASTQ_DIR/${SAMPLE_1_R1}.fastq.gz"
check_file "$FASTQ_DIR/${SAMPLE_1_R2}.fastq.gz"
check_file "$FASTQ_DIR/${SAMPLE_2_R1}.fastq.gz"
check_file "$FASTQ_DIR/${SAMPLE_2_R2}.fastq.gz"
check_file "$FASTQ_DIR/${SAMPLE_3_R1}.fastq.gz"
check_file "$FASTQ_DIR/${SAMPLE_3_R2}.fastq.gz"

# Reference files
check_dir "$GENOME_DIR"
check_file "$GENOME_FASTA"
check_file "$ANNOTATION_GTF"
check_file "$METAMORF_BED"

# Bowtie2 contaminant index files
check_file "${CONTAMINANT_IDX}.1.bt2"
check_file "${CONTAMINANT_IDX}.2.bt2"
check_file "${CONTAMINANT_IDX}.3.bt2"
check_file "${CONTAMINANT_IDX}.4.bt2"
check_file "${CONTAMINANT_IDX}.rev.1.bt2"
check_file "${CONTAMINANT_IDX}.rev.2.bt2"

if [ "$ERRORS" -ne 0 ]; then
    echo
    echo "ERROR: Validation failed. Fix missing files before submitting."
    exit 1
fi

echo
echo "=================================================="
echo "Ribo-seq Pipeline — Two Stage Submission"
echo "Imperial College London HPC cx3"
echo "Project: bbsrc_sva / User: ldeelen"
echo "Date: $(date)"
echo "=================================================="
echo
echo "Stage 1: Upstream processing (parallel)"
echo "Stage 2: Ribo-seq analysis (after Stage 1)"
echo

# ============================================================
# SECTION 4: SUBMIT STAGE 1
# ============================================================

echo "--- Submitting Stage 1 jobs ---"
echo

STAGE1_JOB1=$(qsub \
    -N "STAGE1_${SAMPLE_1_NAME}" \
    -e "$LOG_DIR/stage1_${SAMPLE_1_NAME}_error.log" \
    -o "$LOG_DIR/stage1_${SAMPLE_1_NAME}_output.log" \
    -v ARG1="${SAMPLE_1_R1}",ARG2="${SAMPLE_1_R2}",ARG3="${SAMPLE_1_NAME}" \
    "$STAGE1_SCRIPT"
)
echo "Stage 1 submitted — ${SAMPLE_1_NAME}: $STAGE1_JOB1"

STAGE1_JOB2=$(qsub \
    -N "STAGE1_${SAMPLE_2_NAME}" \
    -e "$LOG_DIR/stage1_${SAMPLE_2_NAME}_error.log" \
    -o "$LOG_DIR/stage1_${SAMPLE_2_NAME}_output.log" \
    -v ARG1="${SAMPLE_2_R1}",ARG2="${SAMPLE_2_R2}",ARG3="${SAMPLE_2_NAME}" \
    "$STAGE1_SCRIPT"
)
echo "Stage 1 submitted — ${SAMPLE_2_NAME}: $STAGE1_JOB2"

STAGE1_JOB3=$(qsub \
    -N "STAGE1_${SAMPLE_3_NAME}" \
    -e "$LOG_DIR/stage1_${SAMPLE_3_NAME}_error.log" \
    -o "$LOG_DIR/stage1_${SAMPLE_3_NAME}_output.log" \
    -v ARG1="${SAMPLE_3_R1}",ARG2="${SAMPLE_3_R2}",ARG3="${SAMPLE_3_NAME}" \
    "$STAGE1_SCRIPT"
)
echo "Stage 1 submitted — ${SAMPLE_3_NAME}: $STAGE1_JOB3"


S1_ID1="$STAGE1_JOB1"
S1_ID2="$STAGE1_JOB2"
S1_ID3="$STAGE1_JOB3"

# ============================================================
# SECTION 5: SUBMIT STAGE 2
# ============================================================

echo
echo "--- Submitting Stage 2 jobs (held until Stage 1 completes) ---"
echo

STAGE2_JOB1=$(qsub \
    -N "STAGE2_${SAMPLE_1_NAME}" \
    -e "$LOG_DIR/stage2_${SAMPLE_1_NAME}_error.log" \
    -o "$LOG_DIR/stage2_${SAMPLE_1_NAME}_output.log" \
    -W "depend=afterok:${S1_ID1}" \
    -v ARG1="${SAMPLE_1_R1}",ARG2="${SAMPLE_1_R2}",ARG3="${SAMPLE_1_NAME}" \
    "$STAGE2_SCRIPT"
)
echo "Stage 2 submitted — ${SAMPLE_1_NAME}: $STAGE2_JOB1"
echo "Waiting for Stage 1 job: $STAGE1_JOB1"

STAGE2_JOB2=$(qsub \
    -N "STAGE2_${SAMPLE_2_NAME}" \
    -e "$LOG_DIR/stage2_${SAMPLE_2_NAME}_error.log" \
    -o "$LOG_DIR/stage2_${SAMPLE_2_NAME}_output.log" \
    -W "depend=afterok:${S1_ID2}" \
    -v ARG1="${SAMPLE_2_R1}",ARG2="${SAMPLE_2_R2}",ARG3="${SAMPLE_2_NAME}" \
    "$STAGE2_SCRIPT"
)
echo "Stage 2 submitted — ${SAMPLE_2_NAME}: $STAGE2_JOB2"
echo "Waiting for Stage 1 job: $STAGE1_JOB2"

STAGE2_JOB3=$(qsub \
    -N "STAGE2_${SAMPLE_3_NAME}" \
    -e "$LOG_DIR/stage2_${SAMPLE_3_NAME}_error.log" \
    -o "$LOG_DIR/stage2_${SAMPLE_3_NAME}_output.log" \
    -W "depend=afterok:${S1_ID3}" \
    -v ARG1="${SAMPLE_3_R1}",ARG2="${SAMPLE_3_R2}",ARG3="${SAMPLE_3_NAME}" \
    "$STAGE2_SCRIPT"
)
echo "Stage 2 submitted — ${SAMPLE_3_NAME}: $STAGE2_JOB3"
echo "Waiting for Stage 1 job: $STAGE1_JOB3"

# ============================================================
# SECTION 6: SUMMARY
# ============================================================

echo
echo "=================================================="
echo "All jobs submitted successfully"
echo "=================================================="
echo
echo "STAGE 1 — Upstream Processing:"
echo "  ${SAMPLE_1_NAME}: $STAGE1_JOB1"
echo "  ${SAMPLE_2_NAME}: $STAGE1_JOB2"
echo "  ${SAMPLE_3_NAME}: $STAGE1_JOB3"
echo
echo "STAGE 2 — Ribo-seq Analysis:"
echo "  ${SAMPLE_1_NAME}: $STAGE2_JOB1 (held)"
echo "  ${SAMPLE_2_NAME}: $STAGE2_JOB2 (held)"
echo "  ${SAMPLE_3_NAME}: $STAGE2_JOB3 (held)"
echo
echo "Monitor jobs:"
echo "  qstat -u ldeelen"
echo
echo "Check one job in detail:"
echo "  qstat -f <job_id>"
echo
echo "QC outputs saved to:"
echo "  $PROJECT/QC"
echo

SUBMISSION_LOG="$LOG_DIR/submission_record_$(date +%Y%m%d_%H%M%S).txt"

cat > "$SUBMISSION_LOG" << EOF
Ribo-seq Pipeline — Two Stage Submission Record
Date: $(date)
User: ldeelen
HPC: Imperial College cx3
Project: bbsrc_sva

Directory structure:
  Project:          $PROJECT
  User dir:         $USER_DIR
  FASTQ dir:        $FASTQ_DIR
  Genome FASTA:     $GENOME_FASTA
  Genome dir:       $GENOME_DIR
  Annotation GTF:   $ANNOTATION_GTF
  Contaminant idx:  $CONTAMINANT_IDX
  MetamORF BED:     $METAMORF_BED
  QC dir:           $PROJECT/QC

Scripts:
  Stage 1:          $STAGE1_SCRIPT
  Stage 2:          $STAGE2_SCRIPT

Samples:
  ${SAMPLE_1_NAME}: R1=${SAMPLE_1_R1} R2=${SAMPLE_1_R2}
  ${SAMPLE_2_NAME}: R1=${SAMPLE_2_R1} R2=${SAMPLE_2_R2}
  ${SAMPLE_3_NAME}: R1=${SAMPLE_3_R1} R2=${SAMPLE_3_R2}

Stage 1 job IDs:
  ${SAMPLE_1_NAME}: $STAGE1_JOB1
  ${SAMPLE_2_NAME}: $STAGE1_JOB2
  ${SAMPLE_3_NAME}: $STAGE1_JOB3

Stage 2 job IDs:
  ${SAMPLE_1_NAME}: $STAGE2_JOB1
  ${SAMPLE_2_NAME}: $STAGE2_JOB2
  ${SAMPLE_3_NAME}: $STAGE2_JOB3
EOF

echo "Submission record saved: $SUBMISSION_LOG"
