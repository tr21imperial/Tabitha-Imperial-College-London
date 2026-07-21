#!/bin/bash
#PBS -l walltime=24:00:00
#PBS -l select=1:ncpus=14:mem=100gb
#PBS -N STAGE1_UPSTREAM

set -euo pipefail

# ============================================================
# CALLING VARIABLES
# ARG1 = R1 FASTQ basename (without .fastq.gz)
# ARG2 = R2 FASTQ basename (without .fastq.gz)
# ARG3 = sample name
# ============================================================

: "${ARG1:?ARG1 not set}"
: "${ARG2:?ARG2 not set}"
: "${ARG3:?ARG3 not set}"
: "${TMPDIR:?TMPDIR not set}"

echo "$ARG1   # R1 FASTQ"
echo "$ARG2   # R2 FASTQ"
echo "$ARG3   # sample name"

set +u
eval "$(~/miniforge3/bin/conda shell.bash hook)"
conda activate NGS
set -u

# ============================================================
# STAGE 1: UPSTREAM PROCESSING
# 1. FastQC on raw reads
# 2. Trim Galore paired-end adapter trimming
# 3. Bowtie2 paired-end contaminant removal
# 4. STAR paired-end genome alignment
# 5. Samtools alignment statistics
# 6. Save BAM to permanent storage
# 7. Write Stage 1 QC summary
# 8. Cleanup TMPDIR
# ============================================================

# ============================================================
# DEFINE PATHS
# ============================================================

PROJECT=/rds/general/project/bbsrc_sva/live/Tabitha
USER_DIR=/rds/general/project/bbsrc_sva/live
FASTQ_DIR=/rds/general/project/bbsrc_sva/live/fastq

GENOME_DIR=/rds/general/project/bbsrc_sva/live/hg38/STAR_index
GENOME_FASTA=/rds/general/project/bbsrc_sva/live/hg38/hg38.fa
ANNOTATION_GTF=/rds/general/project/bbsrc_sva/live/annotations/gencode.v47.annotation.gtf
CONTAMINANT_IDX=$PROJECT/indexes/contaminants/contaminants

ALN_DIR=$PROJECT/aln
QC_DIR=$PROJECT/QC
QC_RAW_DIR=$QC_DIR/raw/$ARG3
QC_TRIMMED_DIR=$QC_DIR/trimmed/$ARG3
QC_ALIGNED_DIR=$QC_DIR/aligned/$ARG3
QC_STATS_DIR=$PROJECT/stats
QC_SUMMARY_DIR=$QC_DIR/summaries

MIN_READS_RAW=10000000
MIN_MAPPING_RATE=70
MAX_RRNA_RATE=80
MIN_RPF_READS=100000

mkdir -p "$ALN_DIR" "$QC_RAW_DIR" "$QC_TRIMMED_DIR" "$QC_ALIGNED_DIR" "$QC_STATS_DIR" "$QC_SUMMARY_DIR"

# ============================================================
# FILE VARIABLES
# ============================================================

RAW_R1="$TMPDIR/$ARG1.fastq.gz"
RAW_R2="$TMPDIR/$ARG2.fastq.gz"
TMPDIR_GENOME_FASTA="$TMPDIR/$(basename "$GENOME_FASTA")"
TMPDIR_STAR_INDEX="$TMPDIR/STAR_index/"
TMPDIR_CONTAMINANT_IDX="$TMPDIR/contaminants/contaminants"
TRIM_REPORT_R1="$QC_STATS_DIR/${ARG1}.fastq.gz_trimming_report.txt"
TRIM_REPORT_R2="$QC_STATS_DIR/${ARG2}.fastq.gz_trimming_report.txt"
STAR_BAM="$TMPDIR/${ARG3}_Aligned.sortedByCoord.out.bam"

# ============================================================
# HEADER
# ============================================================

echo "=================================================="
echo "STAGE 1: Upstream Processing"
echo "Sample: $ARG3"
echo "Date: $(date)"
echo "=================================================="

echo "FASTQ dir:        $FASTQ_DIR"
echo "Genome FASTA:     $GENOME_FASTA"
echo "Genome dir:       $GENOME_DIR"
echo "Annotation GTF:   $ANNOTATION_GTF"
echo "Contaminant idx:  $CONTAMINANT_IDX"
echo "QC dir:           $QC_DIR"
echo

# ============================================================
# COPY FILES TO TMPDIR
# ============================================================

echo "=== Copying files to TMPDIR ==="
cd "$TMPDIR"

cp "$FASTQ_DIR/$ARG1.fastq.gz" "$TMPDIR/"
cp "$FASTQ_DIR/$ARG2.fastq.gz" "$TMPDIR/"
cp "$GENOME_FASTA" "$TMPDIR/"
cp -r "$GENOME_DIR" "$TMPDIR/STAR_index"
cp -r "$PROJECT/indexes/contaminants" "$TMPDIR/"

mkdir -p "$TMPDIR/trimmed"

echo "TMPDIR contents:"
ls -lh "$TMPDIR/"

# ============================================================
# STEP 1: FASTQC ON RAW READS
# ============================================================

echo
echo "=== STEP 1: FastQC on raw reads ==="

fastqc \
  --outdir "$QC_RAW_DIR" \
  --extract \
  "$RAW_R1"

fastqc \
  --outdir "$QC_RAW_DIR" \
  --extract \
  "$RAW_R2"

RAW_TOTAL_READS=$(grep "Total Sequences" "$QC_RAW_DIR/${ARG1}_fastqc/fastqc_data.txt" | awk '{print $3}')
RAW_QUALITY_R1=$(grep ">>Per base sequence quality" "$QC_RAW_DIR/${ARG1}_fastqc/fastqc_data.txt" | cut -f2)
RAW_QUALITY_R2=$(grep ">>Per base sequence quality" "$QC_RAW_DIR/${ARG2}_fastqc/fastqc_data.txt" | cut -f2)

echo "=== Raw QC summary: $ARG3 ==="
echo "Total read pairs: $RAW_TOTAL_READS"
echo "R1 quality status: $RAW_QUALITY_R1"
echo "R2 quality status: $RAW_QUALITY_R2"
echo "Reports saved to: $QC_RAW_DIR"
echo "================================"

if [ "$RAW_TOTAL_READS" -lt "$MIN_READS_RAW" ]; then
    echo "WARNING: Raw reads ($RAW_TOTAL_READS) below minimum ($MIN_READS_RAW)"
fi

# ============================================================
# STEP 2: TRIM GALORE
# ============================================================

echo
echo "=== STEP 2: Trim Galore paired-end trimming ==="

trim_galore \
  --paired \
  --quality 20 \
  --length 20 \
  --dont_gzip \
  --fastqc \
  --fastqc_args "--outdir $QC_TRIMMED_DIR --extract" \
  --output_dir "$TMPDIR/trimmed/" \
  "$RAW_R1" \
  "$RAW_R2"

TRIMMED_R1="$TMPDIR/trimmed/${ARG1}_val_1.fq"
TRIMMED_R2="$TMPDIR/trimmed/${ARG2}_val_2.fq"

if [ ! -f "$TRIMMED_R1" ]; then
    echo "ERROR: Trimmed R1 not found: $TRIMMED_R1"
    ls -lh "$TMPDIR/trimmed/"
    exit 1
fi

if [ ! -f "$TRIMMED_R2" ]; then
    echo "ERROR: Trimmed R2 not found: $TRIMMED_R2"
    ls -lh "$TMPDIR/trimmed/"
    exit 1
fi

cp "$TMPDIR/trimmed/${ARG1}.fastq.gz_trimming_report.txt" "$QC_STATS_DIR/"
cp "$TMPDIR/trimmed/${ARG2}.fastq.gz_trimming_report.txt" "$QC_STATS_DIR/"

DETECTED_ADAPTER_R1=$(grep "Adapter sequence:" "$TRIM_REPORT_R1" | sed -E "s/.*Adapter sequence: '([^']+)'.*/\1/")
DETECTED_ADAPTER_R2=$(grep "Adapter sequence:" "$TRIM_REPORT_R2" | sed -E "s/.*Adapter sequence: '([^']+)'.*/\1/")
READS_AFTER_TRIM=$(grep "Reads written" "$TRIM_REPORT_R1" | awk '{print $5}' | tr -d '(),')
TRIM_RATE=$(echo "scale=2; $READS_AFTER_TRIM / $RAW_TOTAL_READS * 100" | bc)

ADAPTER_STATUS_R1=$(grep ">>Adapter Content" "$QC_TRIMMED_DIR/${ARG1}_val_1_fastqc/fastqc_data.txt" | cut -f2)
ADAPTER_STATUS_R2=$(grep ">>Adapter Content" "$QC_TRIMMED_DIR/${ARG2}_val_2_fastqc/fastqc_data.txt" | cut -f2)

echo "=== Trim Galore summary: $ARG3 ==="
echo "R1 adapter detected: $DETECTED_ADAPTER_R1"
echo "R2 adapter detected: $DETECTED_ADAPTER_R2"
echo "Read pairs after trim: $READS_AFTER_TRIM"
echo "Retention rate: $TRIM_RATE%"
echo "R1 adapter status: $ADAPTER_STATUS_R1"
echo "R2 adapter status: $ADAPTER_STATUS_R2"
echo "=================================="

rm -f "$RAW_R1" "$RAW_R2"
echo "Removed raw FASTQs from TMPDIR."

# ============================================================
# STEP 3: BOWTIE2 CONTAMINANT REMOVAL (PAIRED-END)
# ============================================================

echo
echo "=== STEP 3: Contaminant removal with Bowtie2 (paired-end) ==="

bowtie2 \
  -x "$TMPDIR_CONTAMINANT_IDX" \
  -1 "$TRIMMED_R1" \
  -2 "$TRIMMED_R2" \
  --no-unal \
  --un-conc-gz "$TMPDIR/${ARG3}_clean.gz" \
  -S /dev/null \
  2> "$QC_STATS_DIR/${ARG3}_bowtie2_contaminant.txt"

echo "=== Contaminant removal: $ARG3 ==="
cat "$QC_STATS_DIR/${ARG3}_bowtie2_contaminant.txt"
echo "Saved to: $QC_STATS_DIR"
echo "=================================="

echo "TMPDIR contents after Bowtie2:"
ls -lh "$TMPDIR/"

CLEAN_R1="$TMPDIR/${ARG3}_clean.1.gz"
CLEAN_R2="$TMPDIR/${ARG3}_clean.2.gz"

if [ ! -f "$CLEAN_R1" ]; then
    echo "ERROR: Clean R1 not found: $CLEAN_R1"
    ls -lh "$TMPDIR/"
    exit 1
fi

if [ ! -f "$CLEAN_R2" ]; then
    echo "ERROR: Clean R2 not found: $CLEAN_R2"
    ls -lh "$TMPDIR/"
    exit 1
fi

RRNA_RATE=$(grep "overall alignment rate" "$QC_STATS_DIR/${ARG3}_bowtie2_contaminant.txt" | awk '{print $1}' | tr -d '%')
READS_AFTER_DECONTAM=$(zcat "$CLEAN_R1" | wc -l | awk '{print $1/4}')

echo "rRNA contamination rate: $RRNA_RATE%"
echo "Read pairs after decontamination: $READS_AFTER_DECONTAM"

if (( $(echo "$RRNA_RATE > $MAX_RRNA_RATE" | bc -l) )); then
    echo "WARNING: Very high rRNA rate ($RRNA_RATE%)"
fi

rm -f "$TRIMMED_R1" "$TRIMMED_R2"
rm -rf "$TMPDIR/trimmed/"
echo "Removed trimmed files from TMPDIR."

# ============================================================
# STEP 4: STAR GENOME ALIGNMENT (PAIRED-END)
# ============================================================

echo
echo "=== STEP 4: STAR genome alignment (paired-end) ==="

STAR \
  --genomeDir "$TMPDIR_STAR_INDEX" \
  --readFilesIn "$CLEAN_R1" "$CLEAN_R2" \
  --readFilesCommand zcat \
  --outSAMtype BAM SortedByCoordinate \
  --outSAMattributes NH HI AS nM MD \
  --outFilterMultimapNmax 1 \
  --outFilterMismatchNmax 2 \
  --alignIntronMin 20 \
  --alignIntronMax 1000000 \
  --outTmpDir "$TMPDIR/STAR_tmp_$ARG3/" \
  --outFileNamePrefix "$TMPDIR/${ARG3}_" \
  --sjdbGTFfile "$ANNOTATION_GTF"

if [ ! -f "$STAR_BAM" ]; then
    echo "ERROR: STAR BAM not found: $STAR_BAM"
    ls -lh "$TMPDIR/"
    exit 1
fi

samtools index "$STAR_BAM"

cp "$TMPDIR/${ARG3}_Log.final.out" "$QC_ALIGNED_DIR/${ARG3}_Log.final.out"
cp "$TMPDIR/${ARG3}_Log.out" "$QC_ALIGNED_DIR/${ARG3}_Log.out"

echo "STAR logs saved to: $QC_ALIGNED_DIR"

UNIQUE_MAP=$(grep "Uniquely mapped reads %" "$QC_ALIGNED_DIR/${ARG3}_Log.final.out" | awk '{print $NF}' | tr -d '%')
MULTI_MAP=$(grep "% of reads mapped to multiple loci" "$QC_ALIGNED_DIR/${ARG3}_Log.final.out" | awk '{print $NF}' | tr -d '%')
UNMAPPED_SHORT=$(grep "% of reads unmapped: too short" "$QC_ALIGNED_DIR/${ARG3}_Log.final.out" | awk '{print $NF}' | tr -d '%')
MISMATCH_RATE=$(grep "Mismatch rate per base, %" "$QC_ALIGNED_DIR/${ARG3}_Log.final.out" | awk '{print $NF}' | tr -d '%')

echo "=== STAR alignment QC: $ARG3 ==="
echo "Uniquely mapped: $UNIQUE_MAP%"
echo "Multi-mapped: $MULTI_MAP%"
echo "Unmapped (too short): $UNMAPPED_SHORT%"
echo "Mismatch rate: $MISMATCH_RATE%"
echo "================================"

if (( $(echo "$UNIQUE_MAP < $MIN_MAPPING_RATE" | bc -l) )); then
    echo "WARNING: Mapping rate ($UNIQUE_MAP%) below minimum ($MIN_MAPPING_RATE%)"
fi

rm -f "$CLEAN_R1" "$CLEAN_R2"
rm -rf "$TMPDIR/STAR_index" "$TMPDIR/STAR_tmp_$ARG3" "$TMPDIR/contaminants"
echo "Removed clean FASTQs and indexes from TMPDIR."

# ============================================================
# STEP 5: SAMTOOLS ALIGNMENT STATISTICS
# ============================================================

echo
echo "=== STEP 5: Alignment statistics with samtools ==="

samtools flagstat \
  "$STAR_BAM" \
  > "$QC_STATS_DIR/${ARG3}_alignment_stats.txt"

echo "=== Samtools flagstat: $ARG3 ==="
cat "$QC_STATS_DIR/${ARG3}_alignment_stats.txt"
echo "Saved to: $QC_STATS_DIR"
echo "================================"

FINAL_MAPPED=$(grep " mapped (" "$QC_STATS_DIR/${ARG3}_alignment_stats.txt" | head -1 | awk '{print $1}')

echo "Final mapped alignments: $FINAL_MAPPED"

if [ "$FINAL_MAPPED" -lt "$MIN_RPF_READS" ]; then
    echo "WARNING: Final mapped reads ($FINAL_MAPPED) below minimum ($MIN_RPF_READS)"
fi

# ============================================================
# STEP 6: SAVE BAM TO PERMANENT STORAGE
# ============================================================

echo
echo "=== STEP 6: Saving BAM to permanent storage ==="

cp "$STAR_BAM" "$ALN_DIR/"
cp "${STAR_BAM}.bai" "$ALN_DIR/"

echo "BAM saved to: $ALN_DIR"
echo "  $ALN_DIR/${ARG3}_Aligned.sortedByCoord.out.bam"
echo "  $ALN_DIR/${ARG3}_Aligned.sortedByCoord.out.bam.bai"

# ============================================================
# STEP 7: WRITE STAGE 1 QC SUMMARY
# ============================================================

echo
echo "=== STEP 7: Writing Stage 1 QC summary ==="

QC_SUMMARY="$QC_SUMMARY_DIR/${ARG3}_stage1_qc_summary.txt"

cat > "$QC_SUMMARY" << EOF
========================================
STAGE 1 QC SUMMARY: $ARG3
Date: $(date)
========================================

LIBRARY INFORMATION
Sample name:         $ARG3
Sequencing mode:     Paired-end
Trimming:            Trim Galore (--paired, quality 20, min length 20)

INPUT FILES
Read 1 FASTQ:        $FASTQ_DIR/$ARG1.fastq.gz
Read 2 FASTQ:        $FASTQ_DIR/$ARG2.fastq.gz

REFERENCE FILES
Genome FASTA:        $GENOME_FASTA
Genome directory:    $GENOME_DIR
Annotation GTF:      $ANNOTATION_GTF
Contaminant index:   $CONTAMINANT_IDX

AUTO-DETECTED ADAPTERS
R1 adapter:          $DETECTED_ADAPTER_R1
R2 adapter:          $DETECTED_ADAPTER_R2

RAW READ QC (FastQC)
Total read pairs:    $RAW_TOTAL_READS
R1 quality:          $RAW_QUALITY_R1
R2 quality:          $RAW_QUALITY_R2

TRIMMING (Trim Galore)
Read pairs after trim:   $READS_AFTER_TRIM
Retention rate:          $TRIM_RATE%
R1 adapter status:       $ADAPTER_STATUS_R1
R2 adapter status:       $ADAPTER_STATUS_R2

CONTAMINANT REMOVAL (Bowtie2)
rRNA contamination:      $RRNA_RATE%
Read pairs after removal: $READS_AFTER_DECONTAM

GENOME ALIGNMENT (STAR)
Uniquely mapped:         $UNIQUE_MAP%
Multi-mapped:            $MULTI_MAP%
Unmapped too short:      $UNMAPPED_SHORT%
Mismatch rate:           $MISMATCH_RATE%

FINAL ALIGNMENT (samtools)
Final mapped alignments: $FINAL_MAPPED

OUTPUT FILES
BAM file:            $ALN_DIR/${ARG3}_Aligned.sortedByCoord.out.bam
Raw FastQC:          $QC_RAW_DIR
Trimmed FastQC:      $QC_TRIMMED_DIR
STAR logs:           $QC_ALIGNED_DIR
Stats:               $QC_STATS_DIR

PASS/FAIL ASSESSMENT
EOF

PASS_FAIL=0

if [ "$RAW_TOTAL_READS" -ge "$MIN_READS_RAW" ]; then
    echo "Raw reads: PASS ($RAW_TOTAL_READS)" >> "$QC_SUMMARY"
else
    echo "Raw reads: FAIL ($RAW_TOTAL_READS)" >> "$QC_SUMMARY"
    PASS_FAIL=$((PASS_FAIL + 1))
fi

if (( $(echo "$UNIQUE_MAP >= $MIN_MAPPING_RATE" | bc -l) )); then
    echo "Mapping rate: PASS ($UNIQUE_MAP%)" >> "$QC_SUMMARY"
else
    echo "Mapping rate: FAIL ($UNIQUE_MAP%)" >> "$QC_SUMMARY"
    PASS_FAIL=$((PASS_FAIL + 1))
fi

if (( $(echo "$RRNA_RATE <= $MAX_RRNA_RATE" | bc -l) )); then
    echo "rRNA rate: PASS ($RRNA_RATE%)" >> "$QC_SUMMARY"
else
    echo "rRNA rate: FAIL ($RRNA_RATE%)" >> "$QC_SUMMARY"
    PASS_FAIL=$((PASS_FAIL + 1))
fi

if [ "$FINAL_MAPPED" -ge "$MIN_RPF_READS" ]; then
    echo "Final reads: PASS ($FINAL_MAPPED)" >> "$QC_SUMMARY"
else
    echo "Final reads: FAIL ($FINAL_MAPPED)" >> "$QC_SUMMARY"
    PASS_FAIL=$((PASS_FAIL + 1))
fi

echo >> "$QC_SUMMARY"

if [ "$PASS_FAIL" -eq 0 ]; then
    echo "OVERALL: PASS (4/4 checks passed)" >> "$QC_SUMMARY"
else
    echo "OVERALL: FAIL ($PASS_FAIL/4 checks failed)" >> "$QC_SUMMARY"
fi

echo "========================================" >> "$QC_SUMMARY"

cat "$QC_SUMMARY"
echo
echo "Stage 1 QC summary saved to: $QC_SUMMARY"

# ============================================================
# STEP 8: CLEANUP TMPDIR
# ============================================================

echo
echo "=== STEP 8: Cleaning up TMPDIR ==="

rm -f "$STAR_BAM"
rm -f "${STAR_BAM}.bai"
rm -f "$TMPDIR/${ARG3}_Log.final.out"
rm -f "$TMPDIR/${ARG3}_Log.out"
rm -f "$TMPDIR/${ARG3}_Log.progress.out"
rm -f "$TMPDIR/${ARG3}_SJ.out.tab"
rm -f "$TMPDIR/$(basename "$GENOME_FASTA")"

echo "TMPDIR cleaned."
echo
echo "=================================================="
echo "STAGE 1 COMPLETE: $ARG3"
echo "Date: $(date)"
echo "Stage 2 will start automatically"
echo "=================================================="
