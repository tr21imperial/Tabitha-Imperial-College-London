#STAGE2_RIBOSEQ_ANALYSIS.sh

#!/bin/bash
#PBS -l walltime=48:00:00
#PBS -l select=1:ncpus=14:mem=100gb
#PBS -N STAGE2_RIBOSEQ

set -euo pipefail

# ============================================================
# CALLING VARIABLES
# ARG1 = R1 FASTQ basename (reference only)
# ARG2 = R2 FASTQ basename (reference only)
# ARG3 = sample name
# ============================================================

echo "$ARG1   # R1 FASTQ (reference only)"
echo "$ARG2   # R2 FASTQ (reference only)"
echo "$ARG3   # Sample output name"

module load miniforge/3
eval "$(~/miniforge3/bin/conda shell.bash hook)"
conda activate NGS


# ============================================================
# STAGE 2: RIBO-SEQ ANALYSIS
# 1. Ribo-seQC P-site correction & periodicity QC
# 2. RiboTaper ORF discovery
# 3. BEDTools P-site counting over MetamORF
# 4. Write Stage 2 QC summary
# 5. Run MultiQC when all samples complete
# 6. Verify outputs
# 7. Cleanup TMPDIR
# ============================================================

PROJECT=/rds/general/project/bbsrc_sva/live/Tabitha
USER_DIR=/rds/general/project/bbsrc_sva/live

GENOME_FASTA=/rds/general/project/bbsrc_sva/live/hg38/hg38.fa
ANNOTATION_GTF=/rds/general/project/bbsrc_sva/live/annotations/gencode.v47.annotation.gtf
METAMORF_BED=$PROJECT/MetamORF/MetamORF_library.bed

RIBOTAPER_DIR=$USER_DIR/tools/RiboTaper/scripts
RIBOSEQC_SCRIPT=$USER_DIR/tools/RiboseQC/R/RiboseQC.R

ALN_DIR=$PROJECT/aln
BAM_INPUT=$ALN_DIR/${ARG3}_Aligned.sortedByCoord.out.bam

QC_DIR=$PROJECT/QC
QC_RIBOSEQC_DIR=$PROJECT/riboseqc
QC_STATS_DIR=$PROJECT/stats
QC_SUMMARY_DIR=$QC_DIR/summaries
QC_MULTIQC_DIR=$QC_DIR/multiqc

PSITES_DIR=$PROJECT/psites
RIBOTAPER_OUT=$PROJECT/ribotaper
RIBOTAPER_ANNOT=$PROJECT/ribotaper_annot
BEDTOOLS_DIR=$PROJECT/bedtools

RIBOTAPER_LENGTHS="26,27,28,29,30"
RIBOTAPER_OFFSETS="9,9,12,12,12"
THREADS=14

mkdir -p "$QC_RIBOSEQC_DIR" "$QC_SUMMARY_DIR" "$QC_MULTIQC_DIR"
mkdir -p "$PSITES_DIR" "$RIBOTAPER_OUT" "$RIBOTAPER_ANNOT" "$BEDTOOLS_DIR"

echo "=================================================="
echo "STAGE 2: Ribo-seq Analysis"
echo "Sample: $ARG3"
echo "Date: $(date)"
echo "=================================================="
echo "Input BAM:       $BAM_INPUT"
echo "Genome FASTA:    $GENOME_FASTA"
echo "Annotation GTF:  $ANNOTATION_GTF"
echo "MetamORF BED:    $METAMORF_BED"
echo "Ribo-seQC out:   $QC_RIBOSEQC_DIR"
echo "RiboTaper out:   $RIBOTAPER_OUT"
echo "BEDTools out:    $BEDTOOLS_DIR"
echo

# ============================================================
# VERIFY INPUT BAM EXISTS
# ============================================================

if [ ! -f "$BAM_INPUT" ]; then
    echo "ERROR: Stage 1 BAM not found: $BAM_INPUT"
    echo "Check: $QC_SUMMARY_DIR/${ARG3}_stage1_qc_summary.txt"
    exit 1
fi

echo "Stage 1 BAM found: $BAM_INPUT"

# ============================================================
# COPY FILES TO TMPDIR
# ============================================================

echo
echo "=== Copying files to TMPDIR ==="

cd "$TMPDIR"

cp "$BAM_INPUT" "$TMPDIR/"
cp "$GENOME_FASTA" "$TMPDIR/"
cp "$ANNOTATION_GTF" "$TMPDIR/"
cp "$METAMORF_BED" "$TMPDIR/"

TMPDIR_BAM="$TMPDIR/${ARG3}_Aligned.sortedByCoord.out.bam"
TMPDIR_GENOME_FASTA="$TMPDIR/$(basename "$GENOME_FASTA")"
TMPDIR_ANNOTATION_GTF="$TMPDIR/$(basename "$ANNOTATION_GTF")"
TMPDIR_METAMORF_BED="$TMPDIR/$(basename "$METAMORF_BED")"

samtools index "$TMPDIR_BAM"

if [ -f "$RIBOTAPER_ANNOT/annotation_complete.flag" ]; then
    cp -r "$RIBOTAPER_ANNOT" "$TMPDIR/ribotaper_annot"
fi

mkdir -p "$TMPDIR/ribotaper_out/$ARG3"

echo "TMPDIR contents:"
ls -lh "$TMPDIR/"

# ============================================================
# STEP 1: RIBO-SEQC
# ============================================================

echo
echo "=== STEP 1: Ribo-seQC P-site correction ==="
echo "QC report output: $QC_RIBOSEQC_DIR"

Rscript "$RIBOSEQC_SCRIPT" \
  --annotation_file "$TMPDIR_ANNOTATION_GTF" \
  --annotation_file_type "gtf" \
  --genome_seq "$TMPDIR_GENOME_FASTA" \
  --input_files "$TMPDIR_BAM" \
  --output_dir "$QC_RIBOSEQC_DIR/" \
  --report_file "$QC_RIBOSEQC_DIR/${ARG3}_riboseq_QC" \
  --write_tmp_files TRUE

echo "Ribo-seQC complete."
echo "QC report: $QC_RIBOSEQC_DIR/${ARG3}_riboseq_QC.html"

PSITE_BED_SRC="$QC_RIBOSEQC_DIR/${ARG3}_psites_5.bed"

if [ -f "$PSITE_BED_SRC" ]; then
    cp "$PSITE_BED_SRC" "$PSITES_DIR/${ARG3}_psites_5.bed"
    echo "P-site BED saved to: $PSITES_DIR"
else
    echo "WARNING: P-site BED not found."
    ls -lh "$QC_RIBOSEQC_DIR/"
    exit 1
fi


# ============================================================
# STEP 2: RIBOTAPER ORF DISCOVERY
# ============================================================

echo
echo "=== STEP 2: RiboTaper ORF discovery ==="

if [ ! -f "$RIBOTAPER_ANNOT/annotation_complete.flag" ]; then
    echo "Creating RiboTaper annotation files..."
    echo "This should be done only once."

    "$RIBOTAPER_DIR/create_annotations_files.bash" \
      "$TMPDIR_ANNOTATION_GTF" \
      "$TMPDIR_GENOME_FASTA" \
      false \
      false \
      "$RIBOTAPER_ANNOT/"

    touch "$RIBOTAPER_ANNOT/annotation_complete.flag"
    echo "RiboTaper annotation created."
else
    echo "RiboTaper annotation exists, skipping creation."
fi

if [ ! -d "$TMPDIR/ribotaper_annot" ]; then
    cp -r "$RIBOTAPER_ANNOT" "$TMPDIR/ribotaper_annot"
fi

"$RIBOTAPER_DIR/RiboTaper.sh" \
  "$TMPDIR_BAM" \
  "$TMPDIR_BAM" \
  "$TMPDIR/ribotaper_annot/" \
  "$RIBOTAPER_LENGTHS" \
  "$RIBOTAPER_OFFSETS" \
  "$THREADS" \
  "$TMPDIR/ribotaper_out/$ARG3/"

TOTAL_ORFS=$(wc -l < "$TMPDIR/ribotaper_out/$ARG3/ORFs_final.txt")

echo
echo "=== RiboTaper summary: $ARG3 ==="
echo "Total ORFs discovered: $TOTAL_ORFS"
echo "================================"

cp "$TMPDIR/ribotaper_out/$ARG3/ORFs_final.txt" \
   "$RIBOTAPER_OUT/${ARG3}_ORFs_final.txt"

cp "$TMPDIR/ribotaper_out/$ARG3/ORFs_max.txt" \
   "$RIBOTAPER_OUT/${ARG3}_ORFs_max.txt"

echo "RiboTaper outputs saved to: $RIBOTAPER_OUT"

rm -rf "$TMPDIR/ribotaper_annot/"
echo "Removed RiboTaper annotation from TMPDIR."

rm -f "$TMPDIR_ANNOTATION_GTF" "$TMPDIR_GENOME_FASTA"
echo "Removed annotation and genome from TMPDIR."

# ============================================================
# STEP 3: BEDTOOLS P-SITE COUNTING OVER METAMORF
# ============================================================

echo
echo "=== STEP 3: BEDTools P-site counting ==="

PSITE_BED="$PSITES_DIR/${ARG3}_psites_5.bed"

if [ ! -f "$PSITE_BED" ]; then
    echo "ERROR: P-site BED not found: $PSITE_BED"
    ls -lh "$PSITES_DIR/"
    exit 1
fi

bedtools coverage \
  -a "$TMPDIR_METAMORF_BED" \
  -b "$PSITE_BED" \
  -counts \
  -s \
  > "$BEDTOOLS_DIR/${ARG3}_metamorf_counts.txt"

SORFS_WITH_COVERAGE=$(awk '$NF > 0' "$BEDTOOLS_DIR/${ARG3}_metamorf_counts.txt" | wc -l)
TOTAL_SORFS=$(wc -l < "$BEDTOOLS_DIR/${ARG3}_metamorf_counts.txt")

echo "=== BEDTools MetamORF coverage: $ARG3 ==="
echo "Total MetamORF sORFs: $TOTAL_SORFS"
echo "sORFs with P-site coverage: $SORFS_WITH_COVERAGE"
echo "Saved to: $BEDTOOLS_DIR"
echo "========================================"

rm -f "$TMPDIR_METAMORF_BED"
echo "Removed MetamORF BED from TMPDIR."

# ============================================================
# STEP 4: WRITE STAGE 2 QC SUMMARY
# ============================================================

echo
echo "=== STEP 4: Writing Stage 2 QC summary ==="

STAGE1_SUMMARY="$QC_SUMMARY_DIR/${ARG3}_stage1_qc_summary.txt"

if [ -f "$STAGE1_SUMMARY" ]; then
    RAW_READS=$(grep "Total read pairs:" "$STAGE1_SUMMARY" | awk '{print $NF}')
    MAPPING_RATE=$(grep "Uniquely mapped:" "$STAGE1_SUMMARY" | awk '{print $NF}')
    RRNA=$(grep "rRNA contamination:" "$STAGE1_SUMMARY" | awk '{print $NF}')
    STAGE1_RESULT=$(grep "OVERALL:" "$STAGE1_SUMMARY" | awk '{print $2}')
else
    RAW_READS="See Stage 1 summary"
    MAPPING_RATE="See Stage 1 summary"
    RRNA="See Stage 1 summary"
    STAGE1_RESULT="See Stage 1 summary"
fi

QC_SUMMARY="$QC_SUMMARY_DIR/${ARG3}_stage2_qc_summary.txt"

cat > "$QC_SUMMARY" << EOF
========================================
STAGE 2 QC SUMMARY: $ARG3
Date: $(date)
========================================

SAMPLE INFORMATION
Sample name:            $ARG3
Input BAM:              $BAM_INPUT

REFERENCE FILES
Genome FASTA:           $GENOME_FASTA
Annotation GTF:         $ANNOTATION_GTF
MetamORF BED:           $METAMORF_BED

STAGE 1 RESULTS
Raw read pairs:         $RAW_READS
Mapping rate:           $MAPPING_RATE
rRNA contamination:     $RRNA
Stage 1 overall:        $STAGE1_RESULT

RIBO-seQC
QC report:              $QC_RIBOSEQC_DIR/${ARG3}_riboseq_QC.html
P-site BED:             $PSITES_DIR/${ARG3}_psites_5.bed
Note: check HTML report for periodicity, read-length distribution and frame preference.
Update RIBOTAPER_LENGTHS and RIBOTAPER_OFFSETS in STAGE2_RIBOSEQ_ANALYSIS.sh if needed.

RIBOTAPER ORF DISCOVERY
Total ORFs found:       $TOTAL_ORFS
ORFs_final.txt:         $RIBOTAPER_OUT/${ARG3}_ORFs_final.txt
ORFs_max.txt:           $RIBOTAPER_OUT/${ARG3}_ORFs_max.txt

METAMORF COVERAGE
Total sORFs tested:     $TOTAL_SORFS
sORFs with coverage:    $SORFS_WITH_COVERAGE
Count file:             $BEDTOOLS_DIR/${ARG3}_metamorf_counts.txt

PERMANENT OUTPUT LOCATIONS
Ribo-seQC report:       $QC_RIBOSEQC_DIR
P-site BED:             $PSITES_DIR
RiboTaper ORFs:         $RIBOTAPER_OUT
BEDTools counts:        $BEDTOOLS_DIR
QC summaries:           $QC_SUMMARY_DIR
MultiQC report:         $QC_MULTIQC_DIR

NEXT STEPS
1. Check Ribo-seQC HTML report
2. Update P-site offsets in STAGE2 script if required
3. Run downstream R analysis
========================================
EOF

cat "$QC_SUMMARY"
echo
echo "Stage 2 QC summary saved to: $QC_SUMMARY"

# ============================================================
# STEP 5: RUN MULTIQC WHEN ALL SAMPLES COMPLETE
# ============================================================

echo
echo "=== STEP 5: Checking if MultiQC can run ==="

COMPLETED_STAGE2=$(ls "$QC_SUMMARY_DIR"/*_stage2_qc_summary.txt 2>/dev/null | wc -l)

echo "Stage 2 completed samples: $COMPLETED_STAGE2 / 3"

if [ "$COMPLETED_STAGE2" -ge 3 ]; then
    echo "All samples complete — running MultiQC..."

    multiqc \
      "$QC_DIR/raw/" \
      "$QC_DIR/trimmed/" \
      "$QC_DIR/aligned/" \
      "$QC_STATS_DIR" \
      "$QC_RIBOSEQC_DIR" \
      --outdir "$QC_MULTIQC_DIR" \
      --filename "multiqc_riboseq_report" \
      --title "Ribo-seq sORF Pipeline QC — bbsrc_sva" \
      --comment "Brain organoid Ribo-seq: MT1155, MT1156, MT1157" \
      --force

    echo "MultiQC report saved to:"
    echo "  $QC_MULTIQC_DIR/multiqc_riboseq_report.html"
else
    echo "$COMPLETED_STAGE2/3 samples done."
    echo "MultiQC will run when all samples complete."
fi

# ============================================================
# STEP 6: VERIFY ALL STAGE 2 OUTPUT FILES
# ============================================================

echo
echo "=== STEP 6: Verifying Stage 2 output files ==="

MISSING=0

check_output() {
    if [ ! -f "$1" ]; then
        echo "MISSING: $1"
        MISSING=$((MISSING + 1))
    else
        echo "OK: $(basename "$1")"
    fi
}

echo
echo "Ribo-seQC outputs:"
check_output "$QC_RIBOSEQC_DIR/${ARG3}_riboseq_QC.html"
check_output "$PSITES_DIR/${ARG3}_psites_5.bed"

echo
echo "RiboTaper outputs:"
check_output "$RIBOTAPER_OUT/${ARG3}_ORFs_final.txt"
check_output "$RIBOTAPER_OUT/${ARG3}_ORFs_max.txt"

echo
echo "BEDTools outputs:"
check_output "$BEDTOOLS_DIR/${ARG3}_metamorf_counts.txt"

echo
echo "QC summaries:"
check_output "$QC_SUMMARY_DIR/${ARG3}_stage1_qc_summary.txt"
check_output "$QC_SUMMARY_DIR/${ARG3}_stage2_qc_summary.txt"

if [ "$MISSING" -gt 0 ]; then
    echo "WARNING: $MISSING output file(s) missing!"
else
    echo "All Stage 2 output files verified."
fi

# ============================================================
# STEP 7: CLEAN UP TMPDIR
# ============================================================

echo
echo "=== STEP 7: Cleaning up TMPDIR ==="

rm -f "$TMPDIR/${ARG3}_Aligned.sortedByCoord.out.bam"
rm -f "$TMPDIR/${ARG3}_Aligned.sortedByCoord.out.bam.bai"
rm -f "$TMPDIR/$(basename "$GENOME_FASTA")"
rm -rf "$TMPDIR/ribotaper_out"

echo "TMPDIR cleaned."

echo
echo "=================================================="
echo "STAGE 2 COMPLETE: $ARG3"
echo "Date: $(date)"
echo "=================================================="
echo "All outputs saved to:"
echo "  Ribo-seQC:   $QC_RIBOSEQC_DIR"
echo "  P-sites:     $PSITES_DIR"
echo "  RiboTaper:   $RIBOTAPER_OUT"
echo "  BEDTools:    $BEDTOOLS_DIR"
echo "  QC summary:  $QC_SUMMARY_DIR"
