#BUILD_STAR_INDEX
#!/bin/bash

#PBS -l walltime=24:00:00
#PBS -l select=1:ncpus=32:mem=60gb
#PBS -N BUILD_STAR_INDEX
#PBS -e /rds/general/project/bbsrc_sva/live/Tabitha/logfiles/star_index_error.log
#PBS -o /rds/general/project/bbsrc_sva/live/Tabitha/logfiles/star_index_output.log

# ============================================================
# BUILD STAR GENOME INDEX
# hg38 genome
# Imperial College London HPC cx3
# Project: bbsrc_sva
# User: ldeelen
# ============================================================

module load miniforge/3
eval "$(~/miniforge3/bin/conda shell.bash hook)"
conda activate NGS

# ============================================================
# DEFINE PATHS
# ============================================================

GENOME_FASTA=/rds/general/project/bbsrc_sva/live/hg38/hg38.fa

ANNOTATION_GTF=/rds/general/project/bbsrc_sva/live/annotations/gencode.v47.annotation.gtf

STAR_INDEX_DIR=/rds/general/project/bbsrc_sva/live/hg38/STAR_index


# ============================================================
# CREATE OUTPUT DIRECTORY
# ============================================================

mkdir -p $STAR_INDEX_DIR

echo "=============================================="
echo "  Building STAR Genome Index"
echo "  Genome FASTA: $GENOME_FASTA"
echo "  GTF file:     $ANNOTATION_GTF"
echo "  Output dir:   $STAR_INDEX_DIR"
echo "=============================================="
echo ""

# ============================================================
# VERIFY INPUT FILES EXIST BEFORE STARTING
# ============================================================

echo "=== Verifying input files ==="

if [ ! -f $GENOME_FASTA ]; then
    echo "ERROR: Genome FASTA not found:"
    echo "       $GENOME_FASTA"
    exit 1
else
    FASTA_SIZE=$(ls -lh $GENOME_FASTA | awk '{print $5}')
    echo "  Genome FASTA: FOUND ($FASTA_SIZE)"
fi

if [ ! -f $ANNOTATION_GTF ]; then
    echo "ERROR: GTF annotation file not found:"
    echo "       $ANNOTATION_GTF"
    echo ""
    echo "Checking what is available in that directory..."
    ls -lh /rds/general/project/bbsrc_sva/live/indexes/kallisto/homo_sapiens/
    exit 1
else
    GTF_SIZE=$(ls -lh $ANNOTATION_GTF | awk '{print $5}')
    echo "  Annotation GTF: FOUND ($GTF_SIZE)"
fi

echo "All input files verified."
echo ""

# ============================================================
# BUILD STAR INDEX
# ============================================================
# --runMode genomeGenerate: tells STAR to build an index
# --genomeDir: where to save the index files
# --genomeFastaFiles: input genome FASTA file
# --sjdbGTFfile: annotation GTF for splice junction database
# --sjdbOverhang: read length minus 1
#   For Ribo-seq with max RPF length 30 nt: 30-1 = 29
#   This affects splice junction detection accuracy

echo "=== Building STAR index ==="
echo "This typically takes 1-2 hours..."
echo ""

STAR \
    --runMode genomeGenerate \
    --genomeDir $STAR_INDEX_DIR \
    --genomeFastaFiles $GENOME_FASTA \
    --sjdbGTFfile $ANNOTATION_GTF \
    --sjdbOverhang 29

# Check STAR exit code
if [ $? -ne 0 ]; then
    echo ""
    echo "ERROR: STAR index build failed."
    echo "Check error log:"
    echo "/rds/general/project/bbsrc_sva/live/Tabitha/logfiles/star_index_error.log"
    exit 1
fi

# ============================================================
# VERIFY INDEX WAS BUILT SUCCESSFULLY
# ============================================================

echo ""
echo "=== Verifying STAR index ==="

# Check key index files exist
REQUIRED_FILES=(
    "genomeParameters.txt"
    "Genome"
    "SA"
    "SAindex"
    "chrName.txt"
    "chrLength.txt"
)

MISSING=0
for FILE in "${REQUIRED_FILES[@]}"; do
    if [ -f $STAR_INDEX_DIR/$FILE ]; then
        echo "  OK: $FILE"
    else
        echo "  MISSING: $FILE"
        MISSING=$((MISSING + 1))
    fi
done

echo ""

# Check total index size
INDEX_SIZE=$(du -sh $STAR_INDEX_DIR/ | awk '{print $1}')
echo "Total index size: $INDEX_SIZE"
echo "(Expected: approximately 28 GB)"

echo ""
echo "Index directory contents:"
ls -lh $STAR_INDEX_DIR/

# ============================================================
# FINAL RESULT
# ============================================================

echo ""
if [ $MISSING -eq 0 ]; then
    echo "=============================================="
    echo "  STAR INDEX BUILD COMPLETE"
    echo "  Date: $(date)"
    echo "  Location: $STAR_INDEX_DIR"
    echo "  Size: $INDEX_SIZE"
    echo "=============================================="
    echo ""
    echo "Next steps:"
    echo "  1. Update STAGE1 script GENOME_DIR to:"
    echo "     $STAR_INDEX_DIR"
    echo ""
    echo "  2. Resubmit Stage 1 jobs:"
    echo "     qsub STAGE1_UPSTREAM_PROCESSING.sh"
else
    echo "=============================================="
    echo "  WARNING: $MISSING required files missing"
    echo "  Index may be incomplete"
    echo "  Check error log for details"
    echo "=============================================="
    exit 1
fi