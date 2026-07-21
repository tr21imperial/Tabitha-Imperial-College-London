#Riboseq_contaminants_index_script

#!/bin/bash

#PBS -l walltime=12:00:00
#PBS -l select=1:ncpus=8:mem=32gb
#PBS -N BUILD_CONTAMINANT_INDEX
#PBS -e /rds/general/project/bbsrc_sva/live/Tabitha/logfiles/contaminant_index_error.log
#PBS -o /rds/general/project/bbsrc_sva/live/Tabitha/logfiles/contaminant_index_output.log


# ============================================================
# BUILD CONTAMINANT INDEX FOR RIBO-SEQ PIPELINE
# Imperial College London HPC — cx3
# Project: bbsrc_sva
# User: ldeelen
# ============================================================
# This script downloads all contaminant sequences needed
# for Ribo-seq preprocessing:
#   - Human rRNA (5S, 5.8S, 18S, 28S, 45S)
#   - Human tRNA
#   - Human snRNA
#   - Human snoRNA
# Then builds a bowtie2 index from all sequences combined
#
# HOW TO RUN:
#   qsub BUILD_CONTAMINANT_INDEX.sh
#
# OUTPUT:
#   /rds/general/project/bbsrc_sva/live/Tabitha/
#     indexes/contaminants/
#       sequences/         ← individual FASTA files
#       contaminants.fa    ← combined FASTA
#       contaminants.1.bt2 ← bowtie2 index files
#       contaminants.2.bt2
#       contaminants.3.bt2
#       contaminants.4.bt2
#       contaminants.rev.1.bt2
#       contaminants.rev.2.bt2
#       build_summary.txt  ← build log
# ============================================================


module load miniforge/3
eval "$(~/miniforge3/bin/conda shell.bash hook)"
conda activate NGS

# ============================================================
# DEFINE PATHS
# ============================================================

BASE_DIR=/rds/general/project/bbsrc_sva/live/Tabitha
CONTAM_DIR=$BASE_DIR/indexes/contaminants
SEQ_DIR=$CONTAM_DIR/sequences
LOG_DIR=$BASE_DIR/logfiles


# ============================================================
# STEP 1: DOWNLOAD HUMAN rRNA SEQUENCES
# ============================================================
# rRNA is the dominant contaminant in Ribo-seq libraries
# Typically 60-90% of raw reads are rRNA fragments
# We download all major human rRNA species:
#   5S rRNA:   ~120 nt — transcribed by RNA Pol III
#   5.8S rRNA: ~156 nt — part of 45S pre-rRNA
#   18S rRNA:  ~1869 nt — small ribosomal subunit
#   28S rRNA:  ~5070 nt — large ribosomal subunit
#   45S rRNA:  full pre-rRNA transcript including
#              5.8S, 18S, 28S and spacers

echo "=== STEP 1: Downloading human rRNA sequences ==="

cd $SEQ_DIR
wget -q \
    "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NR_023363.1&rettype=fasta&retmode=text" \
    -O human_5S_rRNA.fa

# Verify download
if [ -s human_5S_rRNA.fa ]; then
    echo "  5S rRNA downloaded: $(grep -c '>' human_5S_rRNA.fa) sequences"
else
    echo "  WARNING: 5S rRNA download may have failed"
fi

wget -q \
    "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NR_003285.3&rettype=fasta&retmode=text" \
    -O human_5_8S_rRNA.fa

if [ -s human_5_8S_rRNA.fa ]; then
    echo "  5.8S rRNA downloaded: $(grep -c '>' human_5_8S_rRNA.fa) sequences"
else
    echo "  WARNING: 5.8S rRNA download may have failed"
fi

wget -q \
    "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NR_003286.4&rettype=fasta&retmode=text" \
    -O human_18S_rRNA.fa

if [ -s human_18S_rRNA.fa ]; then
    echo "  18S rRNA downloaded: $(grep -c '>' human_18S_rRNA.fa) sequences"
else
    echo "  WARNING: 18S rRNA download may have failed"
fi

wget -q \
    "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NR_003287.4&rettype=fasta&retmode=text" \
    -O human_28S_rRNA.fa

if [ -s human_28S_rRNA.fa ]; then
    echo "  28S rRNA downloaded: $(grep -c '>' human_28S_rRNA.fa) sequences"
else
    echo "  WARNING: 28S rRNA download may have failed"
fi

wget -q \
    "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=U13369.1&rettype=fasta&retmode=text" \
    -O human_45S_rRNA.fa

if [ -s human_45S_rRNA.fa ]; then
    echo "  45S pre-rRNA downloaded: $(grep -c '>' human_45S_rRNA.fa) sequences"
else
    echo "  WARNING: 45S pre-rRNA download may have failed"
fi

echo "rRNA downloads complete."
echo ""


# ============================================================
# STEP 2: DOWNLOAD HUMAN tRNA SEQUENCES
# ============================================================
# tRNA fragments are common contaminants in Ribo-seq
# particularly after RNase digestion
# Downloaded from GtRNAdb — the gold standard tRNA database
# Contains all human cytoplasmic tRNA sequences

echo "=== STEP 2: Downloading human tRNA sequences ==="

cd $SEQ_DIR

# Download human tRNA sequences from GtRNAdb
# hg38 (GRCh38) tRNA sequences
echo "Downloading human tRNA sequences from GtRNAdb..."
wget -q \
    http://gtrnadb.ucsc.edu/GtRNAdb2/genomes/eukaryota/Hsapi38/hg38-tRNAs.fa \
    -O human_tRNA.fa

if [ -s human_tRNA.fa ]; then
    echo "  tRNA downloaded: $(grep -c '>' human_tRNA.fa) sequences"
else
    echo "  WARNING: tRNA download may have failed"
    echo "  Trying alternative source..."

    # Alternative: download from UCSC
    wget -q \
        https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.tRNA.fa.gz \
        -O human_tRNA.fa.gz

    if [ -s human_tRNA.fa.gz ]; then
        gunzip human_tRNA.fa.gz
        echo "  tRNA downloaded from UCSC: $(grep -c '>' human_tRNA.fa) sequences"
    else
        echo "  ERROR: tRNA download failed from both sources"
        echo "  Please download manually from:"
        echo "  http://gtrnadb.ucsc.edu/GtRNAdb2/genomes/eukaryota/Hsapi38/"
    fi
fi

echo "tRNA download complete."
echo ""


# ============================================================
# STEP 3: DOWNLOAD HUMAN snRNA SEQUENCES
# ============================================================
# snRNA (small nuclear RNA) fragments can appear in
# Ribo-seq data, particularly U1, U2, U4, U5, U6
# Downloaded from GENCODE annotation

echo "=== STEP 3: Downloading human snRNA sequences ==="

cd $SEQ_DIR

# Download snRNA sequences from NCBI
# U1 snRNA (NR_004430.2)
echo "Downloading U1 snRNA..."
wget -q \
    "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NR_004430.2&rettype=fasta&retmode=text" \
    -O human_U1_snRNA.fa

# U2 snRNA (NR_002716.3)
echo "Downloading U2 snRNA..."
wget -q \
    "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NR_002716.3&rettype=fasta&retmode=text" \
    -O human_U2_snRNA.fa

# U4 snRNA (NR_004380.1)
echo "Downloading U4 snRNA..."
wget -q \
    "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NR_004380.1&rettype=fasta&retmode=text" \
    -O human_U4_snRNA.fa

# U5 snRNA (NR_004381.1)
echo "Downloading U5 snRNA..."
wget -q \
    "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NR_004381.1&rettype=fasta&retmode=text" \
    -O human_U5_snRNA.fa

# U6 snRNA (NR_004394.1)
echo "Downloading U6 snRNA..."
wget -q \
    "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NR_004394.1&rettype=fasta&retmode=text" \
    -O human_U6_snRNA.fa

# Count successful downloads
SNRNA_COUNT=$(grep -c '>' human_U*_snRNA.fa 2>/dev/null \
    | awk -F: '{sum += $2} END {print sum}')
echo "  snRNA sequences downloaded: $SNRNA_COUNT"

echo "snRNA download complete."
echo ""


# ============================================================
# STEP 4: DOWNLOAD HUMAN snoRNA SEQUENCES
# ============================================================
# snoRNA (small nucleolar RNA) fragments are common
# contaminants in Ribo-seq data
# Downloaded from NCBI

echo "=== STEP 4: Downloading human snoRNA sequences ==="

cd $SEQ_DIR

# Download snoRNA sequences from NCBI
# SNORD3A (NR_003085.2) — abundant snoRNA
echo "Downloading SNORD3A snoRNA..."
wget -q \
    "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NR_003085.2&rettype=fasta&retmode=text" \
    -O human_SNORD3A_snoRNA.fa

# SNORA73A (NR_002715.1)
echo "Downloading SNORA73A snoRNA..."
wget -q \
    "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NR_002715.1&rettype=fasta&retmode=text" \
    -O human_SNORA73A_snoRNA.fa

# Download comprehensive snoRNA set from GENCODE
# This provides broader coverage of snoRNA species
echo "Downloading comprehensive snoRNA set from GENCODE..."
wget -q \
    https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.snoRNA_transcripts.fa.gz \
    -O gencode_snoRNA.fa.gz

if [ -s gencode_snoRNA.fa.gz ]; then
    gunzip gencode_snoRNA.fa.gz
    echo "  GENCODE snoRNA downloaded: $(grep -c '>' gencode_snoRNA.fa) sequences"
else
    echo "  WARNING: GENCODE snoRNA download failed"
    echo "  Using individual snoRNA sequences only"
fi

echo "snoRNA download complete."
echo ""


# ============================================================
# STEP 5: VERIFY ALL DOWNLOADS
# ============================================================
# Check all sequence files were downloaded successfully
# before combining them into the contaminant database

echo "=== STEP 5: Verifying all downloads ==="

echo ""
echo "Files in $SEQ_DIR:"
ls -lh $SEQ_DIR/

echo ""
echo "Sequence counts per file:"
for FA in $SEQ_DIR/*.fa; do
    if [ -s $FA ]; then
        COUNT=$(grep -c '>' $FA)
        SIZE=$(ls -lh $FA | awk '{print $5}')
        echo "  $(basename $FA): $COUNT sequences ($SIZE)"
    else
        echo "  $(basename $FA): EMPTY OR MISSING"
    fi
done

echo ""


# ============================================================
# STEP 6: COMBINE ALL SEQUENCES INTO SINGLE FASTA
# ============================================================
# All contaminant sequences are combined into a single
# FASTA file for bowtie2 index building
# Each sequence gets a clear header identifying its type

echo "=== STEP 6: Combining all sequences ==="

cd $CONTAM_DIR

# Combine all downloaded FASTA files
cat \
    $SEQ_DIR/human_5S_rRNA.fa \
    $SEQ_DIR/human_5_8S_rRNA.fa \
    $SEQ_DIR/human_18S_rRNA.fa \
    $SEQ_DIR/human_28S_rRNA.fa \
    $SEQ_DIR/human_45S_rRNA.fa \
    $SEQ_DIR/human_tRNA.fa \
    $SEQ_DIR/human_U1_snRNA.fa \
    $SEQ_DIR/human_U2_snRNA.fa \
    $SEQ_DIR/human_U4_snRNA.fa \
    $SEQ_DIR/human_U5_snRNA.fa \
    $SEQ_DIR/human_U6_snRNA.fa \
    > $CONTAM_DIR/contaminants.fa

# Add snoRNA if download was successful
if [ -s $SEQ_DIR/gencode_snoRNA.fa ]; then
    cat $SEQ_DIR/gencode_snoRNA.fa \
        >> $CONTAM_DIR/contaminants.fa
    echo "  GENCODE snoRNA added to combined file"
fi

if [ -s $SEQ_DIR/human_SNORD3A_snoRNA.fa ]; then
    cat $SEQ_DIR/human_SNORD3A_snoRNA.fa \
        >> $CONTAM_DIR/contaminants.fa
fi

if [ -s $SEQ_DIR/human_SNORA73A_snoRNA.fa ]; then
    cat $SEQ_DIR/human_SNORA73A_snoRNA.fa \
        >> $CONTAM_DIR/contaminants.fa
fi

# Verify combined file
TOTAL_SEQS=$(grep -c '>' $CONTAM_DIR/contaminants.fa)
TOTAL_SIZE=$(ls -lh $CONTAM_DIR/contaminants.fa \
    | awk '{print $5}')

echo "Combined contaminants.fa:"
echo "  Total sequences: $TOTAL_SEQS"
echo "  File size:       $TOTAL_SIZE"

if [ "$TOTAL_SEQS" -lt 10 ]; then
    echo "WARNING: Very few sequences in combined file"
    echo "         Check individual downloads above"
fi

echo "Sequence combination complete."
echo ""


# ============================================================
# STEP 7: BUILD BOWTIE2 INDEX
# ============================================================
# Builds the bowtie2 index from the combined contaminant
# FASTA file. This index is used in Step 3 of the main
# Ribo-seq pipeline to remove contaminant reads.
# The index name "contaminants" produces 6 index files:
#   contaminants.1.bt2
#   contaminants.2.bt2
#   contaminants.3.bt2
#   contaminants.4.bt2
#   contaminants.rev.1.bt2
#   contaminants.rev.2.bt2

echo "=== STEP 7: Building bowtie2 index ==="

cd $CONTAM_DIR

bowtie2-build \
    --threads 8 \
    $CONTAM_DIR/contaminants.fa \
    $CONTAM_DIR/contaminants \
    2>&1 | tee $CONTAM_DIR/bowtie2_build.log

# Verify index was built successfully
echo ""
echo "Checking bowtie2 index files:"
INDEX_FILES=$(ls $CONTAM_DIR/contaminants*.bt2 \
    2>/dev/null | wc -l)

if [ "$INDEX_FILES" -eq 6 ]; then
    echo "  All 6 bowtie2 index files present — SUCCESS"
    ls -lh $CONTAM_DIR/contaminants*.bt2
else
    echo "  ERROR: Expected 6 index files, found $INDEX_FILES"
    echo "  Check bowtie2 build log:"
    echo "  $CONTAM_DIR/bowtie2_build.log"
    exit 1
fi

echo ""
echo "bowtie2 index build complete."
echo ""


# ============================================================
# STEP 8: TEST THE INDEX
# ============================================================
# Quick test to verify the index works correctly
# Uses a known rRNA sequence to confirm alignment

echo "=== STEP 8: Testing bowtie2 index ==="

# Create a short test read from 18S rRNA
# This should align to the contaminant index
head -2 $SEQ_DIR/human_18S_rRNA.fa \
    | tail -1 \
    | cut -c1-50 \
    > /rds/general/project/bbsrc_sva/live/Tabitha/indexes/contaminants/test_read.fa

echo ">test_rRNA_read" > /rds/general/project/bbsrc_sva/live/Tabitha/indexes/contaminants/test_read.fa
echo "TACCTGGTTGATCCTGCCAGTAGCATATGCTTGTCTCAAAGATTAAGCC" \
    >> /rds/general/project/bbsrc_sva/live/Tabitha/indexes/contaminants/test_read.fa

# Run test alignment
bowtie2 \
    -x $CONTAM_DIR/contaminants \
    -f /rds/general/project/bbsrc_sva/live/Tabitha/indexes/contaminants/test_read.fa \
    -S /rds/general/project/bbsrc_sva/live/Tabitha/indexes/contaminants/test_alignment.sam \
    2> /rds/general/project/bbsrc_sva/live/Tabitha/indexes/contaminants/test_stats.txt

# Check alignment result
ALIGNED=$(grep "1 reads; of these:" \
    /rds/general/project/bbsrc_sva/live/Tabitha/indexes/contaminants//test_stats.txt \
    | awk '{print $1}')

ALIGN_RATE=$(grep "overall alignment rate" \
    /rds/general/project/bbsrc_sva/live/Tabitha/indexes/contaminants//test_stats.txt \
    | awk '{print $1}')

echo "Test alignment results:"
echo "  Reads tested:     1"
echo "  Alignment rate:   $ALIGN_RATE"

if [ "$ALIGN_RATE" == "100.00%" ]; then
    echo "  Index test:       PASS"
else
    echo "  Index test:       WARNING"
    echo "  Test read did not align at 100%"
    echo "  This may be acceptable — check manually"
fi

echo ""
