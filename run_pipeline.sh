#!/bin/bash

echo "🚀 Script started at $(date)"

# Extract conda environment
tar -xzf nanopore_gpuenv.tar.gz

# 🛠️ Fix shebangs in Medaka and related binaries
echo "🛠️  Patching hardcoded Python shebangs in environment..."
find nanopore_gpuenv/bin/ -type f -exec sed -i "1s|^#!.*python.*$|#!$PWD/nanopore_gpuenv/bin/python3|" {} +

# Setup environment variables
export PATH="$PWD/nanopore_gpuenv/bin:$PATH"
export LD_LIBRARY_PATH="$PWD/nanopore_gpuenv/lib:$LD_LIBRARY_PATH"
export PYTHONPATH="$PWD/nanopore_gpuenv/lib/python3.9/site-packages:$PYTHONPATH"

# Print key tools
which python
which medaka_consensus

# GPU Check
echo "🔍 Checking for GPU..."
nvidia-smi || echo "⚠️ No nvidia-smi found."

export CUDA_VISIBLE_DEVICES=$(nvidia-smi --query-gpu=index --format=csv,noheader | paste -sd "," -)

# Input arguments
#DIST="$1"
ZIPPED_FASTQ="nora_full_data.fastq.gz"
FASTQ="$nora_full_data.fastq"
REFERENCE="reference.fasta"
#ORF_REGEX='(?P<flank5>CGTAA)(?P<orf>[ATGC]{450,950})(?P<flank3>CGCGG)'
#BARCODE_REGEX='(?P<flankL>GCTGC)(?P<barcode>[ATGC]{16})(?P<flankR>CCATT)'
ORF_REGEX='CTAGT[ATGC]{1000,1200}TAGCT'
BARCODE_REGEX='CCACG[ATGC]{20}TGAGA'
OUTDIR="pipeline_results"

# 🔓 Unzip FASTQ if compressed
if [[ -f "$ZIPPED_FASTQ" ]]; then
    echo "📦 Unzipping $ZIPPED_FASTQ..."
    gunzip -c "$ZIPPED_FASTQ" > "$FASTQ"
fi

mkdir -p "$OUTDIR"

# Step 1: Extract ORFs and barcodes
echo "📥 Extracting barcodes..."
python3 extract_barcodes.py "$FASTQ" "$OUTDIR/parsed.csv" "$ORF_REGEX" "$BARCODE_REGEX"

# Step 2: Chunk and count
echo "🔢 Counting barcodes..."
python3 chunk_and_count.py "$OUTDIR/parsed.csv" "$OUTDIR/chunks" "$OUTDIR/barcode_counts.csv"

# Validate barcode count file
if [[ ! -f "$OUTDIR/barcode_counts.csv" ]]; then
    echo "❌ Error: barcode_counts CSV not created"
    exit 1
fi

# Step 3: Create temp directories for FASTA and consensus
TEMP_FASTA_DIR=$(mktemp -d)
TEMP_CONSENSUS_DIR=$(mktemp -d)

# Step 4: Export FASTA
echo "📦 Exporting barcode FASTAs..."
python3 export_fasta.py "$OUTDIR/barcode_counts.csv" "$TEMP_FASTA_DIR"

# Step 5: Monitor GPU during consensus
MONITOR_LOG="$OUTDIR/gpu_log.log"
(
  echo "timestamp,pid,process_name,gpu_uuid,used_gpu_memory [MiB]"
  while true; do
    nvidia-smi --query-compute-apps=pid,process_name,gpu_uuid,used_gpu_memory --format=csv,noheader,nounits | \
    awk -v timestamp="$(date +%Y-%m-%dT%H:%M:%S)" '{ print timestamp "," $0 }'
    sleep 1
  done
) >> "$MONITOR_LOG" &
MONITOR_PID=$!

# Step 6: Generate consensus
echo "🧬 Generating consensus..."
python3 generate_consensus_medaka.py "$TEMP_FASTA_DIR" "$REFERENCE" "$TEMP_CONSENSUS_DIR"

# Stop GPU monitor
kill $MONITOR_PID

# Step 7: Annotate results
echo "🧾 Annotating consensus..."
python3 consensus_annotation.py "$OUTDIR/barcode_counts.csv" "$TEMP_CONSENSUS_DIR" "$REFERENCE" "$OUTDIR/consensus_annotation.csv"

# Step 8: Cleanup
rm -rf "$TEMP_FASTA_DIR" "$TEMP_CONSENSUS_DIR"

echo "✅ Pipeline finished at $(date). Results in: $OUTDIR"
