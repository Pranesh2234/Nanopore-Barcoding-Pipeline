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
ZIPPED_FASTQ="nora_full_data.fastq.gz"
FASTQ="nora_full_data.fastq"
REFERENCE="reference.fasta"
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

echo $(date)

head -n 2 "$OUTDIR/parsed.csv"

# Step 2: Count barcodes and export FASTA
python3 chunk_and_count.py "$OUTDIR/parsed.csv" "$OUTDIR/chunks" "$OUTDIR/barcode_counts.csv"

# Create temp dirs
TEMP_FASTA_DIR=$(mktemp -d)
TEMP_CONSENSUS_DIR=$(mktemp -d)
CHUNKED_FASTA_BASE=$(mktemp -d)
NUM_GPUS=5

# Step 3: Export FASTA files
python3 export_fasta.py "$OUTDIR/barcode_counts.csv" "$TEMP_FASTA_DIR"

# Step 4: Split FASTA into $NUM_GPUS batches
echo "📦 Splitting FASTA into $NUM_GPUS chunks..."
i=0
for fasta in "$TEMP_FASTA_DIR"/*.fa*; do
    batch_dir="$CHUNKED_FASTA_BASE/batch_$((i % NUM_GPUS))"
    mkdir -p "$batch_dir"
    cp "$fasta" "$batch_dir/"
    ((i++))
done

# Step 5: Parallel Medaka across 5 GPUs
echo "⚙️ Running Medaka in parallel across $NUM_GPUS GPUs..."
PIDS=()
for i in $(seq 0 $((NUM_GPUS - 1))); do
    BATCH_DIR="$CHUNKED_FASTA_BASE/batch_$i"
    export CUDA_VISIBLE_DEVICES="$i"
    python3 generate_consensus_medaka.py "$BATCH_DIR" "$REFERENCE" "$TEMP_CONSENSUS_DIR" > "$OUTDIR/medaka_gpu${i}.log" 2>&1 &
    PIDS+=($!)
done

# Wait for all Medaka jobs
for pid in "${PIDS[@]}"; do
    wait "$pid"
done

# Step 6: Annotate results (final combined CSV)
echo "🧬 Running annotation..."
python3 consensus_annotation.py "$OUTDIR/barcode_counts.csv" "$TEMP_CONSENSUS_DIR" "$REFERENCE" "$OUTDIR/consensus_summary.csv"

# Clean up
rm -rf "$TEMP_FASTA_DIR" "$TEMP_CONSENSUS_DIR" "$CHUNKED_FASTA_BASE"

echo "✅ Pipeline finished at $(date)"
