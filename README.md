# 🔬 Nanopore Barcoding Pipeline

This repository provides a robust and scalable pipeline to process Nanopore sequencing data, extract barcode and ORF sequences, generate consensus sequences using Medaka, and annotate variants with reference-based comparison. This workflow is ideal for experiments involving high-throughput barcode-tagged libraries such as deep mutational scanning, functional screens, or lineage tracing.

---

## 🚀 Features

- Regex-based barcode and ORF detection from noisy Nanopore FASTQ reads
- Dynamic chunking and barcode-level sequence counting
- FASTA export per barcode for consensus calling
- Medaka-based GPU-accelerated consensus generation
- Variant calling using Minimap2, Samtools, and Bcftools
- HTCondor integration for large-scale deployment
- Automatic GPU monitoring and conda environment handling

---

## 📂 Pipeline Overview

```text
[FASTQ.gz] 
   └── extract_barcodes.py
        ↓
   [parsed.csv]
        ↓
   chunk_and_count.py
        ↓
   [barcode_counts.csv]
        ↓
   export_fasta.py
        ↓
   [bc*.fasta]
        ↓
   generate_consensus_medaka.py
        ↓
   [*_consensus.fa]
        ↓
   consensus_annotation.py
        ↓
   [consensus_annotation.csv]
```

---

## 📁 File Descriptions and I/O

### 1. `extract_barcodes.py`

**Purpose:**  
Parses Nanopore FASTQ reads, orients them forward based on a regex-defined ORF, and extracts barcode sequences. Filters low-quality or short reads.

**Input:**  
- FASTQ file (e.g., `nora_full_data.fastq`)  
- ORF regex (e.g., `CTAGT[ATGC]{1000,1200}TAGCT`)  
- Barcode regex (e.g., `CCACG[ATGC]{20}TGAGA`)

**Output:**  
- CSV file (`parsed.csv`) with columns:  
  `Barcode`, `Barcode Quality`, `Sequence`, `Sequence Quality`

---

### 2. `chunk_and_count.py`

**Purpose:**  
Counts occurrences of each barcode and gathers associated sequences. Automatically chunks large input CSVs for memory efficiency.

**Input:**  
- `parsed.csv` from step 1

**Output:**  
- `barcode_counts.csv` with columns:  
  `Barcode`, `Barcode_Count`, `Sequence_Count`, `Sequences`  
- [Optional] `chunks/` directory with chunked versions of input if it exceeds 1M rows

---

### 3. `export_fasta.py`

**Purpose:**  
Creates a FASTA file for each barcode using the sequences listed in `barcode_counts.csv`. Removes barcode prefix if present in the sequence.

**Input:**  
- `barcode_counts.csv`

**Output:**  
- Directory with individual FASTA files: `bc0.fasta`, `bc1.fasta`, etc.

---

### 4. `generate_consensus_medaka.py`

**Purpose:**  
Aligns the sequences in each barcode FASTA to a reference using Minimap2, then calls consensus sequences using Medaka (GPU required).

**Input:**  
- FASTA directory (from step 3)
- Reference FASTA (e.g., `reference.fasta`)

**Output:**  
- `*_consensus.fa` files for each barcode in a new output directory

**Dependencies:**  
- `minimap2`, `samtools`, `medaka` (in conda env)

---

### 5. `consensus_annotation.py`

**Purpose:**  
Aligns consensus sequences to reference and annotates mutations using Samtools and Bcftools.

**Input:**  
- `barcode_counts.csv`  
- Consensus FASTA directory  
- Reference FASTA

**Output:**  
- `consensus_annotation.csv` with columns:  
  `Barcode`, `Barcode_Count`, `Sequence_Count`, `Consensus_Sequence`, `Mutations`

---

## 📦 Conda Environment: `nanopore_gpuenv.tar.gz`

This tarball contains a fully pre-configured conda environment with:

- Python 3.9
- Biopython
- pandas
- Medaka (GPU-enabled)
- Samtools, Bcftools
- Minimap2

**Unpacking and activating:**

```bash
tar -xzf nanopore_gpuenv.tar.gz
export PATH=$PWD/nanopore_gpuenv/bin:$PATH
export LD_LIBRARY_PATH=$PWD/nanopore_gpuenv/lib:$LD_LIBRARY_PATH
export PYTHONPATH=$PWD/nanopore_gpuenv/lib/python3.9/site-packages:$PYTHONPATH
```

The bash wrapper (`run_pipeline.sh`) performs these steps automatically and also patches the Python shebang lines in Medaka scripts to ensure environment compatibility.

---

## 🖥️ HTCondor Integration

The pipeline is optimized for distributed compute clusters with GPU support via HTCondor.

### Submit File: `submit_pipeline.sub`

- Requests 5 GPUs, 64 GB RAM, and 500 GB disk
- Runs `run_pipeline.sh`
- Logs output to `logs/`

### Bash Wrapper: `run_pipeline.sh`

- Unzips environment and FASTQ
- Patches hardcoded Python paths
- Runs all stages sequentially
- Logs GPU usage every second

---

## 🧪 Execution

### Local (for small data or testing)

```bash
bash run_pipeline.sh
```

### HTCondor Cluster

```bash
condor_submit submit_pipeline.sub
```

---

## ⚙️ Regex Customization

You can edit these inside `run_pipeline.sh`:

```bash
ORF_REGEX='CTAGT[ATGC]{1000,1200}TAGCT'
BARCODE_REGEX='CCACG[ATGC]{20}TGAGA'
```

Adapt them based on your barcode/ORF flanking design.

---

## 📈 Output

Final results in:

```
pipeline_results/consensus_annotation.csv
```

This file summarizes per-barcode consensus sequences and variants.

---

## 📊 Logging

- Main logs: `logs/`
- GPU monitoring: `pipeline_results/gpu_log.log`

---

## 📞 Support

For bugs, feature requests, or help with deployment, open an issue or contact:

**Pranesh Kulasekhar**

---

## 📄 License

MIT License
