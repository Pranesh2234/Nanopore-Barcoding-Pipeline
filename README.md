# 🔬 Nanopore Barcoding Pipeline

This repository provides a complete pipeline for extracting barcodes and ORFs from Nanopore FASTQ reads, generating consensus sequences using Medaka, and annotating them with variant information. The pipeline is designed to scale efficiently using HTCondor and supports GPU acceleration for consensus calling.

---

## 🚀 Features

- Extracts barcode and ORF sequences from noisy Nanopore reads using regex.
- Counts barcode occurrences and unique associated sequences.
- Exports per-barcode FASTA files.
- Generates polished consensus sequences with Medaka.
- Annotates mutations relative to a reference using Minimap2, Samtools, and Bcftools.
- Supports large-scale parallel execution via HTCondor and bash automation.
- GPU monitoring and dynamic environment patching.

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
   generate_consensus_medaka.py (uses Medaka + Minimap2)
        ↓
   [*_consensus.fa]
        ↓
   consensus_annotation.py
        ↓
   [consensus_annotation.csv]
```

---

## 📁 Directory Structure

```text
project/
├── extract_barcodes.py
├── chunk_and_count.py
├── export_fasta.py
├── generate_consensus_medaka.py
├── consensus_annotation.py
├── run_pipeline.sh
├── nanopore_gpuenv.tar.gz
├── reference.fasta
├── nora_full_data.fastq.gz
├── pipeline_results/
│   └── consensus_annotation.csv
├── logs/
│   ├── pipeline.out
│   ├── pipeline.err
│   └── pipeline.log
└── submit_pipeline.sub
```

---

## ⚙️ Requirements

The pipeline uses the following tools:

- Python 3.8+
- Biopython
- pandas
- samtools
- minimap2
- bcftools
- medaka (GPU version)

Install dependencies via:

```bash
pip install -r requirements.txt
```

Conda users: Extract and use the provided environment `nanopore_gpuenv.tar.gz`.

---

## 🔧 Setup Instructions

### 1. Clone the repository

```bash
git clone https://github.com/yourusername/nanopore-barcode-pipeline.git
cd nanopore-barcode-pipeline
```

### 2. Set up the Conda environment

```bash
tar -xzf nanopore_gpuenv.tar.gz
export PATH=$PWD/nanopore_gpuenv/bin:$PATH
export LD_LIBRARY_PATH=$PWD/nanopore_gpuenv/lib:$LD_LIBRARY_PATH
export PYTHONPATH=$PWD/nanopore_gpuenv/lib/python3.9/site-packages:$PYTHONPATH
```

### 3. Run pipeline locally (optional)

```bash
bash run_pipeline.sh
```

### 4. Submit to HTCondor

```bash
condor_submit submit_pipeline.sub
```

---

## 📜 File Descriptions

### 🧬 `extract_barcodes.py`

- Parses FASTQ reads.
- Extracts barcode and ORF using regex.
- Saves barcode/ORF info and quality scores to a CSV.

### 🧮 `chunk_and_count.py`

- Counts number of reads per barcode.
- Deduplicates sequences per barcode.
- Splits input if too large (>1M rows).

### 📤 `export_fasta.py`

- Exports barcode-associated sequences to individual FASTA files (`bc0.fasta`, `bc1.fasta`, ...).

### 🧪 `generate_consensus_medaka.py`

- Aligns barcode FASTAs to reference.
- Runs Medaka to generate consensus sequences (GPU accelerated).
- Supports parallel execution using `ProcessPoolExecutor`.

### 🧾 `consensus_annotation.py`

- Aligns consensus to reference using Minimap2.
- Calls variants with Samtools + Bcftools.
- Outputs final CSV with consensus sequences and mutations.

---

## 🖥️ HTCondor Integration

### Submit File: `submit_pipeline.sub`

- Uses `run_pipeline.sh` as the entrypoint.
- Requests high GPU, memory, and disk resources.
- Auto-unpacks environment and input files.

### Bash Wrapper: `run_pipeline.sh`

- Unzips environment and FASTQ.
- Patches hardcoded Medaka shebangs.
- Sets up GPU monitoring (`nvidia-smi` loop).
- Runs all 5 pipeline steps in sequence.
- Cleans up temp files and logs.

---

## 🔍 Regex Parameters

You can customize the regex patterns used for ORF and barcode detection:

```bash
# Inside run_pipeline.sh
ORF_REGEX='CTAGT[ATGC]{1000,1200}TAGCT'
BARCODE_REGEX='CCACG[ATGC]{20}TGAGA'
```

Update these patterns to match your library design.

---

## 📈 Output

Final output file:

```
pipeline_results/consensus_annotation.csv
```

Columns:

- `Barcode`
- `Barcode_Count`
- `Sequence_Count`
- `Consensus_Sequence`
- `Mutations`

---

## 📊 Logs & Debugging

- All stdout/stderr logs go to `logs/` directory.
- GPU memory usage is tracked via `gpu_log.log`.
- Enable `DEBUG = True` in `generate_consensus_medaka.py` for detailed Medaka logs.

---

## 📞 Support

For issues or feature requests, please open an issue on GitHub or contact `Pranesh Kulasekhar`.

---

## 📄 License

This project is released under the MIT License.
