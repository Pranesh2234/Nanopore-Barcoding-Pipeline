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
- Automatic conda environment handling

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
- Reference FASTA which is the wild type nucleotide sequence of the protein whose variants are being studied in the fasta format(e.g., `reference.fasta`)

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

[Click here to download the prebuilt conda environment](https://drive.google.com/file/d/1z4pHmtPRKTFgH8KKrdxAghc7tToUImHS/view?usp=drive_link)


Once downloaded transfer to your working node using `scp command from your local terminal` 

example: `scp nanopore_gpuenv.tar.gz username@domain:path/to/directory`


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

### Submit File: `pipeline.submit`

- Requests 5 GPUs, 64 GB RAM, and 500 GB disk
- Runs `run_pipeline.sh`
- Logs output to `logs/`

### Bash Wrapper: `run_pipeline.sh`

- Unzips environment and FASTQ
- Patches hardcoded Python paths
- Checks GPU avaialabilty
- Runs all stages sequentially
- This bash wrapper acts as the executable that runs all the other codes of the pipeline.
- Make sure that it has execue permission by using the command: `chmod +x run_pipeline.sh`

---

## 🧪 Execution

### Local (for small data or testing)

- Use bash to run locally only for very small datasets that your system can handle
 `bash run_pipeline.sh`

   
### HTCondor Cluster
- Most of the datasets which are large will be run using the HTCondor cluster.
- For submitting a job to HTC, use the following command after creating a submit file. (for example: pipeline.submimt)
```bash
condor_submit pipeline.submit
```

---

## ⚙️ Regex Customization

- The regular expressions that are used to capture the barcodes and ORFs from the sequencing data are specified inside the bash wrapper.
- These exressions need to be changed according to users design
User can edit these inside `run_pipeline.sh`:

```
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
- Make sure that before submit `logs/` directory exists

---
## Full setup instructions

- Download environment from link
  - Transfer to working node using `scp`
- Create submit file: pipeline.submit
  - A submit file is the way to request resources from HTC.
  - The submit file given here is setup to use 5 gpus
  - Make sure to change the names of the fasta and fastq files you want to transfer in the submit file at the following line:
    ```
    transfer input files = ....
    ```
  - The amount of disk usage and memory requested in the given submit file is tailored to a different dataset (24 GB fastq file). Feel free to increase or decrease according to your dataset in the following lines:
    ```
    request_disk = ...
    request_memory = ...
    ```
  - All other lines can be used as they are
     
- Create bash wrapper: run_pipline.sh
   - The bash wrapper displays output as the pipeline processes the different steps in sequential order
   - Make sure to change the regular expressions acoording to your dataset in the following lines:
     ```
     ORF_REGEX=..
     BARCODE_REGEX=.. 
     ```
  - Further change the reference file name if needed:
     ```
     REFERENCE_FASTA=..
     ```
  - A very important section is the assignment of the fastq input file. As certain datasets are very large, a zipped fastq or unzipped one for smaller datasets can be used in the following way:
    **Case 1: Zipped fastq**
    ```
    ZIPPED_FASTQ=zipped.fastq
    FASTQ=
    ```
    ***Case 2: Unzipped fastq**
    ```
    ZIPPED_FASTQ=
    FASTQ=unzipped.fastq
    ```
  The wrapper will unzip and assign the unzipped fastq to `FASTQ=` for further usage.
  - If using more than 5 GPUs make sure to change the value of `NUM_GPUS=5`
  - No other real changes required. Feel free to change the names of the files if you need to.
- Once all the required files are present, recheck if all python scripts and the bash wrapper have permission to execute. If not, use:
    ```
    chmod +x *.py
    chmod +x run_pipeline.sh
    ```
- Ensure `logs/` directory exists.
- To submit the job to a HTCondor node, use:
    ```
    condor_submit pipeline.submit
    ```
- Monitor `logs/pipeline.out` and `logs/pipeline.err` for progress.

---
## Testing the pipeline
- The test data available can be used to check if the pipeline is functioning as intended.
- Regex for test data:
    ```
    ORF_REGEX='CGTAA[ATGC]{450,950}CGCGG'
    BARCODE_REGEX='GCTGC[ATGC]{16}CCATT'
    ```
    
  
