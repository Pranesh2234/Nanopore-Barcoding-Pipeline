#!/usr/bin/env python

import os
import csv
import re
import subprocess
import sys
import tempfile

csv.field_size_limit(100_000_000_000)

MINIMAP2 = "nanopore_gpuenv/bin/minimap2"
SAMTOOLS = "nanopore_gpuenv/bin/samtools"
BCFTOOLS = "nanopore_gpuenv/bin/bcftools"

CHUNK_SIZE = 500000  # Change if needed
CHUNK_THRESHOLD = 500000  # Threshold to start chunked processing

def read_barcode_csv(csv_path):
    with open(csv_path, 'r') as f:
        reader = csv.reader(f)
        header = next(reader)
        return list(reader)

def parse_fasta(fasta_path):
    with open(fasta_path, 'r') as f:
        lines = f.readlines()
        sequence = ''.join(line.strip() for line in lines[1:])
        return sequence

def get_mutations(consensus_path, reference_path):
    with tempfile.TemporaryDirectory() as tmpdir:
        bam_path = os.path.join(tmpdir, "aln.bam")
        vcf_gz_path = os.path.join(tmpdir, "calls.vcf.gz")

        # Step 1: minimap2 + samtools sort
        with open(bam_path, 'wb') as bam_file:
            minimap = subprocess.Popen(
                [MINIMAP, '-a', reference_path, consensus_path],
                stdout=subprocess.PIPE
            )
            sort = subprocess.Popen(
                [SAMTOOLS, 'sort', '-o', bam_path],
                stdin=minimap.stdout,
                stdout=subprocess.PIPE
            )
            minimap.stdout.close()
            sort.communicate()

        # Step 2: samtools index
        subprocess.run([SAMTOOLS, 'index', bam_path], check=True)

        # Step 3: bcftools mpileup + call
        mpileup = subprocess.Popen(
            [BCFTOOLS, 'mpileup', '-f', reference_path, bam_path],
            stdout=subprocess.PIPE
        )
        call = subprocess.Popen(
            [BCFTOOLS, 'call', "--ploidy", "1", '-mv', '-Oz', '-o', vcf_gz_path],
            stdin=mpileup.stdout
        )
        mpileup.stdout.close()
        call.communicate()

        # Step 4: bcftools view
        result = subprocess.run(
            [BCFTOOLS, 'view', vcf_gz_path],
            capture_output=True,
            text=True
        )
        mutations = []
        for line in result.stdout.splitlines():
            if not line.startswith("#"):
                parts = line.split("\t")
                chrom, pos, ref, alt = parts[0], parts[1], parts[3], parts[4]
                mutations.append(f"{chrom}:{pos} {ref}>{alt}")
        return mutations

def write_consensus_summary_csv(output_file, barcodes, consensus_dir, reference_path):
    with open(output_file, 'w') as out_csv:
        writer = csv.writer(out_csv)
        writer.writerow(['Barcode', 'Barcode_Count', 'Sequence_Count', 'Consensus_Sequence', 'Mutations'])
        
        for idx, row in enumerate(barcodes):
            filename = f"bc{idx}_consensus.fa"
            consensus_path = os.path.join(consensus_dir, filename)
            if not os.path.exists(consensus_path):
                continue  # Skip missing consensus files

            barcode, count, seq_count, _ = row
            consensus_seq = parse_fasta(consensus_path)
            mutations = get_mutations(consensus_path, reference_path)
            mutation_text = "; ".join(mutations) if mutations else ""
            writer.writerow([barcode, count, seq_count, consensus_seq, mutation_text])

if __name__ == '__main__':
    barcode_csv = sys.argv[1]
    consensus_path = sys.argv[2]
    reference = sys.argv[3]
    output_csv = sys.argv[4]

    print("🔍 Reading barcodes...")
    barcode_data = read_barcode_csv(barcode_csv)
    total_rows = len(barcode_data)
    print(f"✅ Read {total_rows} barcodes")

    print("📝 Writing summary...")
    if total_rows <= CHUNK_THRESHOLD:
        write_consensus_summary_csv(output_csv, barcode_data, consensus_path, reference)
    else:
        print(f"🔁 Large file detected. Splitting into chunks of {CHUNK_SIZE}")
        for i in range(0, total_rows, CHUNK_SIZE):
            chunk = barcode_data[i:i + CHUNK_SIZE]
            print(f"⚙️ Processing rows {i} to {i + len(chunk)}...")
            write_consensus_summary_csv(
                output_csv,
                chunk,
                consensus_path,
                reference,
                start_index=i,
                append=(i != 0)
            )
    print("✅ Done.")
