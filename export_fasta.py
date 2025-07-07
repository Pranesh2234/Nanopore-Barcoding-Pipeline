#!/usr/bin/env python3
import sys
import os
import csv
from Bio import SeqIO

csv.field_size_limit(100_000_000_000)

def export_fasta(input_csv, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    with open(input_csv) as f:
        reader = csv.DictReader(f)
        for i, row in enumerate(reader):
            barcode = row['Barcode']
            sequences = row['Sequences'].split(';')

            # Trim barcode if present
            trimmed = [s[len(barcode):] if s.startswith(barcode) else s for s in sequences]

            fasta_path = f"{output_dir}/bc{i}.fasta"
            with open(fasta_path, 'w') as out:
                for j, seq in enumerate(trimmed):
                    seq = seq.strip().upper()
                    if not seq:
                        continue  # skip empty sequences
                    out.write(f">seq{j+1}\n{seq}\n")

            # If no valid sequences were written, remove the file
            if os.path.exists(fasta_path) and os.path.getsize(fasta_path) == 0:
                os.remove(fasta_path)

    print(f"[✓] Exported FASTA files to {output_dir}")

if __name__ == "__main__":
    input_csv = sys.argv[1]
    output_dir = sys.argv[2]
    export_fasta(input_csv, output_dir)
