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
            trimmed = [s[len(barcode):] if s.startswith(barcode) else s for s in sequences]
            with open(f"{output_dir}/bc{i}.fasta", 'w') as out:
                for j, seq in enumerate(trimmed):
                    header = f"seq{j+1}" if seq.strip() else f"unknown{j+1}"
                    out.write(f">{header}\n{seq}\n")
    print(f"[✓] Exported FASTA files to {output_dir}")

if __name__ == "__main__":
    input_csv = sys.argv[1]
    output_dir = sys.argv[2]
    export_fasta(input_csv, output_dir)
