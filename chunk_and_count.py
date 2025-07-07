#!/usr/bin/env python3

import os
import sys
import csv
from collections import defaultdict
import pandas as pd

csv.field_size_limit(100_000_000_000)

CHUNK_SIZE = 500000

def chunk_csv(input_csv, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    with open(input_csv, 'r') as infile:
        header = infile.readline()
        count = 0
        out_file = open(f"{output_dir}/chunk_{count}.csv", 'w')
        out_file.write(header)
        for i, line in enumerate(infile, 1):
            if i % CHUNK_SIZE == 0:
                out_file.close()
                count += 1
                out_file = open(f"{output_dir}/chunk_{count}.csv", 'w')
                out_file.write(header)
            out_file.write(line)
        out_file.close()
    print(f"[✓] Chunked input CSV into {count+1} files")

def count_barcodes(input_csv, output_csv):
    stats = defaultdict(lambda: {"count": 0, "sequences": []})
    with open(input_csv) as f:
        reader = csv.reader(f)
        next(reader)
        for row in reader:
            barcode, _, seq, _ = row
            trimmed_seq = seq[len(barcode):] if seq.startswith(barcode) else seq
            stats[barcode]["count"] += 1
            stats[barcode]["sequences"].append(trimmed_seq)

    with open(output_csv, 'w') as f:
        writer = csv.writer(f)
        writer.writerow(["Barcode", "Barcode_Count", "Sequence_Count", "Sequences"])
        for bc, data in stats.items():
            unique = list(set(data["sequences"]))
            writer.writerow([bc, data["count"], len(unique), ";".join(unique)])
    print(f"[✓] Wrote barcode counts to {output_csv}")

if __name__ == "__main__":
    input_csv = sys.argv[1]
    output_dir = sys.argv[2]
    barcode_csv = sys.argv[3]

    # Check the number of rows in the input CSV
    num_rows = sum(1 for _ in open(input_csv)) - 1  # Subtract 1 for the header

    if num_rows >= 1_000_000:
        chunk_csv(input_csv, output_dir)

    count_barcodes(input_csv, barcode_csv)
