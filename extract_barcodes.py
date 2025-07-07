import re
import csv
from Bio import SeqIO
import sys

csv.field_size_limit(100_000_000_000)

def get_phred_scores(quality_string):
    return [ord(char) - 33 for char in quality_string]

def filter_reads_by_length(record, min_length=50):
    return len(record.seq) >= min_length

def orient_forward(record, orf_regex):
    seq_str = str(record.seq)
    rc_seq_str = str(record.seq.reverse_complement())
    phred_quality = record.letter_annotations["phred_quality"]

    match_fwd = re.search(orf_regex, seq_str)
    match_rev = re.search(orf_regex, rc_seq_str)

    if match_fwd and not match_rev:
        return seq_str, phred_quality
    elif match_rev and not match_fwd:
        return rc_seq_str, phred_quality[::-1]
    elif match_fwd and match_rev:
        # If ORF found in both, use the one that appears first in the original
        return seq_str, record.letter_annotations["phred_quality"]
    else:
        return None, None

def extract_barcode_and_quality(sequence, quality, barcode_regex):
    match = re.search(barcode_regex, sequence)
    if match:
        barcode = match.group()
        start, end = match.span()
        barcode_quality = quality[start:end]
        return barcode, barcode_quality
    return None, None

def process_fastq(fastq_file, output_csv, orf_regex, barcode_regex, min_length=250):
    barcodes = []
    orfs = []
    total = 0
    no_orf = 0
    no_barcode = 0
    passed = 0
    with open(output_csv, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["Barcode", "Barcode Quality", "Sequence", "Sequence Quality"])
        
        for record in SeqIO.parse(fastq_file, "fastq"):
            total += 1
            if not filter_reads_by_length(record, min_length):
                continue
            
            oriented_seq, oriented_quality = orient_forward(record, orf_regex)
            if oriented_seq is None:
                no_orf += 1
                continue
            
            barcode, barcode_quality = extract_barcode_and_quality(oriented_seq, oriented_quality, barcode_regex)
            if barcode is None:
                no_barcode += 1
                continue
            orfs.append(oriented_seq)
            barcodes.append(barcode)
            writer.writerow([barcode, barcode_quality, oriented_seq, oriented_quality])
            passed += 1

    print("Total records processed:", total)
    print("ORFs found:", total - no_orf)
    print("ORFs missing:", no_orf)
    print("Barcodes missing (after ORF):", no_barcode)
    print("Successfully written:", passed)
    print("Final ORF count:", len(orfs))
    print("Final Barcode count:", len(barcodes))


if __name__ == "__main__":
    fastq_file = sys.argv[1]
    output_csv = sys.argv[2]
    orf_regex = sys.argv[3].strip()  # Provide as a raw regex string
    barcode_regex = sys.argv[4].strip()  # Provide as a raw regex string
    process_fastq(fastq_file, output_csv, orf_regex, barcode_regex)

