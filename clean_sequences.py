import re
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from concurrent.futures import ThreadPoolExecutor

# Function to clean headers
def clean_header(header):
    header = header.replace(' ', '_')  # Replace spaces with underscores first
    cleaned_header = re.sub(r'[^a-zA-Z0-9_]', '_', header)  # Remove non-standard characters
    return cleaned_header

# Function to clean sequences
def clean_sequence(seq, is_protein):
    if is_protein:
        valid_chars = set("ACDEFGHIKLMNPQRSTVWY")  # All standard amino acids
        replace_char = 'X'
    else:
        valid_chars = set("ACGTN")  # DNA bases
        replace_char = 'N'
    return ''.join([char if char in valid_chars else replace_char for char in seq])

# Function to process a single record
def process_record(record, is_protein):
    cleaned_seq = clean_sequence(str(record.seq), is_protein)
    record.seq = Seq(cleaned_seq)
    return record

# Set up argument parser
parser = argparse.ArgumentParser(description="Clean FASTA sequences and headers.")
parser.add_argument("-i", "--input", required=True, help="Input FASTA file")
parser.add_argument("-o", "--output", required=True, help="Output FASTA file")
parser.add_argument("-t", "--threads", type=int, default=4, help="Number of threads to use (default: 4)")
parser.add_argument("--protein", action='store_true', help="Specify if the input is protein sequences")

args = parser.parse_args()

# Process the FASTA file with multithreading
with open(args.output, "w") as output_handle:
    with ThreadPoolExecutor(max_workers=args.threads) as executor:
        records = list(SeqIO.parse(args.input, "fasta"))

        # First clean the headers
        for record in records:
            cleaned_header = clean_header(record.description)
            record.id = cleaned_header
            record.name = cleaned_header
            record.description = cleaned_header

        # Then clean the sequences in parallel
        cleaned_records = list(executor.map(lambda r: process_record(r, args.protein), records))

        SeqIO.write(cleaned_records, output_handle, "fasta")

print("Sequences and headers cleaned successfully!")

