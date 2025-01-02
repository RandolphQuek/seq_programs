import os
import pandas as pd
from Bio import SeqIO

# Directories
blastn_dir = "path_to_blastn_results"
bed_dir = "path_to_bed_files"
genomes_dir = "path_to_genomes"

# Function to get the length of each contig in a genome file
def get_contig_lengths(genome_file):
    contig_lengths = {}
    for record in SeqIO.parse(genome_file, "fasta"):
        contig_lengths[record.id] = len(record.seq)
    return contig_lengths

# Loop through each BLAST result file and create BED files
for blast_file in os.listdir(blastn_dir):
    if blast_file.endswith("_out_filtered.txt"):
        base_name = blast_file.replace("_out_filtered.txt", "")
        genome_file = os.path.join(genomes_dir, f"{base_name}.fna")
        
        # Check if the genome file exists
        if not os.path.exists(genome_file):
            print(f"Genome file {genome_file} not found. Skipping.")
            continue
        
        # Get the contig lengths
        contig_lengths = get_contig_lengths(genome_file)
        
        # Read the BLAST results into a pandas DataFrame
        blast_df = pd.read_csv(os.path.join(blastn_dir, blast_file), sep='\t', header=None)
        
        # Adjust the coordinates and create BED file
        bed_lines = []
        for _, row in blast_df.iterrows():
            contig = row[1]  # Contig name in column 2
            start = row[8]   # Start position in column 9
            end = row[9]     # End position in column 10
            
            # Validate the start and end positions
            try:
                start = int(start)
                end = int(end)
            except ValueError:
                print(f"Invalid start/end values in line: {row}. Skipping.")
                continue

            # Extend the region by 150 bp on each side (300 bp total)
            extension = 150
            
            # Adjust the start and end positions
            adjusted_start = start - extension
            adjusted_end = end + extension
            
            # Ensure the start and end positions are within the contig boundaries
            adjusted_start = max(1, adjusted_start)  # Cannot be less than 1
            adjusted_end = min(contig_lengths[contig], adjusted_end)  # Cannot be more than the contig length

            # Adjust for 0-based indexing in BED format
            bed_lines.append([contig, adjusted_start - 1, adjusted_end])  # BED is 0-based for start

        # Write the BED file
        bed_df = pd.DataFrame(bed_lines, columns=["contig", "start", "end"])
        bed_df.to_csv(os.path.join(bed_dir, f"{base_name}.bed"), sep='\t', header=False, index=False)

        print(f"Processed {blast_file} and generated BED file.")
