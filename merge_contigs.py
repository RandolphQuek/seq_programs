import argparse
import os
from concurrent.futures import ThreadPoolExecutor

def is_fasta_format(file_path):
    """
    Check if the file at 'file_path' is in FASTA format.
    
    Parameters:
    - file_path: Path to the file to check.
    
    Returns:
    - True if the file is in FASTA format, False otherwise.
    """
    with open(file_path, "r") as f:
        line = f.readline().strip()
        if line.startswith(">"):
            return True
    return False

def clean_sequence(sequence):
    """
    Clean up sequence by removing non-standard characters and returning cleaned sequence.
    
    Parameters:
    - sequence: Sequence string to be cleaned.
    
    Returns:
    - Cleaned sequence string.
    """
    # Define non-standard characters to be removed
    non_standard_chars = set(";\"',.!?()[]{}<>")

    # Remove non-standard characters from sequence
    cleaned_sequence = ''.join([char for char in sequence if char not in non_standard_chars])
    
    return cleaned_sequence

def merge_contigs(input_file, output_file, header="Fasta_header", n_count=50):
    """
    Merge contig sequences from input FASTA file into one long sequence with 'N's between each contig.
    
    Parameters:
    - input_file: Path to the input FASTA file containing contigs.
    - output_file: Path to the output file where the combined sequence will be saved.
    - header: Header for the combined sequence in the output FASTA file. Default is "Fasta_header".
    - n_count: Number of 'N's to insert between merged contigs. Default is 50.
    """
    # Check if input and output files are the same
    if os.path.abspath(input_file) == os.path.abspath(output_file):
        raise ValueError("Error: Input and output files cannot have the same name. Please specify different file names.")

    # Check if input file exists
    if not os.path.isfile(input_file):
        raise FileNotFoundError(f"Error: Input file '{input_file}' not found. Use '-h' or '--help' for more information.")

    # Check if input file is in FASTA format
    if not is_fasta_format(input_file):
        raise ValueError(f"Error: Input file '{input_file}' is not in FASTA format. Use '-h' or '--help' for more information.")

    # Initialize an empty list to store sequences
    sequences = []

    # Open and read the input FASTA file
    with open(input_file, "r") as f:
        lines = f.readlines()
        current_sequence = ""
        for line in lines:
            line = line.strip()
            if line.startswith(">"):
                # If a header line is encountered, append the current sequence (if any) to the sequences list
                if current_sequence:
                    sequences.append(current_sequence)
                    current_sequence = ""
            else:
                # Concatenate the sequence lines
                current_sequence += line

        # Append the last sequence (since the file may end without another header)
        if current_sequence:
            sequences.append(current_sequence)

    # Clean each sequence to remove non-standard characters
    cleaned_sequences = [clean_sequence(seq) for seq in sequences]

    # Initialize the combined sequence with the first cleaned sequence
    combined_sequence = cleaned_sequences[0]

    # Append 'N's and subsequent cleaned sequences
    for seq in cleaned_sequences[1:]:
        combined_sequence += "N" * n_count + seq

    # Write the combined sequence to the output FASTA file with the specified header
    with open(output_file, "w") as out_f:
        out_f.write(f">{header}\n")  # Header for the combined sequence (user-specified)
        # Write the sequence in lines of 80 characters for better readability (optional)
        for i in range(0, len(combined_sequence), 80):
            out_f.write(combined_sequence[i:i+80] + "\n")

def main():
    # Command-line arguments and help message
    parser = argparse.ArgumentParser(description="Merge contigs from a FASTA file into one long sequence.")
    parser.add_argument("--input", "-i", required=True, help="Input FASTA file containing contigs")
    parser.add_argument("--output", "-o", required=True, help="Output file to save the combined sequence")
    parser.add_argument("--header", default="Fasta_header", help="Header for the combined sequence in the output FASTA file (default: 'Fasta_header')")
    parser.add_argument("--threads", type=int, default=1, help="Number of threads to use for merging contigs (default: 1)")
    parser.add_argument("--n-count", type=int, default=50, help="Number of 'N's to insert between merged contigs (default: 50). Use 0 to omit 'N's.")
    
    # Parse command-line arguments
    args = parser.parse_args()

    try:
        # Check if input and output files are the same
        if os.path.abspath(args.input) == os.path.abspath(args.output):
            raise ValueError("Error: Input and output files cannot have the same name. Please specify different file names.")
        
        # Check if input file exists
        if not os.path.isfile(args.input):
            raise FileNotFoundError(f"Error: Input file '{args.input}' not found. Use '-h' or '--help' for more information.")
        
        # Check if input file is in FASTA format
        if not is_fasta_format(args.input):
            raise ValueError(f"Error: Input file '{args.input}' is not in FASTA format. Use '-h' or '--help' for more information.")
        
        # Use ThreadPoolExecutor for multithreading with specified number of threads
        with ThreadPoolExecutor(max_workers=args.threads) as executor:
            # Call merge_contigs function with provided input and output files, header, and n_count
            executor.submit(merge_contigs, args.input, args.output, args.header, args.n_count)
        
        print(f"Successfully merged contigs from '{args.input}' into '{args.output}' with header '{args.header}' using {args.threads} threads.")
    except ValueError as e:
        print(e)
    except FileNotFoundError as e:
        print(e)
        print("Please make sure to specify a valid input FASTA file. Use '-h' or '--help' for more information.")
    except Exception as e:
        print(f"Error: {str(e)}")
        print("An unexpected error occurred. Use '-h' or '--help' for more information.")

if __name__ == "__main__":
    main()

