# seq_programs

Here are some programs that are useful for sequence manipulation.

clean_seqeunces.py: replaces all spaces and special characters in fasta header with '_', and any ambiguous nucleotide or amino acid with N or X respectively

extract_beds.py: to determine the BED coordinates for slicing of fasta, useful for in silico tests.

merge_contigs.py: merges contigs from a fasta file into a single sequence, with the number of Ns that can be specfied
