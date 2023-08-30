import seqcurate as sc
import pandas as pd

# Get the sequences that remain after clustering
seqs_df = sc.get_sequence_df(snakemake.input.fasta)

# Get the sequences we have designated as key sequences
key_df = sc.get_sequence_df(snakemake.input.key_sequences)
