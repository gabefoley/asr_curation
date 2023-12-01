from Bio import AlignIO
import pandas as pd

# Read the multiple sequence alignment from the file
aln = AlignIO.read(snakemake.input.aln, format="fasta")

# Extract the order of sequences from the alignment
msa_order = [record.id for record in aln]

# Read the DataFrame from the CSV file
df = pd.read_csv(snakemake.input.csv)

df.set_index("info", inplace=True)

df_reordered = df.reindex(msa_order)

print(df_reordered.index)

df_reordered.to_csv(snakemake.output.csv, index=False)
