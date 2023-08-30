from Bio import AlignIO
import pandas as pd

# Read the multiple sequence alignment from the file
aln = AlignIO.read(snakemake.input.aln, format="fasta")

# Extract the order of sequences from the alignment
msa_order = [record.id for record in aln]

# Read the DataFrame from the CSV file
df = pd.read_csv(snakemake.input.csv)

df.set_index('info', inplace=True)


# # Create a dictionary to map sequence identifiers to their position in MSA
# msa_order_dict = {seq_id: index for index, seq_id in enumerate(msa_order)}

# print (msa_order_dict)


# print (df.index)

# # Create a custom sorting index based on the MSA order
# sorting_index = df.index.map(lambda x: msa_order_dict.get(x, len(msa_order)))

# print (sorting_index)

# Reorder the DataFrame based on the custom sorting index
df_reordered = df.reindex(msa_order)

print (df_reordered.index)

df_reordered.to_csv(snakemake.output.csv, index=False)
