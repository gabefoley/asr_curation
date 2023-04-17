import seqcurate as sc
import pandas as pd

print("inside add key seqs")
print(snakemake.input.fasta)
print(snakemake.output)

print(snakemake.input.key_sequences)


# Get the sequences that remain after clustering
seqs_df = sc.get_sequence_df(snakemake.input.fasta)

# Get the sequences we have designated as key sequences
key_df = sc.get_sequence_df(snakemake.input.key_sequences)

# Get the subset of sequence info
subset_df = pd.read_csv(snakemake.input.subset_csv)

# Get the full set of sequence info
full_df = pd.read_csv(snakemake.input.full_csv)

key_df_info = key_df["info"].values
print("key df info")
print(key_df_info)

# seqs_df_info = seqs_df["info"].values
# print("seqs df info")
# print(seqs_df_info)

# If we find a sequence in the key sequences that isn't in the sequences that remain we need to add it back in
add_info = [x for x in key_df_info if x not in seqs_df_info]
print("add info")
print(add_info)

add_df = full_df[full_df["info"].isin(add_info)]

concat_df = pd.concat([subset_df, add_df])

concat_df = concat_df.drop_duplicates(subset="info", keep="first")

sc.write_to_fasta(concat_df, snakemake.output.fasta, trim=True)

print("caterpillar")

print(concat_df)


# Write the subset to its own csv file
concat_df.to_csv(snakemake.output.csv, index=False)

# print (halo)

# sc.write_to_fasta(sub_df, snakemake.output.fasta, trim=True)
