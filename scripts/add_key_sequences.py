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
align_df = pd.read_csv(snakemake.input.align_csv)

# Get the full set of sequence info
full_df = pd.read_csv(snakemake.input.full_csv)

if "info" in key_df.columns:
    key_df_info = key_df["info"].values
    # print("key df info")
    # print(key_df_info)

    seqs_df_info = seqs_df["info"].values
    # print("seqs df info")
    # print(seqs_df_info)

    # If we find a sequence in the key sequences that isn't in the sequences that remain we need to add it back in
    add_info = [x for x in key_df_info if x not in seqs_df_info]
    # print("add info")
    # print(add_info)

    add_df = full_df[full_df["info"].isin(add_info)]

    subset_df = align_df[align_df["info"].isin(seqs_df_info)]

    print("len of subset df")

    print(len(subset_df))

    print("len of add df")

    print(len(add_df))

    concat_df = pd.concat([subset_df, add_df])

    print("len of concat df")

    print(len(concat_df))

    concat_df = concat_df.drop_duplicates(subset="info", keep="first")

    print(len(concat_df))

    sc.write_to_fasta(concat_df, snakemake.output.fasta, trim=True)

    # Write the subset to its own csv file
    concat_df.to_csv(snakemake.output.csv, index=False)

else:  # No sequences found in the key sequence file, just write out the subset again
    add_df = full_df[full_df["info"].isin(seqs_df["info"].values)]

    sc.write_to_fasta(seqs_df, snakemake.output.fasta, trim=True)

    add_df.to_csv(snakemake.output.csv, index=False)
# print (halo)

# sc.write_to_fasta(sub_df, snakemake.output.fasta, trim=True)
