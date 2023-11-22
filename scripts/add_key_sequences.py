import seqcurate as sc
import pandas as pd

def main():
    input_fasta = snakemake.input.fasta
    output_fasta = snakemake.output.fasta
    input_key_sequences = snakemake.input.key_sequences
    input_align_csv = snakemake.input.align_csv
    input_full_csv = snakemake.input.full_csv
    output_csv = snakemake.output.csv

    print("inside add key seqs")
    print(input_fasta)
    print(output_fasta)
    print(input_key_sequences)

    # Get the sequences that remain after clustering
    seqs_df = sc.get_sequence_df(input_fasta)

    # Get the sequences we have designated as key sequences
    key_df = sc.get_sequence_df(input_key_sequences)

    # Get the subset of sequence info
    align_df = pd.read_csv(input_align_csv)

    # Get the full set of sequence info
    full_df = pd.read_csv(input_full_csv)

    if "info" in key_df.columns:
        key_df_info = key_df["info"].values
        seqs_df_info = seqs_df["info"].values

        # If we find a sequence in the key sequences that isn't in the sequences that remain, add it back in
        add_info = [x for x in key_df_info if x not in seqs_df_info]
        add_df = full_df[full_df["info"].isin(add_info)]

        subset_df = align_df[align_df["info"].isin(seqs_df_info)]

        print("len of subset df:", len(subset_df))
        print("len of add df:", len(add_df))

        concat_df = pd.concat([subset_df, add_df])

        print("len of concat df:", len(concat_df))

        concat_df = concat_df.drop_duplicates(subset="info", keep="first")
        print("len of deduplicated concat df:", len(concat_df))

        sc.write_to_fasta(concat_df, output_fasta, trim=True)

        # Write the subset to its own CSV file
        concat_df.to_csv(output_csv, index=False)
    else:  # No sequences found in the key sequence file, just write out the subset again
        add_df = full_df[full_df["info"].isin(seqs_df["info"].values)]

        sc.write_to_fasta(seqs_df, output_fasta, trim=True)
        add_df.to_csv(output_csv, index=False)

if __name__ == "__main__":
    main()
