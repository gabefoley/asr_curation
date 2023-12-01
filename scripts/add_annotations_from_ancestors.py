import pandas as pd
import seqcurate as sc
from Bio import AlignIO


def main():
    align_df = pd.read_csv(snakemake.input.csv)

    anc_df = sc.get_sequence_df(snakemake.input.aln, alignment=True, ancestor=True)

    # IF CUSTOMIZING THIS FILE, PLACE CUSTOM CODE HERE

    align_df.reset_index(inplace=True, drop=True)
    anc_df.reset_index(inplace=True, drop=True)

    frames = [align_df, anc_df]
    merge_df = pd.concat(frames)

    merge_df.to_csv(snakemake.output.csv, index=False)


if __name__ == "__main__":
    main()
