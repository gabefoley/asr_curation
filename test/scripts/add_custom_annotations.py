"""Adds custom annotations to an annotation dataframe."""

import pandas as pd


def main():
    annot_df = pd.read_csv(snakemake.input[0])

    # IF CUSTOMIZING THIS FILE, PLACE CUSTOM CODE HERE

    annot_df.to_csv(snakemake.output[0], index=False)


if __name__ == "__main__":
    main()
