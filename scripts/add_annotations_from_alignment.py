import pandas as pd
from Bio import AlignIO
import seqcurate as sc

def main():
    df = pd.read_csv(snakemake.input.csv)

    # Get a sequence DataFrame using seqcurate
    align_df = sc.get_sequence_df(snakemake.input.aln, alignment=True, ancestor=True)

    # Merge the data from the CSV file and the sequence DataFrame based on the "info" column
    merged_df = pd.merge(
        df,
        align_df,
        left_on=["info"],
        right_on=["info"],
        suffixes=["", "_r"],
    )

    # Write the merged DataFrame to a CSV file
    merged_df.to_csv(snakemake.output.csv, index=False)

# Check if this script is being run as the main program
if __name__ == "__main__":
    main()
