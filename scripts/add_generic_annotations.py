import annot_functions as an
import pandas as pd

def main():
    # Load the input CSV file
    annot_df = pd.read_csv(snakemake.input[0])

    # Annotate non-AA sequences
    annot_df = an.annotate_nonAA(annot_df)

    # Annotate AA sequences
    annot_df = an.annotate_AA(annot_df)

    # Annotate Swiss-Prot and TrEMBL sequences
    annot_df = an.annotate_sp_tr(annot_df)

    # Add an additional Length column (temporary)
    annot_df["Length_2"] = annot_df["length"]

    # Save the annotated DataFrame to the output CSV file
    annot_df.to_csv(snakemake.output[0], index=False)

if __name__ == "__main__":
    main()
