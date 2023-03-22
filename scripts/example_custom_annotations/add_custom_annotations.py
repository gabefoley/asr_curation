"""Adds custom annotations to an annotation dataframe."""

import scripts.annot_functions as an
import pandas as pd

print ("Adding custom annotations")

annot_df = pd.read_csv(snakemake.input[0])

# Check if the KARI EC number is in the dataset we're currently looking at
if "ec_1_1_1_86" in snakemake.wildcards.dataset:
    print("Adding KARI specific annotations")
    print("Adding KARI Class")
    annot_df["KARI_Class"] = annot_df.apply(
        lambda row: an.classify_KARI(
            row["feature_count"] if pd.notnull(row["feature_count"]) else ""
        ),
        axis=1,
    )
# Create an image

annot_df.to_csv(snakemake.output[0], index=False)
