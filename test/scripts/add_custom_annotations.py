"""Adds custom annotations to an annotation dataframe."""

import annot_functions as an
import pandas as pd

annot_df = pd.read_csv(snakemake.input[0])

# IF CUSTOMISING THIS FILE PLACE CUSTOM CODE HERE

annot_df.to_csv(snakemake.output[0], index=False)
