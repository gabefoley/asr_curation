import annot_functions as an
import seqcurate as sc
import pandas as pd

annot_df = pd.read_csv(snakemake.input[0])
annot_df.to_csv(snakemake.output[0], index=False)

