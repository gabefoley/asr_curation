''' generic python scripts to add common annotations before alignment of sequences '''

import annot_functions as an
import pandas as pd

annot_df = pd.read_csv(snakemake.input[0])
annot_df = an.annotate_nonAA(annot_df)
annot_df = an.annotate_AA(annot_df)
annot_df = an.annotate_sp_tr(annot_df)
annot_df = an.add_thermo(annot_df, "./additional_data/thermophilic_bacteria.txt")

# Add an additional Length column (this can be removed once we can call more complex, multiple queries on a single column)
annot_df["Length_2"] = annot_df["length"]

annot_df.to_csv(snakemake.output[0], index=False)
