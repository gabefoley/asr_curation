''' generic python scripts to add common annotations before alignment of sequences '''

import annot_functions as an
import seqcurate as sc
import pandas as pd

annot_df = pd.read_csv(snakemake.input[0])
annot_df = an.annotate_nonAA(annot_df)
annot_df = an.annotate_AA(annot_df)
annot_df = an.annotate_sp_tr(annot_df)
# annot_df = an.annotate_fragment(annot_df)

annot_df.to_csv(snakemake.output[0], index=False)
