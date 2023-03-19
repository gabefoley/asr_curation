import annot_functions as an
import pandas as pd
from Bio import AlignIO
import seqcurate as sc

df = pd.read_csv(snakemake.input.csv)

aln = AlignIO.read(snakemake.input.aln, format="fasta")
aln_dict = {seq.name: str(seq.seq) for seq in aln}

align_df = sc.get_sequence_df(snakemake.input.aln, alignment=True, ancestor=True)

merged_df = pd.merge(
        df,
        align_df,
        left_on=["accession"],
        right_on=["accession"],
        suffixes=["", "_r"],
    )

merged_df.to_csv(snakemake.output.csv, index=False)


