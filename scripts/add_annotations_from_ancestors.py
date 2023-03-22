import annot_functions as an
import pandas as pd
import seqcurate as sc
from Bio import AlignIO

align_df = pd.read_csv(snakemake.input.csv)

anc_df = sc.get_sequence_df(snakemake.input.aln, alignment=True, ancestor=True)

aln = AlignIO.read(snakemake.input.aln, format="fasta")

aln_dict = {seq.name: str(seq.seq) for seq in aln}

# IF CUSTOMISING THIS FILE PLACE CUSTOM CODE HERE

align_df.reset_index(inplace=True, drop=True)
anc_df.reset_index(inplace=True, drop=True)

frames = [align_df, anc_df]
merge_df = pd.concat(frames)

merge_df.to_csv(snakemake.output.csv, index=False)
