import scripts.annot_functions as an
import pandas as pd
from Bio import AlignIO
import scripts.seqcurate as sc
import matplotlib.pyplot as plt

df = pd.read_csv(snakemake.input.csv)

aln = AlignIO.read(snakemake.input.aln, format="fasta")
aln_dict = {seq.name: str(seq.seq) for seq in aln}

align_df = sc.get_sequence_df(snakemake.input.aln, alignment=True, ancestor=True)


# Create an image of a fingerprint on an alignment

# Create an image of a domain on an alignment
if "ec_1_1_1_86" in snakemake.wildcards.dataset:
    outpath = f"./workflows/example_workflow/datasets/{snakemake.wildcards.dataset}/subsets/{snakemake.wildcards.subset}/plot"

    for col in [
        "lineage_superkingdom",
        "KARI_Class",
    ]:
        if col in df:
            fig, ax = plt.subplots(figsize=(len(col) / 2, 10))
            chart = df[col].value_counts().plot.barh(title=col, ax=ax)
            plt.savefig(f"{outpath}_{col}.png")

    # print (align_df.columns)

    df = pd.merge(
        df,
        align_df,
        left_on=["accession"],
        right_on=["accession"],
        suffixes=["", "_r"],
    )

    df["combined_domain_bounds"] = df.apply(
        lambda x: an.create_domain_bounds(x.ft_domain), axis=1
    )

    df.set_index("accession")
    combined_boundary_dict = pd.Series(
        df.combined_domain_bounds.values, index=df.accession
    ).to_dict()

    outpath = f"./workflows/example_workflow/datasets/{snakemake.wildcards.dataset}/subsets/{snakemake.wildcards.subset}/domain.html"

    colour_dict = {
        "KARI N-terminal Rossmann": "lawngreen",
        "KARI C-terminal knotted": "mediumpurple",
        "KARI C-terminal knotted 1": "orange",
        "KARI C-terminal knotted 2": "lightblue",
    }

    an.create_annotated_alignment(
        df, combined_boundary_dict, outpath, colour_dict=colour_dict
    )


merged_df = pd.merge(
    df,
    align_df,
    left_on=["accession"],
    right_on=["accession"],
    suffixes=["", "_r"],
)

merged_df.to_csv(snakemake.output.csv, index=False)
