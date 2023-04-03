import pandas as pd
import click
import os


@click.command()
@click.option("--df", help="Dataframe with all annotations")
@click.option("--annot", help="Columns to add to annotation file")
@click.option(
    "--outpath", default="./annotation_cols.txt", help="Outpath for annotation file"
)
def create_annotations(df, annot, outpath):
    if "SNAKEMAKE" in os.environ:
        df = pd.read_csv(snakemake.input.csv)
        annotation_cols = snakemake.params.annotation_cols
        outpath = snakemake.output.tsv

    else:
        df = pd.read_csv(df, dtype="object")

        with open(annot) as annot_file:
            annotation_cols = [line.strip() for line in annot_file]

    # Accession needs to be the first column so if it isn't requested, add it in.
    if annotation_cols[0] != "truncated_info":
        annotation_cols = ["truncated_info"] + annotation_cols

    # Subset the columns
    subset_df = df[[x for x in annotation_cols if x in df.columns]]
    subset_df = subset_df.fillna("None")

    # Write out the annotation columns file
    subset_df.to_csv(outpath, sep="\t", index=False)


if __name__ == "__main__":
    create_annotations()
