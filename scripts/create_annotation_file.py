import pandas as pd


def create_annotations(df, annotation_cols, outpath):
    # Accession needs to be the first column so if it isn't requested, add it in.
    # if annotation_cols[0] != "truncated_info":
    #     annotation_cols = ["truncated_info"] + annotation_cols

    if annotation_cols[0] != "extracted_id":
        annotation_cols = ["extracted_id"] + annotation_cols

    # Subset the columns
    subset_df = df[[x for x in annotation_cols if x in df.columns]]
    subset_df = subset_df.fillna("None")

    # Write out the annotation columns file
    subset_df.to_csv(outpath, sep="\t", index=False)


if __name__ == "__main__":
    df = pd.read_csv(snakemake.input.csv)
    annotation_cols = snakemake.params.annotation_cols
    outpath = snakemake.output.tsv

    create_annotations(df, annotation_cols, outpath)
