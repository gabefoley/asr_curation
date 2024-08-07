""" Python code to fetch annotations from uniprot in batches, use pagination to process results and creating .csv output file """

import pandas as pd
import re
import requests
from requests.adapters import HTTPAdapter, Retry
from configs.uniprot_cols import full_uniprot_cols, reduced_uniprot_cols
import numpy as np
import os

# config parameters
API_URL = "https://rest.uniprot.org/uniprotkb/search?"
retries = Retry(total=5, backoff_factor=0.25, status_forcelist=[500, 502, 503, 504])
session = requests.Session()
session.mount("https://", HTTPAdapter(max_retries=retries))
re_next_link = re.compile(r'<(.+)>; rel="next"')


def get_next_link(headers):
    """get next batch link"""

    if "Link" in headers:
        match = re_next_link.match(headers["Link"])
        if match:
            return match.group(1)


def get_batch(batch_url):
    """get next batch"""

    while batch_url:
        response = session.get(batch_url)
        response.raise_for_status()
        total = response.headers["x-total-results"]
        yield response, total
        batch_url = get_next_link(response.headers)


def split_lineage(x):
    """function to split the lineage strings into seperate columns for use later"""

    lineage_order = [
        "lineage_superkingdom",
        "lineage_kingdom",
        "lineage_subkingdom",
        "lineage_superphylum",
        "lineage_phylum",
        "lineage_subphylum",
        "lineage_superclass",
        "lineage_class",
        "lineage_subclass",
        "lineage_infraclass",
        "lineage_superorder",
        "lineage_order",
        "lineage_suborder",
        "lineage_infraorder",
        "lineage_parvorder",
        "lineage_superfamily",
        "lineage_family",
        "lineage_subfamily",
        "lineage_tribe",
        "lineage_subtribe",
        "lineage_genus",
        "lineage_subgenus",
        "lineage_species_group",
        "lineage_species_group",
        "lineage_species_group",
        "lineage_varietas",
        "lineage_forma",
    ]

    # Split lineage. Return empty string if linage info not present
    if pd.isna(x):
        lin_split = ""
    else:
        lin_split = x.split(",")

    lineage_dict = {}
    lineage_list = []

    for lin in lin_split:
        if "no rank" not in lin and "(" in lin and ")" in lin:
            val = lin.split("(")[0].strip()
            key = "lineage_" + lin.split("(")[1].split(")")[0]
            lineage_dict[key] = val

    # return values in the order
    for l_o in lineage_order:
        try:
            lineage_list.append(lineage_dict[l_o])
        except:
            lineage_list.append("")

    return lineage_list

def split_lineage(lineage):
    # Placeholder function for splitting the lineage string.
    # Replace this with your actual split logic.
    return lineage.split(';') if pd.notnull(lineage) else [None] * 25

def process_results(intermediate_tsv_file):
    """process results in tsv file"""

    annot_df = pd.read_csv(intermediate_tsv_file, sep="\t")

    # replace gaps in column names with '_'
    annot_df.columns = annot_df.columns.str.replace(" ", "_")
    annot_df.columns = annot_df.columns.str.replace("-", "_")
    annot_df.columns = annot_df.columns.str.replace("[()]", "", regex=False)
    print("print the columns")

    print("print the columns")
    for x in annot_df.columns:
        print(x)

    # Separate lineage columns (after the new API changes)
    lineage_columns = [
        "lineage_superkingdom",
        "lineage_kingdom",
        "lineage_subkingdom",
        "lineage_superphylum",
        "lineage_phylum",
        "lineage_subphylum",
        "lineage_superclass",
        "lineage_class",
        "lineage_subclass",
        "lineage_infraclass",
        "lineage_superorder",
        "lineage_order",
        "lineage_suborder",
        "lineage_infraorder",
        "lineage_parvorder",
        "lineage_superfamily",
        "lineage_family",
        "lineage_subfamily",
        "lineage_tribe",
        "lineage_subtribe",
        "lineage_genus",
        "lineage_subgenus",
        "lineage_species_group",
        "lineage_varietas",
        "lineage_forma",
    ]

    # Split lineage and create new columns
    annot_df[lineage_columns] = annot_df["lineage"].apply(split_lineage).apply(pd.Series)

    return annot_df


    # for x in annot_df.columns:
    #     print(x)
    # # seperate lineage columns ( after the new api changes)
    # (
    #     annot_df["lineage_superkingdom"],
    #     annot_df["lineage_kingdom"],
    #     annot_df["lineage_subkingdom"],
    #     annot_df["lineage_superphylum"],
    #     annot_df["lineage_phylum"],
    #     annot_df["lineage_subphylum"],
    #     annot_df["lineage_superclass"],
    #     annot_df["lineage_class"],
    #     annot_df["lineage_subclass"],
    #     annot_df["lineage_infraclass"],
    #     annot_df["lineage_superorder"],
    #     annot_df["lineage_order"],
    #     annot_df["lineage_suborder"],
    #     annot_df["lineage_infraorder"],
    #     annot_df["lineage_parvorder"],
    #     annot_df["lineage_superfamily"],
    #     annot_df["lineage_family"],
    #     annot_df["lineage_subfamily"],
    #     annot_df["lineage_tribe"],
    #     annot_df["lineage_subtribe"],
    #     annot_df["lineage_genus"],
    #     annot_df["lineage_subgenus"],
    #     annot_df["lineage_species_group"],
    #     annot_df["lineage_species_group"],
    #     annot_df["lineage_species_group"],
    #     annot_df["lineage_varietas"],
    #     annot_df["lineage_forma"],
    # ) = zip(*annot_df["lineage"].map(split_lineage))
    #
    # return annot_df


def get_uniprot_annotations(
    uniprot_list_ids,
    uniprot_cols,
    intermediate_tsv_file,
    id_batch_size,
    result_batch_size,
    verbose,
):
    """running uniprot ids batches"""


    splits = np.array_split(
        np.array(uniprot_list_ids), max(1, round(len(uniprot_list_ids) / id_batch_size))
    )
    split_count = 0

    # For each split go to UniProt to retrieve information
    for split in splits:
        split_count += 1
        print ('here is split')
        print (split)
        if verbose:
            print(f"Running stage {split_count} of {len(splits)}")
        split_fetch_success = fetch_annotations_in_batches(
            split, uniprot_cols, intermediate_tsv_file, result_batch_size, verbose
        )

        if split_fetch_success == 1:
            if verbose:
                print(f"Successfully processed {split_count} of {len(splits)}")
        else:
            if verbose:
                print(f"Failed to process {split_count} of {len(splits)}")
            return 0

    return 1


def fetch_annotations_in_batches(
    ids, uniprot_cols, intermediate_tsv_file, result_batch_size, verbose
):
    """process uniprot api results in batches and load it in .tsv file on the disk"""

    cols = ",".join(uniprot_cols)

    batch_url = (
        API_URL
        + "format=tsv&"
        + "query=accession:("
        + " OR ".join(ids)
        + ")&fields="
        + cols
        + "&size="
        + str(result_batch_size)
    )

    progress = 0

    with open(intermediate_tsv_file, "a") as f:
        for batch, total in get_batch(batch_url):
            lines = batch.text.splitlines()
            for line in lines[1:]:
                print(line, file=f)
            progress += len(lines[1:])
            if verbose:
                print(f"{progress} / {total}")

    if int(progress) == int(total):  # All have been fetched
        if verbose:
            print("Created and loaded annotations in .tsv file")
        return 1
    else:
        return 0


def get_uniprot_id_list(df):
    """Function to get list of uniprot ids to make an api call to uniprot"""

    uniprot_df_list = df["UNIPROT"].tolist()
    uniprot_list_ids = []
    for l in uniprot_df_list:
        l_item = (str(l[1:-1])).replace(" ", "")
        for m in l_item.split(","):
            uniprot_list_ids.append(m[1:-1])
    return [x for x in uniprot_list_ids if x]


def uniprot_annotation_lkp(
    input_file,
    output_file,
    uniprot_cols,
    intermediate_tsv_file,
    id_batch_size,
    result_batch_size,
    verbose,
):
    """main function to get uniprot ids annotation"""

    # get all uniprot ids
    validated_df = pd.read_csv(input_file)
    print ('validated df is')
    print (validated_df.head(n=5))
    uniprot_list_ids = get_uniprot_id_list(validated_df)

    print ('ids here are')
    print (uniprot_list_ids)

    # create a new intermediate tsv file with header row for dumping results
    header = "\t".join(uniprot_cols)
    with open(intermediate_tsv_file, "w") as f:
        print(header, file=f)

    # submit all ids to uniprot for getting annotation and load in .tsv file
    results_uniprot_loaded = get_uniprot_annotations(
        uniprot_list_ids,
        uniprot_cols,
        intermediate_tsv_file,
        id_batch_size,
        result_batch_size,
        verbose,
    )

    # process results
    if results_uniprot_loaded:
        results_df = process_results(intermediate_tsv_file)
    else:
        print("Uniprot Results fetched failed.")

    # save in .csv file
    if verbose:
        print("Creating .csv output file")

    print(results_df)

    # We merge on sequence if it exists in both (it won't exist if we started with just a list of IDs)
    if "sequence" in validated_df.columns and "sequence" in results_df.columns:
        merged_df = validated_df.merge(
            results_df, left_on="sequence", right_on="sequence", how="left"
        )
    else:
        merged_df = validated_df.merge(
            results_df, left_on="info", right_on="accession", how="left"
        )

    merged_df = merged_df.drop_duplicates(subset="info", keep="first")
    merged_df.to_csv(output_file, index=False)

    if verbose:
        print("Completed UniProt annotations retrieval\n")

    if os.path.exists(intermediate_tsv_file):
        os.remove(intermediate_tsv_file)


def main():
    input_file = snakemake.input.csv
    output_file = snakemake.output[0]
    intermediate_tsv_file = snakemake.input.csv.split(".")[0] + ".tsv"
    uniprot_col_size = snakemake.params.uniprot_col_size
    verbose = snakemake.params.verbose

    print ('input file is ')
    print (input_file)

    if uniprot_col_size == "full":
        if verbose:
            print("Retrieving all of the UniProt columns")

        uniprot_cols = full_uniprot_cols
    elif uniprot_col_size == "reduced":
        print("Retrieving a reduced set of the UniProt columns")

        uniprot_cols = reduced_uniprot_cols

    if verbose:
        print("UniProt columns to retrieve - ")
        print(uniprot_cols)
        print("\nStarting to retrieve UniProt annotations")

    id_batch_size = 110  # process ids in batches of this parameter
    result_batch_size = 50  # process results in this batches (uses pagination)

    uniprot_annotation_lkp(
        input_file,
        output_file,
        uniprot_cols,
        intermediate_tsv_file,
        id_batch_size,
        result_batch_size,
        verbose,
    )


if __name__ == "__main__":
    main()
