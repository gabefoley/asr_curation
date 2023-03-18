import pandas as pd
import json
from collections import defaultdict
import requests
from ast import literal_eval
from itertools import islice
import numpy as np

original_df = pd.read_csv(snakemake.input[0])
original_df.Genomes_with_identical_proteins = (
    original_df.Genomes_with_identical_proteins.apply(
        lambda x: literal_eval(x) if pd.notnull(x) else "not found"
    )
)

# assembly_dict = pd.Series(original_df.Genomes_with_identical_proteins.str[1:-1].str.split(','),index=original_df.Entry.values).to_dict()

assembly_dict_with_nulls = pd.Series(
    original_df.Genomes_with_identical_proteins.values, index=original_df.Entry.values
).to_dict()

assembly_dict = {}

for k, v in assembly_dict_with_nulls.items():
    if v != "not found":
        assembly_dict[k] = v


def chunks(data, SIZE=800):
    it = iter(data)
    for i in range(0, len(data), SIZE):
        yield {k: data[k] for k in islice(it, SIZE)}


def get_gtdb_dict(assembly_dict):
    gtdb_dict = defaultdict(list)
    count = 0

    for name, entry_list in assembly_dict.items():
        results = []

        print(name)
        print(entry_list)

        for entry in entry_list:

            print(name)
            print(entry)

            response = requests.get(
                f"https://gtdb.ecogenomic.org/api/v1/genome/summary/{entry}"
            )

            if response.status_code == 400:
                results.append("Not found")

            else:
                results.append(response.json())

        gtdb_dict[name] = results

    return gtdb_dict


def check_for_gtdb_taxa_discrepency(gtdb_dict, taxa_rank="gtdb_phylum"):

    gtdb_taxa_discrepency_dict = {}

    for prot_id, gtdb_genomes in gtdb_dict.items():

        print(gtdb_genomes)

        gtdb_genomes = [x for x in gtdb_genomes if x != "Not found"]

        gtdb_genomes_list = [
            gtdb_genome[taxa_rank].split("__")[1] for gtdb_genome in gtdb_genomes
        ]

        gtdb_genomes_set = set(gtdb_genomes_list)
        if len(gtdb_genomes_set) != 1:

            print("Multiple species detected")
            print(gtdb_genomes_set)
            gtdb_taxa_discrepency_dict[prot_id] = gtdb_genomes_list

        else:
            gtdb_taxa_discrepency_dict[prot_id] = "Consistent"

    return gtdb_taxa_discrepency_dict


def check_consistency_of_uniprot_gtdb_taxa(
    df, gtdb_dict, df_rank="Taxonomic_lineage_PHYLUM", gtdb_rank="gtdb_phylum"
):

    uniprot_gtdb_consistency_dict = {}

    for prot_id, gtdb_genomes in gtdb_dict.items():

        gtdb_genomes = [x for x in gtdb_genomes if x != "Not found"]

        print(gtdb_genomes)

        print(len(gtdb_genomes))

        # print (gtdb_genomes[0])

        # print(gtdb_genomes[0][gtdb_rank])

        check_df = df[df["Entry"] == prot_id]

        print("pi")
        print(prot_id)

        if len(gtdb_genomes) > 1:

            print(gtdb_genomes[0])

            print(gtdb_genomes[0][gtdb_rank])

            print("cdf")

            print(check_df)

            print(check_df[df_rank].values)
            print(check_df[df_rank].values[0])

            if check_df[df_rank].values[0] == gtdb_genomes[0][gtdb_rank].split("__")[1]:
                uniprot_gtdb_consistency_dict[prot_id] = "Consistent"

                continue
            else:
                uniprot_gtdb_consistency_dict[prot_id] = "Inconsistent"

        else:
            # print ('diff')
            # print (check_df[df_rank].values[0])
            # print (gtdb_genomes[0][gtdb_rank].split("__")[1])
            uniprot_gtdb_consistency_dict[prot_id] = "No_GTDB_info"

    return uniprot_gtdb_consistency_dict


def get_highest_checkM_completeness(gtdb_dict):

    highest_checkM_dict = {}

    for prot_id, gtdb_genomes in gtdb_dict.items():
        print("check it")
        print(gtdb_genomes)
        gtdb_genomes = [x for x in gtdb_genomes if x != "Not found"]

        if len(gtdb_genomes) > 1:

            highest_checkM_completeness = max(
                [gtdb_genome["checkm_completeness"] for gtdb_genome in gtdb_genomes]
            )
            print(highest_checkM_completeness)
            highest_checkM_dict[prot_id] = highest_checkM_completeness
        else:
            highest_checkM_dict[prot_id] = "No_GTDB_info"

    return highest_checkM_dict


def add_to_df_from_dict(df, dict_to_add, col_name, add_name="Entry"):

    if col_name not in df.columns:
        df[col_name] = np.nan

    df[col_name] = df.apply(
        lambda row: dict_to_add[row[add_name]]
        if row[add_name] in dict_to_add
        else row[col_name],
        axis=1,
    )

    return df


for chunked_dict in chunks(assembly_dict):

    # assembly_dict = {'A0A2V7DB71_9BACT': ['GCA_003220fff655.1'], 'A0A1Q1G3S7_9STAP': ['GCF_001989575.1']}

    gtdb_dict = get_gtdb_dict(chunked_dict)
    print(gtdb_dict)
    print("Check for multiple phyla")
    # gtdb_taxa_discrepency_dict = check_for_gtdb_taxa_discrepency(gtdb_dict)
    print("Check consistency errors")
    ## print (original_df.head(3))

    ## prot_id = 'A0A0R1XS52_9LACO'

    ## check_df = original_df[original_df['Entry'] == prot_id]

    ## print (check_df)

    print(original_df["Entry"])
    # uniprot_gtdb_consistency_dict = check_consistency_of_uniprot_gtdb_taxa(original_df, gtdb_dict)
    print("Check for low CheckM scores")
    # highest_checkM_dict = get_highest_checkM_completeness(gtdb_dict)

    original_df = add_to_df_from_dict(original_df, gtdb_dict, "GTDB_FULL_DICTS")
    # original_df = add_to_df_from_dict(original_df, gtdb_taxa_discrepency_dict, 'GTDB_taxa_discrepency')
    # original_df = add_to_df_from_dict(original_df, uniprot_gtdb_consistency_dict, 'Uniprot_GTDB_consistency')
    # original_df = add_to_df_from_dict(original_df, highest_checkM_dict, 'Highest_checkM_completeness')

# Save the merged annotations to a csv
original_df.to_csv(snakemake.output[0], index=False)
