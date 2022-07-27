import pandas as pd

from ast import literal_eval


def check_for_gtdb_taxa_discrepency(gtdb_dict_list, taxa_rank="gtdb_phylum"):

    gtdb_taxa_discrepency_dict = {}

    gtdb_genomes_list = [
        gtdb_genome[taxa_rank].split("__")[1] for gtdb_genome in gtdb_dict_list
    ]

    gtdb_genomes_list = [x for x in gtdb_genomes_list if len(x) > 0]

    gtdb_genomes_set = set(gtdb_genomes_list)
    if len(gtdb_genomes_set) == 0:
        return None
    if len(gtdb_genomes_set) != 1:

        #         print ('Multiple species detected')
        #         print (gtdb_genomes_set)
        return gtdb_genomes_list

    else:
        return "Consistent"


def check_consistency_of_uniprot_gtdb_taxa(
    entry, gtdb_dict_list, df_rank="Taxonomic_lineage_PHYLUM", gtdb_rank="gtdb_phylum"
):

    alias_dict = {
        "Actinobacteria": "Actinobacteriota",
        "Planctomycetes": "Planctomycetota",
        "Nitrospirae": "Nitrospirota",
        "Spirochaetes": "Spirochaetota",
        "Acidobacteria": "Acidobacteriota",
        "Fibrobacteres": "Fibrobacterota",
        "Verrucomicrobia": "Verrucomicrobiota",
    }

    uniprot_gtdb_consistency_dict = {}

    gtdb_phyla = [x[gtdb_rank].split("__")[1] for x in gtdb_dict_list if len(x) > 0]

    gtdb_phyla = [x for x in gtdb_phyla if len(x) > 0]

    if len(gtdb_phyla) > 0:
        if entry.Taxonomic_lineage_PHYLUM == gtdb_phyla[0] or (
            (entry.Taxonomic_lineage_PHYLUM in alias_dict)
            and (
                alias_dict[entry.Taxonomic_lineage_PHYLUM].strip()
                == gtdb_phyla[0].strip()
            )
        ):
            return "Consistent"
        else:

            return gtdb_phyla
    else:
        return None


def get_highest_checkM_completeness(gtdb_dict_list):

    if len(gtdb_dict_list) > 0:
        highest_checkM_completeness = max(
            [gtdb_genome["checkm_completeness"] for gtdb_genome in gtdb_dict_list]
        )
        #         print (highest_checkM_completeness)
        return highest_checkM_completeness
    else:
        return None


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


def get_gtdb_rank(gtdb_dict_list, gtdb_rank="gtdb_phylum"):
    gtdb_phyla = [x[gtdb_rank].split("__")[1] for x in gtdb_dict_list if len(x) > 0]

    gtdb_phyla = [x for x in gtdb_phyla if len(x) > 0]

    if gtdb_phyla:

        most_freq_phyla = max(gtdb_phyla, key=gtdb_phyla.count)

        return most_freq_phyla
    else:
        return None


df = pd.read_csv(snakemake.input[0])

import numpy as np

first_dict = df.head(10)["GTDB_FULL_DICTS"]

check_m_dict = {}
taxa_discrepency_dict = {}
gtdb_uniprot_consistency_dict = {}
gtdb_rank_dict = {}

for row in df.itertuples():
    #     print (row.Entry)

    # gtdb_dict_list = [literal_eval(x.replace("[", "").replace("]", "") + "}").replace("}}", "}") for x in first_dict[1].split("},")]

    gtdb_dict_list = []

    if not pd.isna(row.GTDB_FULL_DICTS):

        #         print (row.GTDB_FULL_DICTS)

        for x in row.GTDB_FULL_DICTS.replace("'Not found',", "").split("},"):
            #             print()
            #             print (x)

            #     print (x.replace("[", "").replace("]", "") + "}")
            #     print()
            remove_brackets = x.replace("[", "").replace("]", "") + "}"

            remove_brackets = str(remove_brackets.replace("}}", "}"))

            #             print (remove_brackets)

            if remove_brackets.strip() != "'Not found'}":

                gtdb_dict_list.append(str(remove_brackets.replace("}}", "}")))

    #         for x in gtdb_dict_list:
    #             print (x)
    #             print()
    #     print (gtdb_dict_list)

    gtdb_dict_list = [literal_eval(x.strip()) for x in gtdb_dict_list]

    check_m_highest = get_highest_checkM_completeness(gtdb_dict_list)
    taxa_discrepency = check_for_gtdb_taxa_discrepency(gtdb_dict_list)
    gtdb_uniprot_consistency = check_consistency_of_uniprot_gtdb_taxa(
        row, gtdb_dict_list
    )
    gtdb_rank = get_gtdb_rank(gtdb_dict_list)
    check_m_dict[row.Entry] = check_m_highest
    taxa_discrepency_dict[row.Entry] = taxa_discrepency
    gtdb_uniprot_consistency_dict[row.Entry] = gtdb_uniprot_consistency
    gtdb_rank_dict[row.Entry] = gtdb_rank

# print (check_m_dict)
# print (taxa_discrepency_dict)
# print(gtdb_uniprot_consistency_dict)

df = add_to_df_from_dict(df, check_m_dict, "Check_M_Highest")
df = add_to_df_from_dict(df, taxa_discrepency_dict, "Taxa_discrepency")
df = add_to_df_from_dict(df, gtdb_uniprot_consistency_dict, "GTDB_Uniprot_consistency")
df = add_to_df_from_dict(df, gtdb_rank_dict, "GTDB_Rank")

df["Updated_Rank"] = df["GTDB_Rank"].fillna(df["Taxonomic_lineage_PHYLUM"])

# Check the files
check_m_df = df[
    (df["Check_M_Highest"] < 70) & (df["Check_M_Highest"] != "No_GTDB_info")
]

taxa_df = df[
    (df["Taxa_discrepency"] != "Consistent") & (df["Taxa_discrepency"].notnull())
]

uniprot_gtdb_df = df[
    (df["GTDB_Uniprot_consistency"] != "Consistent")
    & (df["GTDB_Uniprot_consistency"].notnull())
]
df.to_csv(snakemake.output[0], index=False)


print(df["GTDB_Uniprot_consistency"])
