import pandas as pd
import urllib.parse
import urllib.request
from io import StringIO

# At the moment this workflow just supports two ID types. In the future if we want to support more we need to change
# this workflow

# Setup the two types of IDs to map to / from
to_from_dict = {"ACC": "P_REFSEQ_AC", "P_REFSEQ_AC": "ACC"}


def get_uniprot_ids(queries, from_id, to_id):
    url = "https://www.uniprot.org/uploadlists/"

    queries = " ".join(queries)

    params = {"from": from_id, "to": to_id, "format": "tab", "query": queries}

    #     print ('getting uniprot')

    data = urllib.parse.urlencode(params)
    data = data.encode("utf-8")
    req = urllib.request.Request(url, data)
    with urllib.request.urlopen(req) as response:
        res = response.read()
    df_fasta = pd.read_csv(StringIO(res.decode("utf-8")), sep="\t")
    df_fasta.columns = ["Entry", "To"]

    df_fasta = df_fasta.groupby(["Entry"]).agg(To=("To", ";".join)).reset_index()

    return df_fasta


def map_df_to_uniprot(chunk_df, try_ids, map_back, from_id, to_id, to_name):
    #     print ('variables')
    #     print(try_ids)
    #     print (map_back)
    #     print (from_id)
    #     print (to_id)
    #     print (to_name)

    # Try mapping the uniprot IDs to UniProt IDs
    up_df = get_uniprot_ids(try_ids, from_id, to_id)

    up_df = up_df.rename(columns={"Entry": map_back})

    # Check if the item returned and create a new column that validates this as being a correct ID

    #     print("this was the up df from here")
    #     print(up_df)
    #     print()

    if not up_df.empty:
        #         print("chunk here is ")

        #         print(chunk_df)
        #         print()

        up_df["join"] = 1
        chunk_df["join"] = 1

        # Joint the two dataframes then drop the join columns

        merged_df = chunk_df.merge(up_df, on="join").drop("join", axis=1)

        del chunk_df["join"]

        #         print ('merged')

        #         print(merged_df)

        merged_df["match"] = merged_df.apply(
            lambda row: row[f"{map_back}_x"].find(row[f"{map_back}_y"]), axis=1
        ).ge(0)

        matched_df = merged_df[merged_df.match]

        matched_df = matched_df.drop(["match", f"{map_back}_y"], axis=1)

        if f"{to_name}_validated" in matched_df.columns:
            matched_df = matched_df.drop(f"{to_name}_validated", axis=1)

        #         print ('matched')
        #         print (matched_df)

        # Rename the Entry column back and specify which ID this has been validated as in the To column
        matched_df = matched_df.rename(
            columns={f"{map_back}_x": map_back, "To": f"{to_name}_validated"}
        )

        #         print ('renamed')
        #         print (matched_df)

        return matched_df

    return chunk_df


def merge_dfs(chunk_df, add_df, map_back):
    # print ('check here')
    # print (chunk_df)

    # print ('****')
    # print (add_df)

    # print ('****')

    merge_cols = [x for x in list(chunk_df.columns) if x != map_back]
    chunk_df = chunk_df.merge(
        add_df,
        how="left",
        on=merge_cols,
    )

    # print (chunk_df)

    if f"{map_back}_x" in chunk_df.columns:
        chunk_df[f"{map_back}"] = chunk_df[f"{map_back}_x"].fillna(
            chunk_df[f"{map_back}_y"]
        )

        chunk_df = chunk_df.drop([f"{map_back}_x", f"{map_back}_y"], axis=1)

    return chunk_df


df = pd.read_csv(snakemake.input[0])


entries = [name for name in df["Entry"].tolist()]

ncbi_headers = ["WP_", "XP_", "YP_", "NP_", "AP_"]

ncbi_ids = [x for x in entries if x[0:3] in ncbi_headers]

df["Try"] = df.apply(
    lambda row: "P_REFSEQ_AC" if row["Entry"] in ncbi_ids else "ID", axis=1
)

# Sort the dataframe so that when we chunk the dataframe most of the calls will be to the same type of ID
sorted_df = df.sort_values(by=["Try"])


# ncbi_ids = [x for x in entries if ]

n = 500
list_df = [sorted_df[i : i + n] for i in range(0, df.shape[0], n)]


# For holding the processed chunks
final_list = []

# For each data base ID, try and get the set of UniProt IDs

for chunk_df in list_df:
    # Map each sequence to a UniProt ID

    for db_id in chunk_df["Try"].unique():
        # Define which column we want to use to search and to map back the search results onto
        map_back = "Entry"

        # IDs we want to map to UniProt identifiers
        try_ids = chunk_df[chunk_df["Try"] == db_id][map_back].to_list()

        # Here we don't know if the ID is UniProt or RefSeq or EMBL so we use our best guess stored in db_id
        add_df = map_df_to_uniprot(chunk_df, try_ids, map_back, db_id, "ID", "UniProt")

        map_back = "UniProt_validated"

        # Some still might not have mapped to UniProt correctly so try them with EMBL_ID as well

        # Maybe they all failed
        if map_back not in add_df.columns:
            missing_ids = try_ids

        # Or maybe just a subset
        else:
            missing_ids = add_df[add_df[map_back].isna()]["Entry"].tolist()

        # If there are missing IDs, try mapping them as EMBL IDs

        if missing_ids:
            add_df = map_df_to_uniprot(
                chunk_df, missing_ids, "Entry", "EMBL", "ID", "UniProt"
            )

        # Add the add_df to the chunk_df

        chunk_df = merge_dfs(chunk_df, add_df, map_back)

    # Now that we have validated UniProt identifiers, lets map them to RefSeq
    try_ids = chunk_df[map_back].dropna().to_list()

    # Here we know that the validated IDs are UniProt so we can hardcode the ID field

    add_df = map_df_to_uniprot(
        chunk_df, try_ids, map_back, "ID", "P_REFSEQ_AC", "RefSeq_Protein"
    )

    chunk_df = merge_dfs(chunk_df, add_df, map_back)

    add_df = map_df_to_uniprot(chunk_df, try_ids, map_back, "ID", "EMBL", "EMBL")

    chunk_df = merge_dfs(chunk_df, add_df, map_back)
    # print ('and now')

    # print (chunk_df)

    # print ("*****")

    final_list.append(chunk_df)

final_df = pd.concat(final_list)
final_df.to_csv(snakemake.output[0], index=False)
