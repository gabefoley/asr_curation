import pandas as pd
import numpy as np
import seqcurate as sc
import warnings
import logging
from collections import defaultdict
import regex as re

logging.captureWarnings(True)
import seaborn as sns
from ast import literal_eval
from sklearn.feature_extraction.text import TfidfVectorizer
from sklearn.metrics.pairwise import cosine_similarity
from collections import Counter
import os
import requests


logging.basicConfig(filename="annotation_issues.log", level=logging.DEBUG)


def exact_match(df, col, match):
    return df[col] == match




# Map between a
def map_interpro_ids(ids, result_dict):
    if pd.isnull(ids):
        return None
    ids_list = ids.split(";")
    summary_list = [result_dict.get(i, "") for i in ids_list if i in result_dict]
    return ";".join(summary_list)


# Function to get OrthoDB names and level names from a dataframe with a column containing OrthoDB IDs
def get_orthodb_names(df, output_dir):
    existing_mapping = output_dir + "/orthodb_mappings.txt"

    print("Fetching OrthoDB names and levels")

    unique_orthodb_ids = get_list_of_unique_orthodb_ids(df)

    # Check for an existing mapping generated from a previous run
    name_mapping = read_to_dict(existing_mapping)

    for orthodb_id in unique_orthodb_ids:
        # If it already exists we don't need to retrieve it again
        if orthodb_id not in name_mapping:
            name, level_name = retrieve_orthodb_data(orthodb_id)
            name_mapping[orthodb_id] = name + "|" + level_name

    # Update the name mapping .txt in the annotations folder to reuse
    write_from_dict(existing_mapping, name_mapping)

    # Apply the add_orthodb_columns function to each row of the DataFrame
    df = df.apply(lambda row: add_orthodb_columns(row, name_mapping), axis=1)

    return df


def get_list_of_unique_orthodb_ids(df):
    # Get all unique OrthoDB IDs
    unique_orthodb_ids = set(x.replace(";", "") for x in df["xref_orthodb"].dropna())
    return list(unique_orthodb_ids)


# Function to query OrthoDB API
def retrieve_orthodb_data(orthodb_id):
    api_url = "https://data.orthodb.org/current/group"
    params = {"id": orthodb_id}

    print(orthodb_id)

    response = requests.get(api_url, params=params)

    if response.status_code == 200:
        print(response)
        data = response.json()
        print(data)
        print(data["data"]["name"])

        name = data["data"].get("name", "Unknown")
        level_name = data["data"].get("level_name", "Unknown")

        print("hello")
        print(name)
        print(level_name)

        return name, level_name
    else:
        return "Unknown", "Unknown"


def add_orthodb_columns(row, orthodb_mapping):
    if pd.isnull(row["xref_orthodb"]):
        row["orthodb_name"] = "Unknown"
        row["orthodb_level_name"] = "Unknown"
        row["orthodb_name_level_name"] = "Unknown Unknown"
        return row

    print("here is ")
    print(row["xref_orthodb"])

    orthodb_ids = row["xref_orthodb"].strip(";").split(";")

    # Use list comprehensions to extract the details
    names = [
        orthodb_mapping.get(oid, "Unknown|Unknown").split("|")[0] for oid in orthodb_ids
    ]
    level_names = [
        orthodb_mapping.get(oid, "Unknown|Unknown").split("|")[1] for oid in orthodb_ids
    ]

    # Using zip to combine the corresponding names and level names
    name_level_names = [f"{n} {ln}" for n, ln in zip(names, level_names)]

    row["orthodb_name"] = ";".join(names)
    row["orthodb_level_name"] = ";".join(level_names)
    row["orthodb_name_level_name"] = ";".join(name_level_names)

    return row


# Load in dictionary file stored as plain text locally
def read_to_dict(filename):
    data_dict = {}

    if not os.path.exists(filename):
        return data_dict

    with open(filename, "r") as file:
        for line in file:
            if ":" in line:
                key, value = line.strip().split(":", 1)
                data_dict[key.strip()] = value.strip()
    return data_dict


def write_from_dict(filename, data_dict):
    existing_data = read_to_dict(filename)

    # Merge dictionaries. If there are overlaps, data_dict will overwrite existing_data
    merged_data = {**existing_data, **data_dict}

    with open(filename, "w+") as file:
        for key, value in merged_data.items():
            file.write(f"{key} : {value}\n")


# Function to get interpro names from a dataframe with a column containing a list of interpro IDs
def get_interpro_names(df, output_dir):
    existing_mapping = output_dir + "/interpro_mappings.txt"

    print("Add the actual InterPro names")

    unique_interpro_ids = get_list_of_unique_interpro_ids(df)

    # Check for an existing mapping generated from a previous run
    name_mapping = read_to_dict(existing_mapping)

    for interpro_id in unique_interpro_ids:
        # If it already exists we don't need to retrieve it again
        if interpro_id not in name_mapping:
            name_mapping[interpro_id] = retrieve_interpro_name(interpro_id)

    # Update the name mapping .txt in the annotations folder to be able to reuse this
    write_from_dict(existing_mapping, name_mapping)
    # Apply the add_interpro_name_column function to each row of the DataFrame
    df = df.apply(lambda row: add_interpro_name_column(row, name_mapping), axis=1)

    return df


def get_list_of_unique_interpro_ids(df):
    # Get all unique InterPro IDs
    unique_interpro_ids = set()
    for entry in df["xref_interpro"]:
        if isinstance(entry, str):
            unique_interpro_ids.update(entry.split(";"))

    unique_interpro_ids = [x for x in unique_interpro_ids if x]

    return unique_interpro_ids


# Function to query InterPro API and retrieve names for multiple IDs
def retrieve_interpro_name(interpro_id):
    api_url = "https://www.ebi.ac.uk/interpro/api/entry/interpro/"

    url = api_url + interpro_id
    response = requests.get(url)
    if response.status_code == 200:
        data = response.json()
        interpro_name = data["metadata"]["name"]["name"]
        return interpro_name
    else:
        return "Unknown"


def add_interpro_name_column(row, interpro_mapping):
    if isinstance(row["xref_interpro"], float):
        row["interpro_name"] = None
    else:
        unique_entries = row["xref_interpro"].split(";")
        row["interpro_name"] = ";".join(
            interpro_mapping.get(ids, "Unknown") for ids in unique_entries
        )
    return row



def process_note_column(df, column):
    pattern = r'note="([^"]+)"'
    df[column + '_notes'] = df[column].str.extractall(pattern)[0].groupby(level=0).apply(lambda x: '|'.join(x))
    binary_matrix = df[column + '_notes'].str.get_dummies()
    result_df = pd.concat([df, binary_matrix.add_prefix(f'{column}||')], axis=1)
    result_df = result_df.drop(columns=[column + '_notes'])
    return result_df



def separate_notes(df):
    note_columns = df.applymap(lambda x: 'note=' in str(x)).any()
    result_columns = note_columns[note_columns].index.tolist()


    # Iterate through each column and process it
    for col in result_columns:
        df = process_note_column(df, col)
        
    return df


def annotate_motif(df, motif):
    # For a dataframe, check if the sequence column contains a certain motif
    df[f"MOTIF_{motif.replace('.', 'x').replace('[', '_').replace(']', '_')}"] = (
        df["sequence"].dropna().str.contains(motif)
    )
    return df


def get_motif_indexes(string, motif):
    # Returns the start and end positions of a motif

    if type(string) == str:
        return [
            [m.start(), m.end()]
            for m in re.finditer(rf"{motif}" "", string, overlapped=True)
        ]
    else:
        return None


def find_value(keyword, possible_words):
    for word in possible_words:
        if word in keyword:
            return word
    return "not found"


def annotate_nonAA(df):
    # Does the sequences have non amino acid characters in it

    non_AA = "B|J|O|U|X|Z"
    df["Non_AA_Character"] = df["sequence"].dropna().str.contains(non_AA)

    return df





def annotate_AA(df):
    # Does the sequences have non amino acid characters in it

    non_AA = "B|J|O|U|X|Z"
    df["AA_Character"] = ~(df["sequence"].dropna().str.contains(non_AA, na=None))

    return df

def get_pos(align_df, seq_id, col_name, pos_type="list"):
    """
    Get positions from a DataFrame column for a specific sequence ID.

    Args:
        align_df (DataFrame): The DataFrame containing alignment data.
        seq_id (str): The ID of the sequence for which positions are to be retrieved.
        col_name (str): The name of the column containing the positions in string format.
        pos_type (str, optional): The type of positions to return, either 'list' or 'range'.
                                  Defaults to 'list'.

    Returns:
        list or range: The positions as either a list or a range, depending on pos_type.
                      Returns None if the positions are not available in the DataFrame.

    Raises:
        NameError: If an invalid pos_type is provided.

    """
    # Query the DataFrame to find the row with the given sequence ID and extract the column value
    pos_str = align_df.query(f"info=='{seq_id}'")[col_name].tolist()[0]

    # Use literal_eval to safely convert the string representation of positions to a Python object
    positions = literal_eval(pos_str)

    if not positions:
        return None

    if pos_type == "list":
        # Convert positions to a list and increment each element by 1
        positions = [x + 1 for x in positions]

    elif pos_type == "range":
        # Convert positions to a range and increment each element by 1
        positions = range(positions[0] + 1, positions[1] + 1)

    else:
        raise NameError("Not a valid position type")

    return positions


# def add_tag_relative_to_seq(align_df, seq_id, tag, unaligned_pos, aln_dict):
#         print(f"Add in {tag}")
#         if seq_id in aln_dict:
#             aligned_seq = aln_dict[seq_id]
#             align_df = track_residues(align_df, seq_id, aligned_seq, tag, *unaligned_pos)
#         return align_df


def track_residues2(align_df, seq_id, aligned_seq, tag, *unaligned_pos):
    # Map the positions to an index in the alignment
    # aligned_pos = get_aligned_positions(
    #     align_df[align_df["info"] == seq_id], aligned_seq, *unaligned_pos
    # )

    aligned_pos = get_aligned_positions(aligned_seq, *unaligned_pos)

    # Add the aligned positions we want to track to the dataframe
    align_df = track_aligned_positions(align_df, seq_id, tag, aligned_pos)

    align_df = get_tracked_content(align_df, tag, *aligned_pos)

    return align_df


# def track_residues(seq_id, aligned_)


def assign_labels(row):
    proteins = row["protein_name"]

    # Create a TF-IDF vectorizer to convert protein names into numerical representations
    vectorizer = TfidfVectorizer()
    proteins_tfidf = vectorizer.fit_transform(proteins)

    # Calculate pairwise cosine similarity between proteins
    cosine_similarities = cosine_similarity(proteins_tfidf, proteins_tfidf)

    # Define a similarity threshold to determine group membership
    similarity_threshold = 0.5

    # Group proteins based on cosine similarity
    groups = []
    assigned = [False] * len(
        proteins
    )  # Track whether a protein has been assigned to a group

    for i in range(len(proteins)):
        if not assigned[i]:
            group = [proteins[i]]  # Start a new group
            assigned[i] = True

            for j in range(i + 1, len(proteins)):
                if (
                    not assigned[j]
                    and cosine_similarities[i, j] >= similarity_threshold
                ):
                    group.append(proteins[j])
                    assigned[j] = True

            groups.append(group)

    # Define labels for each group
    labels = []

    # Iterate over each group
    for group in groups:
        # Concatenate all protein names into a single string
        all_names = " ".join(group)

        # Extract individual keywords from the concatenated string
        keywords = all_names.split()

        # Count the occurrence of each keyword
        keyword_counts = Counter(keywords)

        # Find the most common keywords
        most_common_keywords = keyword_counts.most_common()

        # Concatenate the most common keywords into a string
        label = " ".join(
            keyword[0] for keyword in most_common_keywords[:4]
        )  # Limit to four words

        # Assign the concatenated string as the label
        labels.append(label)

    # Return the labels
    return labels


def get_aligned_pos_and_content(seq_map_from, seq_map_to, *pos_set):
    """"""

    aligned_pos_set = []
    content_set = []

    for pos in pos_set:
        aligned_pos = get_aligned_positions(seq_map_from, *pos)

        content = get_content_at_pos(seq_map_to, *aligned_pos)

        aligned_pos_set.append(aligned_pos)
        content_set.append([content])

    return pd.Series([aligned_pos_set, content_set])


def get_aligned_positions(sequence, *positions):
    """
    Take position/s and return the equivalent positions/s in an alignment - accounting for gap characters added

    Args:
        sequence (str): The aligned sequence to map the position/s to
        positions (int): The unaligned position/s to create a mapping from

    Returns:
        List of the equivalent positions in the aligned sequence
    """

    sequence = "".join(sequence)


    aligned_positions = []

    for pos in positions:
        # Get the current index
        curr_idx = pos

        print(curr_idx)
        print(type(curr_idx))

        # Get offset implied by first position
        offset = sequence[0 : int(curr_idx)].count("-")

        # Get next position based on the first position and the gap offset
        next_idx = curr_idx + offset

        # Search to find the final position

        final_pos = get_final_pos(sequence, pos, curr_idx, next_idx)
        aligned_positions.append(final_pos)

    # Update the indexes
    aligned_positions = [x - 1 for x in aligned_positions]

    return aligned_positions


def get_final_pos(sequence, pos, curr_idx, next_idx):
    # Need to find the position that 1) isn't a gap and 2) takes into account all of the
    # offset implied by previous gaps in the sequence

    # If the content at this index is a gap, proceed to the next actual position
    while curr_idx < len(sequence) and sequence[max(curr_idx - 1, 0)] == "-":
        curr_idx += 1

    # Get the offset implied by the number of gaps in the preceeding sequence

    offset = sequence[0:curr_idx].count("-")

    # Add the offset implied by the gaps in the previous positions
    next_idx = pos + offset

    # If the offset is the same as previous offset we are at our final position
    if curr_idx == next_idx:
        return curr_idx
    else:
        # Keep searching
        return get_final_pos(sequence, pos, next_idx, next_idx)


def get_content_at_pos(seq, *pos):
    # print ('get content at pos')

    # print (seq)
    # print (pos)

    return "".join([seq[int(p)] for p in pos])


# def get_content_from_position(align_df, tag, *pos):
#     align_df[tag] = align_df.apply(
#         lambda row: get_amino_acids(row["Sequence_aligned"], *aligned_pos),
#         axis=1,
#     )
#     return align_df
def track_aligned_positions(align_df, seq_id, tag, aligned_pos):
    seq_index = align_df.index[align_df["info"] == seq_id].tolist()[0]
    # align_df.at[seq_index, f"tracked_{tag}"] = ",".join(str(x) for x in aligned_pos)
    align_df.at[seq_index, f"tracked_{tag}"] = aligned_pos

    return align_df


def get_tracked_content(align_df, tag, *aligned_pos):
    align_df[tag] = align_df.apply(
        lambda row: get_amino_acids(row["Sequence_aligned"], *aligned_pos),
        axis=1,
    )
    return align_df


def add_tag_if_in_fasta(annot_df, filepath, tag):
    seqs = sc.get_entry_ids_from_fasta(filepath)

    annot_df[tag] = annot_df.apply(lambda row: 1 if row["info"] in seqs else 0, axis=1)

    return annot_df


def add_tag_if_sequence_matches(annot_df, filepath, tag):
    seqs = sc.get_sequence_content_from_fasta(filepath)

    annot_df[tag] = annot_df.apply(
        lambda row: True if row["sequence"] in seqs else False, axis=1
    )

    return annot_df


def annotate_sp_tr(df):
    # Is the sequence from SwissProt or TrEMBL
    df.loc[df["info"].str.startswith("sp"), "UniProt_DB"] = "SwissProt"
    df.loc[df["info"].str.startswith("tr"), "UniProt_DB"] = "TrEMBL"

    return df


def add_lab_annotations(annot_df, filepath, seq_col="sequence"):
    lab_df = pd.read_csv(filepath, encoding="utf-8", header=0)

    annot_df.columns = annot_df.columns.str.replace(" ", "")
    lab_df.columns = lab_df.columns.str.replace(" ", "")

    if "info" not in lab_df.columns:
        raise ValueError("Lab annotations are missing info field")

    if "sequence" not in lab_df.columns:
        raise ValueError("Lab annotations are missing sequence field")
    # Get the columns for id / sequence that are different
    df_diff = pd.concat(
        [annot_df[["info", seq_col]], lab_df[["info", "sequence"]]]
    ).drop_duplicates(keep=False)

    # Get the ids from df_diff, if an id is in both then the sequence might be different - need to raise an error

    if not df_diff.empty:
        for val in df_diff["info"].values:
            if val in annot_df["info"].values and val in lab_df["info"].values:
                annot_seq = annot_df.loc[annot_df["info"] == val, seq_col].values[0]
                lab_seq = lab_df.loc[lab_df["info"] == val, "sequence"].values[0]
                if annot_seq != lab_seq:
                    raise ValueError(
                        "Lab annotations contain different sequence to existing annotations"
                    )

    for col in [x for x in lab_df.columns if not x.strip().startswith("lab")]:
        if col not in ["info", "sequence"] and col in annot_df.columns:
            raise ValueError("Duplicate column between lab and existing annotations")

    if len(lab_df.columns) != len(
        lab_df.columns.str.replace(".1$", "").drop_duplicates()
    ):
        raise ValueError("Lab annotations contain multiple identically named columns")

    # existing_lab_annotations =  [x for x in annot_df if x.strip().startswith("lab")]

    # annot_df['id'] = annot_df['id'].astype(object)
    merge_on = ["info", "sequence"]

    # merged_df = pd.concat([annot_df, lab_df], axis=1)

    # merged_df = annot_df.join(lab_df, on='id', how='left', lsuffix='_dup', rsuffix='_dup1')

    merged_df = pd.merge(
        annot_df,
        lab_df,
        how="left",
        left_on=merge_on,
        right_on=merge_on,
        suffixes=["_dup", "_dup1"],
    )

    dup_cols = [i.strip() for i in merged_df.columns if i.endswith("_dup")]

    for dup_col in dup_cols:
        merged_df[dup_col.split("_dup")[0]] = merged_df.apply(
            lambda row: merge_lab_annotation_cells(
                row, row[dup_col], row[dup_col + "1"]
            ),
            axis=1,
        )
        merged_df.drop([dup_col, dup_col + "1"], inplace=True, axis=1)

    return merged_df


def merge_lab_annotation_cells(id, first_cell, second_cell):
    if pd.isnull(first_cell):
        return second_cell

    elif pd.isnull(second_cell):
        return first_cell

    else:
        overwritten = False
        first_list = str(first_cell).split(",")
        second_list = str(second_cell).split(",")
        if first_list:
            for pos, val in enumerate(first_list):
                if pos > len(second_list):
                    warnings.warn(
                        f"Lab annotations have missing values - {id}",
                        UserWarning,
                    )

                if val == second_list[pos]:
                    continue
                else:
                    overwritten = True
                    warnings.warn(
                        f"Lab annotations are overwriting values - {id}",
                        UserWarning,
                    )

            if not overwritten:
                warnings.warn(
                    f"Lab annotations are adding to values - {id}", UserWarning
                )

    return second_cell


def check_terms(check, terms):
    if pd.isnull(check):
        return
    for term in terms:
        if term.lower() in check.lower():
            return True
    return False


def add_thermo(annot_df, filepath):
    thermo = open(filepath)

    thermo_species = [x.strip() for x in thermo.readlines() if len(x) > 1]
    thermo_split = [x.split(" ")[0] for x in thermo_species]
    thermo_species_terms = ["therm", "acid", "sulfur", "methan", "pyro", "lividus"]
    thermo_species_terms_no_methan = ["therm", "acid", "sulfur", "pyro", "lividus"]

    annot_df["thermo_bacteria_split"] = annot_df.apply(
        lambda row: (
            True if row["lineage_genus"].split(" ")[0] in thermo_split else False
        )
        if pd.notnull(row["lineage_genus"])
        else "No genus",
        axis=1,
    )

    annot_df["thermo_bacteria"] = annot_df.apply(
        lambda row: True
        if row["lineage_genus"] in thermo_species
        else False
        if pd.notnull(row["lineage_genus"])
        else "No genus",
        axis=1,
    )

    annot_df["thermo_terms"] = annot_df.apply(
        lambda row: True
        if check_terms(row["lineage_genus"], thermo_species_terms)
        else False
        if pd.notnull(row["lineage_genus"])
        else "No genus",
        axis=1,
    )

    annot_df["lineage_thermo"] = annot_df.apply(
        lambda row: True
        if check_terms(row["lineage"], thermo_species_terms + thermo_species)
        else False
        if pd.notnull(row["lineage"])
        else False,
        axis=1,
    )

    annot_df["lineage_thermo_no_metha"] = annot_df.apply(
        lambda row: True
        if check_terms(row["lineage"], thermo_species + thermo_species_terms_no_methan)
        else False
        if pd.notnull(row["lineage"])
        else False,
        axis=1,
    )
    return annot_df


def create_domain_bounds(domains):
    positions = []

    interval = None

    domain_name = None

    if pd.isnull(domains):
        return positions

    for domain in domains.split(";"):
        if domain.strip().startswith("DOMAIN"):
            pos = domain.split("DOMAIN ")[1].split("..")

            interval = pd.Interval(
                int(pos[0].replace("<", "").replace(">", "")),
                int(pos[1].replace("<", "").replace(">", "")),
            )
        if domain.startswith(" /note="):
            domain_name = domain.split('/note="')[1][
                0:-1
            ]  # Domain name, minus the final quotation

        if interval and domain_name:
            positions.append((domain_name, interval))
            interval = None
            domain_name = None

    return positions


def create_ss_bounds(strand, helix, turn):
    names = ["STRAND", "HELIX", "TURN"]
    positions = defaultdict(list)
    for name, feats in zip(names, [strand, helix, turn]):
        if type(feats) == float and np.isnan(feats):
            pass
        else:
            for feat in feats.split(name):
                if feat:
                    pos = feat.strip().split(";")[0].split("..")
                    interval = pd.Interval(int(pos[0]), int(pos[1]))
                    positions[name].append(interval)
    return positions


def create_combined_ss_bounds(strand, helix, turn):
    names = ["STRAND", "HELIX", "TURN"]
    positions = []
    for name, feats in zip(names, [strand, helix, turn]):
        if type(feats) == float and np.isnan(feats):
            pass
        else:
            for feat in feats.split(name):
                if feat:
                    pos = feat.strip().split(";")[0].split("..")
                    interval = pd.Interval(int(pos[0]), int(pos[1]))
                    positions.append((name, interval))
    return positions


def create_annotated_alignment(df, boundary_dict, outpath, colour_dict=None):
    # If we don't have a supplied colour_dict, lets make one
    if not colour_dict:
        boundary_labels = set(
            [label[0] for entry in boundary_dict.values() for label in entry]
        )
        palette = sns.color_palette(None, len(boundary_labels)).as_hex()

        colour_dict = {col: label for col, label in zip(boundary_labels, palette)}

    # Creating an HTML file
    with open(outpath, "w") as align_html:
        # Get the length needed
        for acc, _bounds in boundary_dict.items():
            aligned_seq = df.loc[df["extracted_id"] == acc]["Sequence_aligned"].values[
                0
            ]
            alignment_length = str(len(aligned_seq) * 10)

        align_html.write(
            '<html>\n<head><style>  #container{    width : 20px;  }  .item{ overflow: hidden; white-space: nowrap; font-family:"Courier New";   float:left;    width: '
            + alignment_length
            + 'px;    height: 20px;    padding 2px;    margin: 0px 2px;  }  .clearfix{    clear: both;  }</style>\n<title> \nOutput Data in an HTML file</title>\n</head> <body>  <div id="container">'
        )

        for acc, bounds in boundary_dict.items():
            if bounds:
                orig_seq = df.loc[df["extracted_id"] == acc]["Sequence_aligned"].values[
                    0
                ]
                formatted_sequence = df.loc[df["extracted_id"] == acc][
                    "Sequence_aligned"
                ].values[0]
                len_offset = 0

                bounds.sort(key=lambda x: x[1])
                for bound in bounds:
                    bound_name = bound[0]
                    pos = bound[1]
                    #                  for domain, pos in list_w_overlaps:
                    gap_offset = 0
                    first_gap_offset = 0
                    second_gap_offset = 0

                    count = 0
                    for aa in orig_seq:
                        if aa == "-":
                            gap_offset += 1
                        else:
                            count += 1
                            if count == pos.left:
                                first_gap_offset = gap_offset

                            if count == pos.right:
                                second_gap_offset = gap_offset
                                break
                            else:
                                continue

                    prev_len = len(formatted_sequence)

                    formatted_sequence = (
                        formatted_sequence[
                            0 : pos.left - 1 + len_offset + first_gap_offset
                        ]
                        + '<span style = "background-color:'
                        + colour_dict[bound_name]
                        + '">'
                        + formatted_sequence[
                            pos.left
                            - 1
                            + len_offset
                            + first_gap_offset : pos.right
                            + len_offset
                            + second_gap_offset
                        ]
                        + "</span>"
                        + formatted_sequence[
                            pos.right + len_offset + second_gap_offset :
                        ]
                    )

                    len_offset = len(formatted_sequence) - len(orig_seq)

            #             formatted_sequence = 'red'
            else:
                formatted_sequence = df.loc[df["extracted_id"] == acc][
                    "Sequence_aligned"
                ].values[0]

            align_html.write(f"<br>>{acc}<br>")
            align_html.write(f'<div class="item">{formatted_sequence}</div>')

        align_html.write("</body></html>")


##### KARI SPECIFIC ####


def classify_KARI(features):
    if "Domain" in features:
        domain_num = features.split("Domain")[1].split(";")[0]
        if "2" in domain_num:
            return "Class_1"
        elif "3" in domain_num:
            return "Class_2"
        else:
            return "Different_domain_number"
    else:
        return "No_domain_info"


def get_amino_acids(seq, *pos):
    return ["".join([seq[int(bp)] for bp in pos])]


def get_binding_pos(id, binding_sites, ligand=None):
    if pd.notnull(binding_sites):
        results = []

        # Split by "BINDING" to get the chunks
        segments = re.split(r"BINDING\s*\d+;", binding_sites)
        binding_positions = re.findall(r"BINDING\s*(\d+);", binding_sites)

        for pos, segment in zip(
            binding_positions, segments[1:]
        ):  # Ignore the first split as it will be empty
            data = {}
            data["binding_pos"] = int(pos) - 1

            # Extract ligand, ligand_id, and evidence
            ligand_match = re.search(r'/ligand="([^"]+)"', segment)
            ligand_id_match = re.search(r'/ligand_id="([^"]+)"', segment)
            evidence_match = re.search(r'/evidence="([^"]+)"', segment)

            if ligand_match:
                data["binding_ligand"] = ligand_match.group(1)
            if ligand_id_match:
                data["binding_ligand_id"] = ligand_id_match.group(1)
            if evidence_match:
                data["binding_evidence"] = evidence_match.group(1)

            # If a specific ligand is provided, filter based on that
            if not ligand or (ligand and "ligand" in data and data["ligand"] == ligand):
                results.append(data)

        return results
    else:
        return []


def create_top_column(df, threshold_percentage, max_threshold_percentage=100):
    new_df = df.copy()

    for column in df.columns:
        # Check if the column name already starts with "TOP_"
        if not column.startswith("TOP_"):
            # Calculate the percentage of the most common value in the column, including NaN
            value_counts = df[column].value_counts(normalize=True, dropna=False)
            if not value_counts.empty:
                most_common_value_percentage = value_counts.max() * 100
                most_common_value = value_counts.idxmax()

                # Check if the most common value is not np.nan and meets the threshold conditions
                if (most_common_value_percentage > threshold_percentage and
                        most_common_value_percentage < max_threshold_percentage and
                        not pd.isna(most_common_value)):
                    new_column_name = f"TOP_{column}||{most_common_value}"
                    new_df[new_column_name] = (df[column] == most_common_value) & (~df[column].isna())

    return new_df.dropna(axis=1, how='all')



def check_sequence_for_loop_length(seq):
    print(seq[0][0].replace("-", ""))
    return len(seq[0][0].replace("-", ""))


def check_sequence_for_acidic(seq):
    print(seq)
    print(seq[0][0])
    if "E" in seq[0][0] or "D" in seq[0][0]:
        return "Acidic_residue"
    else:
        return "No_Acidic_residue"


def check_sequence_for_charged(seq):
    print(seq)
    print(seq[0][0])
    if "R" in seq[0][0] or "K" in seq[0][0]:
        return "Charged_residue"
    else:
        return "No_charged_residue"


def check_binding_for_acidic(seq, bind_pos):
    if bind_pos and bind_pos != "No_binding_positions":
        binding_aa = get_amino_acids(seq, *bind_pos)
        return check_sequence_for_acidic(binding_aa)


def classify_loop_length(bind_pos_set):
    # Offset accounts for the fact that we need to get the total number of positions and also account for the first position in the loop which isn't a binding position

    for bind_pos in bind_pos_set:
        if bind_pos and bind_pos != "No_binding_positions":
            offset = 2

            return bind_pos[-1] - bind_pos[0] + offset
        else:
            return "No_binding_positions"


def check_if_positions_align_with_target(target_pos, seq_pos):
    if target_pos == seq_pos:
        return True


# Function to add a new column to the DataFrame
def add_labels_from_file(df, column_name, file_path):
    with open(file_path, 'r') as file:
        entries = {}
        for line in file:
            parts = line.strip().split('\t')
            if len(parts) >= 2:
                entries[parts[0]] = parts[1]

        print (entries)

        df[column_name] = df['info'].apply(lambda x: entries.get(x, None))

    return df

