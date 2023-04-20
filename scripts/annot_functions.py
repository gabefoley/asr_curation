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

logging.basicConfig(filename="annotation_issues.log", level=logging.DEBUG)


def annotate_motif(df, motif):
    # For a dataframe, check if the sequence column contains a certain motif
    df[f"MOTIF_{motif.replace('.', 'x')}"] = df["sequence"].dropna().str.contains(motif)
    return df


def get_motif_indexes(string, motif):
    # Returns the start and end positions of a motif

    if (type(string) == str):
        return [[m.start(), m.end()] for m in re.finditer(rf"{motif}" "", string, overlapped=True)]
    else:
        return None


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
    if not literal_eval(align_df.query(f"id=='{seq_id}'")[col_name].tolist()[0]):
        return

    pos = literal_eval(align_df.query(f"id=='{seq_id}'")[col_name].tolist()[0])[0]

    if pos_type == "list":
        pos = [x for x in pos]
        pos[0] += 1

    elif pos_type == "range":
        pos = [x + 1 for x in range(pos[0], pos[1])]

    else:
        raise(NameError("Not a valid position type"))

    return pos


# def add_tag_relative_to_seq(align_df, seq_id, tag, unaligned_pos, aln_dict):
#         print(f"Add in {tag}")
#         if seq_id in aln_dict:
#             aligned_seq = aln_dict[seq_id]
#             align_df = track_residues(align_df, seq_id, aligned_seq, tag, *unaligned_pos)
#         return align_df


def track_residues2(align_df, seq_id, aligned_seq, tag, *unaligned_pos):
    # Map the positions to an index in the alignment
    aligned_pos = get_aligned_positions(
        align_df[align_df["id"] == seq_id], aligned_seq, *unaligned_pos
    )

    # Add the aligned positions we want to track to the dataframe
    align_df = track_aligned_positions(align_df, seq_id, tag, aligned_pos)

    align_df = get_tracked_content(align_df, tag, *aligned_pos)

    return align_df


# def track_residues(seq_id, aligned_)

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

        # Get offset implied by first position
        offset = sequence[0:curr_idx].count("-")

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

    return "".join([seq[int(p)] for p in pos])


# def get_content_from_position(align_df, tag, *pos):
#     align_df[tag] = align_df.apply(
#         lambda row: get_amino_acids(row["Sequence_aligned"], *aligned_pos),
#         axis=1,
#     )
#     return align_df
def track_aligned_positions(align_df, seq_id, tag, aligned_pos):
    seq_index = align_df.index[align_df["id"] == seq_id].tolist()[0]
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

    annot_df[tag] = annot_df.apply(
        lambda row: True if row["id"] in seqs else False, axis=1
    )

    return annot_df


def annotate_sp_tr(df):
    # Is the sequence from SwissProt or TrEMBL
    df.loc[df["info"].str.startswith("sp"), "UniProt_DB"] = "SwissProt"
    df.loc[df["info"].str.startswith("tr"), "UniProt_DB"] = "TrEMBL"

    return df


# def get_aligned_positions(entry, sequence, *positions):
#     # print(f"\nSeq name is {entry}")
#     sequence = "".join(sequence)
#
#     aligned_positions = []
#
#     for pos in positions:
#         # Get the current index
#         curr_idx = pos
#
#         # Get offset implied by first position
#         offset = sequence[0:curr_idx].count("-")
#
#         # Get next position based on the first position and the gap offset
#         next_idx = curr_idx + offset
#
#         # Search to find the final position
#
#         final_pos = get_final_pos(sequence, pos, curr_idx, next_idx)
#         aligned_positions.append(final_pos)
#
#     # Update the indexes
#     aligned_positions = [x - 1 for x in aligned_positions]
#
#     return aligned_positions




def add_lab_annotations(annot_df, filepath, seq_col="sequence"):
    lab_df = pd.read_csv(filepath, encoding="utf-8", header=0)

    annot_df.columns = annot_df.columns.str.replace(" ", "")
    lab_df.columns = lab_df.columns.str.replace(" ", "")

    if "id" not in lab_df.columns:
        raise ValueError("Lab annotations are missing id field")

    if "sequence" not in lab_df.columns:
        raise ValueError("Lab annotations are missing sequence field")
    # Get the columns for id / sequence that are different
    df_diff = pd.concat(
        [annot_df[["id", seq_col]], lab_df[["id", "sequence"]]]
    ).drop_duplicates(keep=False)

    # Get the ids from df_diff, if an id is in both then the sequence might be different - need to raise an error

    if not df_diff.empty:
        for val in df_diff["id"].values:
            if (
                val in annot_df["id"].values
                and val in lab_df["id"].values
            ):
                annot_seq = annot_df.loc[annot_df["id"] == val, seq_col].values[
                    0
                ]
                lab_seq = lab_df.loc[lab_df["id"] == val, "sequence"].values[0]
                if annot_seq != lab_seq:
                    raise ValueError(
                        "Lab annotations contain different sequence to existing annotations"
                    )

    for col in [x for x in lab_df.columns if not x.strip().startswith("lab")]:
        if col not in ["id", "sequence"] and col in annot_df.columns:
            raise ValueError("Duplicate column between lab and existing annotations")

    if len(lab_df.columns) != len(
        lab_df.columns.str.replace(".1$", "").drop_duplicates()
    ):
        raise ValueError("Lab annotations contain multiple identically named columns")

    # existing_lab_annotations =  [x for x in annot_df if x.strip().startswith("lab")]

    # annot_df['id'] = annot_df['id'].astype(object)
    merge_on = ["id", "sequence"]

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
            aligned_seq = df.loc[df["extracted_id"] == acc]["Sequence_aligned"].values[0]
            alignment_length = str(len(aligned_seq) * 10)

        align_html.write(
            '<html>\n<head><style>  #container{    width : 20px;  }  .item{ overflow: hidden; white-space: nowrap; font-family:"Courier New";   float:left;    width: '
            + alignment_length
            + 'px;    height: 20px;    padding 2px;    margin: 0px 2px;  }  .clearfix{    clear: both;  }</style>\n<title> \nOutput Data in an HTML file</title>\n</head> <body>  <div id="container">'
        )

        for acc, bounds in boundary_dict.items():
            if bounds:
                orig_seq = df.loc[df["extracted_id"] == acc]["Sequence_aligned"].values[0]
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


def get_binding_pos(id, binding_sites, ligand=None):
    if pd.notnull(binding_sites):
        bp = []

        for site in binding_sites.split(";"):
            if site.strip().startswith("BINDING"):
                found_pos = site.split("BINDING")[1]

            if (
                site.strip().startswith("/ligand=")
                and site.split("/ligand=")[1].startswith('"NADP')
                and ".." not in found_pos
                and int(found_pos) < 100
            ):
                bp.append(int(found_pos))


        return [bp]
    else:
        return []


def get_amino_acids(seq, *pos):

    return ["".join([seq[int(bp)] for bp in pos])]


def check_sequence_for_acidic(seq):
    if "E" in seq or "D" in seq:
        return "Acidic_residue"
    else:
        return "No_Acidic_residue"


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
