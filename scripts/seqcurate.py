import pandas as pd
import random
from Bio import SeqIO, AlignIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq


def get_sequence_df(
    *fasta_paths,
    drop_duplicates=True,
    alignment=False,
    ancestor=False,
    alphabet="ABCDEFGHIJKLMNOPQRSTUVWXYZ-",
):
    seq_list = []
    duplicates = {}

    cols = [
        "info",
        "truncated_info",
        "extracted_id",
        "extracted_name",
        "sequence",
        "original_fasta",
    ]

    if alignment or ancestor:
        print("is alignment")
        cols.append("original_alignment")
        cols.append("Sequence_aligned")

    # if ancestor:
    #     cols.append("Sequence_aligned")

    for fasta_path in fasta_paths:
        # Load FASTA file
        # seqs = sequence.readFastaFile(fasta_path, alpha)

        if alignment:
            seqs = AlignIO.parse(open(fasta_path), format="fasta")

        else:
            seqs = SeqIO.parse(open(fasta_path), format="fasta")

        # Add to annotation file
        for seq in seqs:
            if alignment == False:
                if seq.name in duplicates:
                    print(
                        f"DUPLICATE:{seq.name} is in {duplicates[seq.name]} and {fasta_path}\n"
                    )
                else:
                    duplicates[seq.name] = fasta_path

                curr_seq = [
                    seq.id,
                    seq.id.split(" ")[0],
                    seq.id.split("|")[1]
                    if len(seq.id.split("|")) > 1
                    else seq.id.split(" ")[0],
                    seq.id.split("|")[-1],
                    "".join(str(seq.seq).replace("-", ""))
                    if len(seq.seq) > 0
                    else None,
                    fasta_path,
                ]

                seq_list.append(curr_seq)

            elif alignment:
                for aligned_seq in seq:
                    curr_seq = [
                        aligned_seq.id,
                        aligned_seq.id.split(" ")[0],
                        aligned_seq.id.split("|")[1]
                        if len(aligned_seq.id.split("|")) > 1
                        else aligned_seq.id.split(" ")[0],
                        aligned_seq.id.split("|")[-1],
                        "".join(str(aligned_seq.seq).replace("-", ""))
                        if len(aligned_seq.seq) > 0
                        else None,
                        None,
                        fasta_path,
                        "".join(aligned_seq.seq),
                    ]
                    seq_list.append(curr_seq)

            # if ancestor:
            #     curr_seq.append("".join(aligned_seq.seq))

    df = pd.DataFrame(seq_list, columns=cols)

    if drop_duplicates:
        df = df.drop_duplicates(subset="info", keep="first")

    # Drop the sequence column if there are no sequences (i.e. if we just added a list of identifiers)
    nan_value = float("NaN")

    # df.replace("", nan_value, inplace=True)

    df.dropna(how="all", axis=1, inplace=True)

    return df


def annotate_col_from_dict(df, col, annot_dict, match="accession"):
    df.loc[df[match].isin(annot_dict.keys()), col] = df[match].map(annot_dict)
    return df


def add_col_from_up_dict(df, cols_to_add, up_dict):
    for col in cols_to_add:
        if not col in df:
            df[col] = ""

    for name, annots in up_dict.items():
        for key in annots.keys():
            df.loc[df["accession"].str.contains(name), key] = annots.get(key)

    return df


def annotate_col_from_other_df(
    df, other_df, col, match_from="accession", match_to="accession"
):
    other_dict = dict(zip(other_df[match_from], other_df[col]))

    df = annotate_col_from_dict(df, col, other_dict, match_to)

    return df


def get_subset(df, *cols_dict, include=True):
    for col_dict in cols_dict:
        if type(col_dict) != dict:
            raise TypeError("col_dict must be a dictionary")
        for col, vals in col_dict.items():
            if vals == "*":  # Get all the entries with values that exist
                if include:
                    df = df[~df[col].isna()]
                else:
                    df = df[df[col].isna()]
            else:  # Get all the entries with the specified values
                if include:
                    df = df[df[col].isin(vals)]
                else:
                    df[~df[col].isin(vals)]

    return df


def write_to_fasta(df, outpath, id_column="info", just_ids=False):
    if just_ids:
        print(id_column)
        seq_list = [
            SeqRecord(
                Seq(""),
                id=str(getattr(r, id_column)),
                description=str(getattr(r, id_column)),
            )
            for r in df.itertuples()
        ]
    else:
        seq_list = [
            SeqRecord(
                Seq(r.sequence),
                id=getattr(r, id_column),
                description=getattr(r, id_column),
            )
            for r in df.itertuples()
            if r.sequence
        ]

    SeqIO.write(seq_list, outpath, "fasta")


def randstring(length=10):
    valid_letters = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    return "".join((random.choice(valid_letters) for i in range(length)))


# Function to merge dataframes based on a common column
def add_from_csv(df, add_df, match="info"):
    """
    Merge two DataFrames based on a common column.

    Args:
        df (pd.DataFrame): The main DataFrame.
        add_df (str): Path to the CSV file containing data to be added.
        match (str): The column name to match on.

    Returns:
        pd.DataFrame: Merged DataFrame.
    """
    add_df = pd.read_csv(add_df)
    merged_df = pd.merge(df, add_df, how="left", on=[match], suffixes=["", "_r"])
    return merged_df


# Function to get entry IDs from a FASTA file
def get_entry_ids_from_fasta(fasta_path, alphabet="ABCDEFGHIJKLMNOPQRSTUVWXYZ-"):
    """
    Extract entry IDs from a FASTA file.

    Args:
        fasta_path (str): Path to the FASTA file.
        alphabet (str): Allowed characters in the entry IDs.

    Returns:
        list: List of entry IDs.
    """
    seqs = SeqIO.parse(fasta_path, "fasta")
    return [seq.name for seq in seqs]


# Function to get sequence content from a FASTA file
def get_sequence_content_from_fasta(fasta_path, alphabet="ABCDEFGHIJKLMNOPQRSTUVWXYZ-"):
    """
    Extract sequence content from a FASTA file.

    Args:
        fasta_path (str): Path to the FASTA file.
        alphabet (str): Allowed characters in the sequence.

    Returns:
        list: List of sequences.
    """
    seqs = SeqIO.parse(fasta_path, "fasta")
    return [seq.seq for seq in seqs]
