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

    cols = ["Entry", "Truncated_Info", "Extracted_ID", "Sequence", "Original_FASTA"]

    if alignment:
        cols.append("Sequence_aligned")

    if ancestor:
        cols.append("Sequence_aligned_ancestor")

    for fasta_path in fasta_paths:

        # Load FASTA file
        # seqs = sequence.readFastaFile(fasta_path, alpha)

        if alignment:
            seqs = AlignIO.parse(open(fasta_path), format='fasta')

        else:

            seqs = SeqIO.parse(open(fasta_path), format='fasta')

        # Add to annotation file
        for seq in seqs:

            if alignment==False: 
                if seq.name in duplicates:
                    print(
                        f"DUPLICATE:{seq.name} is in {duplicates[seq.name]} and {fasta_path}\n"
                    )
                else:
                    duplicates[seq.name] = fasta_path

                curr_seq = [
                    seq.id,
                    seq.id.split(" ")[0],
                    seq.id.split("|")[-1],
                    seq.seq.replace("-", "")
                    if len(seq.seq) > 0
                    else None,
                    fasta_path,
                ]

            elif alignment:
      

                for aligned_seq in seq:
                    curr_seq = [
                    aligned_seq.id,
                    aligned_seq.id.split(" ")[0],
                    aligned_seq.id.split("|")[-1],
                    aligned_seq.seq.replace("-", "")
                    if len(aligned_seq.seq) > 0
                    else None,
                    fasta_path,
                ]
                    curr_seq.append(aligned_seq.seq)


            if ancestor:
                curr_seq.append(aligned_seq)

            seq_list.append(curr_seq)


    df = pd.DataFrame(seq_list, columns=cols)

    if drop_duplicates:
        df = df.drop_duplicates(subset="Entry", keep="first")

    # Drop the sequence column if there are no sequences (i.e. if we just added a list of identifiers)
    nan_value = float("NaN")

    df.replace("", nan_value, inplace=True)

    df.dropna(how="all", axis=1, inplace=True)

    return df


def annotate_col_from_dict(df, col, annot_dict, match="Entry"):
    df.loc[df[match].isin(annot_dict.keys()), col] = df[match].map(annot_dict)
    return df


def add_col_from_up_dict(df, cols_to_add, up_dict):
    for col in cols_to_add:
        if not col in df:
            df[col] = ""

    for name, annots in up_dict.items():

        for key in annots.keys():
            df.loc[df["Entry"].str.contains(name), key] = annots.get(key)

    return df


def annotate_col_from_other_df(df, other_df, col, match_from="Entry", match_to="Entry"):
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


def write_to_fasta(df, outpath, trim=False):

    if trim:
        seq_list = [
            SeqRecord(Seq(r.Sequence), id=r.Entry) for r in df.itertuples()
        ]


    else:
        seq_list = [
            SeqRecord(Seq(r.Sequence), id=r.Entry) for r in df.itertuples()
        ]
    # sequence.writeFastaFile(outpath, seq_list)

    SeqIO.write(seq_list, outpath, "fasta")


def randstring(length=10):
    valid_letters = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    return "".join((random.choice(valid_letters) for i in range(length)))


def add_from_csv(df, add_df, match="Entry"):
    add_df = pd.read_csv(add_df)
    merged_df = pd.merge(df, add_df, how="left", on=[match], suffixes=["", "_r"])

    return merged_df


def get_entry_ids_from_fasta(fasta_path, alphabet="ABCDEFGHIJKLMNOPQRSTUVWXYZ-"):

    #TODO: Change to biopython is not working to work

    seqs = SeqIO.read(fasta_path, 'fasta')

    return [seq.name for seq in seqs]
