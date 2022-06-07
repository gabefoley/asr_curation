import pandas as pd
import sequence
import random

all_alpha = sequence.Alphabet('ABCDEFGHIJKLMNOPQRSTUVWXYZ')

def get_sequence_df(*fasta_paths, drop_duplicates=True, alignment=False, ancestor=False, alphabet='ABCDEFGHIJKLMNOPQRSTUVWXYZ-'):
    
    seq_list = []
    duplicates = {}

    alpha = sequence.Alphabet(alphabet)

    cols=['Entry', 'Info', 'Truncated_Info', 'Sequence', 'Original_FASTA']

    if alignment:
        cols.append('Sequence_aligned')

    if ancestor:
        cols.append('Sequence_aligned_ancestor')


    for fasta_path in fasta_paths:

        # Load FASTA file
        seqs = sequence.readFastaFile(fasta_path, alpha)

        # Add to annotation file
        for seq in seqs:
            if seq.name in duplicates:
                print(f'DUPLICATE:{seq.name} is in {duplicates[seq.name]} and {fasta_path}\n')
            else:
                duplicates[seq.name] = fasta_path


            curr_seq = [seq.name, seq.info, seq.info.split(" ")[0], "".join([x for x in seq.sequence if x != "-"]) if len(seq.sequence) > 0 else None, fasta_path]
            
            if alignment:
                aligned_seq = "".join(x for x in seq.sequence)
                curr_seq.append(aligned_seq)
                print ('it is ')
                # print (aligned_seq)
                # print (curr_seq)
                print (len(curr_seq))

            if ancestor:
                curr_seq.append(aligned_seq)


            seq_list.append(curr_seq)

    # print (seq_list)
    print (len(seq_list))
    print (len(seq_list[0]))
    print (cols)
    df = pd.DataFrame(seq_list, columns=cols)

    if drop_duplicates:
        df = df.drop_duplicates(subset='Entry', keep='first')

    # Drop the sequence column if there are no sequences (i.e. if we just added a list of identifiers)
    nan_value = float("NaN")

    df.replace("", nan_value, inplace=True)
      
    df.dropna(how='all', axis=1, inplace=True)



    return df


def annotate_col_from_dict(df, col, annot_dict, match='Entry'):
    df.loc[df[match].isin(annot_dict.keys()), col] = df[match].map(annot_dict)
    return df


def add_col_from_up_dict(df, cols_to_add, up_dict):
    for col in cols_to_add:
        if not col in df:
            df[col] = ""

    for name, annots in up_dict.items():

        for key in annots.keys():
            df.loc[df['Entry'].str.contains(name), key] = annots.get(key)

    return df


def annotate_col_from_other_df(df, other_df, col, match_from='Entry', match_to='Entry'):
    other_dict = dict(zip(other_df[match_from], other_df[col]))

    df = annotate_col_from_dict(df, col, other_dict, match_to)

    return df


def get_subset(df, *cols_dict, include=True):
    for col_dict in cols_dict:
        if type(col_dict) != dict:
            raise TypeError('col_dict must be a dictionary')
        for col, vals in col_dict.items():
            if vals == '*':  # Get all the entries with values that exist

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
        seq_list = [sequence.Sequence(sequence=r.Sequence, name=r.Entry, alphabet = all_alpha) for r in df.itertuples()]
    else:
        seq_list = [sequence.Sequence(sequence=r.Sequence, name=r.Info, alphabet = all_alpha) for r in df.itertuples()]
    
    sequence.writeFastaFile(outpath, seq_list)


def randstring(length=10):
    valid_letters = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    return ''.join((random.choice(valid_letters) for i in range(length)))


def add_from_csv(df, add_df, match='Entry'):
    add_df = pd.read_csv(add_df)
    merged_df = pd.merge(df, add_df, how = 'left', on=[match], suffixes=['', '_r'])

    return merged_df


def get_entry_ids_from_fasta(fasta_path, alphabet='ABCDEFGHIJKLMNOPQRSTUVWXYZ-'):

    alpha = sequence.Alphabet(alphabet)

    # Load FASTA file
    seqs = sequence.readFastaFile(fasta_path, alpha)

    return [seq.name for seq in seqs]
