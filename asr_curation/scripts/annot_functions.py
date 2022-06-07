import pandas as pd


def annotate_motif(df, motif):

    df[f"MOTIF_{motif}"] = df["Sequence"].dropna().str.contains(motif)

    return df


def annotate_nonAA(df):
    # Does the sequences have non amino acid characters in it

    non_AA = "B|J|O|U|X|Z"
    df["Non_AA_Character"] = df["Sequence"].dropna().str.contains(non_AA)

    return df


def annotate_AA(df):
    # Does the sequences have non amino acid characters in it

    # print ('lets check')
    # print (df['Sequence'])
    non_AA = "B|J|O|U|X|Z"
    df["AA_Character"] = ~(df["Sequence"].dropna().str.contains(non_AA, na=None))

    # booleanDictionary = {True: 'TRUE', False: 'FALSE'}
    # df = df.replace(booleanDictionary)

    # print (df['Non_AA_Character'])

    return df


def annotate_sp_tr(df):
    # Is the sequence from SwissProt or TrEMBL
    df.loc[df["Entry"].str.startswith("sp"), "UniProt_DB"] = "SwissProt"
    df.loc[df["Entry"].str.startswith("tr"), "UniProt_DB"] = "TrEMBL"

    return df


def get_final_pos(sequence, pos, curr_idx, next_idx):
    # Need to find the position that 1) isn't a gap and 2) takes into account all of the
    # offset implied by previous gaps in the sequence

    # print (sequence)
    # print (curr_idx)

    # If the content at this index is a gap, proceed to the next actual position
    while curr_idx < len(sequence) and sequence[curr_idx] == "-":
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


def get_aligned_positions(entry, sequence, *positions):
    # print (f'\nSeq name is {entry}')
    sequence = "".join(sequence)
    # print(sequence)
    # print(len(sequence))
    aligned_positions = []

    # print (f'\nSequence is {sequence}')

    # print(positions)

    for pos in positions:
        # print(pos)
        # Get the current index
        curr_idx = pos

        # Get offset implied by first position
        offset = sequence[0:curr_idx].count("-")

        # Get next position based on the first position and the gap offset
        next_idx = curr_idx + offset
        # print(curr_idx)
        # print(offset)
        # print(next_idx)
        # Search to find the final position
        final_pos = get_final_pos(sequence, pos, curr_idx, next_idx)
        aligned_positions.append(final_pos)


    return aligned_positions


##### KARI SPECIFIC ####


def classify_KARI(features):
    # print (features)
    if "Domain" in features:
        domain_num = features.split("Domain")[1].split(";")[0]
        if "2" in domain_num:
            return "Class_I"
        elif "3" in domain_num:
            return "Class_II"
        else:
            return "Different_domain_number"
    else:
        return "No_domain_info"


def get_binding_pos(binding_sites):

    # print (binding_sites)
    if pd.notnull(binding_sites):
        bp = []
        #         print ('bs is')
        #         print (binding_sites)
        for site in binding_sites.split(";"):
            #         print ('***')
            #         print (site)
            if site.strip().startswith("BINDING"):
                found_pos = int(site.split("BINDING")[1]) - 1
            #         if site.strip().startswith('/note='):
            #             print (site)
            #             print(site.split('/note='))
            #             print(site.split('/note=')[1])
            if (
                site.strip().startswith("/note=")
                and site.split("/note=")[1].startswith('"NADP"')
            ):
                bp.append(found_pos)

        # print (bp)
        return bp
    else:
        return []


def get_amino_acids(seq, *pos):
    return "".join([seq[bp] for bp in pos])


def check_sequence_for_acidic(seq):
    if "E" in seq or "D" in seq:
        return "Acidic_residue"
    else:
        return "No_Acidic_residue"


def check_binding_for_acidic(seq, bind_pos):
    if bind_pos and bind_pos != "No_binding_positions":
        binding_aa = get_amino_acids(seq, *bind_pos)
        return check_sequence_for_acidic(binding_aa)


def classify_loop_length(bind_pos):

    # Offset accounts for the fact that we need to get the total number of positions and also account for the first position in the loop which isn't a binding position

    if bind_pos and bind_pos != "No_binding_positions":

        offset = 2
        return bind_pos[-1] - bind_pos[0] + offset
    else:
        return "No_binding_positions"


# ref_df.apply(lambda row : get_binding_aa(row['Sequence'],
#                      get_binding_pos(row['feature(BINDING SITE)'])), axis = 1)
