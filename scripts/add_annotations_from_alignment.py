import annot_functions as an
import pandas as pd
from Bio import AlignIO
import seqcurate as sc

# def get_content_at_position(sequence, *positions)
#     alignment_content = []
#     for pos in positions:


df = pd.read_csv(snakemake.input.csv)

# Get the equivalent positions of the first sequence's binding positions
# aln = sequence.readFastaFile(snakemake.input.aln)

aln = AlignIO.read(snakemake.input.aln, format="fasta")

# print (aln)

# aln_dict = {seq.name: str(seq.seq) for seq in aln}

# align_df = pd.read_csv(snakemake.input.csv)


align_df = sc.get_sequence_df(snakemake.input.aln, alignment=True, ancestor=True)

merged_df = pd.merge(
        df,
        align_df,
        left_on=["accession"],
        right_on=["accession"],
        suffixes=["", "_r"],
    )

merged_df.to_csv(snakemake.output.csv, index=False)




# print (aln_dict)

# if "ec_1_1_1_86" in snakemake.wildcards.dataset:
#     print("Adding KARI specific annotations based on the alignment")

#     # print("Adding dimer / dodecamer key residues")

#     # dimer_df = align_df[align_df["Entry"] == "sp|A6VCE7|ILVC_PSEA7"]

#     # print("The reference binding positions we will use are 256, 288, and 299")

#     # ref_binding_pos = [256, 288, 299]

#     # aln_binding_pos = an.get_aligned_positions(
#     #     align_df[align_df["Entry"] == "sp|A6VCE7|ILVC_PSEA7"],
#     #     aln_dict["sp|A6VCE7|ILVC_PSEA7"],
#     #     *ref_binding_pos,
#     # )

#     # print(
#     #     f"The alignment binding positions for dimer / dodecamr are {' '.join([str(x) for x in aln_binding_pos])}"
#     # )

#     # print("need to minus one")

#     # aln_binding_pos = [x - 1 for x in aln_binding_pos]

#     # align_df["Dimer_alignment_residues"] = align_df.apply(
#     #     lambda row: an.get_amino_acids(
#     #         aln_dict[row["Entry"]], *aln_binding_pos
#     #     ),
#     #     axis=1,
#     # )

#     # print("Adding binding site entries")

#     # print("Reference binding site entry")
#     # print("Set it to E COLI if it is there")
#     # align_df[align_df["Entry"] == "D0WGK0"]
#     # print(align_df["Entry"].iloc[0])
#     # print(align_df["Binding_site"].iloc[0])
#     # print(type(align_df["Binding_site"].iloc[0]))

#     # ecoli = align_df[align_df['Entry'] == 'D0WGK0']

#     # print ('exoli')
#     # print (ecoli)
#     # print (ecoli[['Binding_site']].str)
#     # print ('troubles')

#     # ref_binding_pos = an.get_binding_pos(ecoli[['Binding_site']])

#     # else:
#     # Get the binding positions of the first sequence as a reference
#     ref_binding_pos = an.get_binding_pos(align_df["Binding_site"].iloc[0])

#     print(ref_binding_pos)

#     print(
#         f'The reference binding positions we will use are {" ".join([str(x) for x in ref_binding_pos])}'
#     )

#     print(type(ref_binding_pos[0]))

#     # Find where in the alignment the binding positions are
#     aln_binding_pos = an.get_aligned_positions(
#         align_df["Entry"].iloc[0],
#         aln_dict[align_df["Entry"].iloc[0]],
#         *ref_binding_pos,
#     )

#     print(
#         f'Mapping this positions to alignment columns gives these positions {" ".join([str(x) for x in aln_binding_pos])}'
#     )

#     # align_df['Binding_positions_alignment'] = align_df.apply(
#     #     lambda row: get_aligned_positions(row['Entry'],aln_dict[row['Entry']].sequence, *ref_binding_pos), axis=1)
#     align_df["Binding_positions_alignment_residues"] = align_df.apply(
#         lambda row: an.get_amino_acids(
#             aln_dict[row["Entry"]], *aln_binding_pos
#         ),
#         axis=1,
#     )

#     # align_df['Acidic_Binding_alignment'] = align_df.apply(lambda row: an.check_binding_for_acidic(row['Sequence'],
#     #                                                                                               an.get_binding_pos(
#     #                                                                                                   row[
#     #                                                                                                       'Binding_site']) if pd.notnull(
#     #                                                                                                   row[
#     #                                                                                                       'Binding_site']) else 'No_binding_positions'),
#     #                                                       axis=1)

#     align_df["Acidic_Binding_alignment"] = align_df.apply(
#         lambda row: an.check_sequence_for_acidic(
#             row["Binding_positions_alignment_residues"]
#         ),
#         axis=1,
#     )

# print (align_df.columns)

# align_df.to_csv(snakemake.output[0], index=False)


# frames = [df, align_df]
# merge_df = pd.concat(frames)


# merge_df.to_csv(snakemake.output.csv, index=False)
