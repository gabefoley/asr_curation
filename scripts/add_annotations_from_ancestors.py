import annot_functions as an
import pandas as pd
import seqcurate as sc


align_df = pd.read_csv(snakemake.input.csv)


anc_df = sc.get_sequence_df(snakemake.input.aln, alignment=True, ancestor=True)

merged_df = pd.merge(
        align_df,
        anc_df,
        left_on=["accession"],
        right_on=["accession"],
        suffixes=["", "_r"],
    )

merged_df.to_csv(snakemake.output.csv, index=False)


# if "ec_1_1_1_86" in snakemake.wildcards.dataset:
#     print("Adding KARI specific annotations based on the ancestors")

#     # Get the binding positions of the first sequence as a reference
#     ref_binding_pos = an.get_binding_pos(align_df["Binding_site"].iloc[0])

#     align_df[align_df["Entry"] == "D0WGK0"]

#     print(ref_binding_pos)

#     print(
#         f'The reference binding positions we will use are {" ".join([str(x) for x in ref_binding_pos])}'
#     )

#     print(type(ref_binding_pos[0]))

#     # Find where in the alignment the binding positions are
#     aln_binding_pos = an.get_aligned_positions(
#         align_df["Entry"].iloc[0], align_df["Sequence"].iloc[0], *ref_binding_pos
#     )

#     print(
#         f'Mapping this positions to alignment columns gives these positions {" ".join([str(x) for x in aln_binding_pos])}'
#     )

#     # align_df['Binding_positions_alignment'] = align_df.apply(
#     #     lambda row: get_aligned_positions(row['Entry'],aln_dict[row['Entry']].sequence, *ref_binding_pos), axis=1)

#     aln_binding_pos = [74, 83, 85]

#     anc_df["Binding_positions_alignment_residues_ancestor"] = anc_df.apply(
#         lambda row: an.get_amino_acids(row["Sequence_aligned"], *aln_binding_pos),
#         axis=1,
#     )

#     anc_df["Acidic_Binding_alignment"] = anc_df.apply(
#         lambda row: an.check_sequence_for_acidic(
#             row["Binding_positions_alignment_residues_ancestor"]
#         ),
#         axis=1,
#     )


# frames = [align_df, anc_df]
# merge_df = pd.concat(frames)


# merge_df.to_csv(snakemake.output.csv, index=False)


#     # else:
#     # Get the binding positions of the first sequence as a reference
#     ref_binding_pos = an.get_binding_pos(align_df['Binding_site'].iloc[0])

#     print(ref_binding_pos)

#     print(f'The reference binding positions we will use are {" ".join([str(x) for x in ref_binding_pos])}')

#     print(type(ref_binding_pos[0]))

#     # Find where in the alignment the binding positions are
#     aln_binding_pos = an.get_aligned_positions(align_df['Entry'].iloc[0],aln_dict[align_df['Entry'].iloc[0]].sequence, *ref_binding_pos)

#     print(
#         f'Mapping this positions to alignment columns gives these positions {" ".join([str(x) for x in aln_binding_pos])}')

#     # align_df['Binding_positions_alignment'] = align_df.apply(
#     #     lambda row: get_aligned_positions(row['Entry'],aln_dict[row['Entry']].sequence, *ref_binding_pos), axis=1)
#     align_df['Binding_positions_alignment_residues'] = align_df.apply(
#         lambda row: an.get_amino_acids(aln_dict[row['Entry']].sequence, *aln_binding_pos), axis=1)

#     # align_df['Acidic_Binding_alignment'] = align_df.apply(lambda row: an.check_binding_for_acidic(row['Sequence'],
#     #                                                                                               an.get_binding_pos(
#     #                                                                                                   row[
#     #                                                                                                       'Binding_site']) if pd.notnull(
#     #                                                                                                   row[
#     #                                                                                                       'Binding_site']) else 'No_binding_positions'),
#     #                                                       axis=1)


#     align_df['Acidic_Binding_alignment'] = align_df.apply(lambda row: an.check_sequence_for_acidic(row['Binding_positions_alignment_residues']), axis=1)


# align_df.to_csv(snakemake.output[0], index=False)
