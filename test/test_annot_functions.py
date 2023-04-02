import scripts.annot_functions as an
import pandas as pd
import os
import pytest
import scripts.seqcurate as sc


def test_get_amino_acids():
    pos = an.get_amino_acids("PPGP", 0)
    assert pos == "P"


def test_get_amino_acids_multiple():
    pos = an.get_amino_acids("PPGP", 0, 2)
    assert pos == "PG"


def test_get_binding_pos():
    ft_binding = 'BINDING 123..130; /ligand="NADP(+)"; /ligand_id="ChEBI:CHEBI:58349"; /evidence="ECO:0000250"; BINDING 156..161; /ligand="NADP(+)"; /ligand_id="ChEBI:CHEBI:58349"; /evidence="ECO:0000250"; BINDING 195..199; /ligand="NADP(+)"; /ligand_id="ChEBI:CHEBI:58349"; /evidence="ECO:0000250"; BINDING 309; /ligand="Mg(2+)"; /ligand_id="ChEBI:CHEBI:18420"; /ligand_label="1"; /evidence="ECO:0000255|PROSITE-ProRule:PRU01198"; BINDING 309; /ligand="Mg(2+)"; /ligand_id="ChEBI:CHEBI:18420"; /ligand_label="2"; /evidence="ECO:0000255|PROSITE-ProRule:PRU01198"; BINDING 313; /ligand="Mg(2+)"; /ligand_id="ChEBI:CHEBI:18420"; /ligand_label="1"; /evidence="ECO:0000255|PROSITE-ProRule:PRU01198"; BINDING 486; /ligand="Mg(2+)"; /ligand_id="ChEBI:CHEBI:18420"; /ligand_label="2"; /evidence="ECO:0000255|PROSITE-ProRule:PRU01198"; BINDING 490; /ligand="Mg(2+)"; /ligand_id="ChEBI:CHEBI:18420"; /ligand_label="2"; /evidence="ECO:0000255|PROSITE-ProRule:PRU01198"; BINDING 512; /ligand="substrate"; /evidence="ECO:0000255|PROSITE-ProRule:PRU01198'


def test_add_lab_annotations_correct():
    # Check that the column gets added okay

    print("in here")
    print(os.getcwd())

    annot_df = pd.read_csv("test/files/test_add_lab_annot_df.csv")
    filepath = "test/files/test_add_lab_correct_lab_df.csv"
    annot_df = an.add_lab_annotations(annot_df, filepath)

    assert annot_df.loc[annot_df["accession"] == "P9XVRR", "lab_km"].values[0] == 0.4


def test_add_lab_annotations_duplicate_column_between_annot_and_lab():
    # Should throw an error and not proceed

    annot_df = pd.read_csv("test/files/test_add_lab_annot_df.csv")
    filepath = "test/files/test_add_lab_dups_between_lab_df.csv"
    with pytest.raises(
        ValueError, match="Duplicate column between lab and existing annotations"
    ):
        annot_df = an.add_lab_annotations(annot_df, filepath)


def test_add_lab_annotations_no_accession_column():
    # Should throw an error and not proceed
    annot_df = pd.read_csv("test/files/test_add_lab_annot_df.csv")
    filepath = "test/files/test_add_lab_no_accession_lab_df.csv"

    with pytest.raises(ValueError, match="Lab annotations are missing accession field"):
        annot_df = an.add_lab_annotations(annot_df, filepath)


def test_add_lab_annotations_no_sequence_column():
    # Should throw an error and not proceed
    annot_df = pd.read_csv("test/files/test_add_lab_annot_df.csv")
    filepath = "test/files/test_add_lab_no_sequence_lab_df.csv"

    with pytest.raises(ValueError, match="Lab annotations are missing sequence field"):
        annot_df = an.add_lab_annotations(annot_df, filepath)


def test_add_lab_annotations_duplicate_column_in_lab():
    # Should throw an error and not proceed
    annot_df = pd.read_csv("test/files/test_add_lab_annot_df.csv")
    filepath = "test/files/test_add_lab_dups_within_lab_df.csv"

    with pytest.raises(
        ValueError, match="Lab annotations contain multiple identically named columns"
    ):
        annot_df = an.add_lab_annotations(annot_df, filepath)


def test_add_lab_annotations_overwrite_different_value():
    # Should throw a warning and proceed
    annot_df = pd.read_csv("test/files/test_add_lab_annot_df_with_values_already.csv")
    filepath = "test/files/test_add_lab_overwrite_diff_lab_df.csv"

    with pytest.warns(UserWarning, match=r"^Lab annotations are overwriting values"):
        annot_df = an.add_lab_annotations(annot_df, filepath)


def test_add_lab_annotations_overwrite_add_to_value():
    # Should throw a warning and proceed
    annot_df = pd.read_csv("test/files/test_add_lab_annot_df_with_values_already.csv")
    filepath = "test/files/test_add_lab_overwrite_add_lab_df.csv"

    with pytest.warns(UserWarning, match=r"^Lab annotations are adding to values"):
        annot_df = an.add_lab_annotations(annot_df, filepath)


def test_add_lab_annotations_missing_values():
    # Should throw a warning and proceed
    annot_df = pd.read_csv("test/files/test_add_lab_annot_df_with_values_already.csv")
    filepath = "test/files/test_add_lab_missing_values_lab_df.csv"

    with pytest.warns(UserWarning, match=r"^Lab annotations are overwriting values"):
        annot_df = an.add_lab_annotations(annot_df, filepath)


def test_add_lab_annotations_different_sequence():
    # Should throw an error and not proceed
    annot_df = pd.read_csv("test/files/test_add_lab_annot_df.csv")
    filepath = "test/files/test_add_lab_different_sequence.csv"

    with pytest.raises(
        ValueError,
        match="Lab annotations contain different sequence to existing annotations",
    ):
        annot_df = an.add_lab_annotations(annot_df, filepath)


# TEST TRACK RESIDUES


def test_track_residues():
    unaligned_pos = [89, 296, 297]
    seq_id = "Q5HVD9"
    tag = "activation_critical"

    filepath = "test/files/track_residues.aln"

    align_df = sc.get_sequence_df(filepath, alignment=True, ancestor=False)

    aligned_seq = "----------------------------------------------------------MAITVYY---------------------------------------------------------DKDCDLNLIKS------------------------------KKVAIIGFGSQGHAHAMNLRD------NGVNVTIGLREGS-----VSAVKAKNAGF------EVMSVSEASKIADVIMILAPDEIQADIFNVEIKPNLSEGKAIAFAHGFNIHYG---QIVVPKGVDVIMIAPKAPGHTVRNEFTLGG-----GTPCLIAIH--QDESKNAKNLALSYASAIGGGRTGIIETTFKAETETDLFGEQAVLC-----------------------------------------GGLS----------------------------------------------------------------------------------------------------------------------------ALIQAG----FETLVEAGYEPEMAYFECLHE-MKLIVDLIYQGGIADMRYSISNTAEYGDYITGPK---IITEETKKAMKGVLKDIQNGV---------FAKDFILER-RAG-FARMHAERKNMNDSLIEKTGRNLRAMMPWISAKKLVDK------DKN------"

    align_df = an.track_residues(align_df, seq_id, aligned_seq, tag, *unaligned_pos)

    assert align_df.loc[align_df["accession"] == seq_id, tag].values[0] == "DRR"


# TEST ALIGNMENT ANNOTATION


def test_create_annotated_alignment():
    outpath = "test/files/domains.html"

    # Delete output

    if os.path.isfile(outpath):
        os.remove(outpath)

    align_df = pd.read_csv("test/files/test_annotated_alignment_df.csv")

    align_df["combined_domain_bounds"] = align_df.apply(
        lambda x: an.create_domain_bounds(x.ft_domain), axis=1
    )

    boundary_dict = pd.Series(
        align_df.combined_domain_bounds.values, index=align_df.accession
    ).to_dict()

    # Assert the boundary dict contains a correct key
    assert "Q02138" in boundary_dict.keys()

    an.create_annotated_alignment(align_df, boundary_dict, "test/files/domains.html")

    # Assert that the HTML alignment gets created
    assert os.path.isfile(outpath)




def test_get_motif_indexes():
    seq = "PGPGHGNNG"
    motif = "G.G..G"
    multi_motif = "G.G"
    pos = an.get_motif_indexes(seq, motif)

    assert pos[0] == (3, 9)

    multi_pos = an.get_motif_indexes(seq, multi_motif)

    assert multi_pos[0] == (1, 4)
    assert multi_pos[1] == (3, 6)


def test_unaligned_positions_align():
    # If we have a set of positions derived from an annotation, check that for a given alignment they actually align
    # This tests when they do align
    align_df = pd.read_csv("test/files/annot_functions/test_positions_in_alignment.csv")

    correct_align_ids = ["P05793", "P9WKJ7"]

    incorrect_align_ids = []

    correctly_aligning_df = align_df[align_df["accession"].isin(correct_align_ids)]

    print(correctly_aligning_df)

    # print (align_df)

    target_pos = align_df.query("accession=='P05793'")["Binding_positions_extracted"]

    print(target_pos)

    align_df["aligned_pos"] = align_df.apply(
        lambda x: an.check_if_positions_align_with_target(x.ft_domain), axis=1
    )

def test_indexes_align_in_alignment():
    from ast import literal_eval
    align_df = pd.read_csv("test/files/annot_functions/test_motifs_in_alignment2.csv")

    print (align_df[['MOTIF_GxGxxG_indexes', 'GxGxxG_motif_from_alignment']])



    an.check_if_positions_align_with_target()


