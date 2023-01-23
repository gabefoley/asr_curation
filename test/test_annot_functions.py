import scripts.annot_functions as an
import pandas as pd
import os
import pytest
import scripts.seqcurate as sc

def test_get_amino_acids():

    pos = an.get_amino_acids('PPGP', 0)
    assert pos == 'P'


def test_get_amino_acids_multiple():

    pos = an.get_amino_acids('PPGP', 0, 2)
    assert pos == 'PG'

def test_get_binding_pos():
    ft_binding = 'BINDING 123..130; /ligand="NADP(+)"; /ligand_id="ChEBI:CHEBI:58349"; /evidence="ECO:0000250"; BINDING 156..161; /ligand="NADP(+)"; /ligand_id="ChEBI:CHEBI:58349"; /evidence="ECO:0000250"; BINDING 195..199; /ligand="NADP(+)"; /ligand_id="ChEBI:CHEBI:58349"; /evidence="ECO:0000250"; BINDING 309; /ligand="Mg(2+)"; /ligand_id="ChEBI:CHEBI:18420"; /ligand_label="1"; /evidence="ECO:0000255|PROSITE-ProRule:PRU01198"; BINDING 309; /ligand="Mg(2+)"; /ligand_id="ChEBI:CHEBI:18420"; /ligand_label="2"; /evidence="ECO:0000255|PROSITE-ProRule:PRU01198"; BINDING 313; /ligand="Mg(2+)"; /ligand_id="ChEBI:CHEBI:18420"; /ligand_label="1"; /evidence="ECO:0000255|PROSITE-ProRule:PRU01198"; BINDING 486; /ligand="Mg(2+)"; /ligand_id="ChEBI:CHEBI:18420"; /ligand_label="2"; /evidence="ECO:0000255|PROSITE-ProRule:PRU01198"; BINDING 490; /ligand="Mg(2+)"; /ligand_id="ChEBI:CHEBI:18420"; /ligand_label="2"; /evidence="ECO:0000255|PROSITE-ProRule:PRU01198"; BINDING 512; /ligand="substrate"; /evidence="ECO:0000255|PROSITE-ProRule:PRU01198'


def test_add_lab_annotations_correct():
    # Check that the column gets added okay

    annot_df = pd.read_csv("./files/test_add_lab_annot_df.csv")
    filepath = "./files/test_add_lab_correct_lab_df.csv"
    annot_df = an.add_lab_annotations(annot_df, filepath)

    assert (annot_df.loc[annot_df['accession'] == 'P9XVRR', 'lab_km'].values[0] == 0.4)

def test_add_lab_annotations_duplicate_column_between_annot_and_lab():
    # Should throw an error and not proceed

    annot_df = pd.read_csv("./files/test_add_lab_annot_df.csv")
    filepath = "./files/test_add_lab_dups_between_lab_df.csv"
    with pytest.raises(ValueError, match='Duplicate column between lab and existing annotations'):
        annot_df = an.add_lab_annotations(annot_df, filepath)


def test_add_lab_annotations_no_accession_column():
    # Should throw an error and not proceed
    annot_df = pd.read_csv("./files/test_add_lab_annot_df.csv")
    filepath = "./files/test_add_lab_no_accession_lab_df.csv"

    with pytest.raises(ValueError, match='Lab annotations are missing accession field'):
        annot_df = an.add_lab_annotations(annot_df, filepath)

def test_add_lab_annotations_no_sequence_column():
    # Should throw an error and not proceed
    annot_df = pd.read_csv("./files/test_add_lab_annot_df.csv")
    filepath = "./files/test_add_lab_no_sequence_lab_df.csv"

    with pytest.raises(ValueError, match='Lab annotations are missing sequence field'):
        annot_df = an.add_lab_annotations(annot_df, filepath)

def test_add_lab_annotations_duplicate_column_in_lab():
    # Should throw an error and not proceed
    annot_df = pd.read_csv("./files/test_add_lab_annot_df.csv")
    filepath = "./files/test_add_lab_dups_within_lab_df.csv"

    with pytest.raises(ValueError, match='Lab annotations contain multiple identically named columns'):
        annot_df = an.add_lab_annotations(annot_df, filepath)


def test_add_lab_annotations_overwrite_different_value():
    # Should throw a warning and proceed
    annot_df = pd.read_csv("./files/test_add_lab_annot_df_with_values_already.csv")
    filepath = "./files/test_add_lab_overwrite_diff_lab_df.csv"

    with pytest.warns(UserWarning, match=r'^Lab annotations are overwriting values'):
        annot_df = an.add_lab_annotations(annot_df, filepath)


def test_add_lab_annotations_overwrite_add_to_value():
    # Should throw a warning and proceed
    annot_df = pd.read_csv("./files/test_add_lab_annot_df_with_values_already.csv")
    filepath = "./files/test_add_lab_overwrite_add_lab_df.csv"

    with pytest.warns(UserWarning, match=r'^Lab annotations are adding to values'):

        annot_df = an.add_lab_annotations(annot_df, filepath)

def test_add_lab_annotations_missing_values():
    # Should throw a warning and proceed
    annot_df = pd.read_csv("./files/test_add_lab_annot_df_with_values_already.csv")
    filepath = "./files/test_add_lab_missing_values_lab_df.csv"

    with pytest.warns(UserWarning, match=r'^Lab annotations are overwriting values'):

        annot_df = an.add_lab_annotations(annot_df, filepath)



def test_add_lab_annotations_different_sequence():
    # Should throw an error and not proceed
    annot_df = pd.read_csv("./files/test_add_lab_annot_df.csv")
    filepath = "./files/test_add_lab_different_sequence.csv"

    with pytest.raises(ValueError, match='Lab annotations contain different sequence to existing annotations'):
        annot_df = an.add_lab_annotations(annot_df, filepath)



# TEST TRACK RESIDUES

def test_track_residues():
    unaligned_pos = [89, 296, 297]
    seq_id = 'Q5HVD9'
    tag = 'activation_critical'

    filepath = "files/track_residues.aln"

    align_df = sc.get_sequence_df(filepath, alignment=True, ancestor=False)

    aligned_seq = '----------------------------------------------------------MAITVYY---------------------------------------------------------DKDCDLNLIKS------------------------------KKVAIIGFGSQGHAHAMNLRD------NGVNVTIGLREGS-----VSAVKAKNAGF------EVMSVSEASKIADVIMILAPDEIQADIFNVEIKPNLSEGKAIAFAHGFNIHYG---QIVVPKGVDVIMIAPKAPGHTVRNEFTLGG-----GTPCLIAIH--QDESKNAKNLALSYASAIGGGRTGIIETTFKAETETDLFGEQAVLC-----------------------------------------GGLS----------------------------------------------------------------------------------------------------------------------------ALIQAG----FETLVEAGYEPEMAYFECLHE-MKLIVDLIYQGGIADMRYSISNTAEYGDYITGPK---IITEETKKAMKGVLKDIQNGV---------FAKDFILER-RAG-FARMHAERKNMNDSLIEKTGRNLRAMMPWISAKKLVDK------DKN------'

    align_df = an.track_residues(align_df, seq_id, aligned_seq, tag, *unaligned_pos)

    assert(align_df.loc[align_df['accession'] == seq_id, tag].values[0] == 'DRR')


# TEST ALIGNMENT ANNOTATION

def test_create_annotated_alignment():

    outpath = "files/domains.html"

    # Delete output

    if os.path.isfile(outpath):
        os.remove(outpath)


    align_df = pd.read_csv("files/test_annotated_alignment_df.csv")

    align_df['combined_domain_bounds'] = align_df.apply(lambda x: an.create_domain_bounds(x.ft_domain), axis=1)

    boundary_dict = pd.Series(align_df.combined_domain_bounds.values, index=align_df.accession).to_dict()

    # Assert the boundary dict contains a correct key
    assert ('Q02138' in boundary_dict.keys())

    an.create_annotated_alignment(align_df, boundary_dict, "files/domains.html")

    # Assert that the HTML alignment gets created
    assert (os.path.isfile(outpath))
