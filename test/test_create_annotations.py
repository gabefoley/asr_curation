import scripts.seqcurate as sc
import scripts.create_annotations as ca
import os


def test_create_annotations_main():
    input_file = "test/files/test.fasta"
    outpath = "test/files/create_annotations/output_test_validated_main.csv"

    if os.path.exists(outpath):
        os.remove(outpath)

    ca.create_initial_annotation(input_file, outpath)

    assert os.path.isfile(outpath)


def test_validate_ncbi_and_uniprot_sequences():
    # Tests a FASTA file that contains both NCBI and UniProt headers

    input_file = "test/files/create_annotations/ncbi_and_uniprot.fasta"
    outpath = "test/files/create_annotations/ncbi_and_uniprot.csv"

    if os.path.exists(outpath):
        os.remove(outpath)

    ca.create_initial_annotation(input_file, outpath)

    assert os.path.isfile(outpath)
