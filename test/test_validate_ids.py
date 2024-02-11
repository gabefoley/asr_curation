import os
import scripts.validate_ids as vi


# Tests to check the validation of IDs
def test_validate_ids_main():
    outpath = "test/files/output_test_validated_main.csv"

    if os.path.exists(outpath):
        os.remove(outpath)

    vi.all_ids_lookup("test/files/test_original.csv", outpath, verbose=True)

    assert os.path.isfile(outpath)




# def test_creating_and_validating_diverse_headers():
#     # This test creates a DataFrame that contains full NCBI accession header, full UniProt header, gi header
#     # (as downloaded from Conserved Domain Database alignment) and single NCBI identifier and single UniProt identifier
#
#
#     seq_outpath = "test/files/validate_ids/output_test_diverse_5.csv"
#     validated_outpath = "test/files/output_test_diverse_5_validated.csv"
#
#     if os.path.exists(seq_outpath):
#         os.remove(seq_outpath)
#
#     if os.path.exists(validated_outpath):
#         os.remove(validated_outpath)
#
#     # Create a sequence df
#
#     #


# def test_validating_ncbi_and_uniprot_headers():
#     # This test creates a DataFrame that contains full NCBI accession header, full UniProt header, gi header
#     # (as downloaded from Conserved Domain Database alignment) and single NCBI identifier and single UniProt identifier
#
#     # seq_outpath = "test/files/validate_ids/ncbi_and_uniprot.csv"
#     # validated_outpath = "test/files/output_test_diverse_5_validated.csv"
#     #
#
#     outpath = "test/files/validate_ids/ncbi_and_uniprot_validated.csv"
#
#     if os.path.exists(outpath):
#         os.remove(outpath)
#
#     vi.all_ids_lookup("test/files/validate_ids/ncbi_and_uniprot.csv", outpath)
#
#     assert os.path.isfile(outpath)
