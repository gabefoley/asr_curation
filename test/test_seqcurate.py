import scripts.seqcurate as sc
import pytest
import os
import tempfile
@pytest.fixture
def temporary_fasta_file():
    # Create a temporary FASTA file with test data
    fasta_content = ">seq1\nPPGN\n>seq2\nPPGD"
    fasta_file = tempfile.NamedTemporaryFile(delete=False, mode='w')
    fasta_file.write(fasta_content)
    fasta_file.close()
    yield fasta_file.name

    # Remove the temporary FASTA file
    os.remove(fasta_file.name)

def test_get_sequence_content_from_fasta(temporary_fasta_file):
    # Call the function with the temporary FASTA file
    sequences = sc.get_sequence_content_from_fasta(temporary_fasta_file)

    # Check if the returned sequences match the expected ones
    expected_sequences = ["PPGN", "PPGD"]
    assert sequences == expected_sequences

if __name__ == '__main__':
    pytest.main()

def test_randstring():
    default_rand = sc.randstring()
    five_rand = sc.randstring(5)

    assert len(default_rand) == 10
    assert len(five_rand) == 5


def test_get_sequence_df():
    path = "test/files/test.fasta"
    df = sc.get_sequence_df(path)

    # Test that the columns get added correctly
    assert [x for x in df.columns] == [
        "info",
        "truncated_info",
        "extracted_id",
        "extracted_name",
        "sequence",
        "original_fasta",
    ]

    # Size should be 10 rows x 6 columns = 60
    assert df.size == 60


def test_get_sequence_df_multiple_paths():
    path_1 = "test/files/test_1.fasta"
    path_2 = "test/files/test_2.fasta"
    df = sc.get_sequence_df(path_1, path_2)

    # Test that the columns get added correctly
    assert [x for x in df.columns] == [
        "info",
        "truncated_info",
        "extracted_id",
        "extracted_name",
        "sequence",
        "original_fasta",
    ]

    # Size should be 5 rows x 6 columns = 30
    assert df.size == 30
