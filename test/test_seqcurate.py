import scripts.seqcurate as sc
import pytest
import os
import tempfile
import pandas as pd
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


# Test for add_from_csv function
def test_add_from_csv(tmp_path):
    # Create a sample DataFrame
    df = pd.DataFrame({"info": ["A", "B", "C"], "data": ['1', '2', '3']})

    # Create a sample CSV file
    add_data = pd.DataFrame({"info": ["B", "C"], "additional_data": ["X", "Y"]})
    add_data.to_csv(tmp_path / "add_data.csv", index=False)

    # Call the add_from_csv function with the DataFrame and CSV file path
    merged_df = sc.add_from_csv(df, tmp_path / "add_data.csv", match="info")

    # Define the expected merged DataFrame
    expected_df = pd.DataFrame({
        "info": ["A", "B", "C"],
        "data": ['1', '2', '3'],
        "additional_data": [None, "X", "Y"]
    })

    # Check if the result matches the expected DataFrame
    assert merged_df.equals(expected_df)

# Test for get_entry_ids_from_fasta function
def test_get_entry_ids_from_fasta(tmp_path):
    # Create a sample FASTA file
    fasta_file = tmp_path / "sample.fasta"
    with open(fasta_file, "w") as f:
        f.write(">entry1\nACGT\n>entry2\nTGCA")

    entry_ids = sc.get_entry_ids_from_fasta(fasta_file)

    expected_ids = ["entry1", "entry2"]

    assert entry_ids == expected_ids


# Test for get_sequence_content_from_fasta function
def test_get_sequence_content_from_fasta(tmp_path):
    # Create a sample FASTA file
    fasta_file = tmp_path / "sample.fasta"
    with open(fasta_file, "w") as f:
        f.write(">entry1\nACGT\n>entry2\nTGCA")

    sequences = sc.get_sequence_content_from_fasta(fasta_file)

    expected_sequences = ["ACGT", "TGCA"]

    assert sequences == expected_sequences
