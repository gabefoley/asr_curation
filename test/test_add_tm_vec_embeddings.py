import pytest
import pandas as pd
import numpy as np
import torch
from collections import defaultdict

# Import the functions from the module
from scripts.add_tm_vec_embeddings import calculate_embeddings, process_and_store_embeddings


# Create a pytest fixture for sample data
@pytest.fixture
def sample_data():
    data = {
        'info': ['Seq1', 'Seq2', 'Seq3'],
        'sequence': [
            'MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQANLVLNAIQRWQDMRRSTLSVDAAKTG',
            'MKAILVVLLYTFVTASPLTGLQDQISLTTRSYIAVNNKPTFFDIAVDKIWTQLSKLGSDF',
            'MNKLIYLLVTALLLVTCSSPLGQYQEDYNTYRSITNDTNPTDFFDVLEKIWEQLSKLSSD'
        ]
    }
    df = pd.DataFrame(data)
    return df


@pytest.fixture
def model_config():
    return "/Users/uqgfoley/Dropbox/Code/Python_Workspace/asr_curation_data/workflows/2024_mbl_galaxy/custom_annotations/tm_vec_cath_model_params.json"

@pytest.fixture
def model_checkpoint():
    return "/Users/uqgfoley/Dropbox/Code/Python_Workspace/asr_curation_data/workflows/2024_mbl_galaxy/custom_annotations/tm_vec_cath_model.ckpt"


def test_unique_embeddings(sample_data, model_checkpoint, model_config):
    model_name = 'Prot_T5'

    # Process embeddings
    df_with_embeddings = process_and_store_embeddings(
        sample_data, model_name, model_checkpoint, model_config
    )

    print (df_with_embeddings)

    # Extract embeddings for comparison
    embeddings = df_with_embeddings[f"{model_name} Embed Encoded"].tolist()

    print (embeddings)

    # Ensure embeddings are not empty
    assert len(embeddings[0]) == len(sample_data), "Embeddings not generated for all sequences."

    # Check for uniqueness by converting lists to tuples and using a set
    embedding_tuples = [tuple(embed) for embed in embeddings]
    unique_embeddings = set(embedding_tuples)

    assert len(unique_embeddings) == len(embeddings), "Embeddings are not unique for different sequences."

    # Check for each embedding's uniqueness using numpy allclose with a tolerance
    for i in range(len(embeddings)):
        for j in range(i + 1, len(embeddings)):
            assert not np.allclose(embeddings[i], embeddings[j], atol=1e-6), \
                f"Embeddings for sequence {i} and {j} are too similar."


