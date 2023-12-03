import pandas as pd
import pytest
from scripts.create_itol_files import get_color_dict_and_info_list

@pytest.fixture
def sample_dataframe():
    data = {
        'info': ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H'],
        'column_to_color': ['X', 'Y', 'Z', 'X', 'Y', 'X', 'Z', 'Y']
    }
    return pd.DataFrame(data)

def test_get_color_dict_and_info_list(sample_dataframe):
    col = 'column_to_color'
    color_dict, info_list = get_color_dict_and_info_list(sample_dataframe, col)

    # Check if the color_dict contains the expected keys

    assert list(color_dict.keys()) == ['X', 'Y', 'Z']

    # Check if the info_list contains the expected pairs of (info, color)
    expected_info_list = [
        ('A', 'X'), ('B', 'Y'), ('C', 'Z'),
        ('D', 'X'), ('E', 'Y'), ('F', 'X'),
        ('G', 'Z'), ('H', 'Y')
    ]
    assert info_list == expected_info_list
