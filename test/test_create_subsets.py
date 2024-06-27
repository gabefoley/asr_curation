import os

import pandas as pd


# Tests to check the validation of IDs
def test_create_subset_with_blank():
    annot_df = pd.read_csv(
        "test/files/create_subsets/test_subset_with_missing_values.csv"
    )

    outpath = "test/files/create_subsets/test_subset.csv"

    if os.path.exists(outpath):
        os.remove(outpath)

    os.system(
        f"python ./scripts/validate_ids.py --input_file ./test/files/test_original.csv --output_file {outpath}"
    )


def test_create_subset_with_broken_rule():
    annot_df = pd.read_csv(
        "test/files/create_subsets/test_subset_with_missing_values.csv"
    )

    outpath = "test/files/create_subsets/test_subset.csv"

    if os.path.exists(outpath):
        os.remove(outpath)

    os.system(
        f"python ./scripts/validate_ids.py --input_file ./test/files/test_original.csv --output_file {outpath}"
    )
