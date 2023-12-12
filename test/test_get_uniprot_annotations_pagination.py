import os

import pandas as pd

import scripts.get_uniprot_annotations


def test_split_lineage():
    print(os.getcwd())
    assert True


def test_get_uniprot_id_list():
    data = [
        ["P00893"],
        ["P17597"],
        ["P9WG41"],
        ["P9WG39"],
        ["P07342"],
        ["Q04789"],
        ["Q04524"],
        ["P27696"],
    ]

    # Create the pandas DataFrame with column name is provided explicitly
    df = pd.DataFrame(data, columns=["UNIPROT"])

    print(df)
    print("hello")

    id_list = scripts.get_uniprot_annotations.get_uniprot_id_list(df)

    print(id_list)


# def test_split_lineage_when_no_info():
#     # Check that returning a cell with no lineage still works when splitting the lineage
#     # A0A4Y7BF74 is a deleted entry so has no lineage info
#
#     outpath = "./test.tsv"
#
#     if os.path.exists(outpath):
#         os.remove(outpath)
#
#     uniprot_list_ids = ["A0A873WER9", "A0A4Y7BF74"]
#
#     result_df = pd.read_csv(outpath, sep="\t", header=None)
#
#     assert len(result_df) == len(uniprot_list_ids)
#
#     os.remove(outpath)
