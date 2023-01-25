import scripts.annot_functions as an
import pandas as pd
import os
import pytest
import scripts.seqcurate as sc
import scripts.validate_ids as vi
import subprocess as sp
import shutil


# Tests to check the validation of IDs
def test_validate_ids_main():

    outpath = "./files/output_test_validated_main.csv"

    if os.path.exists(outpath):
        os.remove(outpath)

    vi.all_ids_lookup("./files/test_original.csv", outpath)

    assert (os.path.isfile(outpath))


def test_validate_ids_main_from_command_line():

    outpath = "./files/output_test_validated_command_line.csv"

    if os.path.exists(outpath):
        os.remove(outpath)

    os.system(f"python ../scripts/validate_ids.py --input_file ./files/test_original.csv --output_file {outpath}")

    assert (os.path.isfile(outpath))

def test_validate_ids_call_command_line_directly():

    outpath = "./files/output_test_validated_command_line.csv"

    if os.path.exists(outpath):
        os.remove(outpath)

    vi.all_ids_lookup_cmd(["--input_file", "./files/test_original.csv", "--output_file", outpath], standalone_mode=False)

    assert (os.path.isfile(outpath))


