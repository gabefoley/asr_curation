import asr_curation.scripts.annot_functions as an
import pandas as pd

def test_get_amino_acids():

    pos = an.get_amino_acids('PPGP', 0)
    assert pos == 'P'


def test_get_amino_acids_multiple():

    pos = an.get_amino_acids('PPGP', 0, 2)
    assert pos == 'PG'

