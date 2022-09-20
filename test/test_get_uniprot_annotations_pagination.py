import os

import pandas as pd

import scripts.get_uniprot_annotations_pagination

def test_split_lineage():

    print (os.getcwd())
    assert (True)

def test_get_uniprot_id_list():

    data = [['P00893'],['P17597'],['P9WG41'],['P9WG39'],['P07342'],['Q04789'],['Q04524'],['P27696']]

    # Create the pandas DataFrame with column name is provided explicitly
    df = pd.DataFrame(data, columns=['UNIPROT'])

    print (df)
    print ('hello')


    id_list = scripts.get_uniprot_annotations_pagination.get_uniprot_id_list(df)

    print (id_list)
