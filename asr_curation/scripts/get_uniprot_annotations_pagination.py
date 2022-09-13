''' Python code to fetch annotations from uniprot in batches, use pagination to process results and creating .csv output file '''

import pandas as pd
import urllib.parse
import urllib.request
from io import StringIO
import requests,json
from time import sleep
import re
import time
import json
import zlib
from xml.etree import ElementTree
from urllib.parse import urlparse, parse_qs, urlencode
import requests
from requests.adapters import HTTPAdapter, Retry
from configs.uniprot_cols import full_uniprot_cols
import numpy as np

# config parameters
API_URL = "https://rest.uniprot.org/uniprotkb/search?"
retries = Retry(total=5, backoff_factor=0.25, status_forcelist=[500, 502, 503, 504])
session = requests.Session()
session.mount("https://", HTTPAdapter(max_retries=retries))
re_next_link = re.compile(r'<(.+)>; rel="next"')

def get_next_link(headers):
    ''' get next batch link '''

    if "Link" in headers:
        match = re_next_link.match(headers["Link"])
        if match:
            return match.group(1)

def get_batch(batch_url):
    ''' get next batch '''

    while batch_url:
        response = session.get(batch_url)
        response.raise_for_status()
        total = response.headers["x-total-results"]
        yield response, total
        batch_url = get_next_link(response.headers)

def split_lineage(x):
    ''' function to split the lineage strings into seperate columns for use later '''

    lineage_order = ['lineage_superkingdom','lineage_kingdom',\
                     'lineage_subkingdom','lineage_superphylum','lineage_phylum',\
                     'lineage_subphylum','lineage_superclass','lineage_class',\
                     'lineage_subclass','lineage_infraclass','lineage_superorder',\
                     'lineage_order','lineage_suborder','lineage_infraorder',\
                     'lineage_parvorder','lineage_superfamily','lineage_family',\
                     'lineage_subfamily','lineage_tribe','lineage_subtribe',\
                     'lineage_genus','lineage_subgenus','lineage_species_group',\
                     'lineage_species_group','lineage_species_group','lineage_varietas',\
                     'lineage_forma']

    # split lineage
    lin_split = x.split(',')
    lineage_dict = {}
    lineage_list = []

    for l in lin_split:
        if 'no rank' not in l:
            val = l.split('(')[0].strip()
            key = 'lineage_' + l.split('(')[1].split(')')[0]
            lineage_dict[key] = val

    #print(lineage_dict)
    # return values in the order
    for l_o in lineage_order:
        try:
            lineage_list.append(lineage_dict[l_o])
        except:
            lineage_list.append('')

    return lineage_list


def process_results(intermediate_tsv_file):
    ''' process results in tsv file '''

    annot_df = pd.read_csv(intermediate_tsv_file,sep='\t')

    # replace gaps in column names with '_'
    annot_df.columns = annot_df.columns.str.replace(" ", "_")
    annot_df.columns = annot_df.columns.str.replace("-", "_")
    annot_df.columns = annot_df.columns.str.replace("[()]", "",regex=False)


    # seperate lineage columns ( after the new api changes)
    annot_df['lineage_superkingdom'],annot_df['lineage_kingdom'], \
    annot_df['lineage_subkingdom'],annot_df['lineage_superphylum'],annot_df['lineage_phylum'],\
    annot_df['lineage_subphylum'],annot_df['lineage_superclass'],annot_df['lineage_class'],\
    annot_df['lineage_subclass'],annot_df['lineage_infraclass'],annot_df['lineage_superorder'],\
    annot_df['lineage_order'],annot_df['lineage_suborder'],annot_df['lineage_infraorder'],\
    annot_df['lineage_parvorder'],annot_df['lineage_superfamily'],annot_df['lineage_family'],\
    annot_df['lineage_subfamily'],annot_df['lineage_tribe'],annot_df['lineage_subtribe'],\
    annot_df['lineage_genus'],annot_df['lineage_subgenus'],annot_df['lineage_species_group'],\
    annot_df['lineage_species_group'],annot_df['lineage_species_group'],annot_df['lineage_varietas'],\
    annot_df['lineage_forma'] = zip(*annot_df['lineage'].map(split_lineage))

    return annot_df

def get_uniprot_annotations(uniprot_list_ids,intermediate_tsv_file,id_batch_size,result_batch_size):
    ''' running uniprot ids batches '''

    splits = np.array_split(np.array(uniprot_list_ids), max(1, round(len(uniprot_list_ids) / id_batch_size)))
    split_count = 0

    # For each split go to UniProt to retrieve information
    for split in splits:
        split_count += 1
        print(f"Running stage {split_count} of {len(splits)}")
        split_fetch_success = fetch_annotations_in_batches(split,intermediate_tsv_file,result_batch_size)

        if split_fetch_success == 1:
            print(f"Successfully Processed {split_count} of {len(splits)}")
        else:
            print(f"Failed Processed {split_count} of {len(splits)}")
            return 0

    return 1


def fetch_annotations_in_batches(ids,intermediate_tsv_file,result_batch_size):
    ''' process uniprot api results in batches and load it in .tsv file on the disk'''

    cols = ",".join(full_uniprot_cols)
    batch_url = API_URL + "format=tsv&" + "query=" +  ' OR '.join(ids) + "&fields=" + \
                                                     cols + "&size=" + str(result_batch_size)
    progress = 0

    with open(intermediate_tsv_file, 'a') as f:
        for batch, total in get_batch(batch_url):
            lines = batch.text.splitlines()
            for line in lines[1:]:
                print(line, file=f)
            progress += len(lines[1:])
            print(f'{progress} / {total}')

    if int(progress) ==  int(total): # all have be fetched
        print("Created and loaded annotations in .tsv file")
        return 1
    else:
        return 0

def get_uniprot_id_list(df):
    ''' function to get list of uniprot ids to make an api call to uniprot '''

    uniprot_df_list = df['UNIPROT'].tolist()
    uniprot_list_ids = []
    for l in uniprot_df_list:
        l_item = (str(l[1:-1])).replace(' ','')
        for m in l_item.split(','):
            uniprot_list_ids.append(m[1:-1])
    return uniprot_list_ids


def uniprot_annotation_lkp(input_file,output_file,intermediate_tsv_file,id_batch_size,result_batch_size):
    ''' main function to get uniprot ids annotation '''

    # get all uniprot ids
    validated_df = pd.read_csv(input_file)
    uniprot_list_ids = get_uniprot_id_list(validated_df)
    #print(uniprot_list_ids)

    # create a new intermediate tsv file with header row for dumping results
    header = "\t".join(full_uniprot_cols)
    with open(intermediate_tsv_file, 'w') as f:
        print(header, file=f)

    # submit all ids to uniprot for gettting annotation and load in .tsv file
    results_uniprot_loaded = get_uniprot_annotations(uniprot_list_ids,\
                                                     intermediate_tsv_file,id_batch_size,result_batch_size)

    # process results
    if results_uniprot_loaded:
        results_df = process_results(intermediate_tsv_file)
    else:
        print("Uniprot Results fetched failed.")

    # save in .csv file
    print("Creating .csv output file")
    results_df.to_csv(output_file, index=False)
    print("Completed Uniprot Annotations Job")

def main():

    #input_file  = 'kari_example_ec_1_1_1_86_original_validated.csv'
    #output_file = 'kari_example_ec_1_1_1_86_uniprot.csv'
    #intermediate_tsv_file = 'kari_example_ec_1_1_1_86_temp.tsv'

    input_file  = snakemake.input[0]
    output_file = snakemake.output[0]
    intermediate_tsv_file = snakemake.input[0].split(".")[0] + '.tsv'

    id_batch_size = 100 # process ids in batches of this parameter
    result_batch_size = 50 # process results in this batches (uses pagination)

    print("Starting Uniprot Annotations Job")
    uniprot_annotation_lkp(input_file,output_file,intermediate_tsv_file,id_batch_size,result_batch_size)

if __name__ == "__main__":
    main()
