""" Python script to validate uniprot ids and find the corresponding NCBI,EMBL,Uniparc ids """

import pandas as pd
import requests
import json
import zlib
import re
import time
import math
from xml.etree import ElementTree
from urllib.parse import urlparse, parse_qs, urlencode
from requests.adapters import HTTPAdapter, Retry
from itertools import islice

# Configuration and API parameters
POLLING_INTERVAL = 3
API_URL = "https://rest.uniprot.org"
retries = Retry(total=5, backoff_factor=0.25, status_forcelist=[500, 502, 503, 504])
session = requests.Session()
session.mount("https://", HTTPAdapter(max_retries=retries))

# Dictionary for mapping database names
from_id_db_lookup = {
    "NCBI": "RefSeq_Protein",
    "EMBL": "EMBL-GenBank-DDBJ",
    "EMBL-CDS": "EMBL-GenBank-DDBJ_CDS",
    "UNIPROT-FROM": "UniProtKB_AC-ID",
    "UNIPARC": "UniParc",
}

to_id_db_lookup = {
    "NCBI": "RefSeq_Protein",
    "EMBL": "EMBL-GenBank-DDBJ",
    "EMBL-CDS": "EMBL-GenBank-DDBJ_CDS",
    "UNIPROT": "UniProtKB",
    "UNIPARC": "UniParc",
    "UNIREF50": "UniRef50",
    "UNIREF90": "UniRef90",
}

to_id_db_lookup = {
    "UNIPROT": "UniProtKB",
    "UNIREF50": "UniRef50",
    "UNIREF90": "UniRef90",
}

# Response format for different databases
db_response_format = {
    "UNIPROT": "long",
    "NCBI": "short",
    "UNIPARC": "long",
    "EMBL": "short",
    "EMBL-CDS": "short",
    "UNIREF50": "long",
    "UNIREF90": "long",
}

# Response key to be used in JSON format
db_response_key = {
    "UNIPROT": "primaryAccession",
    "UNIPARC": "uniParcId",
    "UNIREF50": "id",
    "UNIREF90": "id",
}
# db_response_key = {"UNIPROT": "primaryAccession"}


# functions
def submit_id_mapping(from_db, to_db, ids):
    """function to make api call to uniprot"""

    request = requests.post(
        f"{API_URL}/idmapping/run",
        data={"from": from_db, "to": to_db, "ids": ",".join(ids)},
    )
    try:
        request.raise_for_status()
    except:
        return
    return request.json()["jobId"]


def check_id_mapping_results_ready(job_id, verbose):
    """function to check if the api call results are ready"""

    while True:
        request = session.get(f"{API_URL}/idmapping/status/{job_id}")
        request.raise_for_status()
        j = request.json()
        if "jobStatus" in j:
            if j["jobStatus"] == "RUNNING":
                if verbose:
                    print(f"Retrying in {POLLING_INTERVAL}s")
                time.sleep(POLLING_INTERVAL)
            else:
                raise Exception(request["jobStatus"])
        else:
            return bool(j["results"] or j["failedIds"])


def get_id_mapping_results_link(job_id):
    """function to get the mapping results link from uniprot after results are ready"""

    url = f"{API_URL}/idmapping/details/{job_id}"
    request = session.get(url)
    request.raise_for_status()
    return request.json()["redirectURL"]


def decode_results(response, file_format, compressed):
    """function to decode the results based on file format"""

    if compressed:
        decompressed = zlib.decompress(response.content, 16 + zlib.MAX_WBITS)
        if file_format == "json":
            j = json.loads(decompressed.decode("utf-8"))
            return j
        elif file_format == "tsv":
            return [line for line in decompressed.decode("utf-8").split("\n") if line]
        elif file_format == "xlsx":
            return [decompressed]
        elif file_format == "xml":
            return [decompressed.decode("utf-8")]
        else:
            return decompressed.decode("utf-8")
    elif file_format == "json":
        return response.json()
    elif file_format == "tsv":
        return [line for line in response.text.split("\n") if line]
    elif file_format == "xlsx":
        return [response.content]
    elif file_format == "xml":
        return [response.text]
    return response.text


def get_batch(batch_response, file_format, compressed):
    """function to get the batch if the results are too big"""

    batch_url = get_next_link(batch_response.headers)
    while batch_url:
        batch_response = session.get(batch_url)
        batch_response.raise_for_status()
        yield decode_results(batch_response, file_format, compressed)
        batch_url = get_next_link(batch_response.headers)


def get_next_link(headers):
    re_next_link = re.compile(r'<(.+)>; rel="next"')
    if "Link" in headers:
        match = re_next_link.match(headers["Link"])
        if match:
            return match.group(1)


def print_progress_batches(batch_index, size, total):
    n_fetched = min((batch_index + 1) * size, total)
    print(f"Fetched: {n_fetched} / {total}")


def get_xml_namespace(element):
    m = re.match(r"\{(.*)\}", element.tag)
    return m.groups()[0] if m else ""


def merge_xml_results(xml_results):
    merged_root = ElementTree.fromstring(xml_results[0])
    for result in xml_results[1:]:
        root = ElementTree.fromstring(result)
        for child in root.findall("{http://uniprot.org/uniprot}entry"):
            merged_root.insert(-1, child)
    ElementTree.register_namespace("", get_xml_namespace(merged_root[0]))
    return ElementTree.tostring(merged_root, encoding="utf-8", xml_declaration=True)


def get_id_mapping_results_search(url, verbose=True):
    """function to process the results in a batch by batch"""

    parsed = urlparse(url)
    query = parse_qs(parsed.query)
    file_format = query["format"][0] if "format" in query else "json"
    if "size" in query:
        size = int(query["size"][0])
    else:
        size = 500
        query["size"] = size
    compressed = (
        query["compressed"][0].lower() == "true" if "compressed" in query else False
    )
    parsed = parsed._replace(query=urlencode(query, doseq=True))
    url = parsed.geturl()
    request = session.get(url)
    request.raise_for_status()
    results = decode_results(request, file_format, compressed)

    if results["results"]:
        total = int(request.headers["x-total-results"])
        if verbose:
            print_progress_batches(0, size, total)
        for i, batch in enumerate(get_batch(request, file_format, compressed), 1):
            results = combine_batches(results, batch, file_format)
            if verbose:
                print_progress_batches(i, size, total)
        if file_format == "xml":
            return merge_xml_results(results)
    else:
        return  # no results returned. all ids failed
    return results


def combine_batches(all_results, batch_results, file_format):
    """function to combine results from many batches into 1 output json format"""

    if file_format == "json":
        for key in ("results", "failedIds"):
            if key in batch_results and batch_results[key]:
                all_results[key] += batch_results[key]
    elif file_format == "tsv":
        return all_results + batch_results[1:]
    else:
        return all_results + batch_results
    return all_results


def map_ids_in_uniprot(ids, from_db, to_db, verbose=True):
    """map ids to uniprot"""

    # submit job to uniprot
    job_id = submit_id_mapping(from_db, to_db, ids)
    # print('job_id',job_id)

    if job_id != None:  # if the call was made successfully
        # check if ready and get if ready.
        if check_id_mapping_results_ready(job_id, verbose):
            link = get_id_mapping_results_link(job_id)
            results = get_id_mapping_results_search(link, verbose=verbose)
    else:
        return  # this happens if the "from id" is not correct id, so the call was not made
    return results


def try_map_all_ids(ids, from_db_short, to_db_short, verbose=True):
    """function to try guess each database id belongs to and make a call to uniprot"""

    # get uniprot identifier for database
    from_db = from_id_db_lookup[from_db_short]
    to_db = to_id_db_lookup[to_db_short]

    results_this_db = map_ids_in_uniprot(ids, from_db, to_db, verbose=verbose)

    return results_this_db


def merge_results_to_main_data(main_df, add_df, from_db, to_db):
    """merge the result processed output with main master data"""

    # left join on all extracted ids to check if they match on the returned result
    main_df = main_df.merge(add_df, left_on="extracted_id", right_on="from", how="left")

    print("done")
    print(from_db)
    print(to_db)

    print(main_df.head())

    main_df[to_db] = main_df[to_db].str.cat(main_df["grouped_to"], sep=" ", na_rep="")

    # to be done only if it is not uniprot id
    # if from_db != "UNIPROT-FROM":
    #     main_df[from_db] = main_df[from_db].str.cat(main_df["from"], sep=" ", na_rep="")

    if from_db == "EMBL-CDS":
        main_df[from_db] = main_df["grouped_to"].str.cat(
            main_df["from"], sep=" ", na_rep=""
        )
    main_df.drop(columns=["from", "grouped_to"], inplace=True)

    # main_df = main_df.merge(add_df,left_on='Extracted_ID_2', right_on='from', how='left')
    # main_df[to_db]   = main_df[to_db].str.cat(main_df['grouped_to'], sep = " ",na_rep = '')
    # if from_db != 'UNIPROT-FROM':
    #     main_df[from_db] = main_df[from_db].str.cat(main_df['from'], sep = " ",na_rep = '')
    # main_df.drop(columns=['from','grouped_to'],inplace=True)

    # main_df = main_df.merge(add_df,left_on='Extracted_ID_3', right_on='from', how='left')
    # main_df[to_db]   = main_df[to_db].str.cat(main_df['grouped_to'], sep = " ",na_rep = '')
    # if from_db != 'UNIPROT-FROM':
    #     main_df[from_db] = main_df[from_db].str.cat(main_df['from'], sep = " ",na_rep = '')
    # main_df.drop(columns=['from','grouped_to'],inplace=True)

    return main_df


def process_results(data_df, results, from_db, to_db):
    """process the combined output results from uniprot to a grouped cleaner format"""

    # response can be short or long, depending on that need to process seperately
    if db_response_format[to_db] == "short":  # short response recieved
        result_df = pd.DataFrame(results["results"], columns=["from", "to"])
    else:
        id_mapping_list = []
        for d in results["results"]:
            from_id = d["from"]
            to_key = db_response_key[to_db]
            to_id = d["to"][to_key]

            print(to_key)
            print(to_id)
            id_mapping_list.append((from_id, to_id))
        result_df = pd.DataFrame(id_mapping_list, columns=["from", "to"])

    # group the results by from and remove duplicates
    result_df["grouped_to"] = result_df.groupby(["from"])["to"].transform(
        lambda x: " ".join(x)
    )
    result_df = (
        result_df[["from", "grouped_to"]].drop_duplicates().reset_index(drop=True)
    )

    print(result_df)

    # merge the results with the main dataframe
    merged_df = merge_results_to_main_data(data_df, result_df, from_db, to_db)

    return merged_df


def create_output_file(data_df, to_id_lookup, output_file):
    """function is create a clean output file with sequence as the primary key"""

    print("got here")
    print(data_df)

    print(to_id_lookup)

    # combine entry column,id columns by sequence
    data_df["accession_all"] = data_df.groupby(["info"])["info"].transform(
        lambda x: " ".join(x)
    )
    for db in to_id_lookup:
        data_df[db] = data_df.groupby(["info"])[db].transform(
            lambda x: " ".join(x.str.strip())
        )

    # drop duplicates by sequence
    cols_output_file = ["info", "accession_all"] + [db for db in to_id_lookup]
    data_df = data_df.drop_duplicates().reset_index(drop=True)

    # remove duplicates in id each column themselves
    for db in to_id_lookup:
        data_df[db] = data_df[db].apply(lambda x: list(set(x.strip().split(" "))))

    # remove duplicates in the entry column
    data_df["accession_all"] = data_df["accession_all"].apply(
        lambda x: list(set(x.strip().split(" ")))
    )

    # save as csv and return
    data_df.to_csv(output_file, index=False)


def chunk_it(iterable, size):
    it = iter(iterable)
    while True:
        chunk = tuple(islice(it, size))
        if not chunk:
            return
        yield chunk


CHUNK_SIZE = 5000  # Adjust based on your requirements


def all_ids_lookup(
    input_file, output_file, from_id_lookup=None, to_id_lookup=None, verbose=True
):
    """Main function to map input ids to different database specified in the id_lookup list"""

    # Set default from and to id lookups
    if not from_id_lookup:
        from_id_lookup = ["UNIPROT-FROM"]
        # from_id_lookup = ["EMBL-CDS"]

    if not to_id_lookup:
        to_id_lookup = ["UNIPROT"]

    # read data and get ids
    df_data = pd.read_csv(input_file)
    ids = list(df_data["extracted_id"])

    print(ids)

    if verbose:
        print("Validating IDs")
        print(f"Mapping from - {from_id_lookup}")
        print(f"Mapping to - {to_id_lookup}")
        print(f"There are a total of {len(ids)} IDs to map")

    # Create an empty DataFrame to hold concatenated results
    df_final = pd.DataFrame()

    # empty string for each of the lookup id columns
    for db in to_id_lookup:
        df_data[db] = ""

    total_chunks = math.ceil(len(ids) / CHUNK_SIZE)

    # process mapping from uniprot to "to ids list" using chunked IDs
    for chunk_num, chunked_ids in enumerate(chunk_it(ids, CHUNK_SIZE), start=1):
        if verbose:
            print(f"Processing chunk {chunk_num} of {total_chunks}")

        for from_db_short in from_id_lookup:
            if verbose:
                print("Finding the ids from database:", from_db_short)

            for to_db_short in to_id_lookup:
                if from_db_short != to_db_short:
                    if verbose:
                        print("Finding the ids to database:", to_db_short)

                    # call uniprot API with chunked IDs
                    results_uniprot = try_map_all_ids(
                        chunked_ids, from_db_short, to_db_short, verbose
                    )

                    # check if any results returned
                    if results_uniprot:
                        # process the results
                        if verbose:
                            print("OUTCOME :: Results Returned")
                        df_data = process_results(
                            df_data, results_uniprot, from_db_short, to_db_short
                        )
                        if verbose:
                            print("Results processed\n")
                    else:
                        if verbose:
                            print("OUTCOME :: No Results Returned\n")

        # After processing the results for a chunk, concatenate it to the final DataFrame
        df_final = pd.concat(
            [df_final, df_data[df_data["extracted_id"].isin(chunked_ids)]],
            ignore_index=True,
        )

    # create final clean output file - sequence as primary key
    create_output_file(df_final, to_id_lookup, output_file)

    return df_final


if __name__ == "__main__":
    input_file = snakemake.input[0]
    output_file = snakemake.output[0]
    verbose = snakemake.params.verbose
    # all_ids_lookup(input_file, output_file, to_id_lookup=to_id_db_lookup, verbose=verbose)
    all_ids_lookup(
        input_file, output_file, to_id_lookup=to_id_db_lookup, verbose=verbose
    )
