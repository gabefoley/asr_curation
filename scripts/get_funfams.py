import pandas as pd
import os
import requests
import time
import ast
import annot_functions as an

base_url = "https://www.cathdb.info"


# Function to add FunFam information to the DataFrame
def add_funfam_info_to_df(df, funfam_dict):
    for seq_id, info_dict in funfam_dict.items():
        if type(info_dict) != dict:
            info_dict = ast.literal_eval(info_dict)
        if info_dict:
            # If the entry is not None, add the information to the DataFrame
            df.loc[df["info"] == seq_id, "funfam_names"] = info_dict["funfam_names"]
            df.loc[df["info"] == seq_id, "funfam_descriptions"] = info_dict[
                "funfam_descriptions"
            ]
            df.loc[df["info"] == seq_id, "funfam_significance"] = info_dict[
                "funfam_significance"
            ]
        else:
            # If the entry is None, set the columns to None
            df.loc[df["info"] == seq_id, "funfam_names"] = None
            df.loc[df["info"] == seq_id, "funfam_descriptions"] = None
            df.loc[df["info"] == seq_id, "funfam_significance"] = None
    return df


def submit_fasta_sequence(fasta_sequence):
    url = f"{base_url}/search/by_funfhmmer"
    data = {"fasta": fasta_sequence}
    headers = {"Accept": "application/json"}
    response = requests.post(url, data=data, headers=headers)

    if response.status_code == 202:
        return response.json()["task_id"]
    else:
        print(f"Error submitting the FASTA sequence: {response.text}")
        return None


# Function to check the progress of the scan and wait for completion
def check_scan_progress(task_id):
    url = f"{base_url}/search/by_funfhmmer/check/{task_id}"
    headers = {"Accept": "application/json"}

    while True:
        response = requests.get(url, headers=headers)
        if response.status_code == 200:
            data = response.json()["data"]
            status = data.get("status")
            message = data.get("message")

            if status == "done":
                print("Scan completed successfully.")
                return True
            elif status == "running":
                print("Scan is still running. Waiting...")
                time.sleep(5)
            else:
                print(f"Scan status: {status}")
                print(f"Scan message: {message}")
                return False
        else:
            print(f"Error: {response.status_code} - {response.reason}")
            print(response.text)
            return False


# Function to retrieve scan results for a task ID
def retrieve_scan_results(task_id):
    url = f"{base_url}/search/by_funfhmmer/results/{task_id}"
    headers = {"Accept": "application/json"}
    response = requests.get(url, headers=headers)

    if response.status_code == 200:
        results = response.json()
        return results
    else:
        print(f"Error retrieving scan results for Task ID {task_id}: {response.text}")
        return None


# Function to process a single sequence
def process_sequence(seq_id, sequence):
    # Submit the sequence for scanning
    task_id = submit_fasta_sequence(sequence)
    # print(task_id)

    if task_id:
        # Check progress and retrieve results for the current task
        success = check_scan_progress(task_id)
        if success:
            results = retrieve_scan_results(task_id)
            if results:
                # Process the results and return them as a dictionary
                return extract_funfam_info(seq_id, results)

    return None  # Return None if processing failed for the sequence


# Function to extract FunFam information from scan results
def extract_funfam_info(seq_id, results):
    if len(results) > 0:
        funfam_names_list = []
        funfam_descriptions_list = []
        funfam_significance_list = []

        funfam_scan_results = results.get("funfam_scan", {}).get("results", [])

        for result in funfam_scan_results:
            hits = result.get("hits", [])

            for hit in hits:
                match_name = hit.get("match_name", "")
                match_description = hit.get("match_description", "")
                significance = hit.get("significance", "")

                funfam_names_list.append(match_name)
                funfam_descriptions_list.append(match_description)
                funfam_significance_list.append(str(significance))

        funfam_info = {
            "funfam_names": ";".join(funfam_names_list),
            "funfam_descriptions": ";".join(funfam_descriptions_list),
            "funfam_significance": ";".join(funfam_significance_list),
        }

        return funfam_info


# Function to get FunFam information for a DataFrame
def get_funfams(df, output_dir):
    existing_mapping = output_dir + "/funfam_mappings.txt"
    name_mapping = an.read_to_dict(existing_mapping)

    # Process sequences in batches
    results = []
    for index, row in df.iterrows():
        seq_id = row["info"]
        sequence = row["sequence"]

        # print (seq_id)
        # print (name_mapping)

        if seq_id not in name_mapping:
            # print ('not there')
            # print (seq_id)

            seq_funfam_info = process_sequence(seq_id, sequence)
            if seq_funfam_info:
                name_mapping[seq_id] = seq_funfam_info

        an.write_from_dict(existing_mapping, name_mapping)

    # print (name_mapping)

    df = add_funfam_info_to_df(df, name_mapping)

    return df


# Define classification functions
def classify_funfams_specific(row, name_mapping, name_mapping_specific):
    longer_names_specific = list(name_mapping_specific.keys())
    longer_names_general = list(name_mapping.keys())

    description = row['funfam_descriptions']

    if not pd.isnull(description):
        specific_matches = [name_mapping_specific.get(str(longer_name), None) for longer_name in longer_names_specific
                            if longer_name in description]
        if specific_matches:
            return ';'.join(specific_matches)
        else:
            general_matches = [name_mapping.get(str(longer_name), None) for longer_name in longer_names_general if
                               longer_name in description]
            return ';'.join(general_matches) if general_matches else None
    else:
        return None


def classify_funfams_general(row, name_mapping):
    longer_names_general = list(name_mapping.keys())
    description = row['funfam_descriptions']

    # for longer_name in longer_names_general:
    # 	print (longer_name)
    # 	print (type(longer_name))
    # print('description')

    # print(description)

    # print('longer_names_general')
    # print(longer_names_general)

    if not pd.isnull(description):

        general_matches = [name_mapping.get(str(longer_name), None) for longer_name in longer_names_general if
                           longer_name in description]

        return ';'.join(general_matches) if general_matches else None
    else:
        return None


def classify_funfams_specific_top(row, name_mapping, name_mapping_specific):
    longer_names_specific = list(name_mapping_specific.keys())
    longer_names_general = list(name_mapping.keys())

    description = row['funfam_descriptions']

    if not pd.isnull(description):
        first_match = None

        for longer_name in longer_names_specific:
            if longer_name in description:
                first_match = name_mapping_specific.get(longer_name)
                break

        if not first_match:
            for longer_name in longer_names_general:
                if longer_name in description:
                    first_match = name_mapping.get(longer_name)
                    break

        return first_match

    else:
        return None


def classify_funfams_general_top(row, name_mapping):
    longer_names_general = list(name_mapping.keys())
    description = row['funfam_descriptions']

    if not pd.isnull(description):
        first_match = None

        for longer_name in longer_names_general:
            if longer_name in description:
                first_match = name_mapping.get(longer_name)
                break

        return first_match

    else:
        return None