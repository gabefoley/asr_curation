import io
#import sequence
#from sequence import Alphabet
import numpy as np
import pandas as pd
import seqcurate as sc
import os
from configs.uniprot_cols import full_uniprot_cols
import aiohttp
import asyncio
import async_timeout
import time

async def get_uniprot(session, url, params):
    async with session.post(url, data=params) as resp:
        up_resp = await resp.text()
        # print ('got query')
        return up_resp


def format_uniprot_params(ids):
    updated_ids = " or ".join(ids)
    cols = ",".join(full_uniprot_cols)

    params = {"format": "tab", "query": updated_ids, "columns": "id," + cols}
    return params


async def main(split):

    my_timeout = aiohttp.ClientTimeout(
        total=30,  # default value is 5 minutes, set to `None` for unlimited timeout
        sock_connect=30,  # How long to wait before an open socket allowed to connect
        sock_read=30,  # How long to wait with no data being read before timing out
    )

    timeout = 60

    mysession = aiohttp.ClientSession()
    # async with aiohttp.ClientSession(timeout=my_timeout) as session:
    async with mysession as session:

        tasks = []

        chunks = np.array_split(np.array(split), max(1, round(len(split) / 800)))

        for chunk in chunks:
            task_count = 0
            # print (len(chunk))
            params = format_uniprot_params(chunk)
            url = "https://www.uniprot.org/uniprot/"

            task = asyncio.create_task(get_uniprot(session, url, params))
            task.up_ids = chunk
            tasks.append(task)

            # tasks.append(asyncio.ensure_future(get_uniprot(session, url, params)))

        print(
            f"Submitting {len(tasks)} total requests with around {len(chunk)} sequences each to UniProt at this stage"
        )
        try:
            with async_timeout.timeout(timeout):
                up_response = await asyncio.gather(*tasks)

            # for resp in up_response:
            #     count = len(frames) + 1
            #     print (f"Received annotations {count} of {len(chunks)} ")
            #     response = pd.read_csv(io.StringIO(resp),  sep='\t')
            #     frames.append(response)
        except asyncio.TimeoutError:
            pass

        finally:
            for i, task in enumerate(tasks):
                task_count += 1

                if task.done() and not task.cancelled():

                    print(f"Task {task_count} of {len(tasks)} is finished")
                    # print (task.result())
                    response = pd.read_csv(io.StringIO(task.result()), sep="\t")
                    # print (response)
                    frames.append(response)
                    continue
                else:
                    print(f"Task {task_count} failed")
                    failed.append(task.up_ids)

def split_data_for_uniprot(data, chunk_size):
    splits = np.array_split(np.array(names), max(1, round(len(names) / 8000)))
    split_count = 0

    # For each split go to UniProt to retrieve information
    for split in splits:
        split_count += 1
        print(f"Running stage {split_count} of {len(splits)}")
        asyncio.run(main(split))
        time.sleep(5)



# Read in the dataframe
original_df = pd.read_csv(snakemake.input[0])


# Get all the names
names = [name for name in original_df["Extracted_ID"].tolist()]

print(f"Total number of sequences to query UniProt is {len(names)}")

# Setup lists to store completed and failed sequences and start timing
frames = []
failed = []
start_time = time.time()


# Split entries up into chunks and call the main async method to retrieve data from UniProt
split_data_for_uniprot(names, 8000)

# If entries fail then we retry them once
if failed:
    failed_names = [item for sublist in failed for item in sublist]

    print(f"\nWARNING: {len(failed_names)} sequences failed. Let's retry them.")
    split_data_for_uniprot(failed_names, 8000)



# If sequences still fail write them to a file
if failed:
    failed_names = [item for sublist in failed for item in sublist]


    print(f"\nWARNING: {len(failed_names)} sequences failed even after a second attempt.")

    failed_names = [item for sublist in failed for item in sublist]

    # Filepath to a failed directory
    failed_dir = (
        snakemake.input[0].split(snakemake.wildcards.dataset)[0]
        + snakemake.wildcards.dataset
        + "/failed/"
    )

    # Get the full details for the failed seqs
    failed_seqs = original_df[original_df.Entry.isin(failed_names)]

    # Create the failed directory
    if not os.path.isdir(failed_dir):
        os.mkdir(failed_dir)

    # Write all the failed sequences out to the directory
    sc.write_to_fasta(
        failed_seqs,
        failed_dir + snakemake.wildcards.dataset + "_failed.fasta",
    )



# If we have at least some sequences in frames then lets process them
elif frames:
    full_df = pd.concat(frames, axis=0)


    full_df.drop_duplicates(subset="Entry", keep="first")


    merged_df = pd.merge(
        original_df, full_df, left_on=["Extracted_ID"], right_on = ["Entry"], suffixes=["", "_r"]
    )

    # Remove brackets from header names
    merged_df.columns = merged_df.columns.str.replace("[()]", "")

    # Replace spaces in header names with underscores
    merged_df.columns = merged_df.columns.str.replace(" ", "_")

    # Replace dashes in header names with underscores
    merged_df.columns = merged_df.columns.str.replace("-", "_")

    # Save the merged annotations to a csv
    merged_df.to_csv(snakemake.output[0], index=False)
