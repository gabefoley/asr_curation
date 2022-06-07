import io
import sequence
from sequence import Alphabet
import numpy as np
import pandas as pd
import seqcurate as sc
import os
from uniprot_cols import full_uniprot_cols
import aiohttp
import asyncio
import async_timeout
import time




# Allow all letters to be in the sequence
all_letters = Alphabet('ABCDEFGHIJKLMNOPQRSTUVWXYZ')
original_df = pd.read_csv(snakemake.input[0])


names = [name for name in original_df['Entry'].tolist()]

print (f'Total number of sequences to query UniProt is {len(names)}')
# half = round(len(names) / 15)
# print(half)
# names = names[0:half]

frames = []

failed = []
# Split the query to UniProt up into chunks of at most 1000 sequences. Normally we can query a lot more sequences but because we're asking for all the columns we restrict the number of sequences in each request



start_time = time.time()


async def get_uniprot(session, url, params):
    async with session.post(url, data=params) as resp:
        up_resp = await resp.text()
        # print ('got query')
        return up_resp

def format_uniprot_params(ids):
    # print ('Formatting the parameters')
    updated_ids = ' or '.join(ids)
    cols =   ",".join(full_uniprot_cols)


    params = {
        'format': 'tab',
        'query': updated_ids,
        'columns': cols
    }
    return params

async def main(split):

    my_timeout = aiohttp.ClientTimeout(
        total=30, # default value is 5 minutes, set to `None` for unlimited timeout
        sock_connect=30, # How long to wait before an open socket allowed to connect
        sock_read=30 # How long to wait with no data being read before timing out
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
            url = 'https://www.uniprot.org/uniprot/'

            task = asyncio.create_task(get_uniprot(session, url, params))
            task.up_ids = chunk
            tasks.append(task)

            # tasks.append(asyncio.ensure_future(get_uniprot(session, url, params)))

        print (f"Submitting {len(tasks)} total requests with around {len(chunk)} sequences each to UniProt at this stage")
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
                task_count +=1

                if task.done() and not task.cancelled():

                    print (f"Task {task_count} of {len(tasks)} is finished")
                    response = pd.read_csv(io.StringIO(task.result()),  sep='\t')
                    frames.append(response)
                    continue
                else:
                    print(f'Task {task_count} failed')
                    failed.append(task.up_ids)


# print ('len of names is ')
# print (len(names))

# print ('len of splits is ')
# print (len(names) / 8000)

splits = np.array_split(np.array(names), max(1, round(len(names) / 8000)))
split_count = 0
for split in splits:
    split_count +=1
    print (f"Running stage {split_count} of {len(splits)}")
    asyncio.run(main(split))
    time.sleep(5)

# print("--- %s seconds ---" % (time.time() - start_time))






if failed:
    failed_names = [item for sublist in failed for item in sublist]

    print (f"\nWARNING: {len(failed_names)} sequences failed. Let's retry them.")
    splits = np.array_split(np.array(failed_names), max(1, len(names) / 8000))
    failed = []
    split_count = 0
    for split in splits:
        split_count +=1
        print (f"Retrying failed sequences stage {split_count} of {len(splits)}")
        asyncio.run(main(split))
        time.sleep(2)




if failed:
    
    print ('Some stuff still failed')

    failed_names = [item for sublist in failed for item in sublist]


    failed_dir = snakemake.input[0].split(snakemake.wildcards.dataset)[0] + snakemake.wildcards.dataset + "/failed/"

    if not os.path.isdir(failed_dir):
        os.mkdir(failed_dir)

    sc.write_to_fasta(original_df[original_df.Entry.isin(failed_names)], failed_dir + snakemake.wildcards.dataset + "_failed.fasta")



if frames:
    full_df = pd.concat(frames, axis=0)

    #print ('here is the full df')

    #print (full_df)



    full_df.drop_duplicates(subset='Entry', keep='first')

    full_df['Entry'].tolist()

    # full_df.to_csv(snakemake.output[0] + "_full")


    #print ('here is original_df')

    #print (original_df)

    #print ('orig dype')

    #print(original_df.dtypes)

    #print ('full dtype')
    #print(full_df.dtypes)

    merged_df = pd.merge(original_df, full_df, how = 'left', on=['Entry'], suffixes=['', '_r'])

    # Remove brackets from header names
    merged_df.columns = merged_df.columns.str.replace("[()]", "")

    # Replace spaces in header names with underscores
    merged_df.columns = merged_df.columns.str.replace(" ", "_")

    merged_df.columns = merged_df.columns.str.replace("-", "_")


    # merged_df.to_csv(snakemake.output[0] + "_merged")


    # Add the new UniProt annotations to the original annotation file
    # merged_df = pd.merge(original_df, full_df, how = 'left', on = ['Entry']).set_index(original_df.index)

    #print ('here is the merged df')

    #print (merged_df)

    # Save the merged annotations to a csv
    merged_df.to_csv(snakemake.output[0], index=False)
