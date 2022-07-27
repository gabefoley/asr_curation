import io
import sequence
from sequence import Alphabet
import numpy as np
import pandas as pd
import seqcurate as sc
import logger
import os
from uniprot_cols import full_uniprot_cols
import itertools
import requests
import concurrent

from requests.adapters import HTTPAdapter
from requests.packages.urllib3.util.retry import Retry


def requests_retry_session(
    retries=3,
    backoff_factor=0.3,
    status_forcelist=(500, 502, 504),
    session=None,
):
    session = session or requests.Session()
    retry = Retry(
        total=retries,
        read=retries,
        connect=retries,
        backoff_factor=backoff_factor,
        status_forcelist=status_forcelist,
    )
    adapter = HTTPAdapter(max_retries=retry)
    session.mount('http://', adapter)
    session.mount('https://', adapter)
    return session



def getUniProtDict(ids, chunk_id, total_chunks, cols=""):


    # Format the lists of IDs and columns correctly
    updated_ids = ' or '.join(ids)
    cols = ",".join(cols)
    

    url = 'https://www.uniprot.org/uniprot/'
    
#     if chunk_id == 3 and total_chunks > 2:
#         url = 'http://djsdisjds.com'


    params = {
        'format': 'tab',
        'query': updated_ids,
        'columns': "id," + cols
    }
    

    session = requests.Session()

    print (f'Working on {chunk_id} of {total_chunks}')



    try:
        print ('Querying UniProt')
        resp = requests_retry_session().post(url, data=params)
        return (0, resp)

    except requests.exceptions.ConnectionError as e:
        print ('Connection error')
        return (1, ids)
    
    except requests.exceptions.HTTPError as e:
        print ('404 error')
        return (1, ids)
    
    except Exception as e:
        print (f"ERROR - {str(e)} ")
        return (1, ids)

def setup_query(chunks):
    success = 1
    total_chunks = len(chunks)

    with concurrent.futures.ThreadPoolExecutor() as executor:
        # while count < max_retry:
        print ('chunks has length ')
        print (len(chunks))

        res = [executor.submit(getUniProtDict, chunk, chunk_id + 1, len(chunks), full_uniprot_cols) for chunk_id, chunk in enumerate(chunks)]

        try:
            for future in concurrent.futures.as_completed(res, timeout=3):

                print (future)
                print (future.result())


                if future.result()[0] == 1:
                    total_chunks +=1
                    print (f'An ID mapping failed. We are on try {count} of {max_retry}')
                    failed.append(future.result()[1])

                elif future.result()[0] == 0:
                    print (f'Completed {success} of {total_chunks}')
                    success +=1
                    futures.append(future.result()[1])


        except concurrent.futures.TimeoutError as e:
            pass

    return (futures, failed)
        



# Allow all letters to be in the sequence
all_letters = Alphabet('ABCDEFGHIJKLMNOPQRSTUVWXYZ')
original_df = pd.read_csv(snakemake.input[0])


names = [name for name in original_df['Entry'].tolist()]

half = round(len(names) / 2)
print(half)
names = names[0:half]

# Split the query to UniProt up into chunks of at most 1000 sequences. Normally we can query a lot more sequences but because we're asking for all the columns we restrict the number of sequences in each request

chunks = np.array_split(np.array(names), max(1, len(names) / 1000))

frames = []


futures = []
failed = []

count = 1
max_retry = 4

failed_dir = snakemake.input[0].split(snakemake.wildcards.dataset)[0] + snakemake.wildcards.dataset + "/failed/"

if not os.path.isdir(failed_dir):
    os.mkdir(failed_dir)


print ("Retrieving annotations from UniProt. This could take some time.")

complete = True


completed, failed = setup_query(chunks)


print (completed)

# finished_futures = []
# with concurrent.futures.ThreadPoolExecutor() as executor:

#     # Schedule the first N futures.  We don't want to schedule them all
#     # at once, to avoid consuming excessive amounts of memory.
#     futures = [executor.submit(getUniProtDict, chunk, chunk_id + 1, len(chunks), full_uniprot_cols) for chunk_id, chunk in enumerate(chunks)]
    
#     result = False
#     while not result:
#         maybe_futures = concurrent.futures.wait(futures, return_when=concurrent.futures.FIRST_COMPLETED)
#         futures = maybe_futures.not_done

#         for future in maybe_futures.done:
#             print ('Found a maybe')
#             temp = future.result()
#             if not temp:
#                 executor.submit(getUniProtDict, chunk, chunk_id + 1, len(chunks), full_uniprot_cols)
#             else:
#                 result = temp
#                 finished_futures.append(temp)
#                 break



    # while futures:
    #     # Wait for the next future to complete.
    #     done, futures = concurrent.futures.wait(
    #         futures, return_when=concurrent.futures.FIRST_COMPLETED
    #     )

    #     for fut in done:
    #         print(f"The outcome is {fut.result()}")

    #     # Schedule the next set of futures.  We don't want more than N futures
    #     # in the pool at a time, to keep memory consumption down.
    #     for task in itertools.islice(chunks, len(done)):
    #         futures.add(executor.submit(getUniProtDict, chunk, chunk_id + 1, len(chunks), full_uniprot_cols))

# print (len(finished_futures))


# with concurrent.futures.ThreadPoolExecutor() as executor:

#     # Schedule the first N futures.  We don't want to schedule them all
#     # at once, to avoid consuming excessive amounts of memory.
#     futures = {executor.submit(getUniProtDict, chunk, chunk_id + 1, len(chunks), full_uniprot_cols) for chunk_id, chunk in enumerate(chunks)}
    
#     while futures:
#         # Wait for the next future to complete.
#         done, futures = concurrent.futures.wait(
#             futures, return_when=concurrent.futures.FIRST_COMPLETED
#         )

#         for fut in done:
#             print(f"The outcome is {fut.result()}")

#         # Schedule the next set of futures.  We don't want more than N futures
#         # in the pool at a time, to keep memory consumption down.
#         for task in itertools.islice(chunks, len(done)):
#             futures.add(executor.submit(getUniProtDict, chunk, chunk_id + 1, len(chunks), full_uniprot_cols))


# if failed:
    
#     print ('Some stuff failed')

#     failed_names = [item for sublist in failed for item in sublist]

#     print (failed_namesl)

#     print (snakemake.wildcards.dataset)
#     print ()

#     sc.write_to_fasta(original_df[original_df.name.isin(failed_names)])

# df[df.country.isin(countries_to_keep)]
# while not complete or count < max_retry:
#     if failed:
#         print (f"{len(failed)} of the queries failed.")
#         count +=1

#         if count < max_retry:
#             print (f"On try {count} of {max_retry}")
#             setup_query(failed)






if len(futures) == len(chunks):
    print ('All annotations retrieved')
    
else:
    print ('#WARNING# some annotations could not be retrieved')
    print (chunks)

for index in range(len(futures)):
    # print (res[index].result())
    response = pd.read_csv(io.StringIO(futures[index].text),  sep='\t')
    frames.append(response)


# Concatenate all the retrieved UniProt annotations
full_df = pd.concat(frames, axis=0)

print ('here is the full df')

print (full_df)



full_df.drop_duplicates(subset='Entry', keep='first')

full_df['Entry'].tolist()

full_df.to_csv(snakemake.output[0] + "_full")


print ('here is original_df')

print (original_df)

print ('orig dype')

print(original_df.dtypes)

print ('full dtype')
print(full_df.dtypes)

merged_df = pd.merge(original_df, full_df, how = 'left', on=['Entry'])

# Remove brackets from header names
merged_df.columns = merged_df.columns.str.replace("[()]", "")

# Replace spaces in header names with underscores
merged_df.columns = merged_df.columns.str.replace(" ", "_")



# merged_df.to_csv(snakemake.output[0] + "_merged")


# Add the new UniProt annotations to the original annotation file
# merged_df = pd.merge(original_df, full_df, how = 'left', on = ['Entry']).set_index(original_df.index)

print ('here is the merged df')

print (merged_df)

# Save the merged annotations to a csv
merged_df.to_csv(snakemake.output[0], index=False)
