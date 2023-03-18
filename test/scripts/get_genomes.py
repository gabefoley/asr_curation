from collections import defaultdict
import pandas as pd
import numpy as np
import aiohttp
import asyncio
import async_timeout
import urllib.parse
import urllib.request
import time
import io
from Bio import SeqIO, Entrez
import sys


from itertools import islice

Entrez.email = "g.foley@uq.edu.au"


def add_to_df_from_dict(df, dict_to_add, col_name, add_name="Entry"):

    if col_name not in df.columns:
        df[col_name] = np.nan

    df[col_name] = df.apply(
        lambda row: dict_to_add[row[add_name]]
        if row[add_name] in dict_to_add
        else row[col_name],
        axis=1,
    )

    # for k, v in dict_to_add:
    #     df[k]
    # add_df = pd.DataFrame(dict_to_add.items(), columns=[add_name, col_name])

    # print ('merging')
    # print (df)

    # merged_df = df.join(add_df, on=add_name)

    # print (add_df)
    # merged_df = pd.merge(
    #         df, add_df, on=add_name, how='left'
    #     )

    # print (merged_df)
    return df


def chunks(data, SIZE=2):
    it = iter(data)
    for i in range(0, len(data), SIZE):
        yield {k: data[k] for k in islice(it, SIZE)}


async def get_uniprot(session, url, params):
    async with session.post(url, data=params) as resp:
        up_resp = await resp.text()
        # print ('got query')
        return up_resp


def format_uniprot_params(ids):
    # print ('Formatting the parameters')
    updated_ids = " or ".join(ids)

    params = {"format": "xml", "query": updated_ids}
    return params


async def main(split):

    my_timeout = aiohttp.ClientTimeout(
        total=30,  # default value is 5 minutes, set to `None` for unlimited timeout
        sock_connect=30,  # How long to wait before an open socket allowed to connect
        sock_read=30,  # How long to wait with no data being read before timing out
    )

    timeout = 60

    failed = []

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

                    name = ""
                    prot_id = ""

                    print(f"Task {task_count} of {len(tasks)} is finished")

                    for line in task.result().split("\n"):
                        if line.startswith("<name>"):
                            name = line.split("<name>")[1].split("</name")[0]
                        if line.startswith('<property type="protein sequence ID" '):

                            prot_id = line.split('value="')[1].split('"/>')[0]
                            # print ('prot id is ' + prot_id)
                        if len(name) > 1 and len(prot_id) > 1:
                            n2protid[name] = prot_id
                            name = ""
                            prot_id = ""

                    # response = pd.read_csv(io.StringIO(task.result()), sep="\t")
                    # with open('new_file.txt', 'w+') as new_file:
                    #     new_file.write(task.result())
                    # frames.append(response)
                    continue
                else:
                    print(f"Task {task_count} failed")
                    failed.append(task.up_ids)


def get_assembly_dict(id_dict):
    count = 0

    assembly_dict = defaultdict(list)

    for entry_name, gene_ids in id_dict.items():
        if pd.isnull(gene_ids):
            nulls.append(entry_name)

        else:

            for accession in gene_ids.split(";"):
                if len(accession) > 1:
                    added = False
                    handle = Entrez.efetch(db="ipg", id=[accession], rettype="ipg")
                    mapping = Entrez.read(handle)

                    if mapping:
                        print(mapping)
                        print(mapping["IPGReport"])

                        if "ProteinList" not in mapping["IPGReport"].keys():
                            pass

                        else:

                            for prot in mapping["IPGReport"]["ProteinList"]:

                                print("prot")

                                print(prot)

                                for idx, elem in enumerate(prot):

                                    if "CDSList" in prot and len(prot["CDSList"]) > 0:

                                        print("wowzers")
                                        print(prot["CDSList"])

                                        for x in prot["CDSList"]:
                                            print(x)
                                            if "assembly" in x.attributes:
                                                a_id = x.attributes["assembly"]
                                                # If we find a Genbank ID, check if we already have the RefSeq Genome in our list
                                                if "GCA" in a_id:
                                                    if (
                                                        a_id.replace("GCA_", "GCF_")
                                                        not in assembly_dict[entry_name]
                                                    ):
                                                        added = True
                                                        assembly_dict[
                                                            entry_name
                                                        ].append(
                                                            x.attributes["assembly"]
                                                        )
                                                else:
                                                    added = True
                                                    assembly_dict[entry_name].append(
                                                        x.attributes["assembly"]
                                                    )

                            if not added:
                                print("Nothing to add")
                                assembly_dict[entry_name].append("unmappable")

                        count += 1
                        print(count)

    return assembly_dict


original_df = pd.read_csv(snakemake.input[0])
id_dict = {}
nulls = []

# Get the IDs to go to
for row in original_df.itertuples():
    refseq_value = row.Cross_reference_RefSeq
    bacterial_transcript_value = row.EnsemblBacteria_transcript

    # Use the RefSeq value
    if not pd.isnull(refseq_value):
        id_dict[row.Entry] = refseq_value

    # Otherwise use the EnsemblBacteria transcript value
    elif not pd.isnull(bacterial_transcript_value):
        id_dict[row.Entry] = bacterial_transcript_value

    # Otherwise we'll have to grab them
    else:
        nulls.append(row.Entry)

# We can't get the EMBL protein ID directly from UniProt, so we need to manually
# map them if we get a null value

n2protid = {}


print("nulls is ")
print(nulls)

print("length of nulls is ")
print(len(nulls))


# nulls = nulls[0:10]

# nulls = ['A0A6G7EU28_9STAP']

splits = np.array_split(np.array(nulls), max(1, round(len(nulls) / 8000)))
splits = np.array_split(np.array(nulls), 1)

split_count = 0
for split in splits:
    split_count += 1
    print(f"Running stage {split_count} of {len(splits)}")
    asyncio.run(main(split))
    time.sleep(5)

# print (n2protid)

# print ('id dict')
# print (id_dict)


all_ids_dict = {**n2protid, **id_dict}

print("batman")

print("length of nulls")
print(len(nulls))

print("length of nulls to prot id")


print(len(n2protid))


print("lenght of all ids")
print(len(all_ids_dict))

print("length of original df")
print(len(original_df))

print(nulls)

print(n2protid.keys())

for x in nulls:
    if x not in n2protid.keys():
        print(x)

# sys.exit()

print("now we are adding to the IPD_search_IDS field")

# add_df = pd.DataFrame(all_ids_dict.items(), columns=['Extracted_ID', 'IPD_search_IDS'])
# original_df = pd.merge(
#         original_df, add_df
#     )

original_df = add_to_df_from_dict(original_df, all_ids_dict, "IPD_search_IDS")

print(original_df)

print("did it work?")

# original_df['IPD_search_IDS'] = np.where(df['Entry'].astype(float) >= 1000, 'OK', 'NOTOK')


# print (original_df[['Entry', 'IPD_search_IDS']])


# Add back in the specific IDs we will use to search Identical Protein Groups for genomes

# 'IPD_search_IDS'


# original_df["Binding_positions_alignment_residues"] = original_df.apply(
#     lambda row: an.get_amino_acids(
#         aln_dict[row["Entry"]].sequence, *aln_binding_pos
#     ),
#     axis=1,
# )


# Search for genome IDs

counto = 0
for chunked_dict in chunks(all_ids_dict, SIZE=800):
    # print (chunked_dict)
    assembly_dict = get_assembly_dict(chunked_dict)
    counto += 1
    print("here we go joe")
    print(assembly_dict)

    # add_assembly_df = pd.DataFrame(assembly_dict.items(), columns=['Entry', 'Genomes_with_identical_proteins'])

    # print (add_assembly_df)

    # print ('and original')
    # print (original_df)

    # original_df = original_df.merge(add_assembly_df)

    original_df = add_to_df_from_dict(
        original_df, assembly_dict, "Genomes_with_identical_proteins"
    )

    print("here it beeee")
    print(original_df.head(5))
    # original_df.to_csv(snakemake.output[0], index=False)
    # merged_df = pd.merge(
    #         original_df, add_assembly_df
    #     )
    # print (merged_df[['Entry', 'IPD_search_IDS', 'Genomes_with_identical_proteins']])


# Save the merged annotations to a csv
original_df.to_csv(snakemake.output[0], index=False)
