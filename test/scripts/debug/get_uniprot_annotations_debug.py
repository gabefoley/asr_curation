import requests
import io
import sequence
from sequence import Alphabet
import numpy as np
import pandas as pd
import seqcurate as sc
import logger
from uniprot_cols import full_uniprot_cols

# Number of times to retry UniProt database
MAX_RETRIES = 20

# The following is used for debugging purposes to check if a column is returning correctly


def getUniProtDict_OneCol(ids, cols="", db="uniprot", identities=None):
    """

    :param ids: The list of UniProt IDs
    :param cols: The list of UniProt database columns
    :param db: The database to search - uniprot or uniref
    :param identity: The identity to search uniref with
    :return: A dictionary mapping each UniProt ID to another dictionary where the keys are database columns and the
    values are the information stored within those columns
    """

    """
    *** EXAMPLE USAGE ***
    Get a list of UniProt IDs and a list of UniProt columns you're interested in.
    Full list of UniProt column names - https://www.uniprot.org/help/uniprotkb_column_names

    uniprot_names = ['Q9LIR4', 'Q1JUQ1', 'P05791', 'P0ADF6']
    cols = ["lineage(SUPERKINGDOM)", "genes", "lineage(KINGDOM)"]    
    up_dict = getUniProtDict(uniprot_names, cols)

    for record in up_dict:
        print (record, up_dict[record].get("lineage(SUPERKINGDOM)"))

    print()

    for record in up_dict:
        print (record, up_dict[record].get("genes"))


    If a record doesn't have an entry in UniProt for that column it'll just return None

    print (up_dict['Q1JUQ1'])
    print (up_dict['Q1JUQ1']['lineage(KINGDOM)'])


    *** EXAMPLE USAGE FOR UNIREF SEARCHING ***

    up_dict = getUniProtDict(["Q9LIR4", "P99999"], cols=["members"], db="uniref", identities = 1.0)

    You can either pass a list of identities for each UniProt identifier (in which case the list of identities must be
    the same size as the list of identifiers. Or you can just pass a single identity to search Uniref at.
    """

    up_dict = {}

    # Format the lists of IDs and columns correctly
    cols = ",".join(cols)

    if db == "uniprot":
        updated_ids = " or ".join(ids)

    # elif db == 'uniref':
    #
    #     # If we just got a single value for an identity, create a list the same size as the queries with just this value
    #     if type(identities) != list:
    #         identities = [identities] * len(chunk)
    #     elif len(identities) != len(chunk):
    #         raise RuntimeError(
    #             'Either supply a single identity threshold or supply one for each identifier in the list')
    #
    #     # Check that the identity thresholds are valid values
    #     for x in identities:
    #         if x not in [1.0, 0.9, 0.5]:
    #             raise RuntimeError(
    #                 "UniRef threshold values must be either 1.0, 0.9, or 0.5. Supplied value was - " + str(x))
    #
    #     # Add the query syntax around the identifiers
    #     updated_ids = ""
    #     for query, identity in zip(ids, identities):
    #         updated_ids += "(member:" + query + "+AND+identity:" + str(identity) + ")+OR+"
    #
    #     updated_ids = updated_ids[0:-4]
    #

    # else:
    #     raise RuntimeError('Database should be either uniprot or uniref')

    url = "https://www.uniprot.org/" + db + "/"

    for col in cols.split(","):

        print(f"trying col - {col}")

        params = {"format": "tab", "query": updated_ids, "columns": "id," + col}

        #     print (url)
        print(updated_ids)

        data = urllib.parse.urlencode(params).encode("utf-8")
        request = urllib.request.Request(url, data)

        with urllib.request.urlopen(request) as response:

            page = response.read(200000000).decode("utf-8")

            print(f'Length of cols is {len(cols.split(","))}')

            # print(page)

            # For each record we retrieve, split the line by tabs and build up the UniProt dict
            for line in page.split("\n")[1:]:
                #             print (line)
                if line:
                    splitlines = line.split("\t")
                    id_dict = {}
                    pos = 1
                    print(f"splitlines \n")
                    print(splitlines)
                    print(f"THE LENGTH OF THE RETURN VALUE IS {len(splitlines)}")

                    if splitlines and len(splitlines) > 1:
                        # for col in cols.split(","):
                        print(f"col is {col}")
                        print(f"pos is {pos}")
                        id_dict[col] = (
                            None if splitlines[pos] == "" else splitlines[pos]
                        )
                        pos += 1
                    up_dict[splitlines[0]] = id_dict

    print("here")
    print(len(up_dict))
    return up_dict
