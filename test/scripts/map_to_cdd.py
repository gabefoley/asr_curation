import requests
import time
import pandas as pd


def map_to_batch_cdd_search(ids, df):
    # URL to the Batch CD-Search server
    bwrpsb = "https://www.ncbi.nlm.nih.gov/Structure/bwrpsb/bwrpsb.cgi"

    # Parameters
    cddefl = "false"
    qdefl = "false"
    useid1 = "true"
    maxhit = 250
    filter = "true"
    db = "cdd"
    evalue = 0.01
    dmode = "rep"
    clonly = "false"
    tdata = "hits"

    chunk_size = 900

    # Chunk the IDs into groups of three sequence IDs
    chunked_ids = [ids[i : i + chunk_size] for i in range(0, len(ids), chunk_size)]

    # Process each chunk separately
    for chunk in chunked_ids:
        # Submitting the search
        rid = None

        params = {
            "useid1": useid1,
            "maxhit": maxhit,
            "filter": filter,
            "db": db,
            "evalue": evalue,
            "cddefl": cddefl,
            "qdefl": qdefl,
            "dmode": dmode,
            "clonly": clonly,
            "tdata": "hits",
            "queries": "\n".join(chunk),
        }

        response = requests.post(bwrpsb, data=params)
        response.raise_for_status()

        lines = response.text.split("\n")

        for line in lines:
            if line.startswith("#cdsid"):
                rid = line.split()[1]
                print(f"Search with Request-ID {rid} started.")
                break

        if rid is None:
            raise ValueError(
                "Submitting the search failed. Can't retrieve the Request-ID."
            )

        # Checking for completion, wait 5 seconds between checks
        done = False

        while not done:
            time.sleep(5)

            params = {"tdata": "hits", "cdsid": rid}

            response = requests.post(bwrpsb, data=params)
            response.raise_for_status()

            lines = response.text.split("\n")

            for line in lines:
                if line.startswith("#status"):
                    status = int(line.split("\t")[1])
                    if status == 0:
                        done = True
                        print("Search has been completed, retrieving results...")
                        break
                    elif status == 3:
                        message = line.split("\t")[3]
                        print(f"Job is still running. Message: {message}")
                    else:
                        raise ValueError(f"Search status check failed. Status: {line}")

        print(
            "\n===============================================================================\n\n"
        )

        # Retrieve and display results
        results = []
        header_lines = []

        for line in lines:
            if line.startswith("#"):
                header_lines.append(line)
            else:
                results.append(line)

        # Determine the number of columns dynamically from the headers
        num_columns = len(results[1].split("\t"))

        # Create DataFrame from the results
        df_results = pd.DataFrame(
            [row.split("\t")[:num_columns] for row in results if row]
        )
        df_results.columns = df_results.iloc[0]
        #         df_results = df_results[1:]

        df_results.columns = [
            f'cdd_{x.replace(" ", "_").lower()}' for x in df_results.iloc[0]
        ]
        df_results.rename(columns={"cdd_query": "info"}, inplace=True)
        df_results["info"] = df_results["info"].str.replace(">", "")

        print(df_results)

        if df_results.empty:
            raise ValueError("No valid results found. The result may be incomplete.")

        # Extract the 'Query' and relevant columns
        df_results["info"] = df_results["info"].str.split(" - ").str[-1].str.strip()
        df_results["cdd_short_name"] = df_results["cdd_short_name"].str.strip()

        # Merge the retrieved attributes with the original DataFrame
        df = df.merge(
            df_results[
                [
                    "info",
                    "cdd_bitscore",
                    "cdd_accession",
                    "cdd_short_name",
                    "cdd_incomplete",
                    "cdd_superfamily",
                ]
            ],
            on="info",
            how="left",
        )

    return df
