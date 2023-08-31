import seqcurate as sc
import pandas as pd
import numpy as np


def format_subset_val(col, col_val):
    if type(col_val) == str:
        # Check if it is a list of values
        if "," in col_val:
            return f"{col}.isin({[val.strip() for val in col_val.split(',')]})"

        # If it is an equality return it as such
        if col_val.strip()[0] in [">", "<", "="]:
            return f"{col.strip()} {col_val.strip()}"

        # If it is a string encase it in quotation marks
        # return f"{col.strip()} == '{col_val.strip()}'"
        return f"{col.strip()}.str.contains('{col_val.strip()}')"

    else:
        print (col.strip())
        return f"{col.strip()} == {col_val}"


def subset_column_vals(df, col_val_dict, not_col_val_dict, req_col_val_dict):
    qry = " and ".join(
        [format_subset_val(col, col_val) for col, col_val in col_val_dict.items()]
    )

    # sub_df = df.loc[(df[list(col_val_dict)] == pd.Series(col_val_dict)).all(axis=1)]

    print(f"Creating a subset with {qry}")

    # qry = "gene_primary.str.contains('ilvC', na=False)"
    # qry = "Length_2.str.contains('None')"


    print (qry)

    print ('hello')


    sub_df = df.query(qry)

    print ('got here')

    print (sub_df)

    sub_df = sub_df.fillna("None")

    # print(f"Subset length is {len(sub_df)}")

    for col, col_val in not_col_val_dict.items():

        print (col)
        # If it is a list split it up
        for val in col_val.split(","):
            print (val)
            sub_df = sub_df[~sub_df[col].str.contains(val.strip(), na=False)]



    print ('now')
    if req_col_val_dict:
        req_qry = " or ".join(
            [
                format_subset_val(col, col_val)
                for col, col_val in req_col_val_dict.items()
            ]
        )

        print (sub_df)

        req_df = df.query(req_qry)

        frames = [sub_df, req_df]
        merge_df = pd.concat(frames)

        merge_df = merge_df.drop_duplicates(subset="accession")

        return merge_df
    else:
        print ('sub df')
        print (sub_df)
        sub_df = sub_df.drop_duplicates(subset="accession")
        return sub_df


def get_col_val_name(col_val_dict, not_col_val_dict):
    name = "_".join(str(k) + "_" + str(v) for k, v in col_val_dict.items())

    name += "_".join(str(k) + "_" + str(v) for k, v in not_col_val_dict.items())

    return name


def add_val_to_dict(term, add_val, add_dict):
    # If the value is a true or false, turn it into a boolean

    if add_val.lower() in ["true", "false"]:
        add_dict[term] = bool(add_val.lower() == "true")

    else:
        add_dict[term] = add_val

        # # Try and add it as an int
        # try:
        #     add_dict[term] = int(add_val)
        # except ValueError:
        #     add_dict[term] = add_val


# If a subset rule doesn't exist, then just write out the full data set
if snakemake.wildcards.subset == "full_set":
    sc.write_to_fasta(pd.read_csv(snakemake.input.csv), snakemake.output[0])


for line in open(snakemake.input.rules).read().splitlines():
    # Reset the condition that we want to sample from the dataframe
    sample_from = None

    if line and not line.startswith("#"):
        col_val_dict = {}
        not_col_val_dict = {}
        req_col_val_dict = {}

        # Check to see if custom name exists
        name = line.split("=")[0].strip()

        # Make the dictionary
        # col_val_dict = dict(line.split("=")[1].strip())

        dict_def = line.split("=")[1].strip()

        # Wildcard for accepting all columns
        if dict_def != "*":
            for i in dict_def.split("$"):
                term = i.split(":")[0].strip()

                if term:
                    print (i)
                    term_val = i.split(":")[1].strip()

                    print (term_val)

                    # If it is instructions on how to sample from columns create the conditions

                    if term.split(" ")[0].strip().lower() == "sample_from":
                        sample_num = term.split("(")[1].split(")")[0]
                        sample_from = [x.strip() for x in term_val.split(",")]

                    else:
                        # If the value is a negation add it to the negation dictionary
                        if term_val.split(" ")[0].strip().lower() == "not":
                            add_val_to_dict(
                                term, "".join(term_val.split(" ")[1:]), not_col_val_dict
                            )

                        # If the value is required to be there add it to the required dictionary
                        elif term_val.split(" ")[0].strip().lower() == "required":
                            add_val_to_dict(
                                term, term_val.split(" ")[1], req_col_val_dict
                            )

                        # Else add the term value we want
                        else:
                            add_val_to_dict(term, term_val.split(" ")[0], col_val_dict)

        # col_val_dict = {i.split(':')[0].strip(): bool(i.split(':')[1].strip().lower() == 'true') if i.split(':')[1].strip().lower() in ['true', 'false'] else i.split(':')[1].strip() for i in dict_def.split(",")}

        # If no custom name make name based on values of dictionary
        if len(name) < 3:
            name = get_col_val_name(col_val_dict, not_col_val_dict)

        if name == snakemake.wildcards.subset:
            df = pd.read_csv(snakemake.input.csv)

            # df = df.fillna("None")

            # Subset the annotation file
            if dict_def == "*":
                sub_df = df
            else:
                sub_df = subset_column_vals(
                    df, col_val_dict, not_col_val_dict, req_col_val_dict
                )

            if sample_from:
                print(
                    f"We are sampling {sample_num} entries from each of the following columns - {','.join([x for x in sample_from])} "
                )

                frames = []

                for sample in sample_from:
                    frames.append(
                        sub_df.loc[sub_df[sample] == True].sample(int(sample_num))
                    )

                sub_df = pd.concat(frames)

                print("After sampling the length to write out is")
                print(len(sub_df))

                # Remove this from above and just do it here?

                req_qry = " or ".join(
                    [
                        format_subset_val(col, col_val)
                        for col, col_val in req_col_val_dict.items()
                    ]
                )

                if req_qry:
                    req_df = df.query(req_qry)

                    frames = [sub_df, req_df]
                    merge_df = pd.concat(frames)

                    sub_df = merge_df.drop_duplicates(subset="accession")

            # Write the subset to a fasta
            # TODO: Currently this is just working for one subset file (not multiple)

            sub_df = sub_df.replace("None", np.NaN)

            print (sub_df)
            sc.write_to_fasta(sub_df, snakemake.output.fasta, trim=True)


            print ('write to csv')

            # Write the subset to its own csv file
            sub_df.to_csv(snakemake.output.csv, index=False)
