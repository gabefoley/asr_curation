import seqcurate as sc
import pandas as pd



def format_subset_val(col, col_val):

    if type(col_val) == str:
        # Check if it is a list of values
        if "," in col_val:
            return f"{col}.isin({[val.strip() for val in col_val.split(',')]})"

        # If it is an equality return it as such
        if col_val.strip()[0] in [">", "<", "="]:
            return f"{col.strip()} {col_val.strip()}"

        # If it is a string encase it in quotation marks
        #return f"{col.strip()} == '{col_val.strip()}'"
        return f"{col.strip()}.str.contains('{col_val.strip()}')"

    else:
        return f"{col.strip()} == {col_val}"


def subset_column_vals(df, col_val_dict, not_col_val_dict, req_col_val_dict):

    print(list(col_val_dict))

    qry = " and ".join(
        [format_subset_val(col, col_val) for col, col_val in col_val_dict.items()]
    )

    # sub_df = df.loc[(df[list(col_val_dict)] == pd.Series(col_val_dict)).all(axis=1)]

    print(f"Creating a subset with {qry}")

    sub_df = df.query(qry)

    print(f"Subset length is {len(sub_df)}")

    for col,col_val in not_col_val_dict.items():

        #If it is a list split it up
        for val in col_val.split(","):
            print (val)
            sub_df = sub_df[~sub_df[col].str.contains(val.strip(), na=False)]
            print (sub_df.head())


    if req_col_val_dict:
        #print("Something is required")

        #print(req_col_val_dict)

        req_qry = " or ".join(
            [
                format_subset_val(col, col_val)
                for col, col_val in req_col_val_dict.items()
            ]
        )

        req_df = df.query(req_qry)

        #print(len(req_df))

        frames = [sub_df, req_df]
        merge_df = pd.concat(frames)

        #print(merge_df)

        merge_df = merge_df.drop_duplicates(subset="Entry")

        return merge_df
    else:
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

        #print("line is ")
        #print(line)

        col_val_dict = {}
        not_col_val_dict = {}
        req_col_val_dict = {}

        # Check to see if custom name exists
        name = line.split("=")[0].strip()

        #print(len(name))

        #print(line.split("=")[1].strip())

        # Make the dictionary
        # col_val_dict = dict(line.split("=")[1].strip())

        dict_def = line.split("=")[1].strip()

        for i in dict_def.split("$"):

            term = i.split(":")[0].strip()

            if term:

                term_val = i.split(":")[1].strip()

                # If it is instructions on how to sample from columns create the conditions

                if term.split(" ")[0].strip().lower() == "sample_from":
                    sample_num = term.split("(")[1].split(")")[0]
                    sample_from = [x.strip() for x in term_val.split(",")]

                else:

                    print(i)
                    print(term)

                    print(term_val)

                    # If the value is a negation add it to the negation dictionary
                    if term_val.split(" ")[0].strip().lower() == "not":
                        print("was a negation")
                        add_val_to_dict(term, "".join(term_val.split(" ")[1:]), not_col_val_dict)

                    # If the value is required to be there add it to the required dictionary
                    elif term_val.split(" ")[0].strip().lower() == "required":
                        print("required ")
                        print(term_val)
                        print(term_val.split(" ")[1])
                        add_val_to_dict(term, term_val.split(" ")[1], req_col_val_dict)

                    # Else add the term value we want
                    else:
                        add_val_to_dict(term, term_val.split(" ")[0], col_val_dict)

        # col_val_dict = {i.split(':')[0].strip(): bool(i.split(':')[1].strip().lower() == 'true') if i.split(':')[1].strip().lower() in ['true', 'false'] else i.split(':')[1].strip() for i in dict_def.split(",")}

        #print("CVD")
        #print(col_val_dict)

        #print("not CVD")
        #print(not_col_val_dict)

        #print("Required to be there")
        #print(req_col_val_dict)

        # If no custom name make name based on values of dictionary
        if len(name) < 3:
            name = get_col_val_name(col_val_dict, not_col_val_dict)

        if name == snakemake.wildcards.subset:

            print("Name given to this subset is " + name)

            print("wildcards are ")
            print(snakemake.wildcards)

            df = pd.read_csv(snakemake.input.csv)

            # Subset the annotation file
            sub_df = subset_column_vals(
                df, col_val_dict, not_col_val_dict, req_col_val_dict
            )

            # print ('sub df is ')
            # print (sub_df.head(1))

            print("Length to write out is")

            print(len(sub_df))

            if sample_from:
                print(
                    f"We are sampling {sample_num} entries from each of the following columns - {','.join([x for x in sample_from])} "
                )

                frames = []

                for sample in sample_from:
                    print(sample)
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

                    print(merge_df)

                    sub_df = merge_df.drop_duplicates(subset="Entry")

            # Write the subset to a fasta
            # TODO: Currently this is just working for one subset file (not multiple)
            sc.write_to_fasta(sub_df, snakemake.output.fasta, trim=True)

            # Write the subset to its own csv file
            sub_df.to_csv(snakemake.output.csv, index=False)
