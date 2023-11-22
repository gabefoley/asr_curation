import seqcurate as sc


def create_initial_annotation():
    input_file = snakemake.input[0]
    output_file = snakemake.output[0]

    # Create the sequence dataframe and save it
    annotations = sc.get_sequence_df(input_file)
    annotations.to_csv(output_file, index=False)


if __name__ == "__main__":
    create_initial_annotation()
