import seqcurate as sc


def create_initial_annotation(input_file, output_file, verbose=True):
    # Create the sequence dataframe and save it
    if verbose:
        print (f'Creating an initial annotation dataframe from {input_file}')

    annotations = sc.get_sequence_df(input_file)
    annotations.to_csv(output_file, index=False)

if __name__ == "__main__":
    input_file = snakemake.input[0]
    output_file = snakemake.output[0]
    create_initial_annotation(input_file, output_file, verbose=snakemake.params.verbose)
