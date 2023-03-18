import seqcurate as sc
import os
import click


@click.command()
@click.option('--input_file', help='Path to input file')
@click.option('--output_file', help='Path for output file')
def create_initial_annotation_cmd(input_file,output_file):
    create_initial_annotation(input_file,output_file)

def create_initial_annotation(input_file, output_file):

    if 'SNAKEMAKE' in os.environ:
        input_file = snakemake.input[0]
        output_file = snakemake.output[0]


    # Create the sequence dataframe and save it
    annotations = sc.get_sequence_df(input_file)
    annotations.to_csv(output_file, index=False)


# Main function
def main():
    print("Creating the initial annotation file")
    create_initial_annotation_cmd()


if __name__ == "__main__":
    main()