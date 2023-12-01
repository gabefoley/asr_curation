def modify_file(input_file_path, output_file_path):
    """
    Modify an input file as described in the provided code snippet.

    Args:
        input_file_path (str): Path to the input file.
        output_file_path (str): Path to the output file where the modified content will be saved.
    """
    with open(input_file_path, "r") as in_file:
        buf = in_file.readlines()

    with open(output_file_path, "w+") as out_file:
        for line in buf:
            if line.strip() == '"snakemake-job-properties"':
                line += ',\n\t"remove-cell"\n'
            out_file.write(line)

def main():
    input_file_path = snakemake.input[0]  # Replace with the actual input file path
    output_file_path = snakemake.output[0]  # Replace with the actual output file path

    modify_file(input_file_path, output_file_path)

if __name__ == "__main__":
    main()