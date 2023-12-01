from scripts.clean_summary_document import modify_file

def test_modify_file(tmpdir):
    # Create a temporary input file with sample data
    input_file = tmpdir.join("input.txt")
    input_content = 'Some text\n"snakemake-job-properties"\nMore text'
    input_file.write(input_content)

    # Create a temporary output file
    output_file = tmpdir.join("output.txt")

    # Modify the input file using the code snippet
    modify_file(input_file, output_file)

    # Read the modified content from the output file
    with open(output_file, "r") as modified_file:
        modified_content = modified_file.read()

    # Check if the modification was done correctly
    expected_output = 'Some text\n"snakemake-job-properties"\n,\n\t"remove-cell"\nMore text'
    assert modified_content == expected_output
