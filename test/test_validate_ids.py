# import os
# import scripts.validate_ids as vi
# from click.testing import CliRunner
#
# # Tests to check the validation of IDs
# def test_validate_ids_main():
#
#     outpath = "./files/output_test_validated_main.csv"
#
#     if os.path.exists(outpath):
#         os.remove(outpath)
#
#     vi.all_ids_lookup("./files/test_original.csv", outpath)
#
#     assert (os.path.isfile(outpath))
#
#
# def test_validate_ids_main_from_command_line():
#
#     outpath = "./files/output_test_validated_command_line.csv"
#
#     if os.path.exists(outpath):
#         os.remove(outpath)
#
#     os.system(f"python ../scripts/validate_ids.py --input_file ./files/test_original.csv --output_file {outpath}")
#
#     assert (os.path.isfile(outpath))
#
# def test_validate_ids_call_command_line_directly():
#
#     outpath = "./files/output_test_validated_command_line_direct.csv"
#
#     if os.path.exists(outpath):
#         os.remove(outpath)
#
#     runner = CliRunner()
#
#
#     # vi.all_ids_lookup_cmd(["--input_file", "./files/test_original.csv", "--output_file", outpath])
#
#     result = runner.invoke(vi.all_ids_lookup_cmd, ["--input_file", "./files/test_original.csv", "--output_file", outpath])
#
#     assert result.exit_code == 0
#
#     assert (os.path.isfile(outpath))
#
#
# def test_creating_and_validating_diverse_headers():
#     # This test creates a DataFrame that contains full NCBI accession header, full UniProt header, gi header
#     # (as downloaded from Conserved Domain Database alignment) and single NCBI identifier and single UniProt identifier
#
#
#     seq_outpath = "./files/output_test_diverse_5.csv"
#     validated_outpath = "./files/output_test_diverse_5_validated.csv"
#
#     if os.path.exists(seq_outpath):
#         os.remove(seq_outpath)
#
#     if os.path.exists(validated_outpath):
#         os.remove(validated_outpath)
#
#     # Create a sequence df
#
#     #
#
#
