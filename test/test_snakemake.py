# import os
# import shutil
#
#
# def test_snakemake_pipeline():
#
#     # Move the snakefile into the test directory
#     snakefile_dest = "./snakefile"
#
#     scripts_dest = "./scripts"
#     snakemake_output = "./files/test_workflow/example_workflow"
#
#     # Remove snakefile
#     if os.path.exists(snakefile_dest):
#         os.remove(snakefile_dest)
#
#     if os.path.exists(scripts_dest):
#         shutil.rmtree(scripts_dest)
#
#     if os.path.exists(snakemake_output):
#         shutil.rmtree(snakemake_output)
#
#     print ("Copy snakefile and scripts folder into test directory")
#
#     shutil.copyfile("../snakefile", snakefile_dest)
#     shutil.copytree("../scripts", scripts_dest)
#
#     os.system("snakemake --cores 1 --configfile files/config/test_config.yaml")
#
#     assert(os.path.exists("./files/test_workflow/example_workflow/datasets/test_ec_1_1_1_86/subsets/all_kari/grasp_results/GRASP_ancestors.fa"))
#     assert(os.path.exists("./files/test_workflow/example_workflow/datasets/test_ec_1_1_1_86/subsets/eukaryotic_test/grasp_results/GRASP_ancestors.fa"))
#
#
#     # snakemake_output = "./files/test_workflow/example_workflow"
#     # Remove snakefile
#     # if os.path.exists(snakefile_dest):
#     #     os.remove(snakefile_dest)
#     #
#     # shutil.rmtree(scripts_dest)
#     #
#     # shutil.rmtree(snakemake_output)
#
# def test_snakemake_pipeline_with_existing_annotation_file():
#     # This tests the full snakemake pipeline after going online to the databases to build the annotation files - so it should be quicker to run
#     # Move the snakefile into the test directory
#     snakefile_dest = "./snakefile"
#
#     scripts_dest = "./scripts"
#
#     subset_output = "./files/test_workflow/keep_annotations/datasets/test_ec_1_1_1_86/subsets"
#     summary_output = "./files/test_workflow/keep_annotations/datasets/test_ec_1_1_1_86/dataset_summary"
#
#     # Remove snakefile
#     if os.path.exists(snakefile_dest):
#         os.remove(snakefile_dest)
#
#     if os.path.exists(scripts_dest):
#         shutil.rmtree(scripts_dest)
#
#     if os.path.exists(subset_output):
#         shutil.rmtree(subset_output)
#
#     if os.path.exists(summary_output):
#         shutil.rmtree(summary_output)
#
#
#
#     print("Copy snakefile and scripts folder into test directory")
#
#     shutil.copyfile("../snakefile", snakefile_dest)
#     shutil.copytree("../scripts", scripts_dest)
#
#     os.system("snakemake --cores 1 --configfile files/config/test_config_keep_annotations.yaml")
#
#     assert (os.path.exists(
#         "./files/test_workflow/keep_annotations/datasets/test_ec_1_1_1_86/subsets/all_kari/grasp_results/GRASP_ancestors.fa"))
#     assert (os.path.exists(
#         "./files/test_workflow/keep_annotations/datasets/test_ec_1_1_1_86/subsets/eukaryotic_test/grasp_results/GRASP_ancestors.fa"))
#
