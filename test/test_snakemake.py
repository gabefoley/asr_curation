import os
import shutil


# def test_snakemake_pipeline():
#
#     # Move the snakefile into the test directory
#     snakefile_dest = "test/snakefile"
#
#     scripts_dest = "test/scripts"
#     snakemake_output = "test/files/test_workflow/example_workflow"
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
#     shutil.copyfile("snakefile", snakefile_dest)
#     shutil.copytree("scripts", scripts_dest)
#
#     os.system("snakemake --cores 1 --configfile test/files/config/test_config.yaml")
#
#     assert(os.path.exists("test/files/test_workflow/example_workflow/datasets/test_ec_1_1_1_86/subsets/all_kari/1/grasp_results/GRASP_ancestors.fa"))
#     assert(os.path.exists("test/files/test_workflow/example_workflow/datasets/test_ec_1_1_1_86/subsets/eukaryotic_test/1/grasp_results/GRASP_ancestors.fa"))
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

# def test_snakemake_pipeline_with_existing_annotation_file():
#     # This tests the full snakemake pipeline without going online to the databases to build the annotation files - so it should be quicker to run
#     # Move the snakefile into the test directory
#     snakefile_dest = "test/snakefile"
#
#     scripts_dest = "test/scripts"
#
#     subset_output = "test/files/test_workflow/keep_annotations/datasets/test_ec_1_1_1_86/subsets"
#     summary_output = "test/files/test_workflow/keep_annotations/datasets/test_ec_1_1_1_86/dataset_summary"
#
#
#     print (os.getcwd())
#
#     print ('got here')
#
#     print (os.path.exists(snakefile_dest))
#     print ('okay')
#
#     print (snakefile_dest)
#     print (subset_output)
#
#     print ('damn')
#
#
#
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
#     # if os.path.exists(summary_output):
#     #     shutil.rmtree(summary_output)
#
#
#
#     print("Copy snakefile and scripts folder into test directory")
#
#     shutil.copyfile("snakefile", snakefile_dest)
#     shutil.copytree("scripts", scripts_dest)
#
#     os.system("snakemake --cores 1 --configfile test/files/config/test_config_keep_annotations.yaml")
#
#     assert (os.path.exists(
#         "test/files/test_workflow/keep_annotations/datasets/test_ec_1_1_1_86/subsets/all_kari/1/grasp_results/GRASP_ancestors.fa"))
#     assert (os.path.exists(
#         "test/files/test_workflow/keep_annotations/datasets/test_ec_1_1_1_86/subsets/eukaryotic_test/1/grasp_results/GRASP_ancestors.fa"))
#


def test_snakemake_pipeline_with_existing_annotation_file(tmpdir):

    # Use tmpdir provided by pytest to create temporary directories
    test_dir = tmpdir.mkdir("test_workflow")
    snakefile_dest = os.path.join(test_dir, "snakefile")
    scripts_dest = os.path.join(test_dir, "scripts")
    subset_output = os.path.join(test_dir, "subsets")


        # Remove snakefile
        if os.path.exists(snakefile_dest):
            os.remove(snakefile_dest)

        if os.path.exists(scripts_dest):
            shutil.rmtree(scripts_dest)

        if os.path.exists(subset_output):
            shutil.rmtree(subset_output)


    # Copy Snakefile from your Git repository to the test directory
    shutil.copyfile("snakefile", snakefile_dest)

    # Copy the 'scripts' folder if needed
    shutil.copytree("scripts", scripts_dest)

    print ('snakefile dest is')
    print (snakefile_dest)

    print (os.path.exists(snakefile_dest))

    print ('okay')

    shutil.copyfile("snakefile", snakefile_dest)
    shutil.copytree("scripts", scripts_dest)

    os.system("snakemake --cores 1 --configfile test/files/config/test_config_keep_annotations.yaml")

    assert (os.path.exists(
        "test/files/test_workflow/keep_annotations/datasets/test_ec_1_1_1_86/subsets/all_kari/1/grasp_results/GRASP_ancestors.fa"))
    assert (os.path.exists(
        "test/files/test_workflow/keep_annotations/datasets/test_ec_1_1_1_86/subsets/eukaryotic_test/1/grasp_results/GRASP_ancestors.fa"))


    assert os.path.exists(os.path.join(subset_output, "all_kari/1/grasp_results/GRASP_ancestors.fa"))
    assert os.path.exists(os.path.join(subset_output, "eukaryotic_test/1/grasp_results/GRASP_ancestors.fa"))