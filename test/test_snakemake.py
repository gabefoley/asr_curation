import os
import shutil
import pytest

@pytest.fixture
def test_dir(tmpdir):
    # Create a temporary directory for testing
    test_dir = tmpdir.mkdir("test_workflow")

    # Copy the Snakefile and test/files directory to the temporary directory
    snakefile_dest = os.path.join(test_dir, "snakefile")
    scripts_dest = os.path.join(test_dir, "scripts")
    test_dest = os.path.join(test_dir, "test")
    # test_files_src = "test/files"

    shutil.copyfile("snakefile", snakefile_dest)
    shutil.copytree("scripts", scripts_dest)
    shutil.copytree("test", test_dest)

    # shutil.copytree(test_files_src, os.path.join(test_dir, "test_files"))

    os.chdir(test_dir)  # Change working directory to the test directory

    return test_dir

def test_snakemake_pipeline(test_dir):

    print("/n Current workdir")
    print (os.getcwd())
    print(os.listdir())

    os.system("snakemake --cores 1 --configfile test/files/config/test_config.yaml")

    assert os.path.exists(os.path.join(test_dir, "test/files/test_workflow/example_workflow/datasets/test_ec_1_1_1_86/subsets/all_kari/1/grasp_results/GRASP_ancestors.fa"))
